"""
rgi_io.py -- peak-list readers for rgi_msn.py.  Standard library only.

Handles the Thermo FreeStyle "SPECTRUM - MS" CSV export, which is what comes out
of Xcalibur/FreeStyle when you export a spectrum, and which is NOT a peak list:

  * a BOM plus an 8-line metadata header (raw path, scan filter, scans, RT, AV,
    data points, column names), with a variable number of trailing commas;
  * either centroided peaks or, more often, **profile** data -- tens of thousands
    of points spaced ~0.0001-0.003 Da, over half of them exactly zero. Handing
    that to a fragment matcher is meaningless: every theoretical mass finds a
    point, and the match loop takes minutes.

`read_spectrum` returns centroided peaks either way, plus the metadata, and
parses the scan filter so the MS^n precursor path does not have to be typed in:

    FTMS + p ESI Full ms3 1571.6100@cid30.00 1087.4200@cid30.00 [295.0-1098.0]
    -> ms_level 3, precursors [1571.61, 1087.42], activations ['cid','hcd'...],
       polarity '+'
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

FILTER_RE = re.compile(r"\bms(\d*)\b", re.I)
SCANFILTER_RE = re.compile(r"\b(?:ft|it)ms\b|\bfull\s+ms\d?\b", re.I)
PREC_RE = re.compile(r"(\d+\.\d+)@([a-zA-Z]+)([\d.]+)")


@dataclass
class Spectrum:
    peaks: List[Tuple[float, float]] = field(default_factory=list)
    ms_level: Optional[int] = None
    precursors: List[float] = field(default_factory=list)
    activations: List[str] = field(default_factory=list)
    polarity: Optional[str] = None
    scan_filter: str = ""
    raw_file: str = ""
    source_mode: str = ""          # 'centroid' | 'profile'
    n_raw: int = 0
    path: str = ""

    @property
    def mz(self) -> List[float]:
        return [m for m, _ in self.peaks]

    def adduct_hint(self) -> str:
        return "-H" if self.polarity == "-" else "Na"

    def top(self, n: int) -> List[Tuple[float, float]]:
        return sorted(self.peaks, key=lambda p: -p[1])[:n]

    def __repr__(self) -> str:
        return (f"<Spectrum MS{self.ms_level or '?'} {self.polarity or '?'} "
                f"precursors={self.precursors} {len(self.peaks)} peaks "
                f"from {self.n_raw} {self.source_mode} points>")


def parse_filter(line: str) -> dict:
    out: dict = {"scan_filter": line.strip()}
    m = re.search(r"\bms(\d*)\b", line[line.lower().find("full"):] or line, re.I) \
        if "full" in line.lower() else FILTER_RE.search(line)
    if m:
        out["ms_level"] = int(m.group(1)) if m.group(1) else 1
    if " + " in line:
        out["polarity"] = "+"
    elif " - " in line:
        out["polarity"] = "-"
    hits = PREC_RE.findall(line)
    out["precursors"] = [float(h[0]) for h in hits]
    out["activations"] = [h[1].lower() for h in hits]
    return out


def _rows(path: str):
    with open(path, "r", encoding="utf-8-sig", errors="replace", newline="") as fh:
        for row in csv.reader(fh):
            yield [c.strip() for c in row]


def centroid(mz: Sequence[float], inten: Sequence[float],
             rel_threshold: float = 0.0005,
             max_width: int = 12) -> List[Tuple[float, float]]:
    """Intensity-weighted centroid of every local maximum in profile data.

    A profile peak spans several points; its apex is a local maximum and its
    true m/z is the intensity-weighted mean over the points either side, taken
    while the signal is still falling away from the apex. Peaks below
    `rel_threshold` of the base peak are dropped as noise.
    """
    n = len(mz)
    if n == 0:
        return []
    cut = max(inten) * rel_threshold
    out: List[Tuple[float, float]] = []
    for k in range(n):
        v = inten[k]
        if v <= cut:
            continue
        if not (v >= inten[k - 1] if k else True):
            continue
        if not (v >= inten[k + 1] if k + 1 < n else True):
            continue
        lo = k
        while (lo > 0 and k - lo < max_width and inten[lo - 1] > 0
               and inten[lo - 1] <= inten[lo]):
            lo -= 1
        hi = k
        while (hi + 1 < n and hi - k < max_width and inten[hi + 1] > 0
               and inten[hi + 1] <= inten[hi]):
            hi += 1
        wsum = sum(inten[lo:hi + 1])
        if wsum <= 0:
            continue
        cm = sum(mz[i] * inten[i] for i in range(lo, hi + 1)) / wsum
        out.append((cm, v))
    return out


def read_spectrum(path: str, rel_threshold: float = 0.0005,
                  force_mode: Optional[str] = None) -> Spectrum:
    """Read a FreeStyle CSV (or a bare two-column peak list) into a Spectrum."""
    sp = Spectrum(path=path)
    mz: List[float] = []
    inten: List[float] = []
    for row in _rows(path):
        if not row:
            continue
        cells = [c for c in row if c != ""]
        if not cells:
            continue
        head = cells[0]
        if len(cells) >= 2:
            try:
                m = float(cells[0])
                v = float(cells[1])
            except ValueError:
                pass
            else:
                mz.append(m)
                inten.append(v)
                continue
        low = head.lower()
        # Only a real scan filter, never the "SPECTRUM - MS" banner, whose
        # trailing "MS" matches a loose \bms\b and whose dash reads as negative
        # polarity. Require an analyser token or an explicit "Full msN".
        if SCANFILTER_RE.search(low):
            if not sp.scan_filter:
                info = parse_filter(head)
                sp.scan_filter = info["scan_filter"]
                sp.ms_level = info.get("ms_level")
                sp.polarity = info.get("polarity")
                sp.precursors = info.get("precursors", [])
                sp.activations = info.get("activations", [])
        elif low.endswith(".raw") or head[1:3] == ":\\":
            sp.raw_file = head
    sp.n_raw = len(mz)
    if not mz:
        return sp
    nz = sum(1 for v in inten if v > 0)
    gaps = sorted(b - a for a, b in zip(mz, mz[1:]))
    med = gaps[len(gaps) // 2] if gaps else 1.0
    mode = force_mode or ("profile" if (med < 0.01 or nz / len(mz) < 0.6)
                          else "centroid")
    sp.source_mode = mode
    if mode == "profile":
        sp.peaks = centroid(mz, inten, rel_threshold)
    else:
        cut = max(inten) * rel_threshold
        sp.peaks = [(m, v) for m, v in zip(mz, inten) if v > cut]
    return sp


C13 = 1.0033548      # 13C - 12C
NEUTRAL_LOSSES = {               # common losses from per-acylated RG-I ions
    "d3-AcOH": 63.0400,          # CD3COOH, the dominant loss here
    "AcOH": 60.0211,
    "ketene": 42.0106,
    "H2O": 18.0106,
}


def deisotope(peaks, tol_mmu: float = 8.0, max_charge: int = 2,
              max_iso: int = 3):
    """Drop 13C isotopologues, keeping the monoisotopic peak of each cluster.

    In a real Orbitrap spectrum every genuine fragment appears two or three
    times, one per 13C. Left in, they are a third of the strong peaks: they
    dilute intensity-weighted scoring and invite spurious assignments.
    A peak is an isotope if a MORE intense peak sits 1.00335/z below it.
    """
    ps = sorted(peaks)
    keep = []
    tol = tol_mmu / 1000.0
    for i, (m, v) in enumerate(ps):
        iso = False
        for z in range(1, max_charge + 1):
            for k in range(1, max_iso + 1):
                target = m - k * C13 / z
                for mm, vv in ps[max(0, i - 400):i]:
                    if abs(mm - target) <= tol and vv > v:
                        iso = True
                        break
                if iso:
                    break
            if iso:
                break
        if not iso:
            keep.append((m, v))
    return keep


def estimate_offset(peaks, reference_mz, tol_mmu: float = 25.0):
    """Median m/z error of the strongest peaks against known reference masses.

    Instrument calibration drifts; on this data every published fragment sits
    +2 to +4 mmu high. A constant offset is worth removing before matching,
    because otherwise the tolerance has to be widened to cover it and the
    candidate count rises with the square of the window.
    """
    errs = []
    for c in reference_mz:
        near = [(abs(m - c), m) for m, _ in peaks if abs(m - c) * 1000 <= tol_mmu]
        if near:
            errs.append(min(near)[1] - c)
    if not errs:
        return 0.0
    errs.sort()
    n = len(errs)
    return errs[n // 2] if n % 2 else (errs[n // 2 - 1] + errs[n // 2]) / 2


def apply_offset(peaks, offset: float):
    """Shift a peak list by a constant m/z offset (subtracts the bias)."""
    return [(m - offset, v) for m, v in peaks]
