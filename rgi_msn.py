"""
rgi_msn.py -- MS^n interpreter for per-acylated, partially O-acetylated RG-I oligomers.

Chemistry model follows Fernando, Yang, Archer-Hartmann, Shahin, Zhang, Urbanowicz,
Heiss & Azadi, J. Am. Soc. Mass Spectrom. 2026, 37, 1160-1172 (doi:10.1021/jasms.5c00440):

  * RG-I backbone  [->2)-alpha-L-Rha-(1->4)-alpha-D-GalA-(1->]_n
  * lyase digestion leaves a 4,5-unsaturated GalA (dGalA) at the non-reducing end
  * reducing end is Rha
  * every free hydroxyl is blocked by trideuteroacetyl (dAc, +45.0294) or
    propionyl (Prop, +56.0262); native O-acetyl groups (Ac, +42.0106) stay put
  * optional single Gal branch on Rha O-4 (so the Rha becomes ->2,4-linked)

Fragment masses are computed atom-by-atom on an explicit molecular graph, so the
mass of any cross-ring ion automatically reflects which acyl group sits on which
carbon -- that is what distinguishes 3-O- from 2-O- or 4-O-acetylation.

Cleavage conventions
--------------------
Glycosidic bond  C1(donor)-O(acceptor)-Cp(acceptor) carries two cleavable bonds:
  C1-O  ->  B (donor side, -H)  /  Y (acceptor side, +H)
  O-Cp  ->  C (donor side, +H)  /  Z (acceptor side, -H)
Ring bonds are numbered Domon-Costello style:
  0: O5-C1   1: C1-C2   2: C2-C3   3: C3-C4   4: C4-C5   5: C5-O5
A cross-ring cleavage cuts two ring bonds; the piece retaining C1 of that residue
is the X-type fragment, the complementary piece is the A-type fragment.  No
hydrogen is transferred across a ring cleavage (retro-Diels-Alder convention);
optional +/-H variants can be requested at match time.
"""

from __future__ import annotations

import itertools
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

# ---------------------------------------------------------------- constants --

ELEMENTS = {
    "C": 12.0,
    "H": 1.00782503207,
    "D": 2.01410177785,
    "O": 15.9949146196,
    "N": 14.0030740048,
}
ELECTRON = 0.00054857990946
ADDUCTS = {
    "Na": (22.9897692809 - ELECTRON, +1),
    "H": (1.00782503207 - ELECTRON, +1),
    "K": (38.96370668 - ELECTRON, +1),
    "-H": (-1.00782503207 + ELECTRON, -1),
}

# O-substituent groups, rooted at the oxygen that is bonded to the ring carbon.
GROUPS: Dict[str, dict] = {
    "OH": dict(atoms=["O", "H"], bonds=[(0, 1)]),
    "Ac": dict(atoms=["O", "C", "O", "C", "H", "H", "H"],
               bonds=[(0, 1), (1, 2), (1, 3), (3, 4), (3, 5), (3, 6)]),
    "dAc": dict(atoms=["O", "C", "O", "C", "D", "D", "D"],
                bonds=[(0, 1), (1, 2), (1, 3), (3, 4), (3, 5), (3, 6)]),
    "Prop": dict(atoms=["O", "C", "O", "C", "H", "H", "C", "H", "H", "H"],
                 bonds=[(0, 1), (1, 2), (1, 3), (3, 4), (3, 5), (3, 6),
                        (6, 7), (6, 8), (6, 9)]),
    "Me": dict(atoms=["O", "C", "H", "H", "H"],
               bonds=[(0, 1), (1, 2), (1, 3), (1, 4)]),
}
ACYLS = ("Ac", "dAc", "Prop")

# Residue templates.  `slots` = ring positions that can bear an O-substituent.
# `link_pos` = position at which this residue accepts the next residue toward
# the non-reducing end (i.e. the backbone attachment point).
RESIDUES: Dict[str, dict] = {
    "Rha":   dict(slots=(1, 2, 3, 4), link_pos=2, c6="CH3",  unsat=False),
    "GalA":  dict(slots=(1, 2, 3, 4), link_pos=4, c6="COOH", unsat=False),
    "dGalA": dict(slots=(1, 2, 3),    link_pos=4, c6="COOH", unsat=True),
    "Gal":   dict(slots=(1, 2, 3, 4, 6), link_pos=4, c6="CH2OH", unsat=False),
}

# ------------------------------------------------------------------- graph ---


@dataclass(frozen=True)
class Cut:
    """A cleavable bond."""
    kind: str           # 'ring' | 'gly'
    site: object        # ring: residue key; gly: linkage key
    idx: int            # ring: bond 0..5; gly: 0 = C1-O (B/Y), 1 = O-Cp (C/Z)
    a: int
    b: int
    oxy: int = -1       # glycosidic oxygen atom id (gly cuts only)


class Mol:
    def __init__(self) -> None:
        self.el: List[str] = []
        self.adj: Dict[int, Set[int]] = defaultdict(set)
        self.cuts: List[Cut] = []
        self.res_c1: Dict[object, int] = {}      # residue key -> C1 atom id
        self.res_type: Dict[object, str] = {}
        self.res_atoms: Dict[object, Set[int]] = defaultdict(set)
        self.n_backbone: int = 0

    def add(self, element: str) -> int:
        self.el.append(element)
        return len(self.el) - 1

    def bond(self, a: int, b: int) -> None:
        self.adj[a].add(b)
        self.adj[b].add(a)

    def mass(self, atoms: Iterable[int]) -> float:
        return sum(ELEMENTS[self.el[a]] for a in atoms)

    @property
    def all_atoms(self) -> FrozenSet[int]:
        return frozenset(range(len(self.el)))


def _attach_group(mol: Mol, carbon: int, group: str, owner) -> int:
    spec = GROUPS[group]
    ids = [mol.add(e) for e in spec["atoms"]]
    for i, j in spec["bonds"]:
        mol.bond(ids[i], ids[j])
    mol.bond(carbon, ids[0])
    mol.res_atoms[owner].update(ids)
    return ids[0]


# ---------------------------------------------------------------- structure --


@dataclass
class RGIStructure:
    """Declarative description of one derivatized RG-I oligomer isomer.

    backbone : residue types, non-reducing end first, e.g.
               ['dGalA','Rha','GalA','Rha','GalA','Rha']
    acyl     : background derivative on every otherwise-free OH ('dAc'|'Prop'|'Ac'|'OH')
    oac      : {(residue_index, position): 'Ac'} native O-acetyl sites
    branches : {(residue_index, position): 'Gal'} side chains
    branch_acyl / branch_oac : same, keyed by (branch_key, position)
    """
    backbone: List[str]
    acyl: str = "dAc"
    oac: Dict[Tuple[int, int], str] = field(default_factory=dict)
    branches: Dict[Tuple[int, int], str] = field(default_factory=dict)
    branch_oac: Dict[Tuple[object, int], str] = field(default_factory=dict)
    name: str = ""

    # ---- derived -----------------------------------------------------------

    def residue_keys(self) -> List[object]:
        keys: List[object] = list(range(len(self.backbone)))
        keys += [("b", r, p) for (r, p) in sorted(self.branches)]
        return keys

    def free_positions(self) -> List[Tuple[object, int]]:
        """Positions that are acylated (native OAc or blocking group)."""
        out: List[Tuple[object, int]] = []
        n = len(self.backbone)
        for i, rt in enumerate(self.backbone):
            spec = RESIDUES[rt]
            for p in spec["slots"]:
                if p == 1 and i < n - 1:
                    continue                      # C1 donates the glycosidic bond
                if i > 0 and p == spec["link_pos"]:
                    continue                      # accepts the backbone
                if (i, p) in self.branches:
                    continue                      # accepts a side chain
                out.append((i, p))
        for (r, p), bt in sorted(self.branches.items()):
            key = ("b", r, p)
            for q in RESIDUES[bt]["slots"]:
                if q == 1:
                    continue                      # donates to the backbone
                out.append((key, q))
        return out

    def build(self) -> Mol:
        mol = Mol()
        n = len(self.backbone)
        mol.n_backbone = n
        c1_of: Dict[object, int] = {}
        slot_oxy: Dict[Tuple[object, int], int] = {}
        carbons: Dict[object, Dict[int, int]] = {}

        def add_residue(key: object, rtype: str,
                        subs: Dict[int, str], donates: bool) -> None:
            spec = RESIDUES[rtype]
            mol.res_type[key] = rtype
            C = {}
            for i in range(1, 6):
                C[i] = mol.add("C")
                mol.res_atoms[key].add(C[i])
            O5 = mol.add("O")
            mol.res_atoms[key].add(O5)
            ring_pairs = [(0, O5, C[1]), (1, C[1], C[2]), (2, C[2], C[3]),
                          (3, C[3], C[4]), (4, C[4], C[5]), (5, C[5], O5)]
            for idx, a, b in ring_pairs:
                mol.bond(a, b)
                if spec["unsat"] and idx == 4:
                    continue                       # C4=C5 double bond, not cleaved
                mol.cuts.append(Cut("ring", key, idx, a, b))
            # ring hydrogens
            hcount = {1: 1, 2: 1, 3: 1, 4: 1, 5: 1}
            if spec["unsat"]:
                hcount[5] = 0
            for i, nh in hcount.items():
                for _ in range(nh):
                    h = mol.add("H")
                    mol.bond(C[i], h)
                    mol.res_atoms[key].add(h)
            # C6
            c6 = mol.add("C")
            mol.bond(C[5], c6)
            mol.res_atoms[key].add(c6)
            if spec["c6"] == "CH3":
                for _ in range(3):
                    h = mol.add("H")
                    mol.bond(c6, h)
                    mol.res_atoms[key].add(h)
            elif spec["c6"] == "COOH":
                od = mol.add("O")
                mol.bond(c6, od)
                oh = mol.add("O")
                mol.bond(c6, oh)
                h = mol.add("H")
                mol.bond(oh, h)
                mol.res_atoms[key].update([od, oh, h])
            elif spec["c6"] == "CH2OH":
                for _ in range(2):
                    h = mol.add("H")
                    mol.bond(c6, h)
                    mol.res_atoms[key].add(h)
                grp = subs.get(6, self.acyl)
                slot_oxy[(key, 6)] = _attach_group(mol, c6, grp, key)
            # O-substituents on ring carbons
            for p in spec["slots"]:
                if p == 6:
                    continue
                if p == 1 and donates:
                    continue                       # no O1: C1 bonds the acceptor O
                grp = subs.get(p)
                if grp is None:
                    continue                       # position is a link acceptor
                slot_oxy[(key, p)] = _attach_group(mol, C[p], grp, key)
            c1_of[key] = C[1]
            carbons[key] = C
            mol.res_c1[key] = C[1]

        # ---- backbone ------------------------------------------------------
        for i, rt in enumerate(self.backbone):
            spec = RESIDUES[rt]
            subs: Dict[int, str] = {}
            for p in spec["slots"]:
                if p == 1 and i < n - 1:
                    continue
                if i > 0 and p == spec["link_pos"]:
                    continue                       # filled by backbone link
                if (i, p) in self.branches:
                    continue                       # filled by side chain
                subs[p] = self.oac.get((i, p), self.acyl)
            add_residue(i, rt, subs, donates=(i < n - 1))

        # ---- branches ------------------------------------------------------
        for (r, p), bt in sorted(self.branches.items()):
            key = ("b", r, p)
            spec = RESIDUES[bt]
            subs = {q: self.branch_oac.get((key, q), self.acyl)
                    for q in spec["slots"] if q != 1}
            add_residue(key, bt, subs, donates=True)

        # ---- glycosidic linkages ------------------------------------------
        def link(donor_key, acceptor_key, pos, lkey):
            spec = RESIDUES[mol.res_type[acceptor_key]]
            o = mol.add("O")
            mol.res_atoms[acceptor_key].add(o)
            cp = carbons[acceptor_key][pos]
            c1 = c1_of[donor_key]
            mol.bond(c1, o)
            mol.bond(o, cp)
            mol.cuts.append(Cut("gly", lkey, 0, c1, o, oxy=o))
            mol.cuts.append(Cut("gly", lkey, 1, o, cp, oxy=o))

        for i in range(n - 1):
            acc = i + 1
            link(i, acc, RESIDUES[self.backbone[acc]]["link_pos"], i)
        for (r, p), bt in sorted(self.branches.items()):
            link(("b", r, p), r, p, ("b", r, p))

        return mol


# -------------------------------------------------------------- fragmenting --


@dataclass(frozen=True)
class Fragment:
    atoms: FrozenSet[int]
    cuts: FrozenSet[Cut]
    mol: object = field(repr=False, default=None)

    def h_adjust(self) -> int:
        h = 0
        for c in self.cuts:
            if c.kind == "gly":
                h += 1 if c.oxy in self.atoms else -1
        return h

    def neutral_mass(self) -> float:
        return self.mol.mass(self.atoms) + self.h_adjust() * ELEMENTS["H"]

    def mz(self, adduct: str = "Na", charge: int = 1) -> float:
        dm, sign = ADDUCTS[adduct]
        return (self.neutral_mass() + dm * charge) / abs(charge)


def _components(mol: Mol, atoms: FrozenSet[int],
                removed: Set[Tuple[int, int]]) -> List[FrozenSet[int]]:
    seen: Set[int] = set()
    out: List[FrozenSet[int]] = []
    for start in atoms:
        if start in seen:
            continue
        stack = [start]
        comp = set()
        seen.add(start)
        while stack:
            a = stack.pop()
            comp.add(a)
            for b in mol.adj[a]:
                if b not in atoms or b in seen:
                    continue
                if (min(a, b), max(a, b)) in removed:
                    continue
                seen.add(b)
                stack.append(b)
        out.append(frozenset(comp))
    return out


def _events(mol: Mol, atoms: FrozenSet[int],
            cross_ring: bool = True) -> List[Tuple[Cut, ...]]:
    """Single cleavage events available inside `atoms`."""
    inside = [c for c in mol.cuts if c.a in atoms and c.b in atoms]
    ev: List[Tuple[Cut, ...]] = [(c,) for c in inside if c.kind == "gly"]
    if cross_ring:
        by_res: Dict[object, List[Cut]] = defaultdict(list)
        for c in inside:
            if c.kind == "ring":
                by_res[c.site].append(c)
        for cuts in by_res.values():
            ev += [tuple(p) for p in itertools.combinations(cuts, 2)]
    return ev


def fragment(mol: Mol,
             precursor: Optional[Fragment] = None,
             depth: int = 1,
             cross_ring: bool = True,
             min_mass: float = 100.0) -> List[Fragment]:
    """Enumerate fragments produced by up to `depth` sequential cleavage events."""
    if precursor is None:
        precursor = Fragment(mol.all_atoms, frozenset(), mol)
    best: Dict[FrozenSet[int], Fragment] = {}
    frontier = [precursor]
    for _ in range(depth):
        nxt: List[Fragment] = []
        for parent in frontier:
            for ev in _events(mol, parent.atoms, cross_ring):
                removed = {(min(c.a, c.b), max(c.a, c.b)) for c in ev}
                comps = _components(mol, parent.atoms, removed)
                if len(comps) < 2:
                    continue                       # cleavage did not separate
                for comp in comps:
                    bound = frozenset(
                        c for c in set(parent.cuts) | set(ev)
                        if (c.a in comp) != (c.b in comp))
                    f = Fragment(comp, bound, mol)
                    if f.neutral_mass() < min_mass:
                        continue
                    prev = best.get(comp)
                    if prev is None or len(f.cuts) < len(prev.cuts):
                        best[comp] = f
                    nxt.append(f)
        frontier = nxt
    return sorted(best.values(), key=lambda f: f.neutral_mass())


# ---------------------------------------------------------------- labelling --

_SUP = {0: "0", 1: "1", 2: "2", 3: "3", 4: "4", 5: "5"}


def label(frag: Fragment, struct: RGIStructure) -> str:
    mol = frag.mol
    n = len(struct.backbone)
    parts: List[str] = []
    ring_by_res: Dict[object, List[int]] = defaultdict(list)
    for c in sorted(frag.cuts, key=lambda c: (c.kind, str(c.site), c.idx)):
        if c.kind == "gly":
            lk = c.site
            if isinstance(lk, int):
                b_sub, y_sub = lk + 1, n - 1 - lk
            else:
                b_sub, y_sub = "b", "b"
            has_o = c.oxy in frag.atoms
            if c.idx == 0:
                parts.append(f"Y{y_sub}" if has_o else f"B{b_sub}")
            else:
                parts.append(f"C{b_sub}" if has_o else f"Z{y_sub}")
        else:
            ring_by_res[c.site].append(c.idx)
    for res, idxs in sorted(ring_by_res.items(), key=lambda kv: str(kv[0])):
        if len(idxs) != 2:
            continue
        i, j = sorted(idxs)
        rt = mol.res_type[res]
        if isinstance(res, int):
            a_sub, x_sub = res + 1, n - 1 - res
        else:
            a_sub, x_sub = "b", "b"
        is_x = mol.res_c1[res] in frag.atoms
        tag = f"{i},{j}" + ("X" if is_x else "A")
        sub = x_sub if is_x else a_sub
        parts.append(f"{tag}{sub}[{rt}]")
    return " ".join(parts) if parts else "M"


def residue_span(frag: Fragment) -> str:
    mol = frag.mol
    whole, partial = [], []
    for key, atoms in mol.res_atoms.items():
        inter = len(atoms & frag.atoms)
        if inter == len(atoms):
            whole.append(key)
        elif inter:
            partial.append(key)
    def nm(k):
        return str(k) if isinstance(k, int) else "b%s" % (k[1],)
    s = "+".join(sorted(nm(k) for k in whole))
    if partial:
        s += " ~" + "+".join(sorted(nm(k) for k in partial))
    return s


# ----------------------------------------------------------------- matching --


@dataclass
class Match:
    mz_obs: float
    intensity: float
    mz_calc: float
    error_mmu: float
    error_ppm: float
    label: str
    span: str
    formula: str


def formula(frag: Fragment) -> str:
    cnt: Dict[str, int] = defaultdict(int)
    for a in frag.atoms:
        cnt[frag.mol.el[a]] += 1
    cnt["H"] += frag.h_adjust()
    return "".join(f"{e}{cnt[e]}" for e in ("C", "H", "D", "O") if cnt.get(e))


def match_peaks(peaks: Sequence,
                frags: Sequence[Fragment],
                struct: RGIStructure,
                adduct: str = "Na",
                tol_mmu: float = 30.0,
                h_transfer: bool = False) -> List[Match]:
    """Match an observed peak list against theoretical fragments.

    peaks : sequence of m/z, or of (m/z, intensity)
    """
    table: List[Tuple[float, Fragment, int]] = []
    for f in frags:
        base = f.mz(adduct)
        table.append((base, f, 0))
        if h_transfer:
            table.append((base + ELEMENTS["H"], f, +1))
            table.append((base - ELEMENTS["H"], f, -1))
    table.sort(key=lambda t: t[0])
    out: List[Match] = []
    for p in peaks:
        mz, inten = (p if isinstance(p, (tuple, list)) else (p, 0.0))
        best = None
        for mzc, f, dh in table:
            err = (mz - mzc) * 1000.0
            if abs(err) <= tol_mmu and (best is None or abs(err) < abs(best[0])):
                lab = label(f, struct) + ("" if dh == 0 else f" {dh:+d}H")
                best = (err, mzc, lab, residue_span(f), formula(f))
        if best is None:
            out.append(Match(mz, inten, float("nan"), float("nan"),
                             float("nan"), "unassigned", "", ""))
        else:
            err, mzc, lab, span, form = best
            out.append(Match(mz, inten, mzc, err, err / mzc * 1000.0,
                             lab, span, form))
    return out


# ----------------------------------------------------- isomer enumeration ----

# Biologically observed O-acetylation sites in RG-I (Komalavilas & Mort 1989;
# Fernando et al. 2026): Rha 3-O (dominant) and 4-O; GalA 2-O and 3-O.
OAC_SITES = {"Rha": (3, 4), "GalA": (2, 3), "dGalA": (2, 3)}


def enumerate_isomers(backbone: Sequence[str],
                      n_oac: int,
                      acyl: str = "dAc",
                      branches: Optional[Dict[Tuple[int, int], str]] = None,
                      sites: Optional[Dict[str, Tuple[int, ...]]] = None,
                      allow_reducing_c1: bool = False) -> List[RGIStructure]:
    """All positional isomers carrying `n_oac` native O-acetyl groups."""
    branches = dict(branches or {})
    sites = sites or OAC_SITES
    n = len(backbone)
    cand: List[Tuple[int, int]] = []
    for i, rt in enumerate(backbone):
        spec = RESIDUES[rt]
        for p in spec["slots"]:
            if p not in sites.get(rt, ()):
                if not (allow_reducing_c1 and p == 1 and i == n - 1):
                    continue
            if p == 1 and i < n - 1:
                continue
            if i > 0 and p == spec["link_pos"]:
                continue
            if (i, p) in branches:
                continue
            cand.append((i, p))
    out: List[RGIStructure] = []
    for combo in itertools.combinations(cand, n_oac):
        oac = {k: "Ac" for k in combo}
        nm = ",".join(f"{backbone[i]}{i+1}-{p}O" for i, p in combo) or "none"
        out.append(RGIStructure(list(backbone), acyl=acyl, oac=oac,
                                branches=branches, name=nm))
    return out


def precursor_mz(struct: RGIStructure, adduct: str = "Na") -> float:
    mol = struct.build()
    return Fragment(mol.all_atoms, frozenset(), mol).mz(adduct)


def select_precursor(mol: Mol, frags: Sequence[Fragment], mz: float,
                     adduct: str = "Na", tol_mmu: float = 500.0
                     ) -> List[Fragment]:
    """Fragments matching an isolated MS^n precursor m/z (wide isolation window)."""
    hits = [f for f in frags if abs(f.mz(adduct) - mz) * 1000.0 <= tol_mmu]
    return sorted(hits, key=lambda f: abs(f.mz(adduct) - mz))


def assign_all(mz: float,
               frags: Sequence[Fragment],
               struct: RGIStructure,
               adduct: str = "Na",
               tol_mmu: float = 10.0) -> List[Tuple[float, str]]:
    """Every theoretical fragment within tolerance of one peak (not just the best)."""
    out = []
    for f in frags:
        c = f.mz(adduct)
        if abs(mz - c) * 1000.0 <= tol_mmu:
            out.append((c, label(f, struct)))
    return sorted(set(out))


def fragments_for(struct: RGIStructure,
                  precursor_path: Sequence[float] = (),
                  depth: int = 2,
                  adduct: str = "Na",
                  isolation_mmu: float = 1500.0,
                  cross_ring: bool = True) -> Tuple[Mol, List[Fragment]]:
    """Fragments of `struct` after walking down an MS^n precursor selection path."""
    mol = struct.build()
    node = Fragment(mol.all_atoms, frozenset(), mol)
    for mz in precursor_path:
        pool = fragment(mol, precursor=node, depth=1, cross_ring=cross_ring)
        hits = select_precursor(mol, pool, mz, adduct, isolation_mmu)
        if not hits:
            raise ValueError(f"no fragment of {struct.name!r} near precursor m/z {mz}")
        node = hits[0]
    return mol, fragment(mol, precursor=node, depth=depth, cross_ring=cross_ring)


@dataclass
class Interpretation:
    peak_table: List[dict]
    scores: Dict[str, dict]
    ranking: List[Tuple[str, float]]


def interpret(peaks: Sequence,
              structures: Sequence[RGIStructure],
              precursor_path: Sequence[float] = (),
              depth: int = 2,
              adduct: str = "Na",
              tol_mmu: float = 10.0,
              isolation_mmu: float = 1500.0) -> Interpretation:
    """Score candidate positional isomers against an observed peak list.

    An isomer is credited only for peaks it explains with an ion that no other
    candidate can produce ("discriminating" ions); peaks explained by every
    candidate carry no positional information and are reported but not scored.
    """
    frs: Dict[str, List[Fragment]] = {}
    for s in structures:
        _, f = fragments_for(s, precursor_path, depth, adduct, isolation_mmu)
        frs[s.name] = f

    rows: List[dict] = []
    credit: Dict[str, float] = {s.name: 0.0 for s in structures}
    nmatched: Dict[str, int] = {s.name: 0 for s in structures}
    total_int = sum((p[1] if isinstance(p, (tuple, list)) else 1.0) for p in peaks) or 1.0

    for p in peaks:
        mz, inten = (p if isinstance(p, (tuple, list)) else (p, 1.0))
        hits = {s.name: assign_all(mz, frs[s.name], s, adduct, tol_mmu)
                for s in structures}
        explained = [n for n, h in hits.items() if h]
        row = {"mz": mz, "intensity": inten,
               "n_candidates_explaining": len(explained),
               "explained_by": ";".join(explained) or "none",
               "discriminating": len(explained) == 1 and len(structures) > 1,
               "n_assignments": sum(len(h) for h in hits.values())}
        for n, h in hits.items():
            row[f"assign[{n}]"] = " | ".join(f"{m:.4f} {l}" for m, l in h[:4])
        rows.append(row)
        for n in explained:
            nmatched[n] += 1
            credit[n] += inten / len(explained)
    scores = {n: {"discriminating_intensity": credit[n],
                  "discriminating_fraction": credit[n] / total_int,
                  "peaks_explained": nmatched[n],
                  "n_theoretical": len(frs[n])}
              for n in credit}
    ranking = sorted(((n, credit[n]) for n in credit), key=lambda t: -t[1])
    return Interpretation(rows, scores, ranking)


def minimal_mixture(peaks: Sequence,
                    structures: Sequence[RGIStructure],
                    precursor_path: Sequence[float] = (),
                    depth: int = 1,
                    adduct: str = "Na",
                    tol_mmu: float = 10.0,
                    max_components: int = 3,
                    isolation_mmu: float = 1500.0) -> dict:
    """Smallest combination of candidate isomers that explains every peak.

    RG-I digests are routinely mixtures of positional isomers, so the question
    is usually not 'which isomer' but 'which smallest set of isomers'.  All
    equally-small solutions are returned; a tie means the peak list cannot
    separate them (typically 3-O vs 4-O, which needs cross-ring ions).
    """
    cover: Dict[str, Set[int]] = {}
    for s in structures:
        _, frs = fragments_for(s, precursor_path, depth, adduct, isolation_mmu)
        mzs = [f.mz(adduct) for f in frs]
        hit = set()
        for i, p in enumerate(peaks):
            mz = p[0] if isinstance(p, (tuple, list)) else p
            if any(abs(mz - c) * 1000.0 <= tol_mmu for c in mzs):
                hit.add(i)
        cover[s.name] = hit
    need = set(range(len(peaks)))
    names = list(cover)
    for k in range(1, max_components + 1):
        sols = [combo for combo in itertools.combinations(names, k)
                if set().union(*(cover[n] for n in combo)) == need]
        if sols:
            return {"n_components": k, "solutions": sols, "coverage": cover,
                    "unexplained": []}
    best = max(itertools.combinations(names, max_components),
               key=lambda c: len(set().union(*(cover[n] for n in c))))
    got = set().union(*(cover[n] for n in best))
    return {"n_components": None, "solutions": [best], "coverage": cover,
            "unexplained": sorted(need - got)}


def diagnostic_ions(structs: Sequence[RGIStructure],
                    depth: int = 2,
                    adduct: str = "Na",
                    tol_mmu: float = 10.0,
                    precursor: Optional[float] = None,
                    precursor_tol_mmu: float = 1500.0
                    ) -> Dict[str, List[Tuple[float, str]]]:
    """m/z values unique to a single candidate isomer, per isomer."""
    per: Dict[str, List[Tuple[float, str]]] = {}
    for s in structs:
        mol = s.build()
        if precursor is None:
            frs = fragment(mol, depth=depth)
        else:
            root = select_precursor(mol, fragment(mol, depth=1), precursor,
                                    adduct, precursor_tol_mmu)
            frs = []
            for r in root[:3]:
                frs += fragment(mol, precursor=r, depth=depth)
        per[s.name] = [(f.mz(adduct), label(f, s)) for f in frs]
    uniq: Dict[str, List[Tuple[float, str]]] = {}
    for name, lst in per.items():
        others = [m for k, v in per.items() if k != name for m, _ in v]
        keep = []
        for mz, lab in lst:
            if all(abs(mz - o) * 1000.0 > tol_mmu for o in others):
                keep.append((mz, lab))
        uniq[name] = sorted(keep)
    return uniq
