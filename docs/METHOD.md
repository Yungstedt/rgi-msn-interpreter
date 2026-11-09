# `rgi_msn.py` — MS<sup>n</sup> interpreter for acetylated RG-I oligomers



Built for the derivatization scheme in Fernando, Yang, Archer-Hartmann, Shahin, Zhang,

Urbanowicz, Heiss & Azadi, *J. Am. Soc. Mass Spectrom.* **2026**, 37, 1160–1172

(doi:10.1021/jasms.5c00440): per-trideuteroacetylation or per-propionylation of a

partially O-acetylated RG-I oligomer, sodiated positive-mode CID (MS²) → HCD (MS³)

on an Orbitrap.



## Chemistry model



| feature | treatment |

|---|---|

| backbone | `[→2)-α-L-Rha-(1→4)-α-D-GalA-(1→]` , non-reducing end first |

| lyase product | Δ4,5-unsaturated GalA (`dGalA`) at the non-reducing end, Rha at the reducing end |

| blocking group | `dAc` +45.0294 · `Prop` +56.0262 · `Me` +14.0157 |

| native modification | `Ac` +42.0106, locked in place |

| branches | Gal on Rha O-4 (makes the Rha →2,4-linked) |

| acylatable sites | GalA O-2/O-3, Rha O-3/O-4, reducing-Rha anomeric O-1, Gal O-2/3/4/6 — **13 for a DP6 backbone** |



Fragment masses are computed atom-by-atom on an explicit molecular graph rather than

from residue increments, so a cross-ring ion's mass automatically reflects which acyl

group sits on which carbon. That is what distinguishes 3-O- from 2-O- or 4-O-acetylation.



Cleavage conventions: the glycosidic unit C1–O–Cp carries two cleavable bonds, giving

B/Y (C1–O) and C/Z (O–Cp), with the oxygen-retaining side gaining one hydrogen. Ring

bonds are numbered 0: O5–C1, 1: C1–C2, 2: C2–C3, 3: C3–C4, 4: C4–C5, 5: C5–O5; a

cross-ring cleavage cuts two of them, the piece retaining C1 is the X fragment, and no

hydrogen is transferred (retro-Diels–Alder convention; ±H variants available at match

time via `h_transfer=True`).



## Usage



```python

import rgi_msn as R



DP6 = ['dGalA','Rha','GalA','Rha','GalA','Rha']       # non-reducing end first



# one named isomer: 3-O-acetyl on the reducing-end rhamnose

E = R.RGIStructure(DP6, acyl='dAc', oac={(5,3):'Ac'}, name='E')

R.precursor_mz(E)                                      # 1571.62  [M+Na]+



mol, frags = R.fragments_for(E, precursor_path=[1087.43], depth=2)   # MS3 of Y4

R.match_peaks(my_peaks, frags, E, tol_mmu=5)



# candidate enumeration + scoring against an observed peak list

cands = R.enumerate_isomers(DP6, n_oac=1, acyl='dAc')  # 12 positional isomers

res   = R.interpret(my_peaks, cands, depth=1, tol_mmu=5)

res.ranking



# mixtures (the normal case for RG-I digests)

R.minimal_mixture(my_peaks, cands, depth=1, tol_mmu=5, max_components=3)



# which ions can actually tell two candidates apart

R.diagnostic_ions([iso_a, iso_b], depth=2, precursor=1087.43, tol_mmu=5)

```



Peak lists are plain `m/z` values or `(m/z, intensity)` pairs. Set `adduct='-H'` for

negative mode, `acyl='Prop'` for the propionyl derivative.



## Validation



`rgi_validation.csv` — 35 published m/z values (12 precursors across both derivatives

and both polarities, 23 MS² fragments) reproduced, all to the nominal value printed in

the paper.



`rgi_ms3_table2_check.csv` — 15 of the 16 MS³ cross-ring masses in the paper's Table 2

reproduced to ≤1 mmu, with the paper's own Domon–Costello labels recovered

independently (e.g. 401.0961 = `B₅ ³'⁵X_Rha Z₃`; 669.2027 = `B₅ ⁰'³X_Rha Y₄`). The

exception is the entry whose published theoretical value, 477.1757, disagrees with its

own observed value, 477.1486. The model computes 477.1486 (0.03 mmu from the published

observation) and recovers the published label among the fragments at that mass, so the

label stands and only the tabulated theoretical mass is wrong. Recorded in

`rgi_ms3_table2_check.csv` as `mz_model_nearest_observed` /

`label_model_nearest_observed`; `mz_model_at_paper_theoretical` is blank for this row

because nothing in the model sits at 477.1757.



Two published ions carry labels the model does not support: negative-mode m/z 1045 is

Z₄ rather than Y₄, and the propionyl m/z 1175 is Y₄ rather than B₅ (a propionyl B₅ must

exceed the deuteroacetyl B₅ at 1275). The masses themselves are reproduced exactly.



For the branched celery pentasaccharide at m/z 1366, B₃ = 1067.40 requires the Gal

branch on the **inner** rhamnose; with the branch on the reducing-end rhamnose B₃ is

770.26. The paper's Figure 3 text assigns structures J/K/L (which give 1067) to a

reducing-end branch, which the arithmetic does not support.



## Known limits



- `depth=3` on a DP6 generates ~11,800 distinct fragment species. Against the 16 MS³

  peaks of Table 2 at ±30 mmu, a peak matches a median of 4 and up to 45 of them, each

  with its own cleavage path and Domon–Costello label. Those matches collapse onto a

  median of 1 and at most 2 distinct calculated m/z, i.e. the competing interpretations

  are essentially all isobaric and cannot be separated by narrowing the tolerance alone.

  Use ≤5 mmu to cut the count, and report the full candidate list from `assign_all`

  rather than a single assignment. Counts are in `rgi_ms3_table2_check.csv`

  (`n_fragment_species_within_30mmu`, `n_distinct_labels_within_30mmu`,

  `n_distinct_mz_within_30mmu`).

- Acyl migration is not modelled — the derivatization is what prevents it.

- No intensity model: scoring uses supplied intensities if present, otherwise unit weights.

- `dGalA + H₂O = GalA` exactly, so some internal B/Y fragments are isobaric with C/Z

  fragments one residue over (e.g. `C₂` and `B₃Y₅` both at 525.1968 in the DP6 case).


## Reading real instrument data

`rgi_io.read_spectrum(path)` takes the Thermo FreeStyle **SPECTRUM - MS** CSV export
directly. That file is not a peak list: it carries a BOM, an eight-line metadata header
and, usually, *profile* data — tens of thousands of points spaced 0.0001–0.003 Da with
over half of them exactly zero. Handing that to a fragment matcher is meaningless, which
is why v0.1.0 struggled with it.

The reader:

- parses the scan filter, so the MS<sup>n</sup> precursor path, polarity and activation
  come from the file rather than from the operator —
  `FTMS + p ESI Full ms3 1571.6100@cid30.00 1087.4200@hcd30.00` becomes
  MS³, `+`, `[1571.61, 1087.42]`, `['cid','hcd']`;
- detects profile vs centroid data and **centroids** profile peaks by intensity-weighted
  apex fitting (84,889 points → 625 peaks on the MS² of *m/z* 1571);
- `deisotope()` removes ¹³C isotopologues, which are about a third of the strong peaks
  and otherwise dilute intensity-weighted scoring;
- `estimate_offset()` / `apply_offset()` remove a constant calibration bias measured
  against known fragments.

On the published MS² of *m/z* 1571, all nine reported ions are recovered from the raw
profile export, each sitting +2.1 to +4.1 mmu high — a systematic calibration bias, not
scatter. Removing the measured +2.78 mmu offset brings every residual inside ±1.3 mmu.
See [data/rgi_realdata_validation.csv](data/rgi_realdata_validation.csv).

Two matching options matter on real data and are off by default in the library, on in the
application: `losses=rgi_io.NEUTRAL_LOSSES` (d₃-acetic acid 63.040 dominates these
spectra, then acetic acid and water), and `coverage=` on `minimal_mixture`, which asks
for a fraction of the explainable intensity rather than demanding every peak be explained
— a real spectrum also contains the surviving precursor, noise and out-of-model losses,
so the all-or-nothing rule simply fails.
