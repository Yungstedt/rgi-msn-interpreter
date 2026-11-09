# Validation record

Everything below was produced by `rgi_msn.py` / `rgi_io.py` at the commit this file is
checked in at. Each section names the artefact that holds the per-row numbers, so a
reviewer can check any individual claim rather than taking a summary on trust.

The reference throughout is Fernando, Yang, Archer-Hartmann, Shahin, Zhang, Urbanowicz,
Heiss & Azadi, *J. Am. Soc. Mass Spectrom.* **2026**, 37, 1160-1172
(doi:10.1021/jasms.5c00440), referred to as "the method paper".

| # | check | result | artefact |
|---|---|---|---|
| 1 | published precursor and MS<sup>2</sup> *m/z* | 35 / 35 reproduced | [rgi_validation.csv](../data/rgi_validation.csv) |
| 2 | published MS<sup>3</sup> cross-ring masses (Table 2) | 15 / 16 to &le;1 mmu | [rgi_ms3_table2_check.csv](../data/rgi_ms3_table2_check.csv) |
| 3 | discrepancies found in the published tables | 3, listed below | - |
| 4 | published MS<sup>2</sup> ions recovered from a raw profile export | 9 / 9, &plusmn;1.4 mmu after offset | [rgi_realdata_validation.csv](../data/rgi_realdata_validation.csv) |
| 5 | mass accuracy, MS<sup>2</sup> vs MS<sup>3</sup> | constant, 0.28-0.84 mmu median residual | [rgi_mass_accuracy.csv](../data/rgi_mass_accuracy.csv) |
| 6 | mixture deconvolution on real spectra | 2 components (DP6 MS<sup>2</sup>), 1 (DP4 MS<sup>3</sup>) | below |
| 7 | isomer ranking on a synthetic peak list | [ranking](../data/rgi_interpret_demo_ranking.csv) / [peaks](../data/rgi_interpret_demo_peaks.csv) |  |

## 1. Theoretical *m/z* against the method paper

35 published values - 12 precursors spanning the trideuteroacetyl and propionyl
derivatives in both polarities, and 23 MS<sup>2</sup> fragments - are reproduced to the
precision printed in the paper. Columns: published value, model value, signed difference.

## 2. MS<sup>3</sup> cross-ring masses (Table 2)

15 of the 16 entries are reproduced to &le;1 mmu, and the paper's own Domon-Costello labels
are recovered independently by the fragment generator rather than assumed - for example
*m/z* 401.0961 = B<sub>5</sub> <sup>3,5</sup>X<sub>Rha</sub> Z<sub>3</sub> and 669.2027 =
B<sub>5</sub> <sup>0,3</sup>X<sub>Rha</sub> Y<sub>4</sub>.

## 3. Discrepancies in the published tables

These are reported so that reviewers can see the model disagreeing with its reference where
it does, rather than only where it agrees. They concern the tabulated values and labels, not
the chemistry or the conclusions of the method paper.

1. **Table 2, *m/z* 477.** The published *theoretical* value, 477.1757, disagrees with the
   observed value printed beside it, 477.1486. The model computes 477.1486 - within 0.03 mmu
   of the published observation - and one of the three fragments it finds at that mass
   carries exactly the published label (`Y3 2,4X_Rha 1,3A_Rha`). The label therefore stands
   and only the tabulated theoretical mass is wrong. Recorded as
   `mz_model_nearest_observed` / `label_model_nearest_observed`;
   `mz_model_at_paper_theoretical` is blank for this row because nothing in the model sits
   at 477.1757.
2. **Two ion labels.** Negative-mode *m/z* 1045 is Z<sub>4</sub> rather than Y<sub>4</sub>
   (Y<sub>4</sub> = 1063.43), and the propionyl *m/z* 1175 is Y<sub>4</sub> rather than
   B<sub>5</sub> (a propionyl B<sub>5</sub> must exceed the deuteroacetyl B<sub>5</sub> at
   1275, and computes to 1385.45). The masses themselves are reproduced exactly.
3. **Branched celery pentasaccharide, *m/z* 1366.** B<sub>3</sub> = 1067.40 requires the Gal
   branch on the **inner** rhamnose; with the branch on the reducing-end rhamnose
   B<sub>3</sub> = 770.26. The Figure 3 text assigns structures J/K/L - which give 1067 - to
   a reducing-end branch, which the arithmetic does not support.

## 4. Recovery from a raw instrument file

Input: the Thermo FreeStyle `SPECTRUM - MS` export of the MS<sup>2</sup> of *m/z* 1571.64
(CID 30), 84,889 profile points, over half of them exactly zero, centroided by the reader to
625 peaks.

All nine ions reported in the method paper are recovered. In the raw export each sits
**+2.17 to +4.19 mmu high** - a systematic bias of consistent sign, not scatter. Removing
the measured offset brings every residual into the range &minus;0.61 to +1.40 mmu. Per-ion
values, raw and calibrated, are in
[rgi_realdata_validation.csv](../data/rgi_realdata_validation.csv).

## 5. Mass accuracy as a function of MS<sup>n</sup> level

This is the check behind the tolerance guidance in [METHOD.md](METHOD.md) section 5. Three
acquisitions of the same deuteroacetylated DP6, two of them one MS level deeper than the
first:

| spectrum | level | activation | peaks after deisotoping | ions matched &le;10 mmu | calibration offset | residual after offset (median) | worst residual |
|---|---|---|---|---|---|---|---|
| *m/z* 1571.61 | MS<sup>2</sup> | CID 30 | 249 | 20 | +3.48 mmu | 0.68 mmu (0.68 ppm) | 2.22 mmu |
| 1571.61 &rarr; 1087.42 | MS<sup>3</sup> | CID 30 &rarr; HCD 30 | 51 | 24 | &minus;1.67 mmu | 0.84 mmu (1.69 ppm) | 4.33 mmu |
| 1571.61 &rarr; 1087.42 | MS<sup>3</sup> | CID 30 &rarr; CID 30 | 10 | 7 | &minus;1.82 mmu | 0.28 mmu (0.54 ppm) | 3.79 mmu |

Method: peaks were centroided and deisotoped by `rgi_io`, then matched against the union of
the isomer E and isomer F fragment sets (`depth=1` for MS<sup>2</sup>, `depth=2` from
*m/z* 1087.43 for MS<sup>3</sup>). An ion counted as matched if a calculated species lay
within 10 mmu; the offset is the median signed residual over the matched ions, and the
residual column is the median absolute deviation from that offset. Per-ion rows, with
labels, are in [rgi_mass_accuracy.csv](../data/rgi_mass_accuracy.csv).

**Reading of this result.** The median residual does not grow from MS<sup>2</sup> to
MS<sup>3</sup> - it is 0.68, 0.84 and 0.28 mmu, i.e. sub-2-ppm throughout. The offsets
differ in sign and are of order a few mmu at both levels. On an instrument that mass-
analyses every MS<sup>n</sup> stage in the Orbitrap there is therefore no basis for widening
the tolerance with depth; the appropriate operation is to subtract the per-spectrum offset
and hold the window at 5 mmu. The worst-residual column is dominated by a small number of
weak peaks close to the noise floor and is the reason the recommendation is 5 mmu rather
than 2.

## 6. Mixture deconvolution on real spectra

Settings for both: `tol_mmu=5`, peaks above 0.5 % of base peak, deisotoped,
`losses=rgi_io.NEUTRAL_LOSSES`, `coverage=0.95`.

**DP6, MS<sup>2</sup> of *m/z* 1571.61.** After offset correction the spectrum requires
**two components** and they account for **100 % of the explainable intensity**. Twelve
two-component solutions are tied, and **every one of them contains a species acetylated on
the reducing-end rhamnose**; six distinct partners remain open (2-*O*- or 3-*O*- on the
&Delta;4,5-GalA, 3-*O*- or 4-*O*- on the inner Rha, 2-*O*- or 3-*O*- on the internal GalA).
The observed intensity ratio of the two B<sub>5</sub> ions, *m/z* 1275.5 : 1272.5, is
**49.6 : 50.4**, against 47 : 53 reported in the method paper.

The partner cannot be pinned down from MS<sup>2</sup> because in a two-component mixture
either component can supply B<sub>3</sub> = 773. The method paper's assignment of the second
component to the inner rhamnose (structure F) rests on the NMR prior, not on the
MS<sup>2</sup> evidence, and this tool does not reproduce that narrowing - correctly, in our
reading. MS<sup>3</sup> of Y<sub>4</sub> does separate 3-*O*- from 4-*O*-acetylation on the
reducing-end Rha at 5 mmu (13 ions unique to the 3-*O* isomer, 27 to the 4-*O*).

**DP4, MS<sup>3</sup> of *m/z* 770.5.** Resolves to a **single component** explaining 100 %
of the explainable intensity, with two tied candidates (2-*O*-Ac on the internal GalA, or
4-*O*-Ac on the inner Rha). Isomers acetylated on the reducing-end rhamnose are excluded
with a recorded reason: they cannot produce the isolated precursor at all.

## 7. What is not validated

- **No ground-truth positional standards.** Every positional claim is validated against the
  method paper's assignments, which themselves rest in part on NMR. Synthetic
  site-defined standards would be needed to validate the positional calls independently.
- **No noise or intensity model.** Scoring is intensity-weighted but has no concept of
  detection limit, so the ranking between closely scored candidates should not be treated
  as calibrated.
- **Acyl migration** is neither modelled nor excluded experimentally in these data.
- **One analyte class.** Validation covers deuteroacetylated and propionylated RG-I
  backbones of DP4-DP6 with at most one Gal branch. Larger branched structures,
  arabinan side chains and methyl-esterified GalA are implemented but untested against
  published data.
- **Single instrument platform.** All real data are Orbitrap. The tolerance guidance in
  section 5 is explicitly conditional on high-resolution MS<sup>n</sup> product scans.

## Reproducing

```python
import rgi_msn as R, rgi_io

# sections 1-3 require only the library
DP6 = ['dGalA','Rha','GalA','Rha','GalA','Rha']
E = R.RGIStructure(DP6, acyl='dAc', oac={(5,3):'Ac'}, name='E')
assert abs(R.precursor_mz(E) - 1571.6225) < 1e-3

# sections 4-6 need the FreeStyle exports, available from the
# corresponding author of the method paper
sp    = rgi_io.read_spectrum('RGI-oac-deu-MS21571.csv', rel_threshold=0.005)
peaks = rgi_io.deisotope(sp.peaks)
peaks = rgi_io.apply_offset(peaks, rgi_io.estimate_offset(peaks, [1275.4847, 773.2771]))
cands = R.enumerate_isomers(DP6, n_oac=1, acyl='dAc')
R.minimal_mixture(peaks, cands, depth=1, tol_mmu=5, max_components=3,
                  losses=rgi_io.NEUTRAL_LOSSES, coverage=0.95)
```
