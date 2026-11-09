# Running the RG-I MS<sup>n</sup> interpreter



Three ways in, in order of how little setup they need.



## 1. The desktop app (no Python needed)



Double-click **`RGI-MSn-Interpreter.exe`**. Nothing to install — it is a single

self-contained Windows executable; put it anywhere, including a USB stick or a shared

instrument PC. The first launch takes a few seconds while it unpacks itself.



Windows SmartScreen will warn that the publisher is unknown because the file is not

code-signed: *More info → Run anyway*.



### Build provenance



Check that the file you have is the one that was tested:



| | |

|---|---|

| size | 9,391,959 bytes (9.0 MB) |

| SHA-256 | `877e73ccf778e29d0c1b60a6188f4fab9990214e9d5308a907f6f9c41d7a390f` |

| built | PyInstaller `--onefile --windowed`, Windows x64 |

| verified | `--selftest` exits 0: Tcl/Tk load, a Tk window is created, and the model reproduces the DP6 1-OAc precursor at 1571.6225 |



```powershell

(Get-FileHash RGI-MSn-Interpreter.exe -Algorithm SHA256).Hash

```



A different hash means a different build. v0.1.0's `36511d8d...` worked but could not read FreeStyle CSV exports; v0.2.0 adds that reader. In particular, **do not use the original

v0.1.0 binary** (SHA-256 `410b8f33…`): it was missing the Tcl/Tk runtime DLLs and failed

at startup with `ImportError: DLL load failed while importing _tkinter`.



### If it will not start



Run it with the self-test flag — it writes a diagnostic report instead of opening a

window, and exits non-zero on failure:



```powershell

.\RGI-MSn-Interpreter.exe --selftest report.txt

type report.txt

```



The report records whether Tcl/Tk loaded, whether a window could be created, and whether

the fragment engine returns the expected reference mass. Send it along if you need help.



An alternative build is published as **RGI-MSn-Interpreter-win64.zip** — the same

application as a folder rather than a single file. It unpacks everything next to the

executable instead of expanding into the temporary directory at each launch, which is

the more reliable choice on machines where the temp folder is redirected, locked down,

or scanned aggressively. Unzip it and run `RGI-MSn-Interpreter.exe` inside the folder;

keep the `_internal` folder beside it.



### Using the window



**Oligomer** (left panel) — pick a preset or type the backbone yourself,

non-reducing end first, comma separated. `dGalA` is the Δ4,5-unsaturated GalA that lyase

digestion leaves at the non-reducing end; `Rha`, `GalA` and `Gal` are the saturated

residues. Choose the blocking group you actually used (trideuteroacetyl or propionyl),

how many native O-acetyls the precursor carries, and the adduct. For a Gal side chain,

enter the index of the Rha that carries it (0 = non-reducing end); it is placed on O-4,

which makes that Rha →2,4-linked.



**Acquisition** — leave the precursor path blank for MS². For MS³, enter the m/z you

isolated at MS² (e.g. `1087.43`); for MS⁴, enter both, separated by a space. *Cleavage

depth* is how many further cleavage events to allow below that precursor — 1 is right

for a clean MS³, 2 widens the net at the cost of many more candidates. Keep the

tolerance at 5 mmu for Orbitrap data; widen it only for ion-trap MS<sup>n</sup>

detection, and read the degeneracy warning in `README_rgi_msn.md` before you do.



**Observed peak list** — paste `m/z` values one per line, optionally followed by an

intensity, or load a `.csv`/`.txt`/`.tsv` export. Header lines and text are skipped.



Then:



- **Predict fragments** — every theoretical ion for every candidate isomer, with its

  Domon–Costello label, formula and residue span. No peak list required.

- **Interpret peak list** — matches your peaks, reports *how many* candidate isomers can

  explain each peak and *how many distinct assignments* exist for it, flags the ones

  unique to a single isomer, and ranks the candidates underneath the peak rows.

- **Minimal mixture** — the smallest set of isomers that jointly explains every peak.

  Use this one by default on digests: these samples are usually mixtures, and a forced

  single-isomer answer will be wrong. Several equally small solutions means the peak

  list genuinely cannot separate them.

- **Export table to CSV** — whatever is in the results grid.



The strip chart shows predicted ions above the axis (matched ones highlighted) and your

observed peaks below.



## 2. From source, as a GUI



```

python rgi_gui.py

```



Needs only Python 3.9+ with Tkinter, which ships with the standard Windows and macOS

installers. On Debian/Ubuntu: `sudo apt install python3-tk`. No other packages —

`rgi_msn.py` is pure standard library.



## 3. From source, as a library



Put `rgi_msn.py` on your path and script it. This is the route for batch work across

many spectra, or for anything the window does not expose.



```python

import rgi_msn as R



DP6 = ['dGalA','Rha','GalA','Rha','GalA','Rha']



E = R.RGIStructure(DP6, acyl='dAc', oac={(5,3):'Ac'}, name='E')

R.precursor_mz(E)                                    # 1571.6225  [M+Na]+



mol, frags = R.fragments_for(E, precursor_path=[1087.43], depth=2)

for m in R.match_peaks(my_peaks, frags, E, tol_mmu=5):

    print(m.mz_obs, m.label, m.error_ppm)



cands = R.enumerate_isomers(DP6, n_oac=1, acyl='dAc')

R.interpret(my_peaks, cands, depth=1, tol_mmu=5).ranking

R.minimal_mixture(my_peaks, cands, depth=1, tol_mmu=5)

R.assign_all(627.1921, frags, E, tol_mmu=5)          # every competing assignment

```



Residue indices in `oac={(i, pos): 'Ac'}` count from the non-reducing end starting at 0,

so `(5, 3)` is 3-O on the reducing-end Rha of a DP6.



## Rebuilding the executable



```

pip install pyinstaller

pyinstaller --onefile --windowed --name RGI-MSn-Interpreter --add-data "rgi_msn.py;." rgi_gui.py

```



On macOS/Linux the `--add-data` separator is `:` rather than `;`, and the result is an

app bundle or a bare binary rather than an `.exe`. PyInstaller cannot cross-compile:

build on the platform you want to target.


## Loading instrument data

Point **Load...** straight at a Thermo FreeStyle `SPECTRUM - MS` CSV export. The
application reads the metadata header, fills in the MS<sup>n</sup> precursor path and
polarity from the scan filter, centroids profile data, and removes 13C isotopologues.
*Min intensity (% base)* sets the noise floor on load, and the status line reports what
happened -- for example `625 peaks, centroided from 84,889 profile points, 212 isotopes
removed, MS2 + 1571.64`.

Two settings matter on real spectra and are on by default:

- **deisotope on load** -- about a third of the strong peaks in these Orbitrap spectra
  are 13C isotopologues, and leaving them in dilutes intensity-weighted scoring.
- **allow AcOH / d3-AcOH / H2O losses** -- loss of d3-acetic acid (63.040) is one of the
  largest features in a per-deuteroacetylated spectrum; without it those peaks look like
  unexplained noise or get mis-assigned.

Instrument calibration drifts. On the published MS2 of *m/z* 1571 every reported
fragment sits about +2.8 mmu high, so either set the tolerance to 5-8 mmu or remove the
offset first with `rgi_io.estimate_offset()` / `apply_offset()`.
