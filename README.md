# Persistence Module Verification

Code and computational experiments accompanying the manuscript

> **Representations of Type A_n Quivers and Persistence Modules: Interval Decomposition with Computational Validation on Synthetic Point Clouds**

The repository reproduces every figure, table and numerical value reported in Section 5 of the manuscript.

## Contents

- `verification.ipynb`: a self-contained Jupyter notebook. When run, it writes three modules,
  runs all experiments, and saves Figures 1–4 of the manuscript to `figures/` as vector PDF
  and 600 dpi PNG.
- `requirements.txt`: Python dependencies.

The modules written by the notebook are:

| Module | Content |
| --- | --- |
| `qrep.py` | Representations of the quiver A_n over GF(p): the quadratic form and its roots, Algorithm 1 with `track = false`, the rank invariant |
| `qrep_track.py` | Algorithm 1 with `track = true`: the barcode together with bases in normal form (Proposition 3) |
| `tda.py` | Persistence modules of Vietoris–Rips filtrations, built without any persistent-homology library |

## Map from the manuscript to the code

| Item of the manuscript | Code |
| --- | --- |
| Proposition 1, enumeration of roots | `qrep.positive_roots_bruteforce` |
| Algorithm 1, `track = false` | `qrep.decompose` |
| Algorithm 1, `track = true` (Proposition 3) | `qrep_track.decompose_tracked`, `qrep_track.check_normal_form` |
| Route (B), inclusion–exclusion on the rank invariant | `qrep.rank_invariant`, `qrep.intervals_from_rank_invariant` |
| Route (A), construction of the modules | `tda.h0_module`, `tda.h1_module` |
| Discretization map of Section 5 | `tda.discretise_barcode` |
| Random modules of Section 5.1 | `qrep.random_module` |

## Computational environment

The results were produced with Python 3.11.5, NumPy 1.26.4, GUDHI 3.13.0 and Matplotlib.
The notebook was also re-executed with Python 3.13, NumPy 2.5 and GUDHI 3.13.0, with identical results.

## Reproducing the results

```bash
pip install -r requirements.txt
```

Open `verification.ipynb` in Jupyter Notebook or JupyterLab and run all cells from top to bottom
(about three minutes; Experiment 2 takes the longest). The final summary cell lists each check
as PASS or FAIL and ends with `all checks passed: True`.

## Random seeds

Each experiment has its own seeded generator, so the reported values are reproduced exactly
and do not depend on the order in which the experiments are run.

| Experiment | Seed |
| --- | --- |
| Experiment 2, noisy circle | 20260002 |
| Experiment 3, two clusters | 20260003 |
| Random modules of Section 5.1 | 20260007 |
| Test of the tracked reduction (Section 2b of the notebook) | 20261116 |

## License

See `LICENSE`.
