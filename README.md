# Persistence Module Verification

Code and computational experiments accompanying the manuscript

> **Representations of Type A_n Quivers and Persistence Modules: Interval Decomposition with Computational Validation on Synthetic Point Clouds**

The repository reproduces every figure, table and numerical value reported in Section 5 of the manuscript.

## Contents

- `verification.ipynb`: a self-contained Jupyter notebook. When run, it writes the modules
  `qrep.py` (representations of the quiver A_n over GF(p), and the reduction algorithm) and
  `tda.py` (persistence modules of Vietoris–Rips filtrations built without any
  persistent-homology library), runs all experiments, and saves Figures 1–4 of the manuscript
  to `figures/` as vector PDF and 600 dpi PNG.
- `requirements.txt`: Python dependencies.

## Computational environment

The results were produced with Python 3.11.5, NumPy 1.26.4, GUDHI 3.13.0 and Matplotlib.

## Reproducing the results

```bash
pip install -r requirements.txt
```

Open `verification.ipynb` in Jupyter Notebook or JupyterLab and run all cells from top to bottom
(about three minutes; Experiment 2 takes the longest). All random generators are seeded, one per
experiment, so the reported values are reproduced exactly. The final summary cell reports
`all checks passed: True`.
