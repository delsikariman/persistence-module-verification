# Persistence Module Verification

This repository contains the code and computational experiments accompanying the manuscript on persistence modules, quiver representations, and persistent homology.

The repository provides the implementation of the computational procedures and experiments reported in Section 5 of the manuscript.

## Contents

The main files in this repository are:

- `verification.ipynb` — a self-contained Jupyter Notebook containing the
  computational experiments and generating the required `qrep.py` and
  `tda.py` modules during execution.
- `requirements.txt` — Python dependencies required to reproduce the
  computational experiments.

## Computational Environment

The computational experiments were performed using Python 3.11.5,
NumPy 1.26.4, and GUDHI 3.13.0.

The required Python dependencies are listed in `requirements.txt`.

## Reproducing the Experiments

Install the required dependencies using:

```bash
pip install -r requirements.txt

Open `verification.ipynb` in Jupyter Notebook or JupyterLab and run all cells
sequentially from top to bottom.

During execution, the notebook automatically generates the `qrep.py` and
`tda.py` modules required for the subsequent computational experiments.
