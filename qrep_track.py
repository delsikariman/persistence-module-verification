"""
qrep_track.py -- Algorithm 1 with track = true (barcode basis), over GF(p).

Companion to qrep.py.  `decompose_tracked` follows Algorithm 1 of the paper line
by line, including the optional bookkeeping of Proposition 3: in pass i the
unitriangular matrix T is recorded, and the bases B_1, ..., B_i are updated by
the leading principal submatrices T[S_j, S_j].  It returns the barcode together
with bases B_1, ..., B_n in which every structure matrix is in the normal form
of Corollary 1.  `check_normal_form` verifies that claim independently.

Conventions as in qrep.py: phis[i-1] is the matrix of phi_i : V_i -> V_{i+1},
of shape (dim V_{i+1}, dim V_i); intervals [b, d] are 1-based.  Because a list
of matrices cannot encode the dimensions when n = 1 or when a space is zero,
the dimension vector may be passed explicitly through `dims`.
"""
from __future__ import annotations

import numpy as np

__all__ = ["decompose_tracked", "check_normal_form"]


def _inv(a: int, p: int) -> int:
    return pow(int(a) % p, p - 2, p)


def _inv_matrix(A: np.ndarray, p: int) -> np.ndarray:
    d = A.shape[0]
    M = np.concatenate([A % p, np.eye(d, dtype=np.int64)], axis=1)
    row = 0
    for col in range(d):
        piv = next((r for r in range(row, d) if M[r, col]), None)
        if piv is None:
            raise ValueError("matrix is singular")
        M[[row, piv]] = M[[piv, row]]
        M[row] = (M[row] * _inv(M[row, col], p)) % p
        for r in range(d):
            if r != row and M[r, col]:
                M[r] = (M[r] - M[r, col] * M[row]) % p
        row += 1
    return M[:, d:]


def decompose_tracked(phis, p: int = 2, dims=None):
    """Algorithm 1 with track = true.

    Returns (bars, bases, chains):
      bars   -- sorted list of intervals (b, d), with multiplicity;
      bases  -- list of n invertible matrices; the columns of bases[j-1] are
                the vectors of B_j in the input coordinates of V_j;
      chains -- chains[j-1][c] is the identifier of the chain (= summand) to
                which column c of bases[j-1] belongs.
    """
    phis = [np.asarray(M, dtype=np.int64) % p for M in phis]
    if dims is None:
        if not phis:
            raise ValueError("pass dims=[dim V_1] for a module of length 1")
        dims = [phis[0].shape[1]] + [M.shape[0] for M in phis]
    n = len(dims)
    assert len(phis) == n - 1

    # line 1
    bases = [np.eye(dims[0], dtype=np.int64)]
    chains = [list(range(dims[0]))]
    birth = {c: 1 for c in range(dims[0])}          # chain id -> birth index
    next_id = dims[0]
    bars = []

    for i in range(1, n):                           # line 2: pass i
        Phi, m, m_next = phis[i - 1], dims[i - 1], dims[i]
        # line 3: order B_i by nondecreasing birth (stable sort)
        order = sorted(range(m), key=lambda c: birth[chains[i - 1][c]])
        U = bases[i - 1][:, order]
        ids = [chains[i - 1][c] for c in order]
        T = np.eye(m, dtype=np.int64)
        R = []                                      # triples (r, v, a)
        for k in range(m):                          # lines 4-17
            w = (Phi @ U[:, k]) % p
            for (r, v, a) in R:
                if w[r]:
                    c = (w[r] * _inv(v[r], p)) % p
                    w = (w - c * v) % p
                    T[:, k] = (T[:, k] - c * T[:, a]) % p      # line 9
            if not w.any():
                bars.append((birth[ids[k]], i))     # line 13
            else:
                r = int(np.nonzero(w)[0][0])
                R.append((r, w, k))                 # line 15
        # line 18: B_{i+1} = images, completed by unit vectors off the pivots
        pivots = {r for (r, _, _) in R}
        cols, cids = [v for (_, v, _) in R], [ids[a] for (_, _, a) in R]
        for s in range(m_next):
            if s not in pivots:
                e = np.zeros(m_next, dtype=np.int64)
                e[s] = 1
                cols.append(e)
                cids.append(next_id)
                birth[next_id] = i + 1
                next_id += 1
        bases.append(np.stack(cols, axis=1) if cols
                     else np.zeros((m_next, 0), dtype=np.int64))
        chains.append(cids)
        # lines 19-23: materialise the compensating changes of basis (Prop. 3)
        for j in range(1, i + 1):
            s = sum(1 for c in ids if birth[c] <= j)          # |S_j|, an initial segment
            if s == 0:
                continue
            pos = [chains[j - 1].index(c) for c in ids[:s]]   # columns of B_j on these chains
            bases[j - 1][:, pos] = (bases[j - 1][:, pos] @ T[:s, :s]) % p

    for c in chains[n - 1]:                         # line 25
        bars.append((birth[c], n))
    return sorted(bars), bases, chains


def check_normal_form(phis, bases, chains, p: int = 2) -> bool:
    """True iff every B_{i+1}^{-1} Phi_i B_i is a 0-1 matrix with at most one
    nonzero entry per row and per column, linking exactly the columns that lie
    on the same chain."""
    for i, Phi in enumerate(phis):
        Phi = np.asarray(Phi, dtype=np.int64) % p
        Bi, Bn = bases[i], bases[i + 1]
        if Bn.shape[0] == 0 or Bi.shape[1] == 0:
            continue
        X = (_inv_matrix(Bn, p) @ Phi @ Bi) % p
        want = np.array([[int(a == b) for b in chains[i]] for a in chains[i + 1]],
                        dtype=np.int64).reshape(X.shape)
        if not np.array_equal(X, want):
            return False
    return True
