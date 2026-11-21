import sys, numpy as np
from collections import Counter
sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else ".")
import qrep
from qrep_track import decompose_tracked, check_normal_form

rng = np.random.default_rng(20261116)
tot = ok = 0
for p in (2, 3, 5):
    for _ in range(300):
        n = int(rng.integers(2, 8))
        dims = [int(rng.integers(0, 6)) for _ in range(n)]
        phis = [rng.integers(0, p, size=(dims[i+1], dims[i])).astype(np.int64) for i in range(n-1)]
        bars, bases, chains = decompose_tracked(phis, p, dims)
        a = bars == qrep.decompose(phis, p)
        b = bars == qrep.intervals_from_rank_invariant(phis, p)
        c = check_normal_form(phis, bases, chains, p)
        # chains read off the bases reproduce the barcode
        span = {}
        for j, ids in enumerate(chains, 1):
            for cid in ids:
                lo, hi = span.get(cid, (j, j)); span[cid] = (min(lo, j), max(hi, j))
        d = sorted(span.values()) == bars
        e = all(sorted(Counter(ids).values()) in ([], [1]*len(ids)) for ids in chains)
        tot += 1; ok += a and b and c and d and e
print(f"tracked reduction: {ok}/{tot} random modules passed")
# degenerate cases
print("n=1, dim 3      :", decompose_tracked([], 2, [3])[0])
print("zero module n=3 :", decompose_tracked([np.zeros((0,0)), np.zeros((0,0))], 2, [0,0,0])[0])
print("V_2 = 0         :", decompose_tracked([np.zeros((0,2)), np.zeros((3,0))], 2, [2,0,3])[0])
print("Phi = 0         :", decompose_tracked([np.zeros((2,2))], 2, [2,2])[0])
