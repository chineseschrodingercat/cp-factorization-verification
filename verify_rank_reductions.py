"""Finite checks for the supplied-factor reduction and weighted tensor construction.

Requires Python 3. The optional --milp cross-check uses NumPy and SciPy/HiGHS.
Run: python verify_rank_reductions.py --milp
"""
from collections import Counter
from itertools import combinations, permutations, product
from math import comb
import sys


def words(g, q, j):
    for i in combinations(range(g), j):
        for a in product(range(q), repeat=j):
            yield tuple(zip(i, a))


def moments(weighted, d):
    out = Counter()
    for w, mass in weighted:
        for j in range(min(d, len(w)) + 1):
            for u in combinations(w, j):
                out[u] += mass
    return out


def oa_check(d, g, q):
    # The parity-check array works over the additive cyclic group, including q=6.
    assert g == d + 1
    array = [tuple(enumerate(x + (sum(x) % q,))) for x in product(range(q), repeat=d)]
    L = 4*g*q
    full = [(s, L**(d-j)) for j in range(d+1) for s in words(g, q, j)]
    short = [(s, 1) for s in array]
    for j in range(d):
        beta = L**(d-j) - (-1)**(d-j)*comb(g-j-1, d-j)*q**(d-j)
        assert beta > 0
        short.extend((s, beta) for s in words(g, q, j))
    M = moments(full, d)
    assert M == moments(short, d)
    for u in M:
        j = len(u)
        if j == d:
            assert M[u] == 1
            continue
        occupied = {h for h,_ in u}
        extension_mass = sum(M[tuple(sorted(u+((h,a),)))]
                             for h in range(g) if h not in occupied for a in range(q))
        assert extension_mass < M[u]
    assert len(full) == sum(comb(g,j)*q**j for j in range(d+1))
    assert len(short) == q**d + sum(comb(g,j)*q**j for j in range(d))
    print(f'weighted d={d}, g={g}, q={q}: {len(M)} exact moments, {len(full)} -> {len(short)}')


def binary_oa_count(d, g):
    # Normalize the first d columns to lexicographic order. Each additional
    # column must be bijective on each binary coordinate fiber.
    data = list(product(range(2), repeat=d))
    good = []
    for mask in product(range(2), repeat=len(data)):
        if all(len({tuple(row[j] for j in range(d) if j!=i)+(v,)
                    for row,v in zip(data,mask)}) == len(data) for i in range(d)):
            good.append(mask)
    cnt = 0
    for extra in product(good,repeat=g-d):
        rows = [row+tuple(col[n] for col in extra) for n,row in enumerate(data)]
        if all(len({tuple(row[j] for j in I) for row in rows}) == len(data)
               for I in combinations(range(g),d)):
            cnt += 1
    return cnt


def partial_squares(q):
    rows = [set() for _ in range(q)]
    cols = [set() for _ in range(q)]
    def visit(p, clues):
        if p == q*q:
            yield tuple(clues)
            return
        r,c = divmod(p,q)
        yield from visit(p+1,clues)
        for a in range(q):
            if a in rows[r] or a in cols[c]:
                continue
            rows[r].add(a)
            cols[c].add(a)
            clues.append((r,c,a))
            yield from visit(p+1,clues)
            clues.pop()
            rows[r].remove(a)
            cols[c].remove(a)
    yield from visit(0,[])


def latin_squares(q):
    for rows in product(permutations(range(q)),repeat=q):
        if all(len({row[c] for row in rows}) == q for c in range(q)):
            yield frozenset((r,c,a) for r,row in enumerate(rows) for c,a in enumerate(row))


def triple_edges(w):
    return tuple((i,w[i],j,w[j]) for i,j in combinations(range(3),2))


def residual(q,clues):
    edges = {(i,a,j,b) for i,j in combinations(range(3),2)
             for a,b in product(range(q),repeat=2)}
    for clue in clues:
        edges.difference_update(triple_edges(clue))
    return edges


def exact_cover(edges, all_triangles):
    if not edges:
        return ()
    first = min(edges)
    for w,tr in all_triangles:
        if first in tr and tr <= edges:
            answer = exact_cover(edges-tr,all_triangles)
            if answer is not None:
                return (w,)+answer
    return None


def reduction_words(q,edges):
    singles = [((i,a),) for i,a in product(range(3),range(q))]
    supplied = [(((i,a),(j,b)),1) for i,a,j,b in edges]
    supplied += [(u,1) for u in singles] + [((),1)]*(len(edges)//3)
    return singles,supplied


def project_gram(q,weighted):
    n = 3*(q+1)
    out = [[0]*n for _ in range(n)]
    for w,mass in weighted:
        word = dict(w)
        I = [h*(q+1)+1+word.get(h,-1) for h in range(3)]
        for i,j in product(I,repeat=2):
            out[i][j] += mass
    return out


def independent_milp(q, gram):
    # Use every possible one-symbol-per-group column, including the empty word.
    # Constraints are the raw Gram entries. All atom weights are continuous.
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    cols = list(product(range(q+1),repeat=3))
    n = len(gram)
    idxs = [(i,j) for i in range(n) for j in range(i,n)]
    targets = np.array([gram[i][j] for i,j in idxs])
    D = np.zeros((len(idxs),len(cols)))
    for l,col in enumerate(cols):
        active = {i*(q+1)+a for i,a in enumerate(col)}
        D[:,l] = [int(i in active and j in active) for i,j in idxs]
    N = len(cols)
    T = sum(gram[i][i] for i in range(q+1))
    cons = [LinearConstraint(np.hstack((D,np.zeros_like(D))),targets,targets),
            LinearConstraint(np.hstack((np.eye(N),-T*np.eye(N))),-np.inf,0)]
    res = milp(np.r_[np.zeros(N),np.ones(N)],integrality=np.r_[np.zeros(N),np.ones(N)],
               bounds=Bounds(np.zeros(2*N),np.r_[np.full(N,T),np.ones(N)]),
               constraints=cons,options={'mip_rel_gap':0.0})
    assert res.success and res.mip_gap < 1e-7 and abs(res.fun-round(res.fun))<1e-7
    return round(res.fun)


def latin_check(q, milp_enabled):
    complete = list(latin_squares(q))
    tris = [(x,set(triple_edges(x))) for x in product(range(q),repeat=3)]
    yes = no = checked = 0
    for clues in partial_squares(q):
        edges = residual(q,clues)
        ans = exact_cover(edges,tris)
        expected = any(set(clues) <= s for s in complete)
        assert (ans is not None) == expected
        yes += expected
        no += not expected
        special = q==2 or clues in [(),((0,0,0),(0,1,1),(1,2,2))]
        if not special:
            continue
        single,supplied = reduction_words(q,edges)
        supplied_gram = project_gram(q,supplied)
        t = q*q-len(clues)
        assert len(supplied) == 4*t+3*q and len(edges)==3*t
        if expected:
            counts = Counter((h,v) for word in ans for h,v in enumerate(word))
            short = [(tuple(enumerate(w)),1) for w in ans]
            short.extend((s,counts[s[0]]+1) for s in single)
            assert len(short)==t+3*q and project_gram(q,short)==supplied_gram
        if milp_enabled:
            kstar = independent_milp(q,supplied_gram)
            assert kstar>=t+3*q and (kstar==t+3*q)==expected
            checked += 1
    print(f'order {q}: {yes} completable and {no} noncompletable partial squares; '
          f'{checked} raw-Gram mixed-integer optimization checks')


if __name__ == '__main__':
    for case in [(2,3,2),(2,3,6),(3,4,2),(4,5,2),(5,6,2)]:
        oa_check(*case)
    for d in (2,3):
        assert binary_oa_count(d,d+1)==2 and binary_oa_count(d,d+2)==0
    print('binary array equality/nonexistence checks passed at strengths 2 and 3')
    for q in (2,3):
        latin_check(q,'--milp' in sys.argv)
