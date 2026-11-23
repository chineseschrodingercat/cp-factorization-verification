"""Exact finite checks for the completely positive factorization examples."""
from fractions import Fraction
from itertools import combinations, product
from math import comb


def indicator(word, sizes):
    return tuple(int(a == selected)
                 for selected, size in zip(word, sizes) for a in range(size))


def gram(weighted_vectors):
    weighted_vectors = list(weighted_vectors)
    width = len(weighted_vectors[0][0])
    return [[sum(w * v[i] * v[j] for v, w in weighted_vectors)
             for j in range(width)] for i in range(width)]


def rational_rank(matrix):
    a = [list(map(Fraction, row)) for row in matrix]
    count = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(count, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[count], a[pivot] = a[pivot], a[count]
        for i in range(count + 1, len(a)):
            scale = a[i][j] / a[count][j]
            a[i] = [x - scale * y for x, y in zip(a[i], a[count])]
        count += 1
    return count


def deletion_example(g, q, expected_widths, expected_rank):
    # Both test cases use a prime alphabet. The empty symbol has label q.
    sizes = (q + 1,) * g
    empty = (q,) * g
    empty_weight = comb(g - 1, 2) * q**2
    capped = [(empty, empty_weight)]
    for active in combinations(range(g), 2):
        for labels in product(range(q), repeat=2):
            word = list(empty)
            for i, a in zip(active, labels):
                word[i] = a
            capped.append((tuple(word), 1))
    projected = []
    for a, b in product(range(q), repeat=2):
        symbols = (a, b) + tuple((a + j*b) % q for j in range(1, g - 1))
        projected.append((symbols, 1))
    for i, a in product(range(g), range(q)):
        word = list(empty)
        word[i] = a
        projected.append((tuple(word), (g - 2) * q))
    lifted = gram((indicator(s, sizes) + (2 - sum(x != q for x in s),), w)
                  for s, w in capped)
    short = gram((indicator(s, sizes), w) for s, w in projected)
    assert (len(capped), len(projected)) == expected_widths
    assert lifted[-1][-1] == 4 * empty_weight
    assert [row[:-1] for row in lifted[:-1]] == short
    assert rational_rank(lifted) == rational_rank(short) == expected_rank


def pair_margins(table):
    result = {}
    for word, mass in table.items():
        for i, j in combinations(range(3), 2):
            key = (i, j, word[i], word[j])
            result[key] = result.get(key, 0) + mass
    return result


def parity_example():
    words = list(product(range(2), repeat=3))
    for theta, size in [(0, 4), (Fraction(1, 3), 8), (1, 4)]:
        table = {w: theta if sum(w) % 2 == 0 else 1 - theta for w in words}
        assert sum(x > 0 for x in table.values()) == size
        assert set(pair_margins(table).values()) == {1}
    equations = [[int(w[i] == a and w[j] == b) for w in words]
                 for i, j in combinations(range(3), 2)
                 for a, b in product(range(2), repeat=2)]
    assert rational_rank(equations) == 7


def extension_examples():
    # One negative cell cuts the three-bit segment to its odd endpoint.
    words = list(product(range(2), repeat=3))
    rows = {s: indicator(s, (2,) * 3) for s in words}
    base = gram((rows[s], 1) for s in words if sum(s) % 2 == 0)
    lifted = gram((rows[s] + (sum(s) - 1,), 1) for s in words if sum(s) % 2)
    assert [row[:-1] for row in lifted[:-1]] == base
    assert rational_rank(lifted) == rational_rank(base) == 4
    assert [s for s in words if sum(s) - 1 < 0] == [(0, 0, 0)]

    for g in (5, 6):
        words = list(product(range(2), repeat=g))
        even = [s for s in words if sum(s) % 2 == 0]
        odd = [s for s in words if sum(s) % 2]
        rows = {s: indicator(s, (2,) * g) for s in words}
        masks = [()] + [(i,) for i in range(g)] + list(combinations(range(g), 2))
        phi = {s: tuple((-1) ** sum(s[i] for i in mask) for mask in masks)
               for s in words}
        assert rational_rank([phi[e] for e in even]) == len(masks)
        if g == 5:
            assert len(masks) == len(even) == 16
        scale = 2 * len(masks)
        target = gram((rows[e], scale) for e in even)
        for u in odd:
            witness = [(rows[e], scale - sum(a*b for a, b in zip(phi[e], phi[u])))
                       for e in even]
            assert all(weight > 0 for _, weight in witness)
            witness.append((rows[u], len(even)))
            assert gram(witness) == target
        for u, v in combinations(odd, 2):
            i = next(i for i in range(g) if u[i] != v[i])
            a = u[:i] + (v[i],) + u[i+1:]
            b = v[:i] + (u[i],) + v[i+1:]
            assert a in even and b in even
            assert all(x + y == z + w
                       for x, y, z, w in zip(rows[u], rows[v], rows[a], rows[b]))
        if g == 6:
            missing = even[0]
            assert rational_rank([phi[e] for e in even if e != missing]) == len(masks)
            protected = [u for u in odd if sum(a != b for a, b in zip(u, missing)) > 1]
            assert len(protected) == len(odd) - g
    print("One-coordinate segment lift and all single-odd-word witnesses at g=5,6 passed.")


def check_prefix_inequalities():
    words = list(product(range(2), repeat=3))
    for mask in range(1, 1 << len(words)):
        support = {w for i, w in enumerate(words) if mask & (1 << i)}
        rows = {}
        for u in support:
            rows[u] = []
            for i in range(len(u)):
                branches = any(v[:i] == u[:i] and v[i] != u[i] for v in support)
                rows[u].append((1 if branches else -1) * (1 << (len(u) - 1 - i)))
        admitted = {x for x in words if all(
            sum(w for w, a, b in zip(coefficients, x, u) if a != b) >= 0
            for u, coefficients in rows.items())}
        assert admitted == support


def incidence_witnesses(q):
    # Arithmetic modulo q is a field for the tested cases q=2,3.
    n = q*q
    edges = {(q*x + (a*x+b) % q, n+q*a+b)
             for x, a, b in product(range(q), repeat=3)}
    assert len(edges) == q**3
    seed = {frozenset([i]) for i in range(2*n)}
    seed.update(frozenset((i, j)) for i in range(n) for j in range(n, 2*n)
                if (i, j) not in edges)
    for i in range(2*n):
        adjacent = [y if x == i else x for x, y in edges if i in (x, y)]
        seed.update(frozenset((i, a, b)) for a, b in combinations(adjacent, 2))
    if q == 2:
        terms = [()] + [(i,) for i in range(2*n)] + list(combinations(range(2*n), 2))
        moments = [[int(all(i in word for i in term)) for term in terms] for word in seed]
        assert rational_rank(moments) == len(seed)
    for (i, j), (k, l) in combinations(edges, 2):
        if i == k:
            a, b = {i}, {i, j, l}
        elif j == l:
            a, b = {j}, {i, k, j}
        elif (i, l) in edges:
            a, b = {k}, {i, j, l}
        elif (k, j) in edges:
            a, b = {i}, {k, j, l}
        else:
            a, b = {i, l}, {k, j}
        assert frozenset(a) in seed and frozenset(b) in seed
        assert sorted((*a, *b)) == sorted((i, j, k, l))


def compatibility_example(facets, n):
    f = len(facets)
    words = list(product(range(f), range(2), range(2)))
    def marker(i, word):
        a, b, c = word
        return int(i in facets[a]) + b + c - 1
    rows = {w: indicator(w, (f, 2, 2)) + tuple(marker(i, w) for i in range(n))
            for w in words}
    masks = [frozenset(combo) for k in range(n+1) for combo in combinations(range(n), k)]
    support_choices = {j: set() for j in masks}
    target = None
    for corner in range(f):
        table = {}
        for a, b, c in words:
            t = int(a == corner)
            table[a, b, c] = ((t, 1-t), (1-t, f-2+t))[b][c]
        assert sum(table.values()) == f*f and min(table.values()) >= 0
        pair = pair_margins(table)
        # Verify pairwise independence directly from all one-way counts.
        for (i, j, a, b), mass in pair.items():
            ai = sum(v for w, v in table.items() if w[i] == a)
            bj = sum(v for w, v in table.items() if w[j] == b)
            assert mass * f*f == ai * bj
        full = gram((rows[w], weight) for w, weight in table.items())
        if target is None:
            target = full
        assert full == target and all(x >= 0 for row in full for x in row)
        for j in masks:
            if all(marker(i, w) >= 0 for w in words if table[w] for i in j):
                support_choices[j].add(corner)
        for i in range(n):
            bad = sum(weight for w, weight in table.items() if marker(i, w) < 0)
            assert bad == int(i not in facets[corner])
    assert all(support_choices[j] == {a for a, face in enumerate(facets) if j <= face}
               for j in masks)
    base = [row[:f+4] for row in target[:f+4]]
    assert rational_rank(target) == rational_rank(base)


def binary_helly_example(g):
    words = list(product(range(2), repeat=g))
    zero = (0,) * g
    exposed = [w for w in words if sum(w) in (1, 2)]
    outside = set(words) - set(exposed)
    terms = [()] + [(i,) for i in range(g)] + list(combinations(range(g), 2))
    values = {w: [int(all(w[i] for i in term)) for term in terms] for w in words}
    target = [2 * sum(values[w][i] for w in exposed) + values[zero][i]
              for i in range(len(terms))]
    assert len(outside) == 2**g - len(terms) + 1
    for u in outside:
        h = sum(u)
        alpha = Fraction(2, (h-1)*(h-2))
        weights = {}
        for z in exposed:
            active = [i for i, bit in enumerate(z) if bit]
            if len(active) == 1:
                beta = alpha * u[active[0]] * (h-2)
            else:
                beta = -alpha * u[active[0]] * u[active[1]]
            weights[z] = 2 + beta
        weights[u] = alpha
        assert all(w > 0 for w in weights.values())
        assert all(sum(wt*values[word][i] for word, wt in weights.items()) == t
                   for i, t in enumerate(target))
        assert sum(wt*(sum(w)-1)*(sum(w)-2) for w, wt in weights.items()) == 2


def two_bin_minimum(weights):
    # With two bins a table is described by one interval variable per item:
    # x_i lies in [0,2a_i] and all x_i sum to the prescribed bin-bit count.
    # Every vertex has at most one coordinate strictly between its endpoints.
    # Every feasible support contains a vertex support (intersect with its zero face).
    n = len(weights)
    target = sum(weights)
    reference = {(i, j, b): a for i, a in enumerate(weights)
                 for j, b in product(range(2), repeat=2)}
    margins = pair_margins(reference)
    best = 4 * n
    for pivot in range(n):
        others = [i for i in range(n) if i != pivot]
        for signs in product((0, 1), repeat=len(others)):
            mass = [0] * n
            for i, sign in zip(others, signs):
                mass[i] = sign * 2 * weights[i]
            mass[pivot] = target - sum(mass)
            if not 0 <= mass[pivot] <= 2 * weights[pivot]:
                continue
            table = {(i, j, b): x if j != b else 2 * a - x
                     for i, (x, a) in enumerate(zip(mass, weights))
                     for j, b in product(range(2), repeat=2)}
            assert pair_margins(table) == margins
            best = min(best, sum(v > 0 for v in table.values()))
    return best


def three_bins():
    weights = [5, 5, 8, 5, 6, 7, 6, 6, 6]
    h = 3
    supplied = {(i, j, bit): a if bit else (h - 1) * a
                for i, a in enumerate(weights)
                for j in range(h) for bit in (0, 1)}
    witness = {(i, j, bit): h * a * int(bit == int(i // 3 == j))
               for i, a in enumerate(weights)
               for j in range(h) for bit in (0, 1)}
    assert sum(supplied.values()) == h ** 3 * 18
    assert pair_margins(supplied) == pair_margins(witness)
    assert sum(v > 0 for v in witness.values()) == h * len(weights)
    assert all(v > 0 for v in supplied.values())


def rational_square_roots(weight):
    """Return nonnegative rational roots whose squared sum is the weight."""
    from fractions import Fraction
    weight = Fraction(weight)
    assert weight >= 0
    numerator = weight.numerator * weight.denominator
    for bit in range(numerator.bit_length()):
        if (numerator >> bit) & 1:
            root = Fraction(1 << (bit // 2), weight.denominator)
            yield root
            if bit % 2:
                yield root


def supplied_factor_example(B, check_binary=False):
    """Compare both factors to the defining block matrix of Theorem 24."""
    from math import isqrt
    from collections import defaultdict

    n = len(B)
    g = max(max(row) for row in B)
    alpha = isqrt(int(n * g)) + 1
    groups = [tuple(range(n)), tuple(range(n, 2 * n))]
    attached = [(core, alpha) for group in groups for core in combinations(group, 2)]
    attached.extend(((i,), 2 * n * alpha) for i in range(2 * n))
    size = 2 * n + 4 * len(attached)
    target = [[0] * size for _ in range(size)]
    scales = [1] * size
    columns = []

    for group in groups:
        for i, j in product(group, repeat=2):
            target[i][j] = (8 * n ** 2 if i == j else 2) * alpha ** 2
    for i, j in product(range(n), repeat=2):
        target[i][n + j] = target[n + j][i] = B[i][j]

    for pos, (core, gamma) in enumerate(attached):
        start = 2 * n + 4 * pos
        for e in range(4):
            v, w = start + e, start + (e + 1) % 4
            target[v][v] = 2
            target[v][w] = target[w][v] = 1
            for i in core:
                target[i][v] = target[v][i] = gamma
            coeff = alpha * int(e % 2 == 0) if len(core) == 2 else n * alpha
            scales[v] = alpha if len(core) == 2 else n * alpha
            column = {v: 1, w: 1}
            if coeff:
                column.update({i: coeff for i in core})
            columns.append(column)

    diag = (4 * n ** 2 - 2 * n + 2) * alpha ** 2
    weights = {(i, n + j): B[i][j] for i, j in product(range(n), repeat=2)}
    for i in range(n):
        weights[(i,)] = diag - sum(B[i])
        weights[(n + i,)] = diag - sum(row[i] for row in B)
    for support, weight in weights.items():
        for root in rational_square_roots(weight):
            columns.append(dict.fromkeys(support, root))
    accumulated = defaultdict(int)
    for column in columns:
        for i, j in product(column, repeat=2):
            accumulated[i, j] += column[i] * column[j]
    assert all(accumulated[i, j] == target[i][j]
               for i, j in product(range(size), repeat=2))

    if check_binary:
        # Expand every cycle and residual weight into explicit 0/1 columns.
        binary_supports = []
        for column in columns[:4 * len(attached)]:
            levels = {scales[i] * v for i, v in column.items()}
            assert len(levels) == 1
            level = levels.pop()
            binary_supports.extend([tuple(column)] * level ** 2)
        for support, weight in weights.items():
            assert int(weight) == weight
            binary_supports.extend([support] * int(weight))
        counts = defaultdict(int)
        for support in binary_supports:
            assert len(support) <= 4
            for i, j in product(support, repeat=2):
                counts[i, j] += 1
        assert all(counts[i, j] == target[i][j] * scales[i] * scales[j]
                   for i, j in product(range(size), repeat=2))


if __name__ == '__main__':
    import argparse
    argparse.ArgumentParser(description=__doc__).parse_args()
    deletion_example(3, 3, (28, 18), 10)
    deletion_example(5, 5, (251, 50), 26)
    parity_example()
    extension_examples()
    check_prefix_inequalities()
    for q in (2, 3):
        incidence_witnesses(q)
    for facets, n in [([{0}], 1), ([{0}, {1}, {2}], 3),
                      ([{0, 1}, {0, 2}, {1, 2}], 3),
                      ([{0, 1}, {1, 2}, {0, 2}, {3}], 4)]:
        compatibility_example(facets, n)
    for g in (4, 5, 6):
        binary_helly_example(g)
    print("Checked 255 cube subsets, two incidence graphs and four compatibility complexes.")
    print("Binary Helly witnesses at g=4,5,6 passed.")
    assert two_bin_minimum([5, 5, 5, 6, 7, 8]) == 12
    assert two_bin_minimum([5, 5, 5, 7, 7, 7]) == 14
    three_bins()
    supplied_factor_example([[0]], check_binary=True)
    supplied_factor_example([[1, 0], [2, 1]], check_binary=True)
    supplied_factor_example([[0, 3, 1], [2, 1, 0], [0, 4, 2]], check_binary=True)
    supplied_factor_example([[Fraction(3, 7), 0], [Fraction(11, 5), Fraction(1, 9)]])
    print("Supplied rational and sparse binary factors reproduce their target matrices.")
    print("Exact Gram factors and ordinary ranks: 28 -> 18 at rank 10, 251 -> 50 at rank 26.")
    print("Exact checks passed: parity segment;")
    print("two-bin optima 12 and 14; three-bin product-law witness on 27 cells.")
