"""Verify the private-summand ordering convention using exact rational arithmetic."""
from fractions import Fraction
from verify import WORDS, lift_test, model, outer, parity


def private_order(matrix, atoms):
    rows, cols, base_weight, blocks, lift, rays, _ = model(matrix)
    size, width, _, core = lift_test(matrix, atoms)
    assert (size, width, base_weight) == (78, 98, 8)
    first = rows + cols
    raw, normalized = {}, {}
    for block_id, block in enumerate(blocks):
        i = block[0]
        if len(block) == 2:
            offset = -core[i][block[1]] / 2
        else:
            side = range(rows) if i < rows else range(rows, rows + cols)
            offset = (sum(core[i][j] for j in side) - 2 * core[i][i]) / 2
        for word in WORDS:
            ray = rays[block_id, word]
            vector = [sum(value * ray.get(j, 0) for j, value in row.items()) for row in lift]
            weight = base_weight + offset * parity(word)
            trace = weight * sum(x*x for x in vector)
            diagonal = tuple(weight*x*x for x in vector[first:first+8*len(blocks)])
            raw[block_id, word] = diagonal
            normalized[block_id, word] = tuple(value / trace for value in diagonal)
    assert len(set(normalized.values())) == 16 * len(blocks)
    return (sorted(raw, key=raw.get, reverse=True),
            sorted(normalized, key=normalized.get, reverse=True),
            raw[0, (0, 0, 0, 0)][0], raw[0, (0, 0, 0, 1)][0])


if __name__ == "__main__":
    matrix = ((1, 1), (2, 3))
    row_terms = (outer((1, 0), matrix[0]), outer((0, 1), matrix[1]))
    col_terms = (outer((1, 2), (1, 0)), outer((1, 3), (0, 1)))
    by_row = private_order(matrix, row_terms)
    by_col = private_order(matrix, col_terms)
    assert by_row[2:] == (7, 9)
    assert by_col[2:] == (Fraction(101, 12), Fraction(91, 12))
    assert by_row[0] != by_col[0] and by_row[1] == by_col[1]
    print("Both 78-by-78 lifts agree. Trace normalization fixes all 96 private slots.")
