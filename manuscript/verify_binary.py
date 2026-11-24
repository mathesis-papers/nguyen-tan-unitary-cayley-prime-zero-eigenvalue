import json
from collections import Counter, deque
from fractions import Fraction
from math import comb


def cube_neighbors(dimension):
    return [
        [vertex ^ (1 << coordinate) for coordinate in range(dimension)]
        for vertex in range(1 << dimension)
    ]


def connected(adjacency):
    reached = {0}
    frontier = deque([0])
    while frontier:
        for neighbor in adjacency[frontier.popleft()]:
            if neighbor not in reached:
                reached.add(neighbor)
                frontier.append(neighbor)
    return len(reached) == len(adjacency)


def homogeneous_sets(adjacency):
    size = len(adjacency)
    neighbor_masks = [sum(1 << neighbor for neighbor in row) for row in adjacency]
    tested = 0
    modules = []
    for subset in range(1, (1 << size) - 1):
        if subset.bit_count() < 2:
            continue
        tested += 1
        homogeneous = True
        for vertex, neighbors in enumerate(neighbor_masks):
            if subset & (1 << vertex):
                continue
            intersection = neighbors & subset
            if intersection != 0 and intersection != subset:
                homogeneous = False
                break
        if homogeneous:
            modules.append([vertex for vertex in range(size) if subset & (1 << vertex)])
    return tested, modules


def exact_rank(matrix):
    rows = [[Fraction(entry) for entry in row] for row in matrix]
    rank = 0
    for column in range(len(rows[0])):
        pivot = next((index for index in range(rank, len(rows)) if rows[index][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][column]
        rows[rank] = [value / scale for value in rows[rank]]
        for index in range(rank + 1, len(rows)):
            multiple = rows[index][column]
            if multiple:
                rows[index] = [
                    value - multiple * pivot_value
                    for value, pivot_value in zip(rows[index], rows[rank])
                ]
        rank += 1
    return rank


def ring_multiply(left, right, dimension):
    product = 0
    for left_index in range(dimension):
        if not left & (1 << left_index):
            continue
        for right_index in range(dimension):
            if right & (1 << right_index):
                product ^= 1 << ((left_index + right_index) % dimension)
    return product


def ring_check():
    dimension = 4
    size = 1 << dimension
    multiply = [[ring_multiply(left, right, dimension) for right in range(size)] for left in range(size)]
    for left in range(size):
        assert multiply[left][1] == left
        for right in range(size):
            assert multiply[left][right] == multiply[right][left]
            for third in range(size):
                assert multiply[multiply[left][right]][third] == multiply[left][multiply[right][third]]
                assert multiply[left][right ^ third] == multiply[left][right] ^ multiply[left][third]
    units = [left for left in range(size) if 1 in multiply[left]]
    subgroup = [1 << coordinate for coordinate in range(dimension)]
    assert all(value in units for value in subgroup)
    assert all(multiply[left][right] in subgroup for left in subgroup for right in subgroup)
    assert 1 in subgroup
    witnesses = []
    for value in range(1, size):
        coordinate = next(index for index in range(dimension) if value & (1 << index))
        multiplier = 1 << ((-coordinate) % dimension)
        assert multiply[multiplier][value] & 1
        witnesses.append({"element": value, "multiplier": multiplier})
    adjacency = [[right for right in range(size) if left ^ right in subgroup] for left in range(size)]
    assert adjacency == [sorted(row) for row in cube_neighbors(dimension)]
    matrix = [[int(right in adjacency[left]) for right in range(size)] for left in range(size)]
    rank = exact_rank(matrix)
    assert rank == 10
    vector = [(-1) ** ((vertex & 3).bit_count()) for vertex in range(size)]
    assert all(sum(vector[neighbor] for neighbor in row) == 0 for row in adjacency)
    return {
        "ring": "F_2[x]/(x^4-1)",
        "encoding": "bit i is the coefficient of x^i; identity=1; addition=XOR",
        "elements": size,
        "units": units,
        "basis_unit_subgroup": subgroup,
        "associativity_triples": size ** 3,
        "distributivity_triples": size ** 3,
        "frobenius_functional": "coefficient of x^0",
        "frobenius_witnesses": witnesses,
        "adjacency_rank_over_Q": rank,
        "adjacency_nullity": size - rank,
        "zero_eigenvector": vector,
    }


def cube_check(dimension):
    adjacency = cube_neighbors(dimension)
    size = len(adjacency)
    complement = [
        [neighbor for neighbor in range(size) if neighbor != vertex and neighbor not in adjacency[vertex]]
        for vertex in range(size)
    ]
    spectrum = Counter()
    for character_mask in range(size):
        eigenvalue = dimension - 2 * character_mask.bit_count()
        spectrum[eigenvalue] += 1
        values = [(-1) ** ((vertex & character_mask).bit_count()) for vertex in range(size)]
        for vertex, row in enumerate(adjacency):
            assert sum(values[neighbor] for neighbor in row) == eigenvalue * values[vertex]
    assert spectrum == Counter({dimension - 2 * index: comb(dimension, index) for index in range(dimension + 1)})
    result = {
        "dimension": dimension,
        "vertices": size,
        "connected": connected(adjacency),
        "complement_connected": connected(complement),
        "character_eigenvector_equations": size ** 2,
        "eigenvalue_multiplicities": {str(value): count for value, count in sorted(spectrum.items())},
    }
    assert result["connected"]
    assert result["complement_connected"] == (dimension >= 3)
    if dimension <= 4:
        tested, modules = homogeneous_sets(adjacency)
        result.update(nontrivial_subsets_tested=tested, homogeneous_sets=modules)
        if dimension >= 3:
            assert modules == []
        if dimension == 2:
            assert modules == [[1, 2], [0, 3]]
    if dimension == 4:
        assert spectrum[0] == 6
    return result


def main():
    result = {
        "arithmetic": "exact integers and rational Gaussian elimination",
        "cube_checks": [cube_check(dimension) for dimension in range(1, 9)],
        "cyclic_group_algebra_check": ring_check(),
        "scope": "Finite checks for dimensions 1 through 8; exhaustive homogeneous-set checks for dimensions 1 through 4; the general family has the written proof in note.tex.",
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
