import itertools
import json
from collections import Counter
from math import comb, prod

import verify_binary


def hamming_check(alphabet, dimension):
    vertices = list(itertools.product(range(alphabet), repeat=dimension))
    lookup = {vertex: index for index, vertex in enumerate(vertices)}
    adjacency = []
    for vertex in vertices:
        neighbors = []
        for coordinate in range(dimension):
            for symbol in range(alphabet):
                if symbol != vertex[coordinate]:
                    changed = list(vertex)
                    changed[coordinate] = symbol
                    neighbors.append(lookup[tuple(changed)])
        adjacency.append(neighbors)
    spectrum = Counter()
    for factors in itertools.product(range(alphabet), repeat=dimension):
        eigenvalue = dimension * (alphabet - 1) - alphabet * sum(factor != 0 for factor in factors)
        values = [prod(1 if factor == 0 else int(symbol == factor) - int(symbol == 0)
                       for symbol, factor in zip(vertex, factors)) for vertex in vertices]
        assert any(values)
        for index, neighbors in enumerate(adjacency):
            assert sum(values[neighbor] for neighbor in neighbors) == eigenvalue * values[index]
        spectrum[eigenvalue] += 1
    expected = {dimension * (alphabet - 1) - alphabet * index:
                comb(dimension, index) * (alphabet - 1) ** index for index in range(dimension + 1)}
    assert spectrum == Counter(expected)
    universe = set(range(len(vertices)))
    complement = [sorted(universe - set(neighbors) - {index})
                  for index, neighbors in enumerate(adjacency)]
    assert verify_binary.connected(adjacency)
    anti_connected = verify_binary.connected(complement)
    assert anti_connected == (dimension >= 2 and (alphabet, dimension) != (2, 2))
    nullity = spectrum[0]
    assert bool(nullity) == (dimension % alphabet == 0)
    if anti_connected:
        assert 2 * nullity < len(vertices)
    result = {'alphabet': alphabet, 'dimension': dimension, 'vertices': len(vertices),
              'connected': True, 'complement_connected': anti_connected,
              'integer_tensor_eigenvector_equations': len(vertices) ** 2,
              'eigenvalue_multiplicities': {str(value): count for value, count in sorted(spectrum.items())},
              'nullity': nullity, 'nullity_less_than_half': 2 * nullity < len(vertices)}
    if len(vertices) <= 16:
        tested, modules = verify_binary.homogeneous_sets(adjacency)
        result.update(nontrivial_subsets_tested=tested, homogeneous_sets=modules)
        if anti_connected:
            assert not modules
    if (alphabet, dimension) == (3, 3):
        matrix = [[int(column in adjacency[row]) for column in range(len(vertices))]
                  for row in range(len(vertices))]
        rank = verify_binary.exact_rank(matrix)
        assert rank == 15 and nullity == 12
        result['adjacency_rank_over_Q'] = rank
    return result


def cyclic_ternary_ring():
    alphabet, dimension = 3, 3
    vertices = list(itertools.product(range(alphabet), repeat=dimension))
    lookup = {vertex: index for index, vertex in enumerate(vertices)}

    def multiply(left, right):
        coefficients = [0] * dimension
        for first, first_value in enumerate(left):
            for second, second_value in enumerate(right):
                coordinate = (first + second) % dimension
                coefficients[coordinate] = (coefficients[coordinate] + first_value * second_value) % alphabet
        return tuple(coefficients)

    multiplication = [[lookup[multiply(left, right)] for right in vertices] for left in vertices]
    identity = lookup[(1, 0, 0)]
    subgroup = [index for index, vertex in enumerate(vertices) if sum(value != 0 for value in vertex) == 1]
    assert len(subgroup) == 6 and identity in subgroup
    for left in subgroup:
        assert any(multiplication[left][right] == identity for right in subgroup)
        for right in subgroup:
            assert multiplication[left][right] in subgroup
    for left, right, third in itertools.product(range(27), repeat=3):
        assert multiplication[multiplication[left][right]][third] == multiplication[left][multiplication[right][third]]
    witnesses = []
    for index, vertex in enumerate(vertices):
        if not any(vertex):
            continue
        witness = next(multiplier for multiplier in subgroup
                       if vertices[multiplication[index][multiplier]][0] != 0)
        witnesses.append({'element': index, 'multiplier': witness})
    edge_checks = 0
    for left, right in itertools.product(range(27), repeat=2):
        difference = tuple((first - second) % 3 for first, second in zip(vertices[left], vertices[right]))
        assert (lookup[difference] in subgroup) == (sum(first != second for first, second in zip(vertices[left], vertices[right])) == 1)
        edge_checks += 1
    return {'ring': 'F_3[x]/(x^3-1)', 'elements': 27, 'identity_encoding': identity,
            'element_encoding': 'lexicographic coefficient triples (a_0,a_1,a_2)',
            'scalar_basis_unit_subgroup': subgroup, 'associativity_triples': 27 ** 3,
            'frobenius_witnesses': witnesses, 'ring_graph_hamming_agreements': edge_checks}


def cycle_check(order):
    adjacency = [sorted({(vertex - 1) % order, (vertex + 1) % order}) for vertex in range(order)]
    universe = set(range(order))
    complement = [sorted(universe - set(neighbors) - {vertex})
                  for vertex, neighbors in enumerate(adjacency)]
    assert verify_binary.connected(adjacency)
    assert verify_binary.connected(complement) == (order >= 5)
    matrix = [[int(column in neighbors) for column in range(order)] for neighbors in adjacency]
    rank = verify_binary.exact_rank(matrix)
    assert order - rank == (2 if order % 4 == 0 else 0)
    subgroup = {1, order - 1}
    assert all((left * right) % order in subgroup for left, right in itertools.product(subgroup, repeat=2))
    result = {'vertices': order, 'ring': f'Z/{order}Z', 'unit_subgroup': sorted(subgroup),
              'connected': True, 'complement_connected': order >= 5,
              'adjacency_rank_over_Q': rank, 'nullity': order - rank}
    if order % 4 == 0:
        vectors = [[(1, 0, -1, 0)[index % 4] for index in range(order)],
                   [(0, 1, 0, -1)[index % 4] for index in range(order)]]
        for vector in vectors:
            assert all(sum(vector[neighbor] for neighbor in neighbors) == 0 for neighbors in adjacency)
        result['integer_null_basis'] = vectors
    if order <= 16:
        tested, modules = verify_binary.homogeneous_sets(adjacency)
        result.update(nontrivial_subsets_tested=tested, homogeneous_sets=modules)
        if order >= 5:
            assert not modules
    return result


def partitions(elements):
    if not elements:
        yield []
        return
    first, *remaining = elements
    for partition in partitions(remaining):
        yield [(first,), *partition]
        for position in range(len(partition)):
            yield [*partition[:position], (first, *partition[position]), *partition[position + 1:]]


def duality_check(dimension):
    size = 1 << dimension
    character = [[(-1) ** ((group_element & dual_element).bit_count())
                  for group_element in range(size)] for dual_element in range(size)]
    tested = 0
    theories = 0
    equalities = 0
    for partition in partitions(list(range(size))):
        tested += 1
        profiles = {}
        for group_element in range(size):
            profile = tuple(sum(character[dual_element][group_element] for dual_element in block)
                            for block in partition)
            profiles.setdefault(profile, []).append(group_element)
        classes = list(profiles.values())
        if len(classes) != len(partition) or [0] not in classes:
            continue
        theories += 1
        for block, group_class in itertools.product(partition, classes):
            sums = [sum(character[dual_element][group_element] for group_element in group_class)
                    for dual_element in block]
            assert len(set(sums)) == 1
            equalities += len(sums)
    assert theories
    return {'group': f'C_2^{dimension}', 'dual_partitions_enumerated': tested,
            'supercharacter_theories': theories, 'dual_constancy_equalities': equalities}


def main():
    cases = [(3, 1), (3, 2), (3, 3), (3, 4), (4, 2), (4, 3), (4, 4), (5, 2), (5, 3)]
    result = {'arithmetic': 'exact integers and rational Gaussian elimination',
              'cycle_checks': [cycle_check(order) for order in range(3, 33)],
              'binary': {'cube_checks': [verify_binary.cube_check(dimension) for dimension in range(1, 9)],
                         'cyclic_group_algebra_check': verify_binary.ring_check()},
              'hamming_checks': [hamming_check(alphabet, dimension) for alphabet, dimension in cases],
              'ternary_group_algebra': cyclic_ternary_ring(),
              'supercharacter_duality': [duality_check(dimension) for dimension in range(1, 4)],
              'scope': 'Exact cycle adjacency ranks in orders 3 through 32 and exhaustive homogeneous sets through order 16; finite tensor eigenvector, connectivity and nullity checks in the listed Hamming parameters; binary dimensions 1 through 8, exhaustive binary homogeneous sets through dimension 4, ternary cyclic ring and elementary binary character partitions. The general theorems have written proofs in note.tex.'}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
