# Cycle and Hamming answers to the primeness question

The note answers Nguyen and Tân's Question 4.17 in arXiv:2508.10348v1, §4.3, p. 19, and their Question 2.20 in arXiv:2603.21239v1, §2.2.4, p. 9, negatively. The cycle example is the generalized Paley graph in Mináč–Muller–Nguyen–Tân, arXiv:2212.02005v1, Corollary 2.5 and Theorem 3.5. The note gives a self-contained cycle proof and a prime-power Hamming extension of Chen's binary polynomial construction, using the tensor spectrum described by Sander. Its final corollary gives a short proof of Hendrickson duality, answering the first source's Remark 2.2 affirmatively.

## Statements

- For every integer n ≥ 5, the Frobenius ring Z/nZ and unit subgroup {1,−1} give the connected, anti-connected prime cycle C_n. Its nullity is two exactly when 4 divides n, and zero otherwise. The eight-vertex example has rank six.
- For every finite abelian group G of order d and prime power q, the Frobenius group algebra F_q[G] and subgroup {c e_g : c ≠ 0, g ∈ G} realize H(d,q). Its eigenvalues are d(q−1)−qk with multiplicity binomial(d,k)(q−1)^k. For d ≥ 2, except (d,q)=(2,2), it is prime and anti-connected. It is singular exactly when q divides d. In particular, H(3,3) has rank 15 and nullity 12.
- Equal-sized character and group partitions with superclass-constant character sums have the dual constancy property. This is Hendrickson's prior duality theorem, with a Fourier proof included in the note.

## Build

Use Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, use a TeX installation containing `amsmath`, `amssymb`, `amsthm`, `geometry`, `booktabs`, `array`, `lmodern` and `hyperref`:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is embedded in `note.tex`; no separate bibliography or figures are required.

## Exact finite checks

Python 3.10 or later, standard library only:

```sh
python3 verify.py > verification-new.json
cmp verification.json verification-new.json
```

Keep `verify.py` and `verify_binary.py` together. All arithmetic is integral or rational. The output records each parameter and check:

- Cycle orders 3–32: exact adjacency ranks, nullity, unit subgroup and connectivity; exhaustive homogeneous subsets through order 16; explicit integer null bases when the order is divisible by four.
- Binary dimensions 1–8: all character eigenvector equations and connectivity; exhaustive homogeneous subsets through dimension four; the sixteen-element cyclic group algebra and its Frobenius witnesses.
- Hamming parameters (q,d) = (3,1), (3,2), (3,3), (3,4), (4,2), (4,3), (4,4), (5,2), (5,3): exact integer tensor eigenvector equations, spectrum, multiplicities and nullity. Homogeneous subsets are exhaustive for graphs with at most sixteen vertices; H(3,3) has an independent rational rank calculation.
- The 27-element ring F_3[x]/(x³−1): subgroup, associativity, Frobenius witnesses and graph identification.
- Character partitions of C_2, C_2² and C_2³: all dual partitions are enumerated, the superclass partition is determined by character-sum profiles, and every resulting supercharacter theory satisfies dual constancy.

These finite ranges are separate from the general written proofs.

## Formal statement scope

The [pinned Lean declaration](https://github.com/the-omega-institute/trureturing/blob/2eae85555e169007ff5548c46eb63cd2408f877f/D5/S3/ArithUnits/NguyenTanUnitaryCayleyPrimeSpectrum.lean#L669) `D5.S3.ArithUnits.NguyenTanUnitaryCayleyPrimeSpectrum.result` states the negation of the affirmative form of Question 4.17. The supporting [group-algebra lemma](https://github.com/the-omega-institute/trureturing/blob/2eae85555e169007ff5548c46eb63cd2408f877f/D5/S3/ArithUnits/NguyenTanUnitaryCayleyPrimeSpectrum.lean#L631) covers the binary Frobenius ring, connectedness, anti-connectedness, primeness and a zero eigenvalue for even-order finite abelian groups of order at least four. The [cyclic specialization](https://github.com/the-omega-institute/trureturing/blob/2eae85555e169007ff5548c46eb63cd2408f877f/D5/S3/ArithUnits/NguyenTanUnitaryCayleyPrimeSpectrum.lean#L650) supplies the order-four instance used by `result`. The cycle family, complete multiplicities, prime-power extension, polynomial descriptions and duality corollary have written proofs in the note.

## Contents and license

`note.tex` is the mathematical source; `metadata.json` records the statements, authors and license; `coverage.json` gives the claim/source/statement mapping; `verify.py` and `verify_binary.py` are the finite checkers; `verification.json` is their expected output; `MANIFEST.sha256` lists the source-package file hashes.

Authors: Haobo Ma (ChronoAI Pte Ltd, Singapore), then Wenlin Zhang (National University of Singapore; ChronoAI Pte Ltd, Singapore).

© 2026 Haobo Ma and Wenlin Zhang. Licensed under CC BY 4.0: https://creativecommons.org/licenses/by/4.0/.
