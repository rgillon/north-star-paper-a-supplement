# Independent audit of the Fourier conductor classification

## Verdict

**VALID in the stated class.** The original Paper A Fourier-conductor theorem gives a complete, constructive classification of homometric mates among subgroup unions for every positive modulus and every specified subgroup union. The proof does not require a conjectural universal anchoring statement, a bound on the number of primes, a bound on dyadic height, or an unstated binary-feasibility test.

This audit is an independent derivation within the current preparation session, not an external referee report and not a Lean formalization. It supports inclusion of the theorem with its complete proof. It does not certify the separately named generator-lattice compression or factor-cut moments theorem, whose source proofs are being recovered separately.

## Sources examined

- Original source `02-sign-graph.tex` in `/workspace/scratch/95a17a5755a9/manuscript_source_original/North_Star_Collection/A-Subgroup-Unions/`.
- Adjacent `03-boundaries.tex` and `04-small-examples.tex`.
- The public revised manuscript extracted as `/workspace/scratch/060e9f471065/research/original_paper.txt`.
- Root reports that the current manuscript-branch graph source at commit `79350cb7bb6032df7637dcfd3e2bd029e6362a85` is the same graph-source blob.

A replacement, self-contained LaTeX section is in `fourier_graph_complete.tex`. It deliberately calls the new linear-algebra restatement the **conductor-constraint rank** criterion, avoiding conflation with the unavailable generator-lattice theorem.

## Hypotheses and exact conclusion

Fix a positive integer N, with supplied prime factorization, and a specified union A of additive subgroups of Z/NZ. Equality means equality of complete cyclic autocorrelation at every labeled shift. No quotient by translations, automorphisms, or arbitrary permutations of shifts is taken.

An integer radial function is constant on capped-valuation shells. A subgroup union is exactly an increasing binary valuation profile. In physical element-order coordinates the same condition is downward closure. These two opposite-looking monotonicity conventions must not be mixed.

The graph uses only the nonzero conductor Fourier coefficients of A. It joins adjacent conductor vectors s and s+e_i exactly when p_i^(a_i−s_i) does not divide 2F(s). Root conductor 1 and all pure prime-power conductors are anchors, provided their coefficients are nonzero. If u is the number of graph components disjoint from those anchors, there are exactly 2^u subgroup-union mates, including A. Every mate is explicitly reconstructed by negating the chosen component coefficients and applying conductor inversion.

The empty set is handled separately. N=1 and Fourier zeros require no exception to the mathematical argument once the graph is defined correctly.

## Implication audit

### 1. Integer inversion and adjacent congruences

Every integer radial function uniquely expands as a sum of subgroup indicators with integer coefficients alpha_r by Möbius inversion on a product of chains. Its conductor transform is F(s)=sum_{r>=s} W_r alpha_r, with W_r the subgroup order. Thus every adjacent difference is divisible by the appropriate primary factor.

The converse is essential and valid. The mixed finite-difference numerator defining alpha_r is divisible by every p_i^(a_i−r_i), because it can be grouped into adjacent differences along coordinate i. These primary factors are pairwise coprime, so their product divides that numerator. Therefore all coefficients alpha_r are integers. This proves there are no omitted mixed-coordinate integrality conditions.

### 2. Exact graph necessity

A radial binary mate has a real integral conductor transform with the same squared coefficients, so it differs only by signs at nonzero coefficients. Cardinality fixes the root sign. Opposite signs across an adjacent pair are compatible with both congruences if and only if the corresponding primary factor divides 2F(s). Therefore signs must be constant on graph components. This direction does not assume that a nullspace vector is binary.

### 3. Exact graph sufficiency

Flipping any non-root components preserves all adjacent congruences, including pairs with one or two zero entries. Hence the inverse is an integer radial function. Root preservation gives equal total mass; Parseval, with multiplicity phi(c) for conductor c, gives equal squared norm. The identity sum_x b(x)(b(x)−1)=0 then forces every integer value of b to be 0 or 1. This step is valid even if negative values were initially possible. Fourier squares give all labeled autocorrelations, not merely marginal data.

### 4. Prime-power positivity

Fix all coordinates except one primary coordinate. The remaining fiber is empty or a subgroup because the primary subgroup lattice is a chain. Its transform is nonnegative at every character. Summing fibers proves nonnegativity at every pure prime-power conductor. Thus a subgroup-union mate must retain the signs of all nonzero pure-prime coefficients.

### 5. Preserving anchors is sufficient for divisor closure

This is the most substantive issue, and the original proof closes it correctly. Keeping pure-prime coefficients preserves each one-coordinate marginal by Fourier inversion. For a radial binary b and a shift supported on one primary coordinate with valuation s, exact shell counting gives

D_{p,s}(b) = sum_{r<s} w_r M_r(b) + p^(a−s) M_s(b) − C_b(t)
           = 2 sum_y sum_{r>s} w_r b_s(y)(1−b_r(y)).

All terms on the right are nonnegative; all w_r are positive. For an increasing indicator A the defect is zero. The fixed marginals and autocorrelation make it zero for b too, forcing b_s(y)<=b_r(y) whenever s<r. Applying this to each primary coordinate yields an increasing profile and hence a subgroup union. The same-shell term before rearrangement is (p−2)p^(a−s−1); for p=2 it is zero, not negative, and the argument remains valid.

### 6. Exhaustion and exact multiplicity

Each unanchored component may be flipped independently. Different selections change at least one nonzero Fourier coefficient, so their inverses are distinct. Conversely every subgroup-union mate yields precisely one such selection. The number is therefore exactly 2^u, not merely a lower bound.

### 7. A rigorous rank restatement

Let R_A have a column for each nonzero conductor, a row e_s−e_t for each graph edge, and a row e_s for each anchor. Its kernel consists of functions constant on each graph component and zero on each anchored component. Thus nullity(R_A)=u, and the binary part of its kernel is exactly a u-dimensional cube. The sign change epsilon_s=1−2z_s converts each binary kernel point into its mate. This completely justifies 2^nullity for this explicitly defined matrix. No general principle that arbitrary rational nullspaces produce binary mates is used.

## Other statements checked

- Odd-order rigidity for symmetric binary sets: valid by reduction of convolution squares modulo 2 and invertibility of doubling.
- Twice-odd subgroup-union rigidity: the minimal disputed physical-order character has odd transform on the difference; all other odd-weight shells divide that order. The second convolution factor has the same odd parity because the singleton profile vanishes on every divisor of the chosen order. This is a valid contradiction.
- Height-two family with at most one odd prime congruent to 3 modulo 4: the mod-four argument is coherent. At a disagreement, singleton and size-two profiles vanish on all larger valuation shells; any contributing lowered coordinate supplies a factor p−1. If only one prime can supply valuation exactly one modulo 4, minimizing the disputed valuation in that prime rules out the final possible contribution. This is a restricted family theorem, not universal 4M rigidity.
- Distinct subgroup unions cannot be translates: after quotienting by the common period subgroup, each primary coordinate has at least one singleton fiber, making every pure faithful Fourier coefficient strictly positive. A translation taking one subgroup union to another must therefore be trivial in each primary coordinate. The stronger single-subgroup-reflection argument in the source is also valid: an odd-order reflection factor cannot have phase −1, and a dyadic factor of size at least four contradicts the values on a faithful character and its square.
- The N=12 and N=24 examples are radial ambiguities, not nonrigid subgroup-union examples. Their putative alternative subgroup unions fail positivity. The manuscript must retain that scope prominently.

## Independent computational check

Script: `verify_sign_graph_independent.py` (Python standard library only).

Command:

```sh
python3 verify_sign_graph_independent.py 60 > sign_graph_independent_results.json
```

It enumerates every binary physical-order-shell profile, computes the Fourier transform using the integer Ramanujan identity c_d(k)=sum_{j|gcd(d,k)} j mu(d/j), groups by the complete squared spectrum, and compares the exact family sizes with the independently built graph. It also checks that the actual sign masks of all mates are exactly the allowed component masks. Subgroup-union profiles are recognized independently by closure of element orders under divisibility. For every nonempty radial profile at N<=30, the script additionally computes all labeled autocorrelation shifts directly from residues and checks equality of the Fourier and time-domain equivalence classes.

Result: **PASS** for all 60 moduli from 1 through 60; 8,182 radial binary profiles; 465 subgroup unions including empties; 8,122 nonempty radial count checks; 405 nonempty subgroup-union count checks; 980 complete direct time-domain profiles. No contradiction occurred.

These are regression checks of formulas and implementation. They are not the proof of the universal theorem and contain no genuine subgroup-union nonrigid example at these small orders.

## Publication recommendations

1. Retain the full integrality lemma, mass/norm argument, and defect identity. A summary of sign propagation alone would omit the heart of the constructive theorem.
2. State factorization and explicit valuation-profile input when quoting O(omega(N)tau(N)) integer operations. Do not call this polynomial time in the succinct generator-list input, and distinguish arithmetic-operation complexity from bit complexity.
3. If the generator-lattice source remains unavailable, present the audited graph theorem and the conductor-rank corollary under their actual names. Identify the stronger compression as an unverified import, not a failed theorem.
4. State that modulus-level rigidity is the universal quantification of this set-level condition over all subgroup unions. Enumeration of all divisor ideals gives a finite complete algorithm, but it is not a concise arithmetic classification and should not be confused with the separately claimed factor-cut procedure.
5. Fix unresolved bibliographic citations in the older PDF; `[?]` reference fields were visible in the extracted introduction.

## Unresolved proof gaps in the audited theorem

None found. Remaining verification obligations concern the independently unavailable generator-lattice compression, factor-cut moment algorithm, selected deeper arithmetic results, and construction certificates; this audit makes no claims about those components.
