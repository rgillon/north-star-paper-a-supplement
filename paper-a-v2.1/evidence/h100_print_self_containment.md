# H100 construction self-containment and finite-data closure

## Verdict

The complete H100 construction text, together with printed donor arrays, mask arrays, damping list, and the full ordinary-slot ordering, determines every mathematical object by finite exact formulas. No additional research conversation or repository theorem is needed. The remaining verifier calls check finite numerical facts; they are not hidden constructive or mathematical existence arguments.

The printed manuscript should make the following choices explicit. These are clarifications of the existing construction, not changes in hypotheses or its data.

## Required clarifications

### Hidden-state and conductor-mask encodings

For the hidden group of order 5^2 11^2 13*17, a tuple (v1,v2,v3,v4) has exponent ranges (0..2,0..2,0..1,0..1). The corresponding bit position in a FIBERS mask is

v1 + 3 v2 + 9 v3 + 18 v4.

Thus the first coordinate changes fastest. A one bit includes the entire element-order shell. Upper conductor-mask bit i corresponds to the ith prime in the printed upper FIBERS table. Lower conductor-mask bit i corresponds to the ith prime in the printed signed lower-prime list. The ordinary tag ordering instead uses increasing lexicographic tuples with the last coordinate fastest; the different conventions must be stated.

### Lower packet coefficients and uniqueness

The lower 14-by-13 incidence matrix has entries 1_{i in T}. Its first 13 rows, in the order of the printed lower-prime list, have determinant 317. Solve these rows with right-hand side (1−p_i)/p_i. This uniquely defines the 13 coefficients.

The final equation is automatic. Every mask column is orthogonal to the vector (s_i p_i), since its signed prime sum vanishes. The right-hand side is also orthogonal to this vector, since its inner product is sum_i s_i(1−p_i)=0 by the donor balance. The residual of the solved system is supported on its final row, whose multiplier s_14 p_14 is nonzero, so this residual is zero.

The same argument explains the upper system: its column relation is (s_i b_i p_i), its right-hand side is −1/p_i, and their inner product is −sum_i s_i b_i=0. The first15-row minor has determinant193. Thus the full upper system is consistent and uniquely solved without an unstated solver convention.

### Choices in scalar recipes

- Take b0=(-2)^{-1} mod o to be the least nonnegative residue, as the source implementation does.
- Start the three partial upper corrections with B0=max|u0| on the 2^16 activity states. After each correction define the next B to be exactly the right side of the displayed max bound, using the preceding B and the specified gamma. This fixes gamma recursively, not merely relative to some unspecified certified envelope.
- Use nearest rounding with ties to even in the original H100 slope definition, matching the implementation.
- Define denominator as the positive denominator in lowest terms, oddpart(n)=n/2^{v2(n)}, and bitlength(n)=1+floor(log2 n) for positive integers.
- For the ordinary baseline, write J_k=z_{2k}+z_{2k+1} and
  C_l=2 sum_{k<l}J_k+sum_{k=l}^{R+159}2^{l−k}J_k.

## Verification data versus proof dependencies

The following are finite assertions whose exact checking remains part of the computer-assisted theorem:

1. Hidden upsets, their totient masses, the signed mask sums, and the two incidence determinants.
2. Donor packet norms evaluated over finite activity cubes.
3. The eight modular-helper primality certificates and the smaller trial-division claims.
4. Tag-background closure and every ordinary/fine cover relation.
5. Intercept interval nonemptiness, endpoint inequalities, tail bound, height/precision bounds, and the explicit modulus prime map.

Their invariants and exhaustive finite ranges are all stated mathematically. The proof of annihilation, forced boundary values, nesting for every digit choice, grid coverage for every target, and labeled-shift homometry is supplied in prose and equations, not inferred from sampled evaluations.

The bundle includes the actual arrays and exact programs. A referee can therefore reproduce the finite checks without GitHub. If PDF-only arithmetic replay is desired, append the Pocklington factor/witness table available in `h100_independent.json` (the field `pocklington`) so that primality certificates are also printed. This is not needed to make the manuscript-plus-supplement complete, provided the electronic certificate is explicitly included and identified.

## No additional missing oracle

The construction uses no optimization oracle after the data are fixed. The stored ordinary permutation is finite input, the donor coefficients are uniquely determined rational solutions of printed matrices, modular inverses have canonical residues, all floors and rounding rules are exact, and the greedy and reverse parity decoding procedures are deterministic. Historical searches that discovered those arrays are irrelevant to the validity and reproducibility of the fixed construction.
