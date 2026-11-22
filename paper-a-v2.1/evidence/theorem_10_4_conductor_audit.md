# Audit of the alleged Theorem 10.4 conductor gap

Date: 22 November 2026.

## Verdict

**VALID WITH AN EXPLICIT COORDINATE CLARIFICATION.** The particular mathematical gap alleged in the pasted referee report is not present under the definitions in the delivered Version 2. Lemma 10.3 promises a common Fourier zero at a nontrivial character, and its proof explicitly rules out both an empty common zero set and a common zero set consisting only of the trivial character. On a product of three prime-order groups, a nontrivial character necessarily has at least one nontrivial local component. The shell sum at that component is −1.

The all-positive evaluations supplied by the referee can vanish. They occur at the trivial character, which is excluded by Lemma 10.3. No stronger zero-existence theorem or added hypothesis is required to exclude them.

This conclusion concerns the exact local manuscript source and the pinned source note read below. It does not assert that a different draft using differently defined subset labels would be unambiguous.

## Sources examined

- output/North_Star_Paper_A_V2/sections/06_primary_rigidity.tex, including the complete proofs of Lemma 10.1, Corollary 10.2, Lemma 10.3 and Theorem 10.4.
- research/snapshot/proofs/LEVEL1_PRIMARY_AXIS_COMMON_ZEROS.md, including both versions of the common-zero argument and the six-shell rigidity proof.

No repository file was changed by this audit.

## 1. The exact character dictionary

Let K = q_1 q_2 q_3, where the q_i are distinct odd primes, and write

  χ = χ_1 ⊗ χ_2 ⊗ χ_3

under the Chinese remainder identification. Let c be the conductor, meaning the order of χ. Then q_i divides c if and only if χ_i is nontrivial.

For the shell of nonzero elements on the i-th prime-order coordinate, the character sum is

  x_i(c) = Σ_{z≠0} χ_i(z)
         = q_i − 1  if χ_i is trivial,
         = −1       if χ_i is nontrivial.

The second equality follows directly: the complete sum of a nontrivial character of a finite cyclic group is zero, and its value at the origin is one. A shell with physical support S contributes the product of x_i(c) over i in S.

Therefore:

  c > 1  ⇔  at least one χ_i is nontrivial
         ⇔  at least one x_i(c) = −1.

In the manuscript's kernel-exponent convention, h_i is 0 or 1 and

  c = ∏ q_i^(1−h_i).

The trivial residual character is h = (1,1,1); a nontrivial residual character has at least one zero kernel exponent. The all-(q_i−1) evaluation is precisely the trivial residual character.

No x_i itself is zero. The "zero coordinate" terminology in the report should be replaced by "zero kernel exponent" or "nontrivial local character" to avoid mixing these two coordinate systems.

## 2. Nonempty common-zero set versus a nontrivial common zero

These assertions are different. Write Z for the set of residual conductor types at which every section transform vanishes.

- Z nonempty permits Z = {1}.
- A common nontrivial zero means Z contains c > 1.

Lemma 10.3 states the second assertion. Its proof does not stop after proving the first.

The proof establishes that, if the first changing primary section is j, the Fourier transform of T_j(g) vanishes modulo p^(j+1) outside Z. Consequently Z empty forces T_j(g) to be zero modulo p^(j+1), while Z = {1} forces it to be constant modulo p^(j+1). The proof excludes both possibilities separately.

A subset label C≠∅ is not sufficient without its definition. If C indexes nontrivial local character components, C≠∅ is equivalent to c>1. If C indexes trivial local components, it is not. The manuscript should make this dictionary explicit instead of relying on an informal subset description.

## 3. Adversarial reconstruction of Lemma 10.3

Let L be the least disputed dyadic size. At a dyadic character of kernel order L, the normalized partial difference is exactly the signed indicator d_t: larger fibers contribute zero, smaller disputes do not exist, and common fibers cancel. The normalized sums g_t are dyadic rationals.

Let j be the first change in d_t. The primary-axis row

  T_k(u) = Σ_{i=0}^k p^i (u_i−u_{i+1})

vanishes on d below j. At every residual character outside Z, some row T_k(d) is nonzero with k≥j, so homometry forces T_k(g)=0 at that character. Since T_k(g)≡T_j(g) modulo p^(j+1), residual Fourier inversion yields the modular support statement above. Inversion is legitimate because 2K is coprime to p; radial inversion coefficients are integer Ramanujan sums divided by K.

At a disappearing size-L dispute, T_j(g)=p^j. At an appearing dispute, the earlier partial sums are zeros followed by twos, and the possible values are

  p^j, −p^j, p^j−2p^b   (0≤b<j).

All are nonzero modulo p^(j+1), so Z cannot be empty.

Assume Z={1}. Then every d_t has total population zero. If d_0 is nonzero, its two signs force both fibers to be empty at the top residual shell in every section; there T_j(g)=0, contradicting constancy and the changing-shell residue. Hence d_0=...=d_j=0, and d_{j+1} has both signs.

At the residual origin, both fibers are equal and at least L through section j+1. Their normalized partial sums are zeros followed by twos. Thus T_j(g) is zero or −2p^c with 0≤c≤j.

Comparing these values with the changing-shell residues gives a contradiction except possibly when p=3, c=j, and the changing-shell residue is +p^j. In that exception the origin fibers exceed L through j and first reach L at j+1, while every appearing dispute had both earlier fibers equal to L.

A larger dispute cannot occur at or after j+1 because the origin bounds every fiber by L. If a larger dispute occurs earlier at a residual shell, all partial sums there through j are zero: before it the occupied side exceeds L; after it the empty side remains empty; size L is ruled out by d_0=...=d_j=0; smaller disputes are ruled out by minimality. At j+1 its sum is zero or one. Therefore T_j(g) is zero or −p^j, neither congruent to +p^j modulo 3^(j+1). Hence the exceptional case has only one disputed size.

In the one-size case, restrict to the subgroup with dyadic order L. The two sets agree outside it, and shifts inside preserve exterior cosets, so the restricted pair remains homometric. Its d_t are unchanged. Apply the same modular Fourier argument to its recomputed sums. A changing shell has T_j(g)=p^j, whereas either the residual origin or residual top has T_j(g)=0, using the two signs supplied by pure-primary marginal balance. This excludes Z={1}.

Thus Lemma 10.3 does establish a genuinely nontrivial residual character. The p=3 leg is needed, and it is present.

An optional exposition improvement is to say explicitly in the one-size paragraph that the preceding modular Fourier argument is reapplied to the restricted pair; this follows from its proved homometry and unchanged d_t.

## 4. The two alleged cancellation obstructions

### Rank one: residual primes (3,5,7)

The signed polynomial is

  −x_1−x_2+x_3.

At the trivial character x=(2,4,6), this is −2−4+6=0. Its values at the seven nontrivial conductors are:

| Conductor | Fourier value |
|---|---:|
| 3 | 3 |
| 5 | 5 |
| 7 | −7 |
| 15 | 8 |
| 21 | −4 |
| 35 | −2 |
| 105 | 1 |

Thus its only zero conductor is 1.

### Rank two: residual primes (5,7,13)

The signed shell polynomial is

  x_2 x_3 − x_1 x_3 − x_1 x_2.

After division by x_1 x_2 x_3 it is 1/x_1−1/x_2−1/x_3. At the trivial character x=(4,6,12), the stated reciprocal cancellation holds. Its exact shell-polynomial values are:

| Conductor | Fourier value |
|---|---:|
| 5 | 90 |
| 7 | −56 |
| 13 | −26 |
| 35 | −1 |
| 65 | −1 |
| 91 | 9 |
| 455 | −1 |

Again its only zero conductor is 1.

The supplied independent script checks both tables using two routes: the local shell-product formula and literal signed residue polynomials reduced modulo cyclotomic polynomials. It imports no construction or existing verification module. Both tables pass in exact arithmetic. This is a bounded control of the two reported examples, not a proof of the infinite rigidity family.

## 5. Remaining Theorem 10.4 chain

The rest of the printed proof was checked independently:

1. Pure-axis marginals and common-fiber equality imply that each physical shell has one fixed discrepancy sign along the primary axis.
2. The three residual marginals imply equality of the aggregate signed masses of each singleton and its complementary doubleton.
3. Therefore a disputed primary section cannot contain both residual ranks; otherwise all its discrepancies would share one sign, contradicting its zero total mass.
4. The least-size sections cannot all coincide, because each active singleton/doubleton has an active complement and such a complement would create a mixed-rank section.
5. At the common nontrivial character no physical residual shell has zero character sum, so every nonzero least-size section contains at least two shells.
6. If all disputes in such a section had the least size, all three shells would occur with both signs. The printed rank-one parity argument and rank-two reciprocal estimates exclude every nontrivial character. The all-positive character is not in this case.
7. Hence there are exactly two least-size shells of opposite sign and one larger disputed shell.
8. The three global complementary-pair signs have one repeated sign. The larger disputed shell always has that repeated sign. Section mass balance forces every nonzero least-size population to have the opposite sign.
9. Their positive-weighted sum cannot vanish, contradicting the exact-size balance proved in Lemma 10.1.

No use of the §11 enumeration, five-factor computation, or a bounded prime search occurs in this chain. No genuine analytic gap was found in this recheck.

## Files supplied for the revision

- primary_zero_clarification.tex: insertion-ready replacement at the exact invocation point in Theorem 10.4; it preserves theorem numbering.
- primary_zero_referee_response.tex: concise response explaining the conventions and both examples.
- check_primary_conductor_coordinates.py: independent exact controls for the two cited examples.
- primary_conductor_coordinates.json: completed PASS output.

## Limits

This is a same-session mathematical audit and exact finite regression, not an external referee report or formal proof-assistant verification. The parent should reconcile the actual full report if available, since this audit addresses the claim in the pasted summary against the delivered Version 2. The correct revision does not require silently strengthening assumptions or relabeling a numerical search as a proof.
