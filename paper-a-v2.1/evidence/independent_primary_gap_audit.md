# Independent audit of the alleged Theorem 10.4 gap

Date: 22 November 2026.

## Verdict

**The specific reported gap is not present in the printed proof under its stated Fourier convention.** A nontrivial residual character necessarily gives at least one singleton-shell sum equal to −1. The all-positive evaluations quoted in the objection occur at the trivial character only.

A short explicit character-to-kernel dictionary should be inserted, because the section uses physical order exponents, character conductor exponents and kernel exponents. Their differing meanings make a mistaken identification possible. This is an exposition clarification, not a strengthening of the theorem's assumptions or an additional mathematical hypothesis.

This audit independently read the conductor definition in sections/02_fourier_graph.tex, the complete Lemma 10.3 and Theorem 10.4 in sections/06_primary_rigidity.tex, and the recovered LEVEL1_PRIMARY_AXIS_COMMON_ZEROS.md. The paths are within output/North_Star_Paper_A_V2 and research/snapshot/proofs respectively.

It addresses the referee finding pasted into the conversation. It does not claim to have inspected an independently attached full referee report.

## 1. The exact character identity

Let K=q₁q₂q₃, where the qᵢ are distinct odd primes. Write a character on its Chinese-remainder product as χ=χ₁⊗χ₂⊗χ₃ and set T={i:χᵢ is nontrivial}.

The character conductor is c(χ)=∏_{i∈T}qᵢ. Thus χ is nontrivial if and only if c(χ)>1, if and only if T is nonempty.

For the physical singleton shell on coordinate i, xᵢ=∑_{y≠0}overline(χᵢ(y)). If χᵢ is trivial, xᵢ=qᵢ−1=uᵢ. If it is nontrivial, the sum over the whole prime group is zero, and deleting the origin gives xᵢ=−1. Consequently χ≠1 implies at least one xᵢ=−1.

Conversely, xᵢ=uᵢ for all i implies χᵢ=1 for all i, hence χ=1. There is no nontrivial residual character at the all-u evaluation.

The prime-power rows in §10 are indexed by kernel orders. On a residual prime axis, kernel exponent kᵢ=1 means the local character is trivial; kᵢ=0 means it is nontrivial. Therefore χ=1 corresponds to (k₁,k₂,k₃)=(1,1,1), and χ≠1 means at least one kᵢ=0.

If a set C were instead defined as the support of the kernel vector, nontriviality would mean C≠{1,2,3}; the condition C≠∅ would be the wrong translation. The printed lemma does not introduce that alternative C notation: it asserts a **nontrivial conductor** and explicitly excludes the trivial-character-only case in its proof.

None of the xᵢ is zero. A “zero coordinate” here can only mean a zero kernel exponent, not a zero singleton-shell sum.

## 2. Does Lemma 10.3 actually provide a nontrivial character?

Yes. Its conclusion is stronger than mere existence of a common zero, and its proof separately excludes a common-zero set confined to the trivial character.

Let j be the first change in the signed least-size sections dₜ. At any residual conductor outside their common zero set Z, some p-axis Fourier row Tₖ(d) is nonzero. Rows with k<j vanish because d₀=⋯=dⱼ. Therefore the nonzero row has k≥j, homometry forces Tₖ(g)=0, and Tₖ(g)−Tⱼ(g) is divisible by p^(j+1) in the dyadic rational coefficient ring. The residual Fourier transform of Tⱼ(g) vanishes modulo p^(j+1) outside Z.

Inversion on K is legitimate because p∤2K. If Z is empty, Tⱼ(g) is zero modulo p^(j+1). If Z contains only the trivial character, Tⱼ(g) is constant modulo p^(j+1). This constant assertion refers to the actual trivial character, not a nonzero or nonempty kernel label.

At a disappearing or newly appearing least-size dispute, the possible values are pʲ, −pʲ, or pʲ−2pᵇ (0≤b<j). Each is nonzero modulo p^(j+1), which first excludes Z=∅.

Under the remaining hypothesis Z={trivial}, each dₜ has zero population. If d₀≠0, its two signs force both fibers empty at the maximal residual shell, giving Tⱼ(g)=0 there and contradicting constancy. Therefore d₀=⋯=dⱼ=0 and dⱼ₊₁ has both signs. At the residual origin the two fibers are nonempty, equal and at least L through j+1. Their partial sums are zeros followed by twos; Tⱼ(g) is either zero or −2pᶜ with 0≤c≤j.

Comparison with the changing-shell values excludes all cases except p=3 and Tⱼ(g)≡pʲ≡−2pʲ. In that exception, the origin fibers first reach L at state j+1, and every newly disputed shell had both fibers equal to L at all earlier states. A larger earlier dispute would have zero partial sums through j, and sum zero or one at j+1, giving Tⱼ(g)=0 or −pʲ there. Neither is congruent to pʲ modulo 3^(j+1). Thus no larger dispute exists.

The printed one-size argument then applies. Restricting the dyadic coordinate to its subgroup of order L preserves homometry because the two original sets agree on every exterior dyadic coset and shifts inside the subgroup preserve those cosets. Recompute g for this restricted pair at its trivial dyadic character; d remains unchanged. A changing shell gives Tⱼ(g)=pʲ, while either the residual origin or maximal shell gives zero. Constancy is again impossible.

This proves that Z contains an actual nontrivial character. The argument does not stop after showing that Z is merely nonempty.

## 3. The two advertised cancellations

They are genuine arithmetic identities, but each is a trivial-character population cancellation. An independent exact script tests the entire eight-conductor spectrum in each example. It constructs the physical radial function on literal residues and reduces its character polynomial modulo the appropriate cyclotomic polynomial; no North Star verifier is imported.

### Rank one, primes (3,5,7)

Take signs (−1,−1,+1) on physical orders 3,5,7. The polynomial is −x₁−x₂+x₃.

| Character conductor | Fourier value |
|---|---:|
| 1 | 0 |
| 3 | 3 |
| 5 | 5 |
| 7 | −7 |
| 15 | 8 |
| 21 | −4 |
| 35 | −2 |
| 105 | 1 |

The identity −2−4+6=0 occurs only in the first row.

### Rank two, primes (5,7,13)

Take the polynomial x₂x₃−x₁x₃−x₁x₂, whose quotient by x₁x₂x₃ is 1/x₁−1/x₂−1/x₃.

| Character conductor | Fourier value |
|---|---:|
| 1 | 0 |
| 5 | 90 |
| 7 | −56 |
| 13 | −26 |
| 35 | −1 |
| 65 | −1 |
| 91 | 9 |
| 455 | −1 |

The identity 1/4−1/6−1/12=0 occurs only at conductor 1. Neither example can supply the common **nontrivial** zero required by Lemma 10.3.

Program: revision_audit/independent_primary_conductor_check.py.
Completed output: revision_audit/independent_primary_conductor_check.json.
Status: PASS, all 16 conductor-class evaluations checked exactly against literal cyclic residue polynomials.

## 4. Rechecking the subsequent polynomial exclusion

At a nontrivial character, let m be the number of coordinates equal to −1.

For a rank-one polynomial with three coefficients ±1, parity forces m even if the polynomial vanishes. Since m≥1 and m≤3, m=2. If the two coefficients at −1 differ, they cancel and the remaining term is nonzero. If they agree, the third coefficient must have the opposite sign because the section contains both signs; all three evaluated terms then have the same sign. Thus no zero exists.

For a rank-two polynomial, divide by the nonzero product x₁x₂x₃. With m=1, the magnitude-one term cannot cancel the other two, whose reciprocal magnitudes sum to at most 1/2+1/4<1. With m=2, their sum is 0 or has magnitude 2 and cannot cancel a remaining reciprocal of magnitude at most 1/2. With m=3, an odd number of ±1 terms cannot sum to zero.

Thus the proof correctly handles all nontrivial conductors. The omitted m=0 case is exactly the excluded trivial character.

I also reconstructed the final complementary-pair sign argument: every section containing the least disputed size has precisely two least-size shells of opposite signs and one larger shell carrying the globally repeated sign. Section population balance therefore gives the same nonzero sign to every exact-least-size section population, contradicting the global exact-size balance. No further gap was found in that inference.

## 5. Recommended manuscript response

Add the supplied explicit χ/T/kernel dictionary before the first polynomial exclusion. State that the two reported identities occur only at conductor 1, and retain the all-prime theorem if the independent parallel audit agrees.

Do not downgrade a theorem solely because the displayed population identities vanish: those identities concern the very trivial-character case excluded in the preceding lemma. Do not claim these sixteen finite evaluations prove the infinite theorem; they are examples illustrating the general character argument.

The insertion changes no hypotheses and relies on no numerical search, prime restriction, or extra conjecture. No claim of Lean formalization or external journal refereeing is made.
