# Cross-audit of terminal incidence and the modulus moments

Date: 22 November 2026. Same-session collaborator `/root/universal_sources`, shared context. This is an adversarial reconstruction of the displayed arguments, not external peer review or a formal proof-assistant certificate.

## Files inspected

- `audits/terminal_certificate_section.tex`
- `output/North_Star_Paper_A_V2/sections/01_introduction.tex`
- `output/North_Star_Paper_A_V2/sections/03b_modulus_moments.tex`
- `output/North_Star_Paper_A_V2/sections/06_primary_rigidity.tex` (the common-fiber lemma and immediate consequences used by terminal exclusion)
- `output/North_Star_Paper_A_V2/sections/10_consequences.tex`

## Verdict

No substantive mathematical error was identified in the terminal incidence equivalence, its implication for the least disputed class, or the moment theorem's physical-order/conductor/CRT conventions. No new enumeration was run. This note does not independently certify the exhaustive 24-grid implementation or every arithmetic construction referred to by the introduction.

## Terminal incidence reconstruction

For a p^a axis, an inverse Fourier column at kernel exponent h is the Ramanujan character sum at character order p^(a-h). On a point of physical order p^v it equals phi(p^(a-h)) for v<=h, -p^(a-h-1) for v=h+1, and zero above that; for h=a it is identically one inside the grid. Extending values by zero past the top gives, including both boundary cases,

`I_h(v)-I_h(v+1) = p^(a-v) (1_{h=v}-1_{h=v-1})`.

Taking one forward difference in each coordinate gives

`b_v = (K/product p_i^v_i) sum_h M_vh t_h`,

where `t_h = Fourier(g)(h)/K`. Thus the displayed left side `b_v product p_i^v_i/K` has exactly the correct normalization. No factor of K or reciprocal prime power is missing.

Necessity uses only that the prescribed region U is upward closed: all forward neighbors in the difference at v are also prescribed. Conversely a solution defines g with Fourier support in Z. The mixed differences of g-|d| vanish on U. The finite product-of-chains identity recovers a value at v by summing those mixed differences over all points above v; all lie in U. Hence g=|d| on U. This proves sufficiency, rather than merely a necessary row condition.

If Z is empty, the same argument applies with an empty-column incidence matrix. The zero-extension convention supplies the top boundary and avoids any unpinned additive-constant ambiguity.

## Least disputed class reconstruction

The common-fiber lemma states equality of the numbers of positive and negative disputes at each exact dyadic size. For the least size L, its nonzero supports P,Q are therefore balanced and both present. They are mutually incomparable. Each same-sign support is order-convex because on a chain between two size-L disputes the empty side stays empty and the occupied nonincreasing size remains L.

Let D be the downward closure of P union Q. Both D\Up(Q) and D\Up(P) are downsets. If a point in D is above P, it is below a point of P or Q; the latter is impossible by cross-sign incomparability, and the former puts it back in P by order-convexity. This proves their difference is exactly 1_P-1_Q=d.

At the dyadic character with kernel order L, a subgroup fiber contributes its full size if its size is at most L and contributes zero otherwise. Common fibers cancel in the difference, and smaller nonzero disputes do not exist. The normalized difference is therefore d. On the upward closure of its support, one original fiber remains empty and the other is empty or L; a positive smaller fiber would contradict leastness. Thus the normalized sum equals |d| there. Homometry gives disjoint Fourier supports for that sum and d at every odd conductor. This is precisely a terminal interpolant, establishing the exclusion-to-rigidity implication without requiring terminal feasibility to imply a realizable pair.

## Symbolic exclusion logic

The described tests are necessary tests and therefore may safely retain supersets of possible Fourier-zero conductors. Proving infeasibility with such a superset also proves infeasibility for the actual zero set. Nonzero weak-sign coefficient polynomials on nonnegative gaps and odd-constant polynomials on even gaps are nonvanishing for the stated reasons. The cube functional is supported above its chosen disputed shell because the normalized difference of kernel rows v and v-1 has only the physical states v and v+1. A unique least prime exponent in an integer left-null identity gives the stated modular contradiction when its nonzero integer coefficient is smaller in absolute value than that prime. The final displayed polynomial reduces modulo s to q^2 r(p-1)(p+1), which is nonzero for distinct odd primes p<q<r<s. The reported category totals add to the root-survivor totals in all six printed rows.

These checks validate the mathematical implications of the described rejection rules. They do not establish that the implementation enumerated every required pattern or implemented each rejection rule without error; that remains the explicit COMPUTED dependency and the separately logged replay.

## Moment normalization and CRT reconstruction

A character of conductor c may be represented by frequency k=N/c. Summing that character over all physical elements of order d gives c_d(k), so `F_I(c)=sum_{d in I}c_d(N/c)` is correct. A different character of the same conductor only changes the frequency by a unit, which leaves these Ramanujan sums unchanged.

For N=p^a M, c|M and m|M,

`c_{p^t m}(N/(p^j c)) = c_{p^t}(p^(a-j)) c_m(M/c)`.

This follows from coprime multiplicativity and unit invariance: M/c is a unit on the p-axis and p^(a-j) is a unit modulo m. Therefore the transfer recurrence's exponent contributions use the correct conductor order j and physical order t; no CRT inverse-unit adjustment is missing.

The inverse formula `1_A(N/d)=N^(-1) sum_{c|N}F_A(c)c_c(N/d)` groups the ordinary inverse transform by conductor. The order-one case N/d=N denotes zero modulo N and is valid. Squaring only after adding all section vectors preserves the mixed cross terms. The stated coefficient moments count actual ideal pairs and the bit-complexity claim is explicitly exponential in tau(N), with the factorization supplied.

## Introduction and consequences

The inspected text distinguishes set-level classification from modulus-level decision, and exhaustive decision from a short arithmetic characterization. It does not identify the independent moment construction with the unrecovered source note. Closure under multiples follows from subgroup-union preimages and the exact autocorrelation scaling; h(M) monotonicity and existence of divisibility-minimal nonrigid divisors follow. The restriction corollary correctly expresses each minimal missing-order condition as a union of subgroups with pairwise coprime prime-power indices and chooses the first nonrigid intersection. The numerical height bounds and the existence of the 99-prime construction rely on their separate construction audits and were not re-established here.
