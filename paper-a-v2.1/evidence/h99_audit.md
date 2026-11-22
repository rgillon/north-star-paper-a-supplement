# Audit of the H100 rank refinement and the 99-odd-prime construction

## Verdict and scope

**VALID, as computer-assisted theorems with their stated finite data.** The rank-dependent decoder rigorously improves the H100 dyadic height from 3446 to 3360 at the identical odd part. The 19-damping-block construction with an additional 65537 spectator guard rigorously gives a nonrigid modulus

N99 = 2^3408 * 65537 * M100 / (2101327^2 * 4415573059603),

with exactly 99 distinct odd prime divisors and 8,017 bits. No claim of minimality of this modulus, odd-prime count, or dyadic exponent follows.

This assessment depends on the complete H100 donor/interpolation proof already audited by the parent and the exact H100 donor arithmetic already independently recomputed. The current audit independently checks the changed mathematics and changed finite arithmetic. It is not an external referee report or a Lean proof.

## Sources

- `proofs/LEVEL1_99_ODD_PRIME_PAIR.md`.
- `proofs/LEVEL1_H100_RANK_PRECISION.md`.
- `experiments/ppi_99_prime_pair.py`, `ppi_99_slot_order.json`, and `ppi_100_rank_precision.py`.
- The complete self-contained H100 manuscript source `manuscript/paper-a/05-construction.tex`.
- The original H100 certificate data and previously independently verified donor report `audits/h100_independent.json`.

## Main mathematical issue: why truncation is legitimate

Let the rank level be l and offset b=R−l. The profile proof bounds the denominator exponent of every requested mass divided by the fine unit M0 by an absolute A, before division by the level scale lambda=2^(−b). Multiplication by 2^b removes b binary denominator places, so the required grid is 2^(−(A−b)) Z, as long as A−b is nonnegative. Thus only D_l=A−b fine positions are needed.

The baseline must be truncated too: use T_l=sum_{j<=D_l}u_j2^(−j), not the old fixed A-term baseline. Both source refinements do this. The endpoint interval checks use this actual truncated baseline.

At a retained position j the corrected subgroup depth is b+j+epsilon_j, with 0<=epsilon_j<=7. Since j<=A−b, the depth never exceeds A+7. This gives 3360 for H100 and 3408 for H99 without discarding precision required by any actual free-order value.

## Nesting audit

Within a level, a prefix of the fine ordering is a downset because all strict divisors occur earlier. Control indices divide along divisibility, and dyadic corrections become smaller toward multiples, as in H100.

Across an increase of l by one, the retained prefix grows by one. Every old fine slot's new minimum equals its previous maximum; the newly included slot was previously empty. This is exactly the all-choice containment needed for arbitrary independent choices at successive ranks. No monotonicity of a greedy digit representation is being assumed.

The ordinary slots retain their factor-of-two divisor gap and linear extension ordering. H99's new fourth guard does not create cross-comparisons with fine slots: every ordinary tag misses a guard (so is not above the fine tag), while the certificate puts every strict divisor of the fine tag in the background (so no ordinary tag is below it).

Equal-rank comparable odd-core orders retain identical interpolation inputs. Deterministic decoding therefore gives identical sections there, as required when the rank does not increase strictly.

## Precision and denominator audit

The H99 slope s is a dyadic multiple of the odd clearing factor E. The intercept b is dyadic even though m=27/10: the common +s/m and −s/m terms cancel in the sum of the extremal endpoints defining its midpoint. It follows that Kc and Kdelta are dyadic.

For the interpolated profile, multiplication by K first cancels the parameter denominator K. Each donor occurs at most once per monomial, so the odd part of the remaining denominator divides E. The numerator of the final blend clears E. It is Kw (and the physical target Kw/2), rather than necessarily w itself, that has the claimed dyadic precision. This wording matters.

The modified absolute precision is 3401. The dominant contribution is the paired anchor exponent eu+el=967+1986=2953, together with the 400 blend bits, 19 damping midpoint bits, and the stated slack/normalization allowance. The independent ordinary bound is only 466. Dividing by the level scale improves these absolute bounds by 84−l, which justifies the D_l range 3317 through 3401.

The endpoints of the H99 target interval need not lie on the dyadic grid because the radius s/m can have denominator 27. This is harmless: the endpoints certify interval coverage; the actual profile masses satisfy the independently proved dyadic grid. The source verifier correctly avoids treating a nondyadic endpoint as a physical target.

## Interpolation audit

Removing one damping factor replaces the common error factor 2^(−20) by 2^(−19) and reduces R from 86 to 84. All support-convexity, sign, forced boundary, and equal-rank-invariance arguments are factorwise and remain valid. The proof uses many remaining balanced factors, so deleting this particular final damping factor does not invalidate an annihilation argument relying on an independent balanced coordinate.

The 400-bit blend preserves annihilation and boundary values and clears odd denominators as above. With the prescribed s and b, exact rational arithmetic gives total interpolation error divided by delta less than 339/1000, hence strictly below 10/27. The reconstructed decimal is approximately 0.3387485931919912. Decimal evaluation is not used for certification.

## Independent finite replay

A new implementation was written as `audits/h99_independent_audit.py`. It imports no construction or verification modules. It uses the stored slot lists and H100 manifest as data, plus the previously independently audited exact donor rational bounds. It implements the totient masses, tag-poset checks, normalization, fine ordering, control indices, interval intersection, physical-height tests, rank-dependent prefix, and reverse-parity decoding directly with standard-library integer/Fraction arithmetic.

The independent implementation reconstructs both refinements:

| Quantity | H100 rank refinement | H99 |
|---|---:|---:|
| Rank bound | 86 | 84 |
| Levels | 87 | 85 |
| Fine absolute precision | 3353 | 3401 |
| Maximum dyadic height | 3360 | 3408 |
| Background tags | 121210 | 245736 |
| Available fine slots | 3840 | 3840 |
| Endpoint inequalities | 174 | 170 |
| Independent exact grid reconstructions | 435 | 425 |

H99 normalized constants independently reconstructed:

- K = 28228205074382983540537992554159963315796843436605343386029.
- n = 195266343792994032798472735290025075808081607128255154605.
- H0 = 1322414346673510530563674176921864782267552196370366464.
- M0 = 4194304.
- Ordinary denominator bound = 466.
- Flat-helper precision = 2003.
- Positive intercept-interval width divided by s is approximately 0.009679379299510998.

The output is `audits/h99_independent_results.json`, including exact rationals and SHA-256 hashes of its inputs. It reports PASS. The grid reconstructions include the first and last actual grid points in each admissible interval and neighboring/interior grid points, avoiding the nondyadic-endpoint issue.

The sample grid replays are regression checks. Universality over every actual free order follows from the symbolic precision bound, affine endpoint containment, greedy residual induction, and exact parity decoder, as reproduced in `audits/h99_section.tex`.

## Primality, modulus, and minimality

The inherited H100 certificate covers all retained core and spectator primes. Independent trial division proves 65537 prime. Deleting the final pair and adjoining 65537 gives exactly 99 distinct odd prime divisors, without changing any other prime exponent. Exact integer arithmetic gives bit length 8017 for N99.

The existence of some divisibility-minimal nonrigid modulus with at most 99 odd primes follows simply by taking a minimal nonrigid divisor of N99. This does not identify that divisor or rule out 98 or fewer odd primes. H99 itself is not a divisor of H100 because it introduces the new prime 65537; H100's height-3360 refinement is a proper divisor of H100.

## Publication recommendation

Include both refinements as **PROVED (computer-assisted)** with their complete mathematical proofs and the included exact finite certificates. Preserve the original H100 interpolation method in full; the refinements can refer internally to that established machinery. Include the 504-slot H99 ordering and both verification scripts in the referee package. Do not present the endpoint checks or the example grid decodings as a proof by themselves.

No unresolved mathematical gap was found in these refinements, conditional on the separately completed H100 proof audit and its arithmetic certificates.
