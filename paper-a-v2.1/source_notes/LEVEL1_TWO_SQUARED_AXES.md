# Two squared prime axes and the six-factor exception region

**Theorem (exact symbolic certificate).** For distinct odd primes
$p,q,r,s$ the core $K=p^2q^2rs$ satisfies terminal Fourier exclusion
for all nonzero equal-population downset differences. Consequently
$2^eK$ is rigid for every $e\ge0$. Prime-power localization allows
additional sufficiently large primes to arbitrary powers using
$B=36(36!)2^4 K$ in the reciprocal-sum bound.
[Prime-power extension theorem](LEVEL1_TERMINAL_PRIME_POWER_EXTENSION.md)

Use the [ordered-grid certificate](LEVEL1_FIVE_FACTORS_WITH_MULTIPLICITY.md)
on the 36-cell grid of shape $(2,2,1,1)$ in each prime ordering.
The number of cells changes the storage of signed differences but not
the Fourier or incidence identities. In increasing cell order, encode a
downset with membership bits $b_i$ by $\sum b_i3^i$. The absolute
difference of two such keys determines the signed pattern up to overall
orientation by uniqueness of balanced ternary representation. Subtract
only after rejecting comparable downsets, because zero population
requires both signs. To recover both masks, add $(3^{36}-1)/2$ and read
the ordinary ternary digits, where two denotes positive, zero negative,
and one neutral. All values fit signed 64-bit integers. The byte
lookup for polynomial coefficients uses five bytes per mask and its
signed intermediate sums have magnitude at most $36\cdot13^6<2^{31}$.

The same cube, incidence and prime-factor tests as in the
[cubed-axis proof](LEVEL1_CUBED_FOUR_AXIS_CORE.md) leave just one
ordered pattern. To explain the last exclusion explicitly, let
$p<q<r<s$ with $p,q$ the two squared primes. The positive and negative
masks are 18258070720 and 9296846848, respectively, when the four order
exponents are listed lexicographically starting at bit zero. Summing
the signed shell populations gives
$$
 p^2q^2(r-s)-p^2rs+p^2s+q^2rs-q^2r=0.
$$
Modulo $s$ this requires $s\mid q^2 r(p^2-1)$. All the factors
$q,r,p-1,p+1$ are positive and smaller than $s$, an impossibility.
The verifier implements this last check as reduction of the root
polynomial modulo the largest prime, admitting the smaller primes and
both their predecessors and successors as units.

Run `python experiments/level1_grid_gap_certificate.py` from the
repository root with Python 3.12.14, NumPy 2.3.5 and SymPy 1.14.0. Each row below starts with 11,194,933 distinct signed
patterns, up to orientation, and covers all choices of the primes. The
second column counts possible zeros at the population row and the
others count successive exclusions.
[Verifier](../experiments/level1_grid_gap_certificate.py)

| Exponents in increasing prime order | Root survivors | Cube | Incidence | Root factor | Row factors | Largest prime |
|---|---:|---:|---:|---:|---:|---:|
| 1,1,2,2 | 2915359 | 2914533 | 817 | 9 | 0 | 0 |
| 1,2,1,2 | 3743527 | 3742267 | 1248 | 12 | 0 | 0 |
| 1,2,2,1 | 4326542 | 4324839 | 1690 | 13 | 0 | 0 |
| 2,1,1,2 | 3772758 | 3771046 | 1688 | 24 | 0 | 0 |
| 2,1,2,1 | 4305351 | 4303010 | 2312 | 28 | 1 | 0 |
| 2,2,1,1 | 5166787 | 5163679 | 3072 | 31 | 4 | 1 |

The same-session replay completed all 24 ordered grids, including 18
from earlier certificates. Every output had an empty residual list,
and in every row the exclusions summed to the number of population
survivors. The six new rows matched the table entry by entry.

## Finiteness for all six-factor odd parts

Only finitely many odd parts with six prime factors **counted with
multiplicity** can admit a genuine pair at any dyadic height. The
finite region is effective. Indeed, order the axes increasingly.
At any proper prefix the product $J$ of their full prime powers has
$\Omega(J)\le5$, so it satisfies the [five-factor terminal theorem](LEVEL1_FIVE_FACTORS_WITH_MULTIPLICITY.md).
Let $n$ be the number of remaining axes. With
$B=\tau(J)\tau(J)!2^{\omega(J)}J$, localization implies that the
least remaining prime in a nonrigid example is at most
$$
 n\max\{8J4^n,\;2^{n+2}\max(J,B)\}.
$$
Otherwise the sum of the remaining prime reciprocals is strictly
below both required thresholds. This applies also to the empty prefix,
where $J=1$ and there are no equal-population differences. At every
step there are only finitely many prime choices and at most six exponent
choices, yielding a finite tree of odd parts. For each fixed odd part
at most one height can give a divisibility-minimal nonrigid modulus,
by comparability along the dyadic chain.

With the [cubed four-axis theorem](LEVEL1_CUBED_FOUR_AXIS_CORE.md)
and the lower-factor results, every order having at most six odd factors
and at most four distinct odd primes is rigid. At exactly six factors
the only unresolved shapes are $(2,1,1,1,1)$ and $(1,1,1,1,1,1)$,
both in the finite region above. The bounds do not decide which, if
any, of its candidates are actually nonrigid.

For both remaining shapes the least disputed size must have a common
nontrivial Fourier zero on the complement of every primary axis.
The same argument also proves all-height rigidity on the four-axis
family with exponents $(a,1,1,1)$ for every positive $a$.
[Primary-axis proof](LEVEL1_PRIMARY_AXIS_COMMON_ZEROS.md)

