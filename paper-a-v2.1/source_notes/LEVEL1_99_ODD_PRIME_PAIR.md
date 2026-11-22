# Nonrigidity with 99 distinct odd prime divisors

Let $M$ be the odd integer in the [H100 manifest](../experiments/ppi_100_certificate.json).
Put $p=2101327$, the last entry of the H100 damping list, and $q=p(p-1)+1$.
There exist distinct, nontranslated divisor-generated subsets of the same
additive cyclic group of order
$$
N=2^{3408}\frac{65537M}{p^2q}
$$
with equal autocorrelation at every labeled shift. In particular there
exists a divisibility-minimal nonrigid modulus with **at most 99 distinct
odd prime divisors**. Its identity is not determined here.

## Odd trade and normalized target masses

Use [H100](PPI_100_PRIME_PAIR.md) with its original upper donor, lower
donor, modular helpers, and first **19** damping blocks. The signed
difference $e$ on this odd core is the product of the signed local
shell functions, zero off their supports. It has nonempty order-convex
support and agrees in sign at comparable supported orders, by the
factorwise H100 proof. These assertions hold for any number of damping
factors. The rank bound is now $R=46+2(19)=84$.

We use the interpolation formula and donor witnesses of H100 section 1.
Its error estimate for variable damping count is derived in
[PPI-UP1](PPI_UNBALANCED_PAIR.md) and
[the mixed-slot interpolation](PPI_MIXED_SLOTS_PAIR.md).
In their notation the rational profile annihilates $e$, has the prescribed
forced boundary values, and differs from $c-\delta r$ on free orders by at most
$$
2^{-19}\left[(1+U)(3+L)/4+(1+L)\eta+
\tfrac32R\delta(1+U)(1+L)\right].
$$
Blend it with the original paired dyadic anchor for these 19 blocks using
$\lambda=E\lfloor2^{400}/E\rfloor/2^{400}$ on the rational profile, where
$E$ is the donors' odd clearing factor. Both terms have the same boundary
values and annihilate $e$. The sum of the error above and the blend error
$(1-\lambda)[5+3(1+U')(1+L')]$ is less than $10\delta/27$ for the parameters
below. The [exact verifier](../experiments/ppi_99_prime_pair.py) checks this
inequality using the original donor bounds, with $0<c-R\delta<c<2$.
The blended profile $w$ satisfies $e*w=0$. The numerator of $\lambda$
clears the odd denominator, so $Kw$ is dyadic.

## The spectator and the ordinary slots

Append the prime 65537 to the ordered H100 tag tuple with exponent one.
Its full exponent is a fourth guard in addition to 3, 7, and 113.
Use the same variable kernel of size $B$ and control kernel of size
$C=34891$. In the distinguished tag $g$ all four guards are full and
the other coordinates are zero. Its mass is the former guard mass
times $65536$, so its odd normalization remains $h=83349$ and the new
fine unit is $M_0=64\cdot65536$. Write $Q$ for the enlarged tag quotient
order and $K=QB/h$. As before, normalized masses are physical cardinalities
divided by $hC2^{3408}$; the full spectator has normalized mass $K$.

The [ordinary-slot input](../experiments/ppi_99_slot_order.json) lists
504 distinct zero-based tags in the lexicographic convention of H100,
with 65537 last. Pair them successively. The extra guard permits a final
supply of tags missing 65537 while full on the original three guards.
Fill every tag missing some guard and absent from the list as the full
background. Let $n$ be its normalized mass. All other tags outside $g$
and the list are empty. The verifier checks all immediate divisor relations
of every background and ordinary tag, and checks that every proper subtag
of $g$ is background. Every listed tag misses a guard. Thus none is above
$g$, and none is below it. Membership in this fixed list is part of the
construction; no optimization oracle is required.

For slot $i$ let $W_i$ be its tag mass times $B/h$. Set $H_0$ to the
floor of $0.99$ times the least of the first 174 weights. Choose the unique
power of two $\beta_i$ with $H_0/2<W_i\beta_i\le H_0$.
Use groups numbered 0 through $R+159$, leaving later listed slots empty.
At level $l=R-r$ take the H100 frozen and unfrozen heights: groups below
$l$ use $2\beta_i$; group $k\ge l$ uses $\beta_i2^{l-k}$ with optional
doubling. All heights fit the dyadic group as checked exactly. An ordinary
multiple has weight at least twice that of its strict ordinary divisor,
so its base height is at most half as large. The slot order and group
scaling preserve nestedness for every digit choice, both within a level
and across consecutive levels.

Let $s$ be the nearest multiple of $E/2^{192}$ to
$130352847K/10^{12}$, ties to even. Let $C_l$ be the ordinary minimum and
$J_l$ the sum of the two unscaled normalized digit weights at group $l$.
Put $m=27/10$ and $\lambda_l=2^{l-R}$. Define $b$ as the midpoint of the
intersection, over all levels, of
$$
[C_l+5M_0\lambda_l/4-ls+s/m,\quad
 C_l+5M_0\lambda_l/4-ls+J_l+H_0-s/m].
$$
The verifier proves this interval has positive width. The $s/m$ terms
cancel in the midpoint of the extremal endpoints, so $b$ is dyadic.
Take $\delta=2s/K$, $c=2(n+b+Rs)/K$. The interpolation estimate gives
the requested mass $Kw/2$ in
$[n+b+ls-s/m,n+b+ls+s/m]$ at each free order.

## Grid and exact realization

Use the H100 fine-shell ordering and control subgroup indices above $g$.
Multiplying the guard mass by a power of two leaves the odd unit weights
$u_j$ and control indices unchanged. Apply the
[rank-dependent decoder](LEVEL1_H100_RANK_PRECISION.md) with absolute
precision $A=3401$: at level $l$ use exactly the first $D_l=A-(R-l)$
slots. Set $T_l=\sum_{j\le D_l}u_j2^{-j}$. Each slot has nominal height
$2^{l-R-j}$ times its old correction $M_0/2^{v_2(V_j)}$ and optional
doubling. The largest relative exponent correction is seven, so all
depths fit $3408$. Prefix inclusion, increasing base heights with $l$,
and the unchanged control divisibilities ensure the required nesting.

For specificity the verifier uses the H100 section 4 precision estimate.
For a dyadic rational $x$ write $d(x)$ for its reduced denominator
exponent. The profile and blend bound before dividing by $\lambda_l$ is
$$
400+19+\max\{d(Kc)+d_y+f+3,
\max(d(s),d(K\delta))+d_y,\ e_u+e_l\}+6+
\operatorname{bitlength}(M_0),
$$
where $f$ is the rounded flat-helper precision and $d_y,e_u,e_l$ are
the original donor precision bounds. This follows from the same product
expansion: each stationary damping midpoint costs one bit, each donor
occurs at most once per monomial, and the final blend costs at most 400
bits. Take the maximum with
$\operatorname{bitlength}(K)+159+R+5+\operatorname{bitlength}(M_0)$
for ordinary contributions. Evaluation gives 3401. Dividing by
$\lambda_l$ clears $R-l$ bits, giving exactly the required grid bound.

After subtracting $n+C_l+M_0\lambda_l(T_l+1/4)$, both interval endpoints
lie in $[0,J_l+H_0]$ at every level, by exact rational evaluation.
Subtraction is affine, so every interior target passes. Greedy selection
of the two ordinary digits per active group leaves at most $H_0$ after
the first group, and at most $H_0/2^j$ after depth $j$, because the two
digits sum to more than that bound and each is at most that bound.
After the final depth $159+R-l$, adding back the reserve and dividing
by $M_0\lambda_l$ puts the residual in $[1/4,3/4]$ on the grid
$2^{-D_l}\mathbb Z$. The odd-weight reverse parity decoder is exact
there: all selected tail weights have normalized sum less than $1/4$,
and the first eight unit weights encode the remaining prefix integer.
The target is therefore realized physically on every free order, with
deterministic equal choices on comparable orders of equal rank.

## Actual sets and labeled shifts

Use identical decoded sections at free orders. At signed support give
the two sides opposite full/empty spectator sections according to the
sign. Below support outside it give both full sections, and above it give
both empty sections. The forced choices and the nested decoded sections
are divisor downsets, hence actual subgroup unions. Their indicator
difference is $e(x)$ independently of the spectator point; their summed
indicator over the fiber is $|S|w(x)$ for spectator $S$. Evenness of both
indicators gives at every individual shift $(z,y)$
$$
C_A(z,y)-C_B(z,y)=|S|(e*w)(z)=0 .
$$
The population profiles differ on the signed support. Both are even
on an odd group, so if translation by $t$ identified them, $2t$ would
be a period and hence so would $t$, a contradiction. Finally the odd
prime count is $48+38+13=99$, all distinct. Among the finitely many
nonrigid divisors of $N$, any divisibility-minimal one is also globally
divisibility-minimal.

Run `python3 experiments/ppi_99_prime_pair.py` with assertions enabled.
The checker verifies the changed finite inequalities and reuses the original
donor calculations. The slot order and approximate trial slope were selected
numerically; only the stored order, rational slope rule and exact checks
are used in the proof. This result does not determine a minimal divisor
or prove an exclusion at 98 odd primes.

