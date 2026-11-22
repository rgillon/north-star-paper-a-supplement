# A nonrigid proper divisor of the H100 modulus

Let $M$ be exactly the odd number factored in the
[H100 manifest](../experiments/ppi_100_certificate.json).
There exist distinct, nontranslated subgroup unions in
$\mathbb Z_{2^{3360}M}$ with equal autocorrelation at every labeled shift.
In particular the H100 modulus $2^{3446}M$ is not divisibility-minimal.
For this fixed odd part, the first nonrigid dyadic height is at most 3360.

## Construction and proof

Use the odd core, odd spectator, signed difference $e$, normalized
population profile $w$, ordinary slots and fine-slot order of
[H100, sections 1–4](PPI_100_PRIME_PAIR.md). Those sections prove that
$e*w=0$, that the forced full/empty sections are compatible, and that the
requested free-section masses are in the prescribed intervals. All odd
prime powers are exactly as in that construction. Change the dyadic
factor to $2^{3360}$ and use the following rank-dependent decoder.

Write $R=86$, $A=3353$, $m=7/2$, and $M_0=64$ (the fine mass unit).
At level $l=R-r$, put $b_l=R-l$, $\lambda_l=2^{-b_l}$ and
$D_l=A-b_l$. Keep only the first $D_l$ fine slots and leave the
others empty. On those slots use the H100 nominal heights
$2^{-b_l-j}$ and optional doubling, with the same control subgroup
and the same correction $M_0/2^{v_2(V_j)}$. Fill ordinary and
background slots as in H100; perform the H100 greedy algorithm with
the new baseline. In particular all ranks use the same target profile.

**Physical realization.** Every used fine slot has corrected dyadic
depth at most
$$
 b_l+j+v_2(V_j)-6\le b_l+D_l+7=3360.
$$
Ordinary heights are legal at this precision as verified below.
Truncation preserves closure within a level, since every strict divisor
on the fine tag comes earlier in the slot order. Across increasing levels
$D_l$ increases. A surviving fine slot's new minimum height equals its old
maximum, and a newly introduced slot was previously empty. This gives all
the required containments between different ranks. Equal-rank comparable
orders have identical targets by H100 and receive identical choices by
deterministic decoding.

**Precision and range for every target.** In normalized mass units write
$X=Kw/2$ for a requested mass. The precision estimate of H100 section 4
bounds the denominator exponent of $X/M_0$ by $A$ *before division by*
$\lambda_l$. The same bound covers the background and every ordinary
contribution divided by $M_0$ (the explicit ordinary bound there is
$\operatorname{bitlength}(K)+159+R+5+\operatorname{bitlength}(M_0)$).
Dividing these dyadic quantities by $\lambda_l$ multiplies by $2^{b_l}$,
so their grid is now $2^{-D_l}\mathbb Z$. This is the reason the last
$b_l$ slots can be discarded. Their contributions to the baseline must
also be removed.

More precisely let $u_j$ be the H100 odd fine weights and
$T_l=\sum_{j=1}^{D_l}u_j2^{-j}$. Subtract the minimum mass (background,
ordinary minimum and $M_0\lambda_l T_l$), then reserve $M_0\lambda_l/4$.
For both endpoints of the full H100 target interval at each level the
result lies in $[0,J_l+H_0]$, as checked with rational arithmetic by the
[rank-precision verifier](../experiments/ppi_100_rank_precision.py).
The subtraction is affine, so the check covers the whole interval.
The H100 greedy bound applies without change to this starting interval:
the first pair leaves at most $H_0$, and each subsequent pair halves
the guaranteed upper bound through the final depth $159+R-l$.
After adding back the reserve and dividing by $M_0\lambda_l$, the
required fine increment belongs to $[1/4,3/4]$. The grid argument above,
$T_l\in2^{-D_l}\mathbb Z$, and $D_l\ge8$ place it on the required grid.

Decode the increment in reverse order as in H100 with $D_l$ digits.
Every tail weight is odd. Its selected normalized mass is at most the
full tail mass, strictly less than $1/4$ by the H100 tail bound.
Consequently after tail decoding the prefix integer lies between zero
and $2^8-1$, and the first eight unit weights complete the exact mass.

**Homometry.** On signed support use opposite full/empty sections, on
the strictly lower and upper forced regions use the common full and empty
sections respectively, and elsewhere use the common decoded sections.
These binary sets are downsets by the inclusions proved above and in H100.
For their indicators $f,g$ the difference on each spectator fiber is
the constant $e(x)$ and the fiber sum of $f+g$ is $|S|w(x)$ for the
new spectator $S$. Both are even indicators. Thus at every shift $(z,y)$,
$$
C_f(z,y)-C_g(z,y)=((f-g)*(f+g))(z,y)=|S|(e*w)(z)=0 .
$$
Their population profiles on the odd core differ on signed support. If
one were a translate of the other by $t$, their evenness would make
$2t$ a period, hence also $t$ by odd order, contradicting the different
profiles. This completes the construction.

**Reproduction.** Run with assertions enabled:

~~~sh
python3 experiments/ppi_100_prime_pair.py
python3 experiments/ppi_100_rank_precision.py
~~~

The first command reran the original arithmetic certificate, including
435 original mass checks. The second checked the changed interval bounds
at all 87 levels, recovered 261 target masses exactly, and obtained
3360 as the maximum permitted fine depth. The proof for interior targets
uses the interval and grid arguments, not sampling. This is a modified
construction with the same odd part; minimality of the new modulus
and the full prime-factor classification remain open.

