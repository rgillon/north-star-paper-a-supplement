# Common zeros on a prime-power axis and rigidity with three squarefree residual primes

For any distinct odd primes $p,q,r,s$, every divisor-generated subset
of the cyclic group of order $2^e p^a qrs$ is determined, among
divisor-generated subsets, by its labeled autocorrelation. Both
$e\ge0$ and $a\ge1$ are arbitrary. The proof uses a common nontrivial
Fourier zero for the least disputed dyadic size on a changing primary
axis and the signed marginals of the three residual primes. We use
the common-fiber theorem DF1, which says that when both dyadic fibers
over an odd point are nonempty they have equal size.
[Common-fiber theorem](https://github.com/rgillon/north-star-math/blob/0930c856256add528c967cd1b35c742f88b4cf95/proofs/DYADIC_FIBER_SIZES_2026_11_17.md)

## Common zero for a single disputed size

**Theorem.** Suppose two distinct homometric divisor ideals have one
disputed dyadic fiber size $L$ on $2^e Kp^a$, where $K$ is odd, $p$
is odd prime and $p\nmid K$. Let $d_t$ be the signed dispute indicator
on $K$ at physical $p$-order $p^t$, $0\le t\le a$. If these section
functions are not all identical, then at some **nontrivial** conductor
of $K$ every one of $\widehat d_0,\ldots,\widehat d_a$ vanishes over
the integers. There is no lower bound on $p$ relative to $K$ and no
upper bound on $a$ or $e$.

First restrict to the subgroup with dyadic order $L$. Outside it the
sets agree, and shifts in the subgroup preserve each exterior coset,
so the two restricted autocorrelations agree. Disputes persist by DF1.
At the trivial dyadic character divide the two partial profiles by $L$.
Their difference at section $t$ is $d_t$ and their sum $g_t$ has
values in $\mathbb Z[1/2]$. To see that each $d_t$ has total zero, observe that the pure
$p$-axis Fourier transforms of the restricted ideals are nonnegative
sums of subgroup character sums, hence agree with sign, and inversion
recovers the section populations. Each nonzero $d_t$ therefore has
both signs.

Let $j<a$ be the first change, $d_j\ne d_{j+1}$, so $d_0=\cdots=d_j$.
For a section vector $u$ set $u_{a+1}=0$ and define
$$
 P_k(u)=\sum_{i=0}^k p^i(u_i-u_{i+1}).
$$
These are exactly its $p$-Fourier rows at kernel order $p^k$, by
expansion in subgroup indicators. At each residual character outside
the common zero set of all the $d_t$ there must be a nonzero difference
row, by inversion on the $p$ axis. Every row below $j$ is zero, so
some $k\ge j$ has nonzero difference. At that measurement the sum
row is zero. All sum rows with $k\ge j$ are congruent to $P_j(g)$
modulo $p^{j+1}$. Hence the Fourier transform of $P_j(g)$ on $K$
vanishes modulo $p^{j+1}$ off the common zero set. Reduction takes
place in $\mathbb Z[1/2]$. On radial functions inversion uses integer
Ramanujan sums divided by $K$, which is invertible modulo $p^{j+1}$.

If every nontrivial residual conductor lay outside the common zero set,
$P_j(g)$ would be constant modulo
$p^{j+1}$. At a residual shell where $d_j\ne d_{j+1}$, opposite signs
are impossible by nesting. If the dispute disappears, all the first
$j+1$ sections there have sum one, followed by zero. If it first
appears, the preceding sections have sum two, followed by one, because both
earlier fibers coincide and contain the full restricted fiber. In
either case $P_j(g)=p^j$ there. If $d_{j+1}$ is nonzero, its two
dispute signs force both restricted fibers to be full at the residual
origin in that section and all preceding ones, giving $P_j(g)=0$.
If $d_{j+1}=0$, then $d_0=d_j$ has both signs, making both fibers
empty at the full residual order in every section by nesting. Again
$P_j(g)=0$. The two values differ modulo $p^{j+1}$, proving the claim.

## Least class with arbitrarily many larger disputes

Let $L$ now be the *least* global disputed size. If its signed sections
$d_t$ change along the primary axis, their transforms still have a
common **nontrivial** zero for every odd $p$, even when larger disputes
occur. The proof first reduces to one special residue at $p=3$ and
then uses nesting to rule it out.

Use the dyadic character with kernel size $L$ on the original group. The normalized difference is
$d$, and the congruence for $P_j(g)$ outside the common zero set
holds exactly as before. At a disappearing dispute all the earlier
sums are one and the next sum is zero, because above it the empty side stays
empty and no smaller disputed fiber is possible. Thus $P_j(g)=p^j$.
At a newly appearing dispute the next sum is one. Earlier the occupied
side contains size $L$, but a larger fiber contributes zero to this
dyadic character. While its size exceeds $L$ the opposite fiber
is either equally large or empty by DF1. If its size equals $L$ before
the first change, both fibers there equal $L$ by the definition of $j$.
Therefore the preceding sums are a sequence of zeros followed by twos.
Using
$$
 P_j(g)=\sum_{t=0}^j\varphi(p^t)g_t-p^jg_{j+1}
$$
gives $p^j$ if all preceding sums are two, $-p^j$ if all are zero,
and $p^j-2p^b$ if the twos start at state $b+1$ with $0\le b<j$.
All these values are nonzero modulo $p^{j+1}$. Hence an empty common
zero set is impossible by inversion on the residual group.

Suppose the only common zero is trivial. All signed sections of this
class then have zero population. If $d_0$ is nonzero, its two signs
make both fibers empty at the full residual shell from the first state
onward. This gives $P_j(g)=0$ there and contradicts constancy.
Otherwise $d_0=\cdots=d_j=0$ and $d_{j+1}$ has both signs. At the
residual origin, both sides have common size at least $L$ through
state $j+1$. Their normalized sums through that state are again
zeros followed by twos. The origin value of $P_j(g)$ is consequently
zero or $-2p^c$ with $0\le c\le j$.

These origin residues cannot agree with the changing-shell residues
when $p\ge5$. For the mixed value $p^j-2p^b$ with $b<j$, comparison
with $-2p^c$ first forces $b=c$ by valuation and then leaves $p^j$.
For $\pm p^j$ comparison forces $c=j$ and then requires
$\pm1=-2\pmod p$. Only $p=3$ and the positive sign work. In this
exception the origin fibers exceed $L$ through $j$ and first attain
$L$ at $j+1$, while at *each* newly disputed shell both fibers are
exactly $L$ at every earlier state. It remains to exclude this case.

In that exceptional case no dispute larger than $L$ can occur. The
origin at $j+1$ caps all fiber sizes there and later by $L$. Suppose
instead that an earlier section has a larger dispute at a residual
shell $s$. All still earlier occupied fibers on its occupied side
are larger than $L$, and by DF1 the opposite side there is either empty
or equally large. Their sums at kernel $L$ are zero. From the larger
dispute onward its empty side stays empty, and up to state $j$ the
other side either remains larger than $L$ or vanishes, since no size-$L$
dispute occurs there and smaller disputes are impossible. Thus every
preceding sum through $j$ is zero at $s$. At $j+1$ the empty side is
still empty and the sum is zero or one. This gives $P_j(g)=0$ or
$-p^j$, neither congruent to the required $p^j$ modulo $p^{j+1}$
for $p=3$. Hence the supposed exception would be a one-size pair,
already excluded under the trivial-only hypothesis by the first theorem.

**Six-factor consequence.** Suppose the odd part has six prime factors
counted with multiplicity. In a genuine pair the least class must
have an exact common nontrivial section zero on the complement of
*each* odd prime-power axis. In particular, this applies to both the
two sections on a simple axis and the three sections on a squared axis.

A single odd axis already has terminal exclusion. For two or more
axes, the signed shell array cannot be a tensor product across any
partition into two nonempty groups of axes. A nonzero rank-one array
with entries in $\{0,1,-1\}$ factors, up to overall sign, into two
actual nonzero sections. Each is a signed downset difference since
the full least class is one and restricting downsets to a section
preserves downward closure. Each block has at most five factors
counted with multiplicity. A balanced factor fails terminal
interpolation there, and an unbalanced factor pairs nontrivially with
its own prescribed absolute value and hence also fails. The product
theorem now contradicts feasibility of the global least class supplied
by the homometric pair. In particular the sections on any one axis
cannot all be identical, and the common-zero theorem applies.
[Five-factor exclusion](LEVEL1_FIVE_FACTORS_WITH_MULTIPLICITY.md) and
[terminal product theorem](LEVEL1_TERMINAL_PRODUCT_CRITERION.md)

## Rigidity for an arbitrary prime power and three squarefree primes

Suppose there is a distinct homometric pair on $2^e p^a qrs$. Index
the three residual primes by $1,2,3$ and write $S$ for a subset of
these indices. On shell $(t,S)$ let $F_t(S)$ be the difference of
the two full dyadic fiber sizes. Its sign is fixed as $t$ varies
where it is nonzero. Indeed, at a dispute one side is empty by
common-fiber rigidity (DF1), and that side stays empty at higher
$p$-orders. The same observation forbids opposite signs on
comparable residual shells within any one section.

Put $u_i=q_i-1$ for the three prime predecessors and set
$$
 W_S=\Big(\prod_{i\in S}u_i\Big)
             \sum_{t=0}^a\varphi(p^t)F_t(S).
$$
The fixed sign implies that $W_S=0$ precisely when $S$ is never
disputed. The pure $p$-axis measurements give zero total mass
difference at each $t$. Indeed the pure-axis transforms of either
ideal are nonnegative, as sums of transforms of subgroup fibers or
empty fibers on that axis. Their squares agree, so inversion recovers
the signed marginals. Every disputed section therefore has both signs.
Neither the bottom nor the top shell of the residual triple can be
disputed, since each is comparable with every shell.

The pure residual-axis measurements likewise give zero signed mass
over the points with coordinate $i$ nonzero, and the total signed mass
is zero. Thus $\sum_{S\ni i}W_S=0$ for each $i$ and $\sum_S W_S=0$.
Only rank-one and rank-two shells contribute. Summing the three
marginals and subtracting the total yields $\sum_{|S|=2}W_S=0$.
For each $i$ the individual marginal now gives
$$
                  W_{\{i\}}=W_{\{1,2,3\}\setminus\{i\}}.
$$
In particular each singleton and its complementary doubleton carry
one fixed sign whenever either is disputed.

At any given $t$, every disputed singleton has the same sign as every
disputed doubleton present there. They are either comparable or
complementary, and both cases have just been proved. If both ranks
occurred, all disputes would have the same sign, contrary to balance
at $t$. Hence every disputed section lies in a single rank and has at
most three disputed shells.

Let $L$ be the least disputed size. Its signed sections cannot all
coincide. Otherwise some shell of that class occurs at every $t$,
while its complementary shell must also be disputed somewhere. This
contradicts the single-rank conclusion. Hence all size-$L$ sections
have a common nontrivial Fourier zero, by the primary-axis theorem.
At a residual conductor each shell contributes a product of nonzero
entries, one factor $x_i\in\{-1,u_i\}$ for each $i\in S$.
There must therefore be at least two size-$L$ shells in each section
where this size occurs.

In each such section there must also be a larger dispute. To see
this directly, suppose every dispute there has size $L$. The signed
shell sum vanishes at $(u_1,u_2,u_3)$ by total mass balance in that
section. Two opposite terms cannot balance, since the three physical
shell weights in either rank are distinct. So all three terms occur,
with both signs. In rank one the polynomial is $\sum_i\epsilon_i x_i$
with $\epsilon_i\in\{-1,1\}$. At a nontrivial zero the number of
coordinates evaluated at $-1$ must be even, by parity, so exactly two
are used. If their signs differ they cancel and the third term is
nonzero. If their signs agree, the third sign is opposite, and all
three evaluated terms have the same sign. In rank two divide by the
nonzero $x_1x_2x_3$ to obtain $\sum_i\epsilon_i/x_i$. With a single
$-1$, the other two reciprocal magnitudes sum to at most $1/2+1/4<1$.
With two, their signed sum is zero or has magnitude two, and the
remaining reciprocal has magnitude in $(0,1/2]$. With three, the
sum is odd. This exhausts all nontrivial conductors, contradicting
the common zero. Thus the section has exactly two size-$L$ shells
and one larger disputed shell.

The two size-$L$ shells must have opposite signs. In rank one their
two-term Fourier polynomial is $\pm x_i\pm x_j$, and in rank two it
has the same form times a common nonzero $x_k$. A same-sign pair
cannot vanish because $x_i+x_j\ne0$ for $x_i,x_j$ equal to $-1$ or
positive even integers. All three complementary shell pairs are now
active. Of their three fixed signs, one occurs twice. In every section
with a size-$L$ dispute the two size-$L$ shells have opposite signs,
so the single larger dispute must have this repeated sign, independent
of the section and of its rank. Balance of the full dyadic mass at
that $t$ forces the signed population of its size-$L$ shells to have
the opposite sign. Consequently these nonzero populations all have
the same sign. Their sum weighted by the positive $\varphi(p^t)$
cannot vanish. Yet at the dyadic character of kernel order $L$, the
normalized odd difference $d$ and sum $g$ have disjoint Fourier
supports by homometry. On the support of $d$ the sum equals one,
so Parseval gives $0=\langle d,g\rangle=\sum_x d(x)$, the total signed
population of the least class. This contradiction proves rigidity
for both arbitrary tower heights.

