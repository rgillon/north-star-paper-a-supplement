# Recovering dyadic fiber sizes from autocorrelation

Robert Gillon’s North Star project · 17 November 2026 · proof candidates submitted for independent review

Let $H_d=\{x\in\mathbb Z_N:d\mid x\}$ for $d\mid N$. For a union $A$ of these subgroups write $C_A(t)=|A\cap(A-t)|$. Throughout, $A$ and $B$ are actual unions on the same labeled cyclic group. The claim of rigidity means literal equality follows from $C_A=C_B$ at all shifts.

**Theorem.** Write $N=2^eM$ with $M$ odd. At any $e\ge0$, rigidity holds in each of the following cases.

1. $\Omega(M)\le5$, counting odd prime factors with multiplicity.
2. $\Omega(M)=6$ and $\omega(M)\le4$, where $\omega$ counts distinct odd primes.
3. $M=p^a q^b r$ for distinct odd primes and arbitrary positive $a,b$. There is also no restriction on exponents when $M$ has at most two distinct primes.

The first two statements have finite symbolic proofs. The third has an elementary proof. The single-size five-factor step uses [FM1–FM2](FIVE_THRESHOLD_MARGINALS_2026_11_17.md). The common-fiber argument below reconstructs the prior [CU6 lemma](CODEX_SHARED_FIBER_RIGIDITY.md). All remaining steps, including the six-factor single-size check, are given here.

## DF1. What a single-prime Fourier measurement recovers

Use the Chinese remainder coordinates. An odd shell $s=(s_1,\ldots,s_n)$ means exact additive orders $p_i^{s_i}$, with $0\le s_i\le a_i=v_{p_i}(M)$. It contains $w_s=\prod_i\varphi(p_i^{s_i})$ points. At every odd point the dyadic fiber of $A$ is empty or the unique subgroup of order $2^j$. Its size $h_A(s)$ is nonincreasing in the coordinatewise order. Its nonempty projection $I$ is a downset. Use $h_B,J$ for the competitor.

For conductor exponent $c$ in a $p^a$-axis, the real, unnormalized Fourier row on shell $u$ is
$$
 R_{p,a}(c,u)=
 \begin{cases}
  \varphi(p^u)&u\le a-c,\\
  -p^{u-1}&u=a-c+1,\\
  0&u>a-c+1 .
 \end{cases}
$$
In particular its value at $u=0$ is always 1. This follows by subtracting adjacent subgroup character sums. The product of these rows over coordinates is the full odd transform. Since each set is invariant under negation, homometry gives equality of squares of the real Fourier coefficients.

A Fourier coefficient on a pure prime-power axis, including 2, is nonnegative. In each fixed complementary fiber the sum is the size of its subgroup or zero, and the empty fiber contributes zero. We can therefore recover its signed value. Fourier inversion on a pure odd axis gives the population in every fixed physical state along that axis. Equivalently, the total signed dyadic population in each fixed shell on that axis is zero.

For an odd conductor $c$ let $a_j(c)$ be the odd Fourier transform of the indicator that $h_A=2^j$. Define $b_j$ similarly. A dyadic character with kernel order $2^k$ has full coefficient $P_k=\sum_{j=0}^k2^ja_j$ at $c$. Define $P'_k$ for $B$. The two integers $D_k=P_k-P'_k$ and $S_k=P_k+P'_k$ obey $D_k S_k=0$. If $a_j-b_j$ is odd, the increments of both integers at step $j$ are $2^j$ times odd numbers. Exactly one incoming component is zero and the other is nonzero, since otherwise both outgoing components would be nonzero. The incoming zero becomes a number of valuation exactly $j$. It stays nonzero through every subsequent step, whose increments are divisible by $2^{j+1}$. The other component must stay zero. A second odd increment is impossible. Use $(D_{-1},S_{-1})=(0,0)$ at the start.

It follows that the fibers coincide wherever both are nonempty. If they differed, choose a minimal such odd shell $s$. Each smaller shell below it belongs to both projections and has identical fiber size. In axis $i$ choose $c_i=a_i-s_i+1$ when $s_i>0$ and $c_i=0$ otherwise. Modulo two, the product row selects only corners having each coordinate either 0 or $s_i$. Both unequal size labels at $s$ would have odd transform difference, which the preceding argument excludes.

Consequently $d=1_I-1_J$ signs every disputed shell. Write $H_s$ for its nonempty dyadic size. Opposite signs never occur at comparable shells. On comparable disputed shells the size is nonincreasing. We use these three necessary equations.

* For each nonempty distinct size $H$, $\sum_{s:H_s=H}w_sd(s)=0$. At the odd root both partial dyadic sums are nonnegative. Subtract consecutive equal partial sums and cancel the common fibers.
* For each fixed shell value $s_i=t$, $\sum_{s:s_i=t}w_s H_s d(s)=0$. This is the pure odd-axis population equality above. Each such slice contains both signs or neither. The projections agree on every one-axis shell, because its presence is equivalent to positive population in that slice.
* At each odd conductor, at most one size class has odd signed transform. In particular, where the transform of $d$ is even, the intersection of *each* size class with the parity row has even cardinality.

All sums over disputed shells ignore common cells.

## DF2. Three odd coordinates, with one exponent equal to one

With two odd axes, fix a state in the first. The second-axis projections of the two profiles are initial segments of a chain. On their intersection the fibers agree by DF1. The two resulting sections are nested. The population equality in the chosen first state forces their equality. One or zero odd axes follows as well.

Now write shells $(i,j,k)$ with $k=0$ or $1$. Assume some discrepancy and take the globally least $j$ at which it occurs. A discrepancy at $(i,j,0)$ would need an opposite discrepancy at the same $i$. Such a point cannot have $k=0$ by comparability along $j$, and at $k=1$ it must have $j'<j$ to be incomparable. This contradicts minimality. Thus every discrepancy at the selected $j$ has $k=1$. Those points form a chain in $i$, so all have the same sign, contrary to their zero signed population in that $j$-slice. This proves the third assertion of the theorem.

## DF3. A Fourier test for the smallest disputed size

Let $H$ be the least size that occurs at a disputed shell, and let $z$ be the set of cells with this size among the disputes. Take a dyadic character with kernel order $H$. Its odd-space profile is $f_A=h_A1_{h_A\le H}$ for $A$ and analogously for $B$. Hence $f_A-f_B=H d 1_z$. At every odd conductor where $\widehat{d1_z}\ne0$, the sum $g=f_A+f_B$ has Fourier coefficient zero.

Choose a maximal element $s$ in $z$. Then $g(s)=H$, and $g$ is zero strictly above $s$. Indeed the empty side stays empty above it. On the other side any nonempty upper fiber would itself be disputed, of size at most $H$. Minimality of $H$ and maximality of $s$ exclude this. For $s_i>0$, the difference of the two rows with conductors $a_i-s_i$ and $a_i-s_i+1$ vanishes below order $p_i^{s_i}$, and equals $p_i^{s_i}$ at that shell. If $s_i=0$, use one common conductor of choice in that coordinate. Its value at zero is 1. Alternating the rows in this product of conductor pairs therefore yields
$$
 \sum_x g(x)\prod_i R_i^{\rm difference\ or\ fixed}(x_i)
       = H\prod_{i:s_i>0}p_i^{s_i}\ne0 .
$$
Here $x$ runs over shells and the row factors already include their multiplicities. If every row of the product cube is known to have nonzero transform of $d1_z$, this is a contradiction. Odd parity suffices to mark such a row. It also suffices if a positive-coefficient polynomial certificate excludes simultaneous vanishing of that row and the root of $d1_z$.

## DF4. A finite symbolic proof through five factors

For one disputed size, the size cancels from all pure odd-axis measurements. The odd projections are downsets with identical cardinality and pure marginals. For $\Omega(M)\le5$ they coincide by FM1–FM2 applied on the odd group, and DF1 recovers the entire sets.

Suppose several sizes occur. DF2 handles every grid of at most three distinct axes through five factors. Four axes require shape $(1,1,1,1)$ or $(2,1,1,1)$, up to permutation. Five axes require five binary coordinates. The following procedure treats all distinct prime values and dyadic exponents at once.

1. Enumerate downsets by visiting cells in a linear extension. A cell can be added exactly when all its immediate predecessors are present. Group downsets by their one-axis shells. Compute every pair's two signed difference masks, remove pairs failing both-sign presence in a coordinate state, and quotient by overall sign and permutations of equal-bound axes.
2. For each such support search the linear system over $\mathbb F_2$ from the even parity rows in DF1. Enumerate its entire kernel by elimination. A proper candidate for a size class and its complement both need both signs. A test on free positive cell masses can preclude some entire supports. It uses only the signed state equalities of DF1.
3. Restore **every** equal-bound axis permutation before imposing numerical order on the primes. Write $p_i-1=i+x_1+\cdots+x_i$ for $1\le i\le n$. The $x_i$ are strictly positive for ordered distinct odd primes. The physical weight for a positive exponent $u$ is $(p_i-1)p_i^{u-1}$. For each candidate class both its signed root and the whole signed root must vanish. Eliminate a candidate only when an integer linear combination has one weak coefficient sign with at least one nonzero coefficient.
4. Enumerate every partition of the support by the surviving proper classes. A partition requires the signed sum of the **actual weighted masses** $w_s H_s$ to vanish in each state and in each class. The latter follows by multiplying that class’s root by its common size. Treating these masses as independent positive numbers gives a valid affine separation test. For each remaining possible smallest class apply DF3. For the rest, take every order of the other size classes and exclude orders violating the monotonicity in DF1.
5. Normalize the smallest size to 1 and write the next sizes as $1+y_1,1+y_1+y_2,\ldots$ with all $y_i>0$. These ratios allow every possible assignment of dyadic sizes. Expand the unweighted class equations and weighted coordinate-state equations in $x,y$. Common positive factors in an equation may be divided out. If an integer linear combination, possibly including positive monomial multiples of these equations, has nonnegative coefficients and at least one positive coefficient, it cannot vanish. This certifies the exclusion for unbounded primes and unbounded height. For the six-factor extension we also use the elementary arithmetic check in DF5.

The C++ implementation uses integer arithmetic for the grid, shell masses and first two-root exclusions. For the later affine and polynomial tests a linear program proposes multipliers. After rational approximation, the program recomputes every coefficient using arbitrary-precision integers. An unsuccessful proposal retains the candidate. Integer coefficient bounds in the C++ product operations follow by evaluating the nonnegative cell products at all $x_i=1$. Each $p_i=2i+1$ at that point, and the sum of weights of all grid cells is $\prod_i(2i+1)^{a_i}$. This is at most $103{,}950$ over the listed ordered shapes. Individual signed coefficients are bounded by this number, triple-score coefficients by three times it, and two-row minors by twice its square. The last operation uses signed 64-bit arithmetic. No numerical optimization outcome alone eliminates a case.

| Ordered grid | Downsets | Labeled signed differences after sign-presence tests | Proper size partitions after root screens | Surviving assignments |
| --- | ---: | ---: | ---: | ---: |
| $1111$ | 168 | 24 | 0 | 0 |
| $2111,\ 1211,\ 1121,\ 1112$ | 887 each | 372 each | 75, 70, 11, 5 | 0 |
| $11111$ | 7,581 | 321,367 | 4,911 | 0 |

The last line first reduces the supports to 3,065 permutation classes. Of the 4,911 partitions, 3,636 fail the affine test. Among the others, DF3 excludes 2,390 choices of smallest class. A coefficient certificate excludes 172 remaining ordered assignments, and monotonicity excludes two. The units of these counts differ because one partition can admit several choices of smallest class.

## DF5. Six odd factors on at most four axes

With three odd axes and no exponent one, the only new grid is $(2,2,2)$. DF1’s sign-presence and axis tests leave three differences up to orientation. In each there is a slice with one coordinate zero and exactly two opposite cells having the other orders $(1,2),(2,1)$. Their odd-shell masses have ratio $p/q$ for distinct odd primes. Their dyadic sizes have ratio a power of two. Equal populations in this slice are impossible. This handles every size assignment at once.

With four axes the new exponent patterns are $(3,1,1,1)$ and $(2,2,1,1)$, with four and six placements in prime order. For several disputed sizes, the DF4 procedure excludes every assignment in the four placements of the first shape. For the second shape, twelve ordered assignments across the six placements survive its polynomial and Fourier tests. Each has a size class with two cells on one side and one on the other. After canceling identical coordinates the class root is
$$
 (p-1)+p(p-1)=r-1 .
$$
This would give $r=p^2$ for distinct primes. The verifier detects this relation from the shell coordinates of each remaining assignment.

For one size, its common value cancels from every state equation. For each labeled sign support form the root equation, the equations on exact single-axis states, and the equations on single-axis upper intervals obtained by adding states. Try sums of up to three of these polynomials with signs and repetition. A one-sign nonzero polynomial excludes the case for every positive gap vector. The first shape leaves zero. The second shape leaves ten labeled patterns across all six orderings. Eight force the same $p^2=r$ identity on one state slice. The other two have identical rectangular structure after relabeling.

Here is the exact Fourier contradiction for that structure. The two doubled coordinates are $p,q$ and the singles are $r,t$. Group the axes into $(p,r)$ and $(q,t)$. The four disputed cells are the Cartesian product of the two order choices $\{p^2,r\}$ with $\{q^2,t\}$, with opposite signs on the two diagonals. Coordinates omitted from an order choice are zero in shell notation. All four cells are maximal in the disputed set, so $g=f_A+f_B$ from DF3 equals $H$ there and vanishes at every strictly larger cell.

In the first block use either the row of conductor $p$ or that of conductor $r$, trivial on the other coordinate. Modulo two, each row selects exactly one of the two disputed block cells. The analogous two rows in the second block give four product rows. Each is odd on the signed difference, so all four must vanish on $g$. Subtract the two rows in each block and tensor the differences. In the first block the difference vanishes away from the upper closure of the two disputed cells. At $p^2,r$ its values are $-p^2,r$. In the second block they are $-q^2,t$. The resulting alternating sum on $g$ equals $H(r-p^2)(t-q^2)\ne0$. This excludes the last two patterns and completes the proof.

The [companion programs and exact counts](../research/collaboration/all-height-five-factors-2026-11-17/) reproduce the finite parts. At six factors the five-axis pattern $(2,1,1,1,1)$ and the six-axis squarefree pattern are outside this enumeration. The latter is constrained further in the [maximal-shell argument](MAXIMAL_SHELL_LOCALIZATION_2026_11_17.md).

