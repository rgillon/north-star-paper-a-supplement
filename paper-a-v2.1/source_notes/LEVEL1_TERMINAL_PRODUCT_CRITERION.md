# Product criterion for terminal Fourier interpolation

For a nonzero real radial shell function $d$ on an odd cyclic group,
write $U_d$ for the upward closure of its support in the divisor lattice.
Call its interpolation feasible if a real radial $g$ exists with
$g=|d|$ on $U_d$ and $\widehat g$ zero wherever $\widehat d$ is nonzero.
This definition requires neither population balance nor an ideal pair.

**Product theorem.** For coprime odd orders and nonzero radial functions
$d_1,d_2$, interpolation for $d_1\otimes d_2$ is feasible if and only
if interpolation for at least one factor is feasible.

Indeed $U_{d_1\otimes d_2}=U_{d_1}\times U_{d_2}$ and the nonzero
Fourier support is also the product of the two nonzero supports. A
solution for the first factor tensored with $|d_2|$ is therefore a
solution, and the roles may be exchanged. Conversely, in the finite
real radial function space use the inner product summing over individual
group elements. If factor $i$ is infeasible, the restriction of $|d_i|$
to $U_{d_i}$ lies outside the image under restriction of the radial
Fourier space supported on the zeros of $\widehat d_i$. Orthogonal
projection gives a real radial functional $w_i$, extended by zero off
$U_{d_i}$, annihilating that space and pairing nontrivially with
$|d_i|$. Fourier orthogonality gives $\operatorname{supp}\widehat w_i
\subseteq\operatorname{supp}\widehat d_i$ (all spaces are radial and
the zero set is a union of conductor shells). If both factors fail,
$w_1\otimes w_2$ annihilates every allowed product-group interpolant
by Fourier orthogonality, but pairs nontrivially with the prescription
by multiplication of the two pairings. This proves necessity.

**Application to actual mates.** Let the odd part $K$ be squarefree.
In any distinct homometric divisor-ideal pair on $2^eK$, let $L$ be the
least disputed dyadic fiber size and $d$ its signed shell indicator.
Then the multiaffine polynomial $D(X)=\sum_s d(s)\prod_{i\in s}X_i$
has an irreducible factor using at least five odd axes. More precisely,
one irreducible factor is itself a feasible signed downset difference
with zero population at the actual prime predecessors. This conclusion
concerns the least single size even when many sizes are disputed.

To prove it, take the dyadic character of kernel order $L$. Normalize
the partial transforms by $L$. Each physical fiber contributes its size
divided by $L$ if the size fits inside the kernel, and zero otherwise.
The signed difference is $d$ by the common-fiber theorem DF1 and the
choice of $L$. The sum $g$ has Fourier support in the zeros of
$\widehat d$ by equality of squared real transforms at each labeled
character. Above a dispute of size $L$ the empty side stays empty,
and any nonempty fiber on the other side has size exactly $L$ by
nesting and minimality of $L$. Thus $g=|d|$ on $U_d$.

The sign supports of $d$ are mutually incomparable and individually
order-convex: along a chain between same-sign size-$L$ disputes the
occupied side stays at size $L$ and the empty side stays empty. Hence
$d$ is a difference of two downsets. Explicitly, in the downward closure
of the union of its two sign supports, remove the upward closure of the
negative support for the first downset and of the positive for the second.
In that downward closure each removal deletes exactly its own sign.

Factor the polynomial over the rationals. Degrees in each variable add,
so the factors use disjoint axes. By Gauss's lemma take primitive integer
factors, with the overall sign absorbed in one factor. Every nonzero
coefficient is a single product and has magnitude one, so each factor
has coefficients in $\{0,1,-1\}$. Include an unused axis as the local
shell vector $(1,0)$. The coefficient array is now a tensor product.
Each factor is, up to sign, a downset difference by restricting the two
actual downsets to supported fixed shell states in the complementary
coordinates. Repeated use of the product theorem requires at least one
feasible factor. An unbalanced sign factor $f$ is infeasible: $f$ itself
is a supported dual functional since $\langle f,|f|\rangle=\sum_x f(x)
\ne0$. This also rules out unused axes and constants. Every balanced
downset difference on at most four squarefree axes is infeasible by the
[small-core terminal theorems](LEVEL1_PRIME_POWER_TERMINAL_GRIDS.md).
The feasible irreducible factor therefore uses at least five axes.

In particular, on a fixed five-prime core any feasible nonzero balanced
signed downset difference uses all five variables and is irreducible.
For the large-prime tower test this eliminates every missing-axis or
reducible pattern on the five-prime core, even for multiple disputed
sizes: the [tower proof](LEVEL1_FIVE_CORE_PRIME_TOWERS.md) gives the
prescription up to an error at most $2/p$ on the last disputed section.
For fixed $K$ the finitely many infeasible patterns have positive
distance in sup norm from their respective restricted Fourier spaces.
The Fourier support inclusion in that proof holds for $p-1>K$.
Thus beyond a bound depending only on $K$ any surviving least-size
pattern on that section is irreducible in five variables, uniformly in
both tower heights. The product criterion alone does not decide feasibility
of these irreducible patterns or their realization by mates.

