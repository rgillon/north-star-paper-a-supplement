# Rigidity proof audit for Paper A Version 2

Audit date: 22 November 2026. This is a mathematical source audit within the
current ChatGPT manuscript-preparation session. It is not an external referee
report or a Lean formal verification. The computational replay uses the
repository implementation; it is not an independent implementation.

## Source versions

The following files were fetched in full from `rgillon/north-star-math` at
commit `20dd528d5fdd4bf7a16cc5ce05f9cd25ae155c25`:

- `proofs/LEVEL1_PRIMARY_AXIS_COMMON_ZEROS.md`, blob
  `a079243b20ed8d4b7d628b408e606f04c0dbf446`.
- `proofs/LEVEL1_TWO_SQUARED_AXES.md`, blob
  `c79bf9a9888ae379122797dc5656b5b8f0196762`.
- `proofs/LEVEL1_TERMINAL_PRODUCT_CRITERION.md`, blob
  `dafe266435d087a8bdc1d96ae8a8510a27aa5b8a`.
- `proofs/LEVEL1_FIVE_FACTORS_WITH_MULTIPLICITY.md`, blob
  `00f2db30a104409e1ccd3a7e4e7ac8071d197e67`.
- `proofs/LEVEL1_PRIME_POWER_TERMINAL_GRIDS.md`.
- `proofs/LEVEL1_CUBED_FOUR_AXIS_CORE.md`.
- `proofs/LEVEL1_UNIVERSAL_FIVE_PRIME_CORE.md`.
- `experiments/level1_grid_gap_certificate.py` and its six transitive
  Python dependencies, preserved in `research/rigidity_sources/experiments`.

The common-fiber source was fetched from its expressly cited version,
commit `0930c856256add528c967cd1b35c742f88b4cf95`:

- `proofs/DYADIC_FIBER_SIZES_2026_11_17.md`, blob
  `249065c67c2fd59e0e91207250d5fe120d398a65`.

## Verdicts

| Claim | Verdict | Exact scope |
|---|---|---|
| Nonnegative pure-prime-axis Fourier coefficients | PROVED | Every divisor ideal in a cyclic group, all exponents. |
| Common nonempty dyadic fibers agree for a homometric pair (DF1) | PROVED | Every odd part and every dyadic height. |
| Balance of each exact disputed size and each pure odd-axis marginal | PROVED | Consequences of the signed pure-axis Fourier values. |
| Two distinct odd axes with arbitrary exponents | PROVED | Every `2^e p^a q^b`, allowing omitted axes. |
| Three axes, one of prime order | PROVED | Every `2^e p^a q^b r`, with distinct odd primes, arbitrary `a,b>=1`. |
| Common nontrivial residual zero for a changing least-size class | PROVED | Every odd primary prime including 3, every primary exponent, any coprime odd residual group. |
| Four axes with only one nonsquarefree axis | PROVED | Every `2^e p^a q r s`, distinct odd primes, arbitrary `e>=0,a>=1`. No computational premise. |
| Terminal exclusion on `p^2 q^2 r s` | COMPUTED, with exact mathematical reduction | All distinct odd prime values, not a bounded numerical list. Full replay status is recorded separately below. |
| Every four-axis prime-power family | OPEN beyond stated families | Neither the preceding theorem nor the certificate establishes all exponent patterns. |
| Terminal feasibility implies an actual homometric pair | NOT CLAIMED | The terminal criterion used here is a necessary condition for a mate, sufficient only for exclusion. |

No genuine proof gap was found in the independently checked proofs of DF1,
DF2, or the `2^e p^a qrs` theorem. The complete reconstructed proofs are in
`audits/rigidity_section.tex` and use no inaccessible conversation premise.

## Detailed adversarial checks

### DF1 common fibers

At a fixed odd conductor, the partial dyadic difference and sum satisfy
`D_k S_k=0`. If one exact-size coefficient difference is odd, both
increments have 2-adic valuation exactly its size index. One incoming
zero component becomes nonzero at that valuation forever, so a second odd
coefficient difference is impossible. This does not assume nonnegativity
of mixed Fourier coefficients.

For a minimal shell with unequal nonempty fibers, choose the odd Fourier
row with kernel exponent one below each positive physical order exponent,
and the trivial row at zero coordinates. Modulo two the product row
selects precisely the corners of that shell. All smaller corners have
equal common fibers. The two unequal sizes therefore give two odd
coefficient differences at that one conductor, a contradiction.

The pure-axis positivity argument is separate and elementary: every
section along a prime-power coordinate is empty or a subgroup, whose
character sum is zero or its nonnegative order. It supplies exact signed
marginals and exact-size population balance, not positivity at general
mixed conductors.

### Common-zero lemma

The relevant Fourier row is indexed by kernel exponent:

`T_k(u)=sum_{i<=k} p^i (u_i-u_{i+1})`, with `u_{a+1}=0`.

When `j` is the first change in the signed sections, all difference rows
below `j` vanish. At each residual conductor outside their common zero
set, some row `k>=j` is nonzero, forcing its sum row to vanish. Therefore
the transform of `T_j(g)` vanishes modulo `p^(j+1)` off the common zero
set. Inversion is legitimate in `Z[1/2]` because `p` is odd and does not
divide the residual group order. This is where coprimality is required.

At a changing shell the possible residues are `p^j`, `-p^j`, and
`p^j-2p^b` for `0<=b<j`; none is zero modulo `p^(j+1)`. An empty common
zero set is impossible. If its only zero were trivial, the transform
argument would make `T_j(g)` constant modulo `p^(j+1)`. Balance and
nesting give an origin or top-shell residue of zero or `-2p^c`.
Valuation comparison excludes every match except `p=3`, `c=j`, and the
positive changing-shell residue. The proof does not simply discard 3.
In that exceptional situation, a larger dispute in any earlier section
would give residue zero or `-p^j`, contradicting the required constant.
The origin bounds every later fiber by the least size. Thus there are no
larger disputes, and the separately proved one-size argument excludes
the exception.

In the one-size argument, restriction to the subgroup of that dyadic
order is valid: the original sets agree on every exterior coset, and
shifts in the subgroup preserve those cosets. Removing their identical
autocorrelation contributions gives an actual homometric restricted pair.
This restriction is not asserted for a general multiple-size pair.

### Final four-axis argument

The pure `p`-axis marginals balance each three-prime residual section.
Aggregate pure residual-axis marginals equate the signed masses of each
singleton shell and its complementary doubleton. A sign cannot change
along the primary axis, so these aggregate equalities really give fixed
signs, not just a potentially cancelling total.

Comparability plus complementary signs then place every disputed section
in a single rank of the residual three-cube. A changing least-size class
has a common nontrivial residual zero. A single shell cannot vanish there;
three equally sized disputed shells cannot vanish there either, by the
explicit rank-one parity and rank-two reciprocal arguments. Thus every
least-size section has exactly two opposite-sign least-size shells and
one larger shell. The larger shell always has the globally repeated sign
of the three complementary pairs. Each least-size section therefore has
population of the same nonzero sign, contradicting global exact-size
balance. All inequalities use only that distinct odd prime predecessors
are distinct positive even integers, with the smallest two at least 2,4.

The reconstructed proof replaces the source's final Parseval invocation
by the already established exact-size population equality. This is an
explicit simplification, not an added assumption or a repair of a gap.

### Computational certificate audit

The grid program enumerates all downward-closed masks in a linear
extension, retains every two-sided difference, and quotients only by
overall sign. For 36 cells it uses balanced ternary keys. Uniqueness of
balanced ternary proves lossless decoding, and `3^36-1 < 2^63` proves
the encoding and offset fit signed 64-bit integers.

For increasing odd primes it substitutes
`p_i=3+2i+g_0+...+g_i` (zero-based indices), with nonnegative even gaps.
Every strictly increasing odd prime tuple has this representation.
Its one-sign/nonzero-constant and odd-constant tests are exact sufficient
tests for nonvanishing on this full unbounded domain; they may retain
impossible zeros but never reject a possible zero by approximation.

The cube test is valid even at nonmaximal disputed cells: on every
positive physical coordinate, a normalized adjacent kernel-row difference
has values +1 at that order, -1 at the next, and zero elsewhere. The
product functional is supported above the disputed shell, where the sum
profile is prescribed to equal the absolute difference. Its provably
nonzero value forces at least one actual Fourier zero in the conductor
cube. Treating all merely possible zeros as columns is conservative.

The incidence matrix is the tensor product of adjacent-difference
matrices. Integer-scaled null vectors give necessary polynomial identities.
A unique term of minimal exponent on an axis contradicts reduction
modulo that axis prime when its nonzero coefficient has absolute value
below the least possible prime. All other prime powers are units.

The factor test only discards a row when every nonconstant rational
factor is proved nonzero. Its reductions use nonzero powers of two,
other prime values, and smaller prime predecessors. The final largest-
prime test may additionally use smaller-prime successors; these are
strictly smaller than the largest odd prime. An unused formal variable
can be interpreted as an auxiliary larger odd prime, so it introduces
no unsupported arithmetic assumption.

Finite integer-array operations are bounded: each shell polynomial has
absolute coefficient sum at most `13^6`, and at most 36 shells contribute,
so every signed sum is bounded by `36*13^6 < 2^31`. The masks are
disjoint, so subtracting positive and negative byte tables does not double
the 36-shell bound. Nullspaces, factorization, and modular factor tests
use arbitrary-precision integers or exact rationals.

## Computational replay and reproduction

Runtime versions: Python 3.12.14; NumPy 2.3.5; SymPy 1.14.0. This exactly
matches the versions specified by the two-squared-axis source.

An independent check of the final exceptional pattern, without importing
the exhaustive verifier, passed in `audits/check_last_squared_pattern.py`.
It reconstructs both downsets from the two published signed masks,
verifies their difference exactly, expands the physical population
polynomial, and checks its residue `q^2*r*(p-1)*(p+1)` modulo the largest
prime. The structured result is `audits/check_last_squared_pattern.json`.
Its scope is one symbolic pattern, not the full exhaustive enumeration.

The smallest six-factor control `(4,1,1)` was separately run successfully:
5,706 signed patterns, 3,643 root survivors, 3,546 cube exclusions,
97 incidence exclusions, zero residual patterns, 0.314 seconds.

A complete replay of all 24 ordered grids **completed successfully** in
this session, with exit code zero and no residual patterns. Total wall
time was 1,039.802 seconds (17 minutes 19.8 seconds). Its authoritative
per-row counts, residual lists, runtime and `complete: true` flag are in
`audits/grid_replay.json`; the console record is `audits/grid_replay.log`.
A separate post-run check verified that the exclusions sum to the root
survivors in every one of the 24 rows, and that all six two-squared-axis
rows agree entry by entry with the source table. That validation is
recorded in `audits/grid_replay_validation.json`.

From a checkout pinned to the source commit, the unmodified command is:

```sh
python experiments/level1_grid_gap_certificate.py
```

The delivered standalone wrapper regenerates the same receipt schema,
captures the runtime versions and elapsed time, and optionally writes an
atomic per-row checkpoint:

```sh
python replay_grid_certificate.py --experiments experiments \
  --checkpoint grid_replay_checkpoint.json > grid_replay.json
```

It emits per-row progress on standard error and the final JSON document
on standard output. It asserts that all 24 shapes run, each residual list
is empty, and the exclusion counts balance. The wrapper is preserved as
`audits/replay_grid_certificate.py`; its help and syntax were checked.
The completed run used the same shape loop and the same unmodified
`certify` function through the original session wrapper. The additional
count-balance assertions were performed on the completed receipt, so the
published wrapper does not imply an unperformed second exhaustive run.

From the delivered archive (which must include the named Python files),
use its preserved `experiments` directory. Every printed residual list
must be empty; a nonempty list raises an exception. The six new ordered
shapes are `(1,1,2,2)`, `(1,2,1,2)`, `(1,2,2,1)`, `(2,1,1,2)`,
`(2,1,2,1)`, `(2,2,1,1)`. Each starts from 11,194,933 distinct signed
patterns. The other 18 rows cover the previously stated grids. Do not
misdescribe these as 24 tested numerical prime tuples: they are 24
exponent/order grids, each using symbolic exclusions for all prime values.

## Dependencies and exclusions

The analytic family theorem needs only cyclic Fourier inversion,
integer shell rows, pure-axis positivity, common fibers, marginal
balance, and the changing-axis common-zero lemma. It does not require
the five-factor terminal theorem, the terminal product theorem,
large-prime localization, or an enumerated grid.

The source's separate six-factor consequence invokes the product theorem
and the five-factor terminal exclusion. That consequence is outside the
analytic family proof. Its computational and mathematical dependencies
must remain visible if it is included elsewhere.

The terminal product theorem itself has a valid finite-dimensional
orthogonal-separation argument: infeasibility in each factor supplies
two supported dual functionals whose tensor product separates the product
prescription; feasibility in either factor supplies a tensor product
interpolant. This equivalence concerns terminal interpolation, not actual
binary divisor ideals. The more extensive localization and finiteness
claims in the linked sources were not fully audited here and should not
be promoted solely on the strength of this report.

## Highest-priority referee targets

1. Verify the common-fiber parity argument and its corner-selection row.
2. Verify that the primary common-zero lemma handles every `p=3` case.
3. Verify the passage from aggregate complementary-shell masses to fixed
   signs; its validity depends on common fibers and monotonicity.
4. Check every exclusion step in the exact grid program, independently
   reimplementing the enumeration or the linear incidence test if possible.
5. Keep the necessity-only terminal test separate from the constructive
   universal mate-classification theorems elsewhere in Paper A.

No unresolved analytic gap was found in the scope explicitly marked PROVED.
Remaining independent verification is needed for the computational
implementation and for any wider linked theorem not audited above.
