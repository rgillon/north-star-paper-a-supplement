# Universal claims audit and defensible replacement

Date: 22 November 2026. Actor: same-session AI collaborator `/root/universal_sources`, shared parent context. This is not independent human refereeing or proof-assistant verification.

## Verdict on the two named source claims

**PENDING SOURCE RECOVERY**, not REFUTED and not VERIFIED:

- `LEVEL1_GENERATOR_LATTICE_RANK.md`
- `LEVEL1_FACTOR_CUT_MOMENTS.md`

The attached four-page Version 2 explicitly says that these complete sources were unavailable to its authoring environment. It supplies only proposed claims and proof obligations. We did not recover the complete files or a precise formula for their claimed matrix/polynomial. Accordingly no verdict on the correctness of their actual proofs is possible, and the claimed 191/2741/73 verification counts have not been reproduced.

### Recovery performed

- Root inspected main and multiple recent branch trees, including `20dd528d5fdd4bf7a16cc5ce05f9cd25ae155c25`, `2ddc150f363da5ae410ed4ac1de573e7cbd7bb80`, `a0f857e436c29eabddd5d3ed312d6886f510dadf`, `e1592ec11eb750de4bdc1a764812be59e13aa653`, and `f804d4f1fc26a995cf640477fa9b79e7d71206f9`. The requested filenames were absent from those inspected trees.
- This collaborator queried all returned matching branch refs, the complete returned recent PR inventory, and 100 commits on the spectator-obstructions branch. Its head was `20dd528d5fdd4bf7a16cc5ce05f9cd25ae155c25`; its history names many new arithmetic proofs but no generator-lattice-rank or factor-cut-moments source. PR #49 is the latest returned PR, with head `79350cb7bb6032df7637dcfd3e2bd029e6362a85`.
- GitHub commit searches for `generator lattice rank`, `factor cut moments`, `moments`, and `lattice` did not recover either theorem. Search nonresults are not proof that a file has never existed.
- Exact-name Library searches returned no sources. Broader equivalent searches recovered the same provisional Version 2 and historical manuscripts, but no complete new source. A Library inventory of files created since 21 November showed no separate Version 2 source archive, rank proof, or cut-moment proof.
- Personal Context returned unavailable. No claims about absent conversations are made.

The correct manuscript treatment is to state this recovery limitation explicitly, omit claims specific to the unavailable formulations, and present only complete proofs supplied in the finished manuscript. Source absence must not be called a genuine mathematical gap or counterexample.

## Independent universal fallback

`universal_fallback.tex` provides complete source-independent proofs of:

1. A cut coefficient criterion for modulus-level rigidity. For every physical order shell d, partition all divisor ideals by membership of d, encode complete squared Fourier signatures in nonnegative-coefficient polynomials, and take the coefficientwise inner product of the two polynomials. Positivity is equivalent to an actual binary, downward-closed homometric pair separated by that cut. All cuts zero is equivalent to rigidity of the modulus.
2. An exact prime-factor transfer recurrence generating those polynomials from nested ideal slices. It retains all mixed Fourier measurements and postpones squaring until all linear section contributions have been added. It includes constructive witness extraction and an explicit exponential complexity bound in tau(N).
3. A graph-incidence rank corollary conditional only on the complete sign-graph theorem proved elsewhere in the paper. The rational nullspace is spanned by indicators of unanchored components; that special basis proves the 2^nullity count.

These are **not** presented as recovered proofs or exact restatements of the unavailable source theorems. The incidence rank matrix is not asserted to be the missing compact generator-lattice matrix. The transfer recurrence is not asserted to avoid all implicit ideal enumeration or antichain-state growth. No polynomial bit-time claim for a fixed number of generators is made.

## Adversarial proof checks

- The sum-of-squares expression is a *difference* of squared coefficient norms. Nonnegativity is justified by its independent expansion as four times a sum of products of nonnegative counts. A difference of arbitrary squared norms need not be nonnegative.
- The cut partition is by physical element order. Subgroup order/index conventions must not be silently interchanged.
- Complete squared Fourier signatures, not populations or only pure-prime measurements, define the polynomial exponents. Thus equality includes every mixed Fourier equation.
- The coefficient-square operation is additive and not multiplicative. Applying it to a factor before the transfer completes would be incorrect because it would omit cross terms.
- The map from nested flags to ideals is bijective. The recurrence proves both directions and counts every flag once; it imposes binary membership and downward closure before spectral comparison.
- Positive moments retain witnesses from both cut classes. They cannot select the same ideal because the cut predicates differ.
- The constructive inverse formula is used only on recorded full Fourier vectors generated by actual ideals. An arbitrary square-root/sign choice is not asserted to have a binary inverse.
- N=1, M=1, the empty ideal, the full ideal, zero Fourier entries, and arbitrary prime powers are included.
- The complexity statement is exponential in tau(N), with factorization supplied. It is not a polynomial-time arithmetic classification of all N.
- The rank corollary requires an incidence matrix with anchor rows. Rational nullity alone does not count binary vectors for general matrices.

## Recommended final wording

“The two proposed universal source notes named in the working draft could not be recovered at the inspected repository revisions. This version therefore does not certify their particular formulations. It proves a complete sign-graph classification for individual targets and a separate explicit coefficient-moment decision criterion for each modulus. The latter includes an exact transfer recurrence and constructive witnesses, with exponential worst-case complexity. A concise arithmetic characterization of rigid moduli remains open.”

## Remaining independent checks

The mathematical proofs of the independent replacement are complete at the written scope. They should receive external adversarial review. The foundational sign-graph theorem is audited by another same-session collaborator; this report does not certify that dependency by repetition. The missing named source notes and their reported regression counts remain unverified until their actual source bytes are supplied.

## Replayed exact controls

Command (Python standard library only, assertions enabled):

```sh
python3 verify_cut_moments.py > cut_moment_results.json
```

The delivered script and saved output cover exactly `N=1,...,72` and `84,90,96,108,120,126,132,140,144,150,168,180,210`: **85 moduli, 1,567 actual divisor ideals, 503 cut moments, and 2,248 transfer-recurrence comparisons**, all PASS. Every factor choice p^a exactly dividing each modulus, every order-shell cut and both cut sides are checked. The recurrence's full unsquared spectral dictionaries and recovered ideal witnesses are compared against a separate recursive downward-ideal enumeration. Literal residue-pair counting independently produces C_A at every shift, and its exact Ramanujan transform is compared entrywise with the squared signature. The complete Fourier inverse recovers every original order-shell indicator. The cut counts and difference-of-squared-norm identities are checked against literal-autocorrelation collision classes.

All subgroup-union moments in this bounded control are zero. To exercise positivity and witness extraction, a separate explicitly labeled positive control uses **all 64 radial binary subsets of Z_12**, a larger class. Its moments in order-shell order `(1,2,3,4,6,12)` are `(24,24,22,8,22,8)`. One extracted pair has order sets `{1}` and `{2}`, corresponding to `{0}` and `{6}`; the latter is not a subgroup union. This is an implementation control only and is never claimed to be a Level 1 counterexample.

These are new counts for the source-independent fallback. They do not reproduce or validate the missing source's reported 73 comparisons.
