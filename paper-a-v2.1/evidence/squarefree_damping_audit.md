# Squarefree damping ingredient audit

Date: 22 November 2026. Same-session collaborator audit, not external human or formal verification.

## Inspected source and scope

Inspected `research/snapshot/proofs/LEVEL1_SQUAREFREE_STAR_PAIR.md`, pinned by root's source snapshot at `20dd528d5fdd4bf7a16cc5ce05f9cd25ae155c25`. The source asserts a full squarefree pair with 196 odd primes and dyadic height 3214. That assertion depends on the weighted-fiber construction WF1, two partial upper extensions UP2, prime lists and rational parameters from `experiments/level1_squarefree_star_pair.py`, fresh spectator tags, a squarefree odd-part product B, decoder precision, and the exact all-rank window checks. This local audit does not establish those separate dependencies and does not promote the full construction.

## Defensible complete lemmas

`squarefree_damping.tex` supplies complete self-contained proofs of the following local method.

1. A balanced signed singleton antichain in a squarefree odd group is the difference of two actual divisor ideals of equal cardinality. Population balance is the exact equation sum epsilon_i(p_i-1)=0. The examples +3,+5,-7 and +5,+13,-7,-11 satisfy it.
2. Its downward and upward closures D,U satisfy DU=|delta|. The constant one-half annihilates delta under unnormalized convolution and approximates each of D,U,DU with sup error exactly one-half.
3. The local rank rho=1-D+U is increasing in physical divisibility order, has values 0,1,2, and preserves D,U,DU along comparable points of equal rank. The sum rank on independent blocks has the same all-choice property for the entire indicator tuple.
4. The tensor error-cancellation formula W=product b_i - product(b_i-q_i) is an exact annihilator whenever each delta_i*q_i=0. Expanding proves this at every character, not only at the trivial character. For half-error blocks the sup error is exactly 2^-t and W lies on a dyadic grid. The same argument multiplies a fixed extra approximation error eta by 2^-t.
5. Deleting s error factors changes the remaining pure product bound from 2^-t to 2^{-(t-s)}, with an explicit multiplicative bound for all replacement factors. This is the rigorous local meaning of 'a removed damping coordinate costs a factor of two.'

These claims were reconstructed algebraically from the definitions; no finite computation is substituted for their proofs.

## Important limits

- Equal population is not homometry. The displayed small seeds are ingredients, not Level 1 counterexamples.
- A tensor annihilator need not be positive, monotone or binary. Its binary section realization still needs a separate argument.
- A rank-dependent multiplication does not automatically preserve the annihilator identity. The appendix supplies only the exact rank compatibility and error bounds it proves.
- The fixed-core extension with eta requires its own annihilator; error size alone is insufficient.
- The group remains cyclic and squarefree under tensoring only when all participating odd primes are distinct across blocks.
- No arbitrary-excluded-prime existence statement, 196-prime pair, 230-prime minimality alternative, or full spectator realization is asserted in this appendix.

## Introductory scope check

Inspected `output/North_Star_Paper_A_V2/sections/01_introduction.tex`. Its description of the independent modulus procedure and its explicit source-recovery limitation are properly scoped. It correctly distinguishes exhaustive completeness from efficiency and an arithmetic characterization. One optional wording refinement is to replace 'two stronger source results' by 'two proposed source results': the original statements are incompletely recovered, so their precise comparative strength is not audited. No other overstatement of the fallback was found.

The broader arithmetic and numerical-construction claims in the introduction depend on other assigned audits. This report neither approves nor rejects them without those proofs.
