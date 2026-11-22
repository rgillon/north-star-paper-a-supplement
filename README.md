# North Star — Paper A

## Current edition: Version 2.1, 22 November 2026

**Reconstructing Unions of Subgroups from Cyclic Autocorrelation: Exact mate classification and arithmetic obstructions**, Robert Gillon.

This referee revision makes the conductor convention in Theorem 10.4 explicit, corrects the verification paths and target-count descriptions, and includes all separately written audit programs. The theorem retains its original hypotheses; the nontrivial-character condition already follows from Lemma 10.3.

| Deliverable | File |
|---|---|
| Complete manuscript, 50 pages | [PDF](paper-a-v2.1/North_Star_Paper_A_V2_1.pdf) |
| Editable single-file manuscript | [LaTeX](paper-a-v2.1/North_Star_Paper_A_V2_1.tex) |
| Proof dependencies and exact verification commands, 10 pages | [Referee support](paper-a-v2.1/North_Star_Paper_A_V2_1_Referee_Support.pdf) |
| Response to supplied referee comments, 4 pages | [Referee response](paper-a-v2.1/North_Star_Paper_A_V2_1_Referee_Response.pdf) |
| Changes, deferred claims and open questions, 4 pages | [Changes and verification status](paper-a-v2.1/North_Star_Paper_A_V2_1_Changes_and_Gaps.pdf) |
| Complete source, evidence and verification package | [Download ZIP](North_Star_Paper_A_V2_1_Review_Package.zip) |

The [versioned directory](paper-a-v2.1/) also exposes the modular LaTeX, bibliography, verifier sources, fixed data and recorded outputs for direct inspection. The package is self-contained and requires no access to the private research repository or research conversations.

## Mathematical scope

For a specified subgroup union, the Fourier-conductor graph classifies every subgroup-union mate. Its conductor-grid rank formulation is exact. A separate exhaustive coefficient-moment procedure decides a specified modulus and constructs witnesses. This is a complete algorithmic classification at those scopes; a concise arithmetic classification of all rigid moduli remains open.

The primary-axis theorem is analytic and uses no enumeration. The two-squared-axis family remains explicitly **COMPUTED**, with its complete symbolic-grid certificate. The H100 and H99 constructions are proved with identified fixed-data exact certificates. The two stronger generator-compressed source formulations named in the preliminary outline were not recovered and are not certified by this edition.

## Verify

From the extracted archive root, follow [README.md](paper-a-v2.1/README.md) and the referee guide. Assertions must remain enabled. All documented commands use `python3`. The package specifies Python 3.12.14, NumPy 2.3.5 and SymPy 1.14.0; construction comparisons use exact integers and rationals.

Key independently written checks are [cut moments](paper-a-v2.1/verification/verify_cut_moments.py), [changed H100/H99 data](paper-a-v2.1/verification/h99_independent_audit.py), and [literal conductor controls](paper-a-v2.1/verification/independent_primary_conductor_check.py). The full [evidence directory](paper-a-v2.1/evidence/) and [checksum manifest](paper-a-v2.1/SHA256SUMS.txt) identify the supplied results and exact bytes. The [revision manifest](paper-a-v2.1/evidence/referee_revision/revision_manifest.json) records the checked scope.

The original H99 program performs 255 reconstructions. The separate changed-data audit performs 425, comprising 170 endpoints and 255 strict interior targets. Its H100 refinement performs 435, comprising 174 endpoints and 261 strict interior targets. The current cut-moment regression covers 85 moduli, 1,567 ideals, 503 cuts and 2,248 transfer comparisons; these counts do not certify the missing outline's different 73-comparison claim.

## Historical edition

The earlier [November 21 public supplement](North_Star_Paper_A_Public_Supplement_2026-11-21.zip) is retained for provenance. Use the Version 2.1 package above for the current manuscript and verification instructions.
