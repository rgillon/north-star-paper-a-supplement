# Reconciliation of the reproducibility objections

## Scope and input identity

This audit used a fresh extraction of the delivered archive
output/North_Star_Paper_A_V2_Review_Package.zip, SHA-256
88d8caa73419bfd461a6109ba2150dc086b1cd403296f3130a122091c99796a0.

The archive has 99 entries. All 98 files covered by its SHA256SUMS.txt
matched their recorded hashes immediately after extraction. The working
copy contains no Git checkout. The following statements concern this
specific delivered package, not an uninspected repository branch or a
possibly different referee download. The full referee attachment was not
available to this subtask; the pasted summary was the objection source.

## Findings

| Referee objection | Finding in the delivered package | Required action |
|---|---|---|
| Corollary 17.2's separate height-3360 audit does not exist | It exists as verification/h99_independent_audit.py, a combined H100/H99 changed-data audit. Its H100 section freshly reproduced 174 endpoint checks, 435 exact grid reconstructions, precision 3353 and height 3360. The manuscript calls it a “separate supplementary audit” without naming the file. | Name the exact script and output in Corollary 17.2; give its two-argument invocation. |
| H99 claims 425 targets, but the named verifier has 255 | Two different programs are involved. The original experiments/ppi_99_prime_pair.py visits three interval targets per level and completes exactly 255 reconstructions (verified by tracing the unmodified program). The separate verification/h99_independent_audit.py chooses five actual dyadic-grid targets per level and freshly reproduces exactly 425 reconstructions. For these fixed parameters the five comprise two endpoints and three strict interior points, so the word “interior” in the manuscript is inaccurate. | Attribute 425 explicitly to the separate script. Say “actual dyadic-grid targets” rather than leave the implementation unidentified. State 255 for the original program, 425 for the separate implementation, and replace “interior” by “actual dyadic-grid”. |
| Evidence path is wrong | The documented explicit command using experiments and evidence/h100_independent.json succeeds, and that evidence file is in the package. However, the separate script's **no-argument defaults are genuinely stale**: it looks for research/snapshot/experiments and for verification/h100_independent.json. Bare invocation fails. | Change defaults to ROOT.parent/'experiments' and ROOT.parent/'evidence/h100_independent.json'; retain compatibility with the documented positional arguments. |
| Use python3, not python | Correct: Section 11 contains python experiments/level1_grid_gap_certificate.py. The guide and README already consistently use python3. | Change that manuscript command to python3. |
| Cut-moments row is unverifiable | The package contains verification/verify_cut_moments.py and evidence/cut_moment_results.json. A fresh standard-library-only run reproduced exactly the archived JSON: 85 moduli, 1,567 ideals, 503 membership cuts, 2,248 transfer comparisons. Its separate positive control is radial at N=12, explicitly not a subgroup-union counterexample. | Name the script/output directly beside the table row. Make these files accessible wherever the supplement is hosted. Do not count these bounded regressions as the proof of the universal theorem. |
| Supplement deleted at branch tip | Not a local-archive finding. The 99-file archive is complete; branch-tip accessibility is a separate publication defect to check and fix. | Root agent is checking repository hosting and should ensure the revised package is available through stable, tested links. |

## Fresh runs

Executed from the fresh extraction's bundle root with assertions enabled:

~~~sh
python3 verification/verify_cut_moments.py
python3 experiments/ppi_h100_independent_audit.py experiments
python3 experiments/ppi_100_rank_precision.py
python3 experiments/ppi_99_prime_pair.py
python3 verification/h99_independent_audit.py experiments /absolute/path/to/fresh/h100_arithmetic.stdout
~~~

All five runs exited 0. The last command consumed the freshly regenerated
H100 donor report, not the archived report. The exact mathematical fields
H100_rank and H99 were identical to the archived changed-data report;
the fresh cut-moment result was identical to the archived JSON as an object.

A separate line-count trace of the unmodified original H99 program verified 255 candidate visits, zero nondyadic skips and 255 completed exact reconstructions. Exact Fraction checks of the separate audit parameters establish that its H99 sample consists of 170 endpoints and 255 strictly interior targets; the H100 sample consists of 174 endpoints and 261 strictly interior targets. These are stored in original_h99_target_counts.json and grid_target_categories.json.

The new H99 report gives the error-ratio decimal
0.3387485931919912; the script compares the underlying rational to
339/1000 < 10/27, rather than relying on that decimal.

The archive and evidence were not edited. Outputs and command receipts
are under revision_audit/repro/fresh/ and
revision_audit/repro/fresh_run_receipt.json.

## Distinguish the unrecovered source counts

The printed Appendix B regression row (85 moduli, 1,567 ideals, 503 membership cuts and 2,248 transfers) is reproducible. The separate provisional-outline count of 73 cut-moment comparisons was tied to LEVEL1_FACTOR_CUT_MOMENTS.md, which was not recovered and remains unverified. Nothing in this replay certifies that original formulation or its 73-count claim. The new 503-cut regression must not be represented as a replay of the missing original 73 comparisons.

## Wording ready for the revision

Corollary 17.2, final sentence:

> The separate implementation verification/h99_independent_audit.py,
> run with the included experiment data and a freshly checked H100 donor
> report, verifies all 174 changed endpoint inequalities and the maximum
> depth 3360, and reconstructs 435 actual dyadic-grid targets. Its recorded
> output is evidence/h99_independent_results.json, under H100_rank.

Theorem 17.3, finite-check paragraph:

> These checks comprise 170 endpoint inequalities over 85 levels.
> The separate implementation verification/h99_independent_audit.py
> imports no construction-verifier module and uses the previously checked
> H100 donor rationals as data. It reconstructs five actual dyadic-grid
> targets at each level (both endpoints and three interior points), 425 in total, and reproduces the changed precision
> and error bounds. Its output is
> evidence/h99_independent_results.json, under H99.
> The proof for every actual target is the interval, denominator and parity
> argument above; the finite targets are supplementary checks.

Appendix B, cut-moment row or following paragraph:

> Run python3 verification/verify_cut_moments.py from the bundle root.
> The JSON output reproduces evidence/cut_moment_results.json: 85 moduli,
> 1,567 ideals, 503 membership cuts and 2,248 transfer comparisons.
> These are bounded formula checks; Theorem 7.1 and Proposition 7.2 supply
> the general proof.

## Limits

These are fresh replays of the programs included in the preparation
package. They establish that the claimed separate artifacts exist and
that their documented outputs reproduce. They are not external refereeing
or an independent reimplementation of every mathematical formula.
The 24-grid certificate was not rerun here because the user reports a new
external complete replay and the local task concerned distinct artifact
and command objections.

## Tested default-path correction

A scratch copy containing only the proposed two-line default-path change was run without arguments. It exited 0 and produced a JSON object exactly equal to the archived h99_independent_results.json. The test source is verification/h99_independent_audit_fixed_paths.py in the extracted scratch copy; output is fresh/h99_fixed_defaults.stdout. The original delivered program and archive were left unchanged.
