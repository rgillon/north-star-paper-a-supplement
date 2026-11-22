# North Star — Paper A, Version 2.1

**22 November 2026. Referee-revised audited-source edition.**

## Start here

- North_Star_Paper_A_V2_1.pdf: the complete manuscript.
- North_Star_Paper_A_V2_1.tex: a single editable LaTeX file with bibliography and finite defining arrays embedded. It builds without the modular section files.
- North_Star_Paper_A_V2_1_Referee_Support.pdf: proof dependencies, exact commands, adversarial targets, and scope of the evidence.
- North_Star_Paper_A_V2_1_Changes_and_Gaps.pdf: concise changes, deferred claims, and open questions.
- main.tex, sections/, appendices/, references.bib: modular manuscript sources.
- referee_support.tex, changes_and_gaps.tex, referee_response.tex: editable companion sources.
- North_Star_Paper_A_V2_1_Referee_Response.pdf: point-by-point response, including the conductor calculation and exact cancellation controls.
- experiments/: pinned construction and symbolic-certificate programs, fixed arrays, and all local import dependencies.
- verification/: separately written audit programs and a JSON-reporting wrapper for the grid certificate.
- evidence/: exact completed outputs, audit notes, source provenance, and computational validation records.
- source_notes/: recovered research inputs, included for comparison. These notes may contain broader claims than the manuscript certifies. The manuscript and its explicit status distinctions determine the scope of this edition.

## Mathematical scope

The complete Fourier-conductor theorem classifies all subgroup-union mates of a specified set. The printed rank theorem is on the explicit conductor grid. The separately proved coefficient-moment criterion decides a modulus by exhaustive enumeration and constructs witnesses.

The exact files LEVEL1_GENERATOR_LATTICE_RANK.md and LEVEL1_FACTOR_CUT_MOMENTS.md were not recovered from accessible source revisions or supplied materials. Their proposed compression, fixed-generator polynomial bit-time and avoidance-of-antichain-backtracking claims remain unverified. The results printed here are disclosed, independently proved formulations; they are not certification of those missing sources.

The manuscript includes analytic rigidity families, the explicitly COMPUTED two-squared-axis family, the original 100-prime construction, its height-3360 refinement, and a 99-prime construction at height 3408. The local squarefree damping method is retained without importing the full 196-prime construction.

This revision responds to the supplied Grok report and the subsequent “Holds with fixes” report summary. The included preparation audits are distinct from those reports. No complete Lean formalization is claimed.

## Referee corrections

Theorem 10.4 is retained with its original hypotheses. A nontrivial residual character has conductor c > 1; at every prime dividing c the nonzero-shell sum is -1. All-positive shell sums occur only at conductor 1, and Lemma 10.3 excludes a common zero set consisting only of that conductor. The new exact conductor controls check the two reported arithmetic cancellations, which vanish only at conductor 1. They are illustrative checks, not the proof of the infinite rigidity family.

The combined H100/H99 changed-data checker is verification/h99_independent_audit.py. Its default paths now resolve the packaged experiments and evidence/h100_independent.json. The explicit command below uses a freshly checked donor report. The manuscript now names both audit files and separates their target counts. The cut-moment row is reproduced by verification/verify_cut_moments.py and evidence/cut_moment_results.json; the older outline’s separate 73-comparison claim remains unverified.

## Reproduce the calculations

Run from this directory with assertions enabled. Do not use python -O or set PYTHONOPTIMIZE. The preparation used Python 3.12.14, NumPy 2.3.5 and SymPy 1.14.0. A conventional optional environment is:

~~~sh
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r requirements.txt
mkdir -p replay
~~~

These commands operate on included data and do not require GitHub access.

~~~sh
python3 verification/verify_sign_graph_independent.py 60 > replay/sign_graph.json
python3 verification/verify_cut_moments.py > replay/cut_moments.json
python3 verification/independent_primary_conductor_check.py > replay/independent_primary_conductor_check.json
python3 verification/check_primary_conductor_coordinates.py > replay/primary_conductor_coordinates.json

python3 experiments/ppi_100_prime_pair.py > replay/h100_original.log
python3 experiments/ppi_h100_independent_audit.py experiments > replay/h100_independent.json
python3 experiments/ppi_h100_decoder_audit.py experiments > replay/h100_decoder.json

python3 experiments/ppi_100_rank_precision.py > replay/h100_rank_precision.log
python3 experiments/ppi_99_prime_pair.py > replay/h99_original.log
python3 verification/h99_independent_audit.py experiments replay/h100_independent.json > replay/h99_independent_results.json

python3 verification/check_last_squared_pattern.py > replay/last_squared_pattern.json
python3 verification/replay_grid_certificate.py --experiments experiments --checkpoint replay/grid_checkpoint.json > replay/grid_replay.json
~~~

The last-pattern script also writes a JSON file beside itself. The original H100 construction may regenerate its fixed certificate. Run in a working copy to preserve every supplied byte. None of these documented entry points requires a local Git checkout; historical helper routines in some source modules are not needed.

For the original plain-text grid output instead of the JSON wrapper:

~~~sh
python3 experiments/level1_grid_gap_certificate.py > replay/grid_replay.log
~~~

The wrapper calls the same certify function on the same 24 ordered grids, checks that exclusions balance root survivors, and emits per-grid progress to standard error. A second full run through the original entry point is unnecessary when reproducing the wrapper result.

Expected key results:

| Check | Expected result |
|---|---|
| Sign graph | PASS; 60 moduli; 8,182 radial profiles; 465 subgroup unions; 980 complete direct time-domain profiles |
| Cut moments | PASS; 85 moduli; 1,567 ideals; 503 moments; 2,248 transfer comparisons |
| H100 arithmetic | PASS; 100 distinct odd primes; height 3446; 8,123 bits; exact error ratio below 2/7 |
| H100 decoder | PASS; 2,175 target reconstructions and all stated cover/height checks |
| H100 rank program | 87 levels; 261 reconstructions; maximum height 3360 |
| Original H99 program | 85 levels; 255 exact reconstructions; height 3408 |
| Changed-data audit | PASS; H100: 435 targets = 174 endpoints + 261 strict interior targets; H99: 425 targets = 170 endpoints + 255 strict interior targets; 99 odd primes, height 3408, 8,017 bits |
| Conductor controls | PASS; all eight conductors for each cited triple; zero-conductor list [1] in both cases |
| Terminal certificate | Complete=true; 24 rows; every residual list empty; six 36-cell rows each start with 11,194,933 signed patterns |

The complete grid replay took approximately 17 minutes 20 seconds in this preparation environment. Timing is descriptive, not a promised runtime. A finite regression is not the proof of a universal theorem. The printed proofs establish the graph converse and all-target decoder coverage.

Compare exact integers, rational numerators and denominators, maps, and completion conditions. Timing, Python build strings, and hashes of newly generated report bytes can differ. The source input hashes identify the fixed data.

## Build the documents

With a normal LaTeX installation containing the declared packages:

~~~sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error referee_support.tex
pdflatex -interaction=nonstopmode -halt-on-error referee_support.tex
pdflatex -interaction=nonstopmode -halt-on-error changes_and_gaps.tex
pdflatex -interaction=nonstopmode -halt-on-error referee_response.tex
pdflatex -interaction=nonstopmode -halt-on-error referee_response.tex
~~~

The single-file manuscript can instead be built with two passes of:

~~~sh
pdflatex -interaction=nonstopmode -halt-on-error North_Star_Paper_A_V2_1.tex
~~~

The single-file version needs neither BibTeX nor external section files. The modular sources are provided for convenient revision. The companion guide identifies claims by titles and manuscript numbers.

## Provenance and integrity

evidence/source_provenance.json records every bundled upstream program/data/note, its full commit, Git blob, delivered SHA-256 and byte count. Some text transfers appended one terminal newline; the provenance records this explicitly and verifies that removing it recovers the pinned Git blob. No executable statement is changed by that transfer detail.

Primary pins:

- Referee manuscript: 79350cb7bb6032df7637dcfd3e2bd029e6362a85
- Arithmetic and refinement snapshot: 20dd528d5fdd4bf7a16cc5ce05f9cd25ae155c25
- Common-fiber source: 0930c856256add528c967cd1b35c742f88b4cf95
- Original H100 lineage: 09d77232d96df1e16ef52917ddf46ccff90deb0b

SHA256SUMS.txt covers the delivered files except itself. A standard check is:

~~~sh
sha256sum -c SHA256SUMS.txt
~~~

The review package is self-contained. Links to historical repository revisions provide attribution and provenance, not missing principal proofs.

## Public edition

The versioned public directory is https://github.com/rgillon/north-star-paper-a-supplement/tree/main/paper-a-v2.1. Its top-level repository index links the current PDF, editable source, referee documents and complete archive. The older November 21 archive is retained as historical material. All verification files and their inputs are also directly browsable.
