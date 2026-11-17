# Core mutation audit

The source suite passes 1,224 assertions in 131 test blocks, with no failures,
errors, or warnings. Validation used R 4.6.1 and muttest 0.3.0. The first
mutation table records the earlier core-decision audit; the upstream extension
is reported separately below.

## Corrections

- Relative cut-site overlap now takes its origin before buffer expansion.
- Primer-overlap filtering recognizes relative events starting at zero.
- Normalization counts each control read/event once across paired strands.
  Aggregated inputs without read IDs remain supported. The frequency comparison
  remains strictly `> min_freq`.
- Consensus compares the original clipped-read sets before updating either set,
  so two clipped reads do not cause strand-biased disagreement resolution.
- Off-target features and safety-cap weights now aggregate once per read.

Each correction has a regression that failed on the previous implementation.
The paired FASTQ fixture checks exact events and counts for controls,
normalization, primer-dimer removal, both HDR modes, threshold equality, and fresh/resumed execution.
Cached filtered/normalized event tables and summaries require regeneration.

## Mutation results

K means assertion failures without test errors; S means survived; E means test
errors without assertion failures; M means both. Errors are not assertion kills.
The original tests ran against an isolated pre-change copy. Final survivors were
also checked against the existing regression files.

| Function | Before K/S/E/M | After K/S/E/M |
| --- | --- | --- |
| `amplicanFilter` | 1/8/1/0 | 7/2/1/0 |
| `findEOP` | 26/11/6/0 | 34/3/6/0 |
| `findPD` | 1/4/0/0 | 5/0/0/0 |
| `findLQR` | 21/14/0/0 | 34/1/0/2 |
| `goodBaseQuality` | 2/7/2/0 | 4/2/5/0 |
| `goodAvgQuality` | 2/7/2/0 | 4/2/5/0 |
| `alphabetQuality` | 0/6/2/0 | 1/2/5/0 |
| `amplicanNormalize` | 5/9/9/0 | 10/3/10/0 |
| `amplicanOverlap` | 5/10/5/0 | 8/2/10/3 |
| `amplicanConsensus` | 24/19/7/0 | 37/2/8/3 |
| `getHits` | 3/2/0/0 | 3/2/0/0 |
| `is_hdr` | 20/8/8/2 | 25/2/9/2 |
| `is_hdr_strict` | 22/8/5/4 | 27/5/5/2 |

Before: 298 mutants, 132 K, 113 S, 47 E, 6 M.
After: 303 mutants, 199 K, 28 S, 64 E, 12 M.
Corrected code introduced five additional mutation locations, so these are not
identical mutant sets. These counts do not estimate a probability of correctness.
All 11 explicit fault replays failed assertions, including the earlier extraction
and frameshift replays and all five corrections above.

The high-editing fixture still requires zero removals. With read-level aggregation
it can bypass clustering; any safety-cap warning is validated if emitted.
Controlled boundary cases require the exact cap warning.

The 28 remaining survivors are reviewed individually in
[mutation-review.csv](mutation-review.csv). They change batch partitioning,
permitted hit ordering, empty-set work, ignored arguments, or mathematically
unchanged comparisons/aggregations. Their recorded input-domain assumptions are
part of the review; the mutants are not silently excluded from the counts.

Raw reports, source snapshots, per-mutant timings, regression rechecks, and
session information are under `/tmp/amplican-core-audit/`, including `before.csv`,
`after.csv`, and `known-faults.csv`. The final classification is in
`after.csv`'s `final_status`; original counters and separate regression rechecks
remain available. These temporary files may be cleared by an
environment reset; this summary and the survivor review remain in the repository.

## Repeating the audit

Build/install the checkout first so `src/` contains its compiled library. Install
`pkgload`, `testthat >= 3.2.0`, and `muttest` as development dependencies. From the
package root:

```sh
Rscript tools/mutation-pilot.R /path/to/reports
Rscript tools/mutation-pilot.R /path/to/reports amplicanNormalize,amplicanFilter
```

Without target selection the runner includes the two original pilot functions
and all 13 additional functions. It applies one mutation in memory, restores
namespace/export bindings afterward, validates the unmodified baseline, and
rechecks survivors against existing regressions. CSV reports retain failures,
errors, warnings, timings, and the actual source changes. No score threshold is
used. `muttest` is not a runtime dependency of amplican.

The assessment covers core event decisions rather than primer assignment,
alignment optimization, or statistical calibration of the clustering heuristic.
`findEOP()` still infers relative coordinates from a non-positive start; use
absolute coordinates for positive-only event tables whose coordinate mode would
otherwise be ambiguous.

## Upstream audit

This extension covers FASTQ assignment, orientation and coordinates,
alignment-to-event metadata, duplicate counts, and pipeline execution.
It preserves assignment to every matching amplicon, selection of the best primer
alignment, and resume only with unchanged inputs and parameters.

Two implementation defects were reproduced and corrected:

- Reverse-only extraction supplied a missing forward-primer position to
  `getEventInfo()`. It now uses the same origin of 1 as an absent position column.
- In-memory alignment assembly converted entirely absent unassigned-read tables
  into a data frame without the required Barcode column. It now retains NULL for
  that empty result. Fully filtered and empty FASTQ inputs return valid objects.

Primer-matching behavior is unchanged. Its documentation now describes the best
overlap alignment rather than promising the earliest valid match.

The fixtures specify expected coordinates, replacement bases, weighted counts,
HDR flags and summaries. They distinguish barcode-level unique sequence pairs
from weighted per-experiment counts, which may overlap across amplicons.
They cover both directions, all four FASTQ modes, primer mismatch/overlap
boundaries, an N-containing read, truncated primers, differing reverse mates,
quality failures on either mate, empty experiments, and paired seeded sampling.

Pipeline checks compare batch sizes 1, 7 and 1000, reordered reads/configuration,
fresh execution, complete and partial resume, and the legacy alignment-object
path. They verify both exported normalization defaults, strict-consensus exclusion,
and off-target filter wiring. Only the wiring check uses a filter stub; it does
not claim to calibrate the clustering algorithm. The parallel check requires two
distinct worker PIDs and compares alignment, extraction and final numerical results.

Generated mutations target `locate_pr_start`, `get_seq`, `get_left_primer`,
`get_right_primer`, `upperGroups`, `flipRanges`, `amplicanMap`, `getEventInfo`, and
the actual `extractEvents` method body. Larger functions `makeAlignment`,
`amplicanAlign`, and `amplicanPipe` use explicit fault seeds. The runner restores
temporary namespace/export bindings and recreates both exported pipeline closures
when testing their factory. It records execution counts and verifies changes to
function bodies or argument defaults before running a mutant.

The baseline snapshot includes the previously accepted, uncommitted fixes.
Existing plotting tests are included where they exercise these helpers; no plot
implementation changed. Raw baseline records, legacy-test replays, final reports,
source snapshots and session information are under `/tmp/amplican-upstream-audit/`.
Those temporary artifacts may be cleared by an environment reset.

### Upstream mutation results

| Function | Before K/S/E/M | After K/S/E/M |
| --- | --- | --- |
| `locate_pr_start` | 8/4/2/0 | 9/0/4/1 |
| `get_seq` | 2/4/15/0 | 7/0/9/5 |
| `get_left_primer` | 0/6/3/0 | 6/0/3/0 |
| `get_right_primer` | 0/6/3/0 | 6/0/3/0 |
| `upperGroups` | 0/0/4/1 | 0/0/4/1 |
| `flipRanges` | 5/20/10/0 | 21/0/8/6 |
| `amplicanMap` | 4/7/2/0 | 5/2/6/0 |
| `getEventInfo` | 10/3/5/1 | 8/3/7/1 |
| `extractEvents` | 12/18/18/0 | 23/6/20/0 |

Before: 173 mutants, 41 K, 68 S, 62 E, 2 M.
After: 174 mutants, 85 K, 11 S, 64 E, 14 M.
The offset fix adds one mutation location, so these sets are not identical.
Before counts include 74 replays against the existing plotting/helper tests;
`before.csv` preserves both initial and replay statuses. The development batches
wrote their complete CSVs before a trailing parser error caused by an in-flight
runner edit. The final current-code audit ran from a frozen runner and exited
successfully. Baseline suite validation and legacy replays also completed.

All 14 explicit upstream fault seeds fail assertions without test errors.
The configuration-order seed initially survived, then failed after adding
explicit object-order assertions. `after.csv` retains its initial result and
recheck; `seeded-faults.csv` lists all final seed outcomes. New artificial faults
are labeled `seeded_fault`, separately from discovered implementation defects.

The 11 survivors are reviewed in `mutation-review.csv`: two ignored ID arguments,
two positive-reference-length comparisons, one scalar-offset validation guard,
four redundant strand-count guards, and two internal-container/name changes.
The vector-offset validation guard is contract-limited rather than universally
equivalent. Error-only mutations remain separate from assertion kills.
Fast cases run first; only survivors continue to existing regressions. An earlier
test error can therefore move a case from K to E. This is not a full-suite
mutation score.

The final suite passes 1,224 assertions in 131 blocks, with zero failures,
errors or unexpected warnings. The parallel check verified two distinct worker
PIDs. Static guards also verified that all 538 registered generated mutations
across the complete runner change executable bodies or argument defaults; only
174 new generated mutations were executed in this extension.

Repeat the extension from the package root with the existing development
dependencies and a local worker environment that permits communication:

```sh
Rscript tools/mutation-pilot.R /tmp/upstream-helpers locate_pr_start,get_seq,get_left_primer,get_right_primer,upperGroups,flipRanges,amplicanMap,getEventInfo,extractEvents
Rscript tools/mutation-pilot.R /tmp/upstream-seeds makeAlignment,amplicanAlign,amplicanPipe
```

This audit does not establish biological accuracy on independent real data,
calibrate off-target discrimination, or validate the external alignment library.
The main pipeline uses scalar amplicon offsets; mixed-validity vector-offset
validation remains outside this assessment. Changed inputs or parameters require
a fresh results folder. A complete pipeline run with no extracted events still
stops by design; empty experiments within a populated run are covered. Existing unused `uBarcode` and `config_order` locals in
`amplicanAlign()` were left in place.
