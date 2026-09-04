# Changes from pRFv2

This document summarizes the pRFv3 work relative to the current `pRFv2`
directory in `gru`. It is written as an implementation and commit-review
record; [README.md](README.md) is the user-facing guide.

pRFv3 retains the original source attribution from pRFv2. The pRFv3
refactoring, new components, and documentation were developed by Austin Kuo
with Codex.

## Summary

pRFv3 keeps the established pRF analysis interface and standard-model fitting
objectives while reorganizing expensive work so it is performed once per scan,
worker, unique stimulus, or exact parameter vector instead of once per voxel.
It also adds bounded parallel execution, aggregate progress and QC reporting,
safer stimulus handling, command-line save conflict handling, faster ROI
graphics, and fixes for several analysis and MATLAB R2025b compatibility bugs.

The main compatibility target is numerically equivalent fit output for the
standard model paths. Gaussian construction and stimulus projection now use
mathematically equivalent reordered double-precision operations, so pRFv3 is
not bit-for-bit identical to pRFv2. `gaussian-DoG-CSS`, `divNorm`, and the
optional `fmincon` solver contain additional intentional numerical or
validation changes and should be regarded as experimental.

## Performance changes

### Prepared scan context

- Scan-invariant setup is assembled once in `pRFFit` and passed through
  `pRFFitPrepared` for direct voxel fitting.
- Repeated argument parsing, view queries, stimulus preparation, model setup,
  and prefit-bank construction are removed from the per-voxel path.
- Fixed HRF and stimulus-delta calculations are retained in the prepared
  context.
- Only the current numeric voxel-time-series block is supplied to the fitting
  loop. The prepared context is kept resident on parallel workers when MATLAB
  supports `parallel.pool.Constant`.

### Adaptive voxel blocks

- ROI time series are loaded in blocks instead of retaining the entire scan's
  voxel data at once.
- The default block size balances worker utilization against an estimated
  512-MiB working-set limit and scales with the resolved worker count.
- An internal `voxelBlockSize` override remains available for controlled
  testing and unusual memory constraints.

### Fast Gaussian construction and stimulus projection

- `pRFPrepareStimulusProjection` detects exactly identical stimulus frames and
  repeated frame runs.
- `pRFProjectStimulus` evaluates unique content once and expands the response
  back to the original time sequence.
- `pRFMakeGaussian` detects exact Cartesian `ndgrid` and `meshgrid` layouts and
  constructs each RF as the outer product of one-dimensional Gaussian
  components. Arbitrary grids and unsupported types retain the general 2-D
  expression.
- Dense double stimulus frames are reshaped into a pixel-by-frame matrix and
  projected through MATLAB's optimized matrix multiplication. This reorders
  floating-point operations relative to pRFv2 but remains double precision.
- Model files share the projection implementation instead of independently
  multiplying every repeated raster by every RF.
- All center-surround model paths use the same Gaussian helper for their
  secondary RF.

### Bounded exact caches

- Prefit model-response banks are cached by the exact stimulus and relevant
  model settings. The cache is bounded by entry count and memory.
- Repeated nonlinear-objective evaluations use exact parameter keys local to a
  fit. The memo is bounded and is destroyed when that fit completes.
- Worker-local projection results use exact identity/value keys and bounded
  storage.
- Cache telemetry remains available for debugging, but it is not saved as a
  primary analysis result.

### Worker management

- `numWorkers` is configured in the main pRF parameters rather than a separate
  worker-count dialog.
- `numWorkers = 1` is true serial execution: neither prefit nor voxel fitting
  inspects or starts a parallel pool.
- Values greater than one request exactly that many workers. An existing pool
  of the correct size is reused; a mismatched pool is replaced.
- Failure to start or initialize a pool emits a warning and falls back to the
  serial path.
- The same resolved setting controls prefit construction and voxel fitting.
- `mlrRunPRF` uses the same resolver and can request all locally available
  workers with `doParallel = -1`.

## Progress, verbose output, and QC

- Parallel workers report compact completion events to a client-side progress
  reporter rather than each streaming ordinary console output.
- Progress lines retain the original pRF presentation while adding completed
  percentage, elapsed time, and estimated time remaining.
- `progressEvery` controls reporting frequency; its automatic setting uses
  approximately one-percent intervals with a minimum of 500 voxels.
- `verbose` retains detailed per-voxel fit output but defaults off.
- `diagnosticVoxelCount` prints a reproducible stratified-random sample in the
  same format as verbose output without affecting optimizer randomness.
- End-of-run checks warn about identical sampled inputs, identical finite fits,
  or a lack of valid sampled results.

## Stimulus reconstruction and orientation

- A newly reconstructed stimulus is converted to the canonical pRF
  world-coordinate orientation inside `pRFGetStimImageFromStimfile` before it
  is returned or saved.
- The main fitting path no longer performs an implicit extra flip. With both
  flip checkboxes unchecked, it performs zero spatial flips on the helper's
  output.
- `xFlipStimulus` and `yFlipStimulus` remain as explicit corrections for a
  stimulus that was actually displayed in the wrong orientation.
- Manual helper output and internally generated output now follow the same
  orientation convention.
- Reconstruction supports stimulus types 1--5, handles the capitalization
  variants of `mglDoubleBars.m`, validates task indexing, and consistently
  includes type 5 in applicable preprocessing.
- `recomputeStimImage` bypasses an embedded image; `saveStimImage` saves the
  reconstructed canonical structure back as `pRFStimImage`.

### Migration warning

pRFv2 embedded `pRFStimImage` before its later implicit orientation correction
and did not save a convention-version marker. An image embedded by pRFv2 must
therefore be reconstructed and resaved once before it is reused by pRFv3. Use
`recomputeStimImage` with `saveStimImage`, or use `saveStimImageOnly`.

## Averaged-scan stimulus reuse

- A multi-constituent scan is identified generically from `averageTSeries`
  provenance (`params.scanList`), not from a subject, session, group name,
  filename, or hard-coded scan number.
- Existing embedded stimuli are reusable only when `im`, `x`, `y`, and `t`
  match exactly.
- Forced reconstruction compares the supported task program, full geometry,
  task/stimulus state, randomized variables, derived trial-volume sequence,
  `volTrigRatio`, and screen geometry before taking the one-construction path.
- If all constituents match, pRFv3 reconstructs once and, when requested,
  writes that same canonical structure to every constituent stimfile.
- Missing, unsupported, or unequal metadata falls back to reconstructing each
  constituent and applying the legacy exact-image comparison.
- Shifted/reversed averages, nested averages, and averages of concatenations
  remain unsupported.

## Stimulus-only workflow

- The parameter GUI includes a read-only `pRFStimImageExists` checkbox for
  `dispStimScan`. It checks MAT-file metadata only and is checked only when
  every linked constituent stimfile already contains the saved variable.
  Changing the selected group or display scan refreshes the indicator.
- `saveStimImageOnly` reconstructs and saves the selected stimulus data, then
  exits before voxel loading, pool creation, fitting, or analysis creation.
- It forces recomputation and saving for that one action.
- When parameters are obtained with
  `[v, params] = pRF(v, [], 'justGetParams=1')`, the one-shot save fields are
  cleared after success so reusing `params` performs a normal fit rather than
  saving the stimulus a second time.

## Analysis saving and continuation

- Analysis-name conflicts use a Command Window prompt rather than a modal
  merge/overwrite/copy dialog.
- The temporary mrTools overwrite policy is restored after saving.
- Continuing an analysis uses the analysis name selected in the GUI instead of
  a stale embedded `params.saveName`; for example, continuing `pRF_test` no
  longer merges into `pRF`.
- Save and merge handling is centralized in `pRFSaveAnalysis`.

## Plotting and interaction

- Aggregate ROI outlines are batched into a small number of patch objects and
  use lower-resolution circles while retaining individual r2 shading and the
  drawing order that places higher-r2 outlines on top.
- RF-center and eccentricity/size point loops are replaced by scatter objects.
- Both aggregate spatial plots use equal x/y physical scaling.
- The first title says `pRF coverages`; ROI titles show the plotted voxel count
  before the r2 cutoff, and underscores are rendered literally.
- Shift-click preserves the original single-voxel behavior and suppresses the
  aggregate ROI plot.
- Clicking a voxel without a saved fit asks for confirmation before stimulus
  reconstruction and fitting.
- A scoped compatibility guard handles modifier-key events that legacy
  `getpts` cannot parse under MATLAB R2025b.
- The interactive stimulus preview clamps its initial and mouse-selected time
  points to the available fit/stimulus range.

## Model and solver changes

- Standard model files call the shared stimulus-projection implementation.
- `gaussian-DoG-CSS` and `divNorm` remain available, with explicit bounded-
  solver recommendations and validation for widths, signs, and nonlinear
  constants.
- Their neural nonlinearities are applied before HRF convolution. This is an
  intentional correction and can change results relative to pRFv2.
- The `gaussian-divs` maximum-bound typo is corrected.
- `fmincon` with SQP is available as an optional experimental solver when the
  Optimization Toolbox is installed.
- `nelder-mead`, `levenberg-marquardt`, and `nelder-mead-bnd` remain available.

## Additional correctness and robustness fixes

- The junk-frame check uses the actual number of entries rather than the
  length of a logical comparison.
- Canonical-model recovery searches successful finite voxel fits and uses the
  same scan, coordinates, parameters, and prepared stimulus as the analysis.
- Unsupported restrictions and model/canonical selections now raise explicit
  errors instead of entering debug breakpoints or failing later.
- Older saved parameter structures receive defaults for new pRFv3 fields.
- Scan output assignment and selected-analysis naming are kept consistent
  across fitting and merging.
- GUI help, callback labels, model headers, and stale helper comments were
  corrected during the cleanup pass.

## Plugin activation

- `pRFv3Plugin` exposes pRFv3 as a distinct `mlrPlugin` choice and activates
  its runtime directories at the front of the path.
- If pRFv2 and pRFv3 are both selected through the plugin callback, pRFv3
  removes pRFv2 from `selectedPlugins` because the versions export the same
  public function names.
- `pRFv3ActivatePath` discovers the plugin root from any location inside the
  tree, places the complete current pRFv3 tree first, refreshes MATLAB's path
  cache, and verifies all cross-file runtime names resolve inside pRFv3.
- Exact helper subdirectories are not encoded in the verification manifest;
  moved functions are discovered from the current tree and duplicate names
  inside pRFv3 are rejected.
- `pRFv3Startup` provides an opt-in startup hook for scripts that use pRF
  without first opening an mrLoadRet view.
- Direct `pRF`/`pRFGUI` use refreshes the getpts compatibility installation,
  replacing listener and keypress callbacks retained from an older helper
  version or file location before the parameter dialog opens.

## Removed legacy components

pRFv3 intentionally omits the following pRFv2 components:

- The obsolete controller and split-analysis stack (`pRFController`,
  `pRFSplit`, `pRFRunSplits`, and `pRFMergeSplits`).
- The `machines` profiles, including Sherlock, GRU-lab, and local batch-job
  helpers.
- The passwordless-SSH test helper.
- The unused `fminsearchcon` copy and model template.
- The unimplemented cross-validation GUI option.

The cleanup audit found no remaining whole `.m` file that is clearly orphaned.
Developer-facing cache statistics and progress summaries were retained as
diagnostic interfaces.

## Validation performed during development

- The legacy arithmetic and optimized kernels were compared on two independent
  deterministic samples of 100 left-hemisphere voxels from the s0129 session.
  Nonlinear voxel fitting was 3.37x and 4.00x faster, while context/prefit
  preparation in the newer comparison was 4.57x faster.
- In the newer `pRF_lV1` comparison, all 100 fits had the same valid/NaN and
  optimizer-exit status. Ninety-four parameter vectors were bit-identical;
  the median parameter and absolute r2 differences were zero. The maximum
  parameter difference was `0.00863` degrees (5.7% of one stimulus pixel), the
  maximum absolute r2 difference was `6.47e-5`, and the worst prediction
  correlation was `0.99999127`.
- Neither 100-voxel comparison changed pass/fail classification at r2 cutoffs
  0.05, 0.10, 0.20, 0.30, 0.50, or 0.70. The nonidentical cases were flat,
  sub-pixel or low-r2 fits rather than meaningfully different predictions.
- Direct tests on the actual 192-by-108 stimulus found RF differences no larger
  than `2.22e-16` and projected-response differences on the order of `1e-13`.
- Rectilinear `ndgrid` and `meshgrid` layouts, duplicate-frame prepared banks,
  ordinary movies, output-frame padding, all secondary-Gaussian models, and
  general-path fallbacks were exercised.
- The averaged-stimulus fast path was exercised with matching constituent
  stimfiles: reconstruction occurred once and the saved payloads matched
  exactly.
- A deliberately changed constituent exercised the fallback path and produced
  independently reconstructed, distinct stimuli.
- MATLAB Code Analyzer was run across pRFv3; the cleanup introduced no parse
  errors. Remaining reports are primarily inherited style and compatibility
  warnings.
- Path smoke tests from a neutral working directory resolved shared pRF entry
  points to pRFv3.

These checks are targeted regression tests, not a substitute for validating a
representative full-session analysis before release.

## Follow-up items to review before upstreaming

These were identified by the final cleanup audit and intentionally left
unchanged because they affect behavior or merit a separate focused change:

1. Validate averaged-scan `shiftList` and `reverseList` before any
   `saveStimImage` write. The current reuse shortcut can reach its save call
   before the later unsupported-transform check.
2. Consider adding an orientation/convention marker to newly saved
   `pRFStimImage` structures so unversioned pRFv2 payloads can be detected
   automatically rather than handled only by migration documentation.
3. Make `pRFv3Startup` persist the same pRFv2/pRFv3 mutual-exclusion decision
   as `pRFv3Plugin`, preferably through one shared helper.
4. Correct inherited `mlrRunPRF` issues in a separate patch: three error paths
   call `results.errorString(...)` rather than assign it, stimulus height is
   read from the first dimension instead of the second, and temporary
   directories are removed through unquoted shell `rm -rf` commands.
5. Decide whether the plugin callback's uninstall branch should remain. The
   current mrTools `mlrPlugin` workflow does not appear to call it, and its
   activation-before-removal sequence is misleading.

## Suggested commit message

```text
Add pRFv3 with accelerated fitting and updated workflows

- prepare scan state once and reuse stimulus/prefit/objective work
- use separable Gaussian construction and matrix stimulus projection
- add bounded worker, memory, progress, and diagnostic-voxel handling
- make reconstructed stimulus orientation and saved images consistent
- reuse verified equivalent stimuli across averaged scans
- improve analysis continuation, ROI plotting, and R2025b interaction
- retain experimental DoG-CSS/divNorm models with clearer validation
- remove obsolete split, machine, and Sherlock infrastructure
```
