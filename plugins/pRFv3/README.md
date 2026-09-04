# pRFv3

pRFv3 is a refactor of the mrTools pRFv2 plugin focused on faster fitting,
bounded memory use, clearer progress reporting, and safer stimulus and analysis
workflows. The standard pRF models continue to use the established fitting
objectives and double-precision arithmetic. Gaussian construction and stimulus
projection use mathematically equivalent reordered operations, so results are
not expected to be bit-for-bit identical to pRFv2.

The original authorship recorded in the pRFv2 source files is retained.
pRFv3-specific refactoring and additions were developed by Austin Kuo with
Codex.

See [CHANGES_FROM_PRFV2.md](CHANGES_FROM_PRFV2.md) for a detailed,
commit-oriented comparison with pRFv2.

## Compatibility goals

- The standard Gaussian, Gaussian-HDR, CSS, difference-of-Gaussians, and
  divisive-Gaussian paths are intended to remain numerically equivalent to
  pRFv2, but not bit-for-bit identical.
- Exact duplicate detection and exact cache keys are used. Approximate model
  responses are not substituted for legacy calculations.
- Rectilinear Gaussian RFs are constructed from separable one-dimensional
  components, and dense double stimuli are projected using MATLAB's optimized
  matrix multiplication. These change floating-point reduction order by a few
  last-place bits and can occasionally change an optimizer's path on a nearly
  flat objective.
- `gaussian-DoG-CSS`, `divNorm`, and the optional `fmincon` solver are explicit
  exceptions. Their validation and numerical behavior were intentionally
  changed and they should be treated as experimental.

## Requirements

- mrTools/mrLoadRet and its normal pRF dependencies.
- MATLAB Parallel Computing Toolbox for runs with `numWorkers > 1`.
- MATLAB Optimization Toolbox only when selecting `fmincon`.
- An MGL installation capable of reconstructing the supported retinotopy
  stimulus files.

Parallel support is optional. If a requested pool cannot be started, pRFv3
warns and continues serially.

## Installation and activation

Place the `gru` repository on the MATLAB path and use `mlrPlugin` to select
`pRFv3`. Do not select `pRFv2` and `pRFv3` together: both plugins intentionally
provide functions named `pRF`, `pRFFit`, `pRFGUI`, and `pRFPlot`.

When an mrLoadRet view installs pRFv3, `pRFv3Plugin` places the pRFv3 runtime
directories first, verifies the active implementations, and removes pRFv2
from `selectedPlugins` if both were selected.

Activation discovers the plugin root from the files currently present and
adds the complete pRFv3 tree. Internal helpers may therefore be reorganized
without encoding their subdirectory in the path verifier. Direct calls to
`pRF` also bootstrap this tree before opening the parameter GUI.

Users who call pRF functions without first opening mrLoadRet can add this to
their personal or site `startup.m`, after mrTools and `gru` have been added to
the path:

```matlab
pRFv3Startup;
```

`pRFv3Startup` activates pRFv3 only when it is present in the saved
`selectedPlugins` preference. Merely placing the helper on the MATLAB path
does not execute it.

MATLAB's Current Folder has precedence over ordinary path entries. Do not
start MATLAB with the Current Folder set inside `plugins/pRFv2` when trying to
activate pRFv3; change to a neutral session or project directory first.

## Running an analysis

From mrLoadRet, choose **Analysis > pRF Analysis (v3)**. The analysis produces
the same primary overlays as pRFv2:

- `r2`
- `polarAngle`
- `eccentricity`
- `rfHalfWidth`

The programmatic parameter workflow remains available:

```matlab
v = newView;
[v, params] = pRF(v, [], 'justGetParams=1');
v = pRF(v, params);
```

Older saved parameter structures are reconciled with defaults for newly added
fields before fitting.

## Important parameters

| Parameter | Behavior |
| --- | --- |
| `numWorkers` | Positive integer used by prefit construction and voxel fitting. `1` is true serial execution. Values greater than `1` request exactly that pool size. |
| `quickPrefit` | Uses the smaller legacy prefit grid. It is faster but may provide a less favorable nonlinear starting point. |
| `prefitOnly` | Returns the best prefit-grid result without nonlinear optimization. |
| `verbose` | Prints every voxel fit and detailed setup output. It defaults to off because worker output is expensive at scale. |
| `progressEvery` | Completed voxels between aggregate progress lines. `0` selects approximately 1% intervals with a minimum interval of 500 voxels; `-1` disables progress lines. |
| `diagnosticVoxelCount` | Number of reproducible, stratified-random voxel fits printed in verbose format for quality control. `0` disables sampling. |
| `diagnosticSeed` | Private seed used only for diagnostic-voxel selection; it does not alter MATLAB's global random state or fit results. |
| `saveStimImage` | Saves the canonical reconstructed structure as `pRFStimImage` in the linked stimulus file. |
| `pRFStimImageExists` | Read-only GUI indicator for `dispStimScan`; checked only when every linked constituent stimfile already contains `pRFStimImage`. It is derived display state and is not saved in the returned analysis parameters. |
| `recomputeStimImage` | Ignores an embedded `pRFStimImage` and reconstructs it from the recorded task information. |
| `saveStimImageOnly` | One-shot action that reconstructs and saves selected stimuli, then exits without loading voxel time series, starting a pool, or creating an analysis. |
| `xFlipStimulus`, `yFlipStimulus` | Explicit additional corrections. An unchecked box performs no additional spatial flip. |
| `timeShiftStimulus` | Circularly shifts the reconstructed movie by the requested number of volumes. |
| `stimImageDiffTolerance` | Percentage of differing frames tolerated by the legacy averaged-scan image check before asking whether to continue. |

## Parallel execution and memory use

The GUI reports the local profile's available worker count without starting a
pool. At run time:

- `numWorkers = 1` uses ordinary `for` loops and does not inspect, start, or
  close a parallel pool.
- `numWorkers = N`, for `N > 1`, reuses an existing pool of exactly `N`
  workers or replaces a mismatched pool.
- Pool startup or worker path-activation failure falls back to serial fitting.
- The same resolved setting controls both prefit construction and voxel fits;
  prefit cannot start an independent pool.

Voxel time series are loaded in adaptive blocks. The default targets enough
voxels to keep workers occupied while capping the estimated block working set
at roughly 512 MiB. Only the numeric time-series block is sent through the
voxel loop, while scan-invariant stimulus, model, and prefit state is prepared
once and kept resident on workers when supported.

The standalone `mlrRunPRF` interface uses the same policy: `doParallel=-1`
requests all locally available workers, `0` or `1` is serial, and an integer
greater than `1` requests that exact count.

## Fast double-precision kernel

For the ordinary dense-double stimulus path, pRFv3 uses two related numerical
optimizations:

- A Gaussian on an exact Cartesian `ndgrid` or `meshgrid` is formed as the
  outer product of its one-dimensional x and y components. Arbitrary grids and
  unsupported numeric types use the general two-dimensional expression.
- Unique stimulus frames are reshaped into a pixel-by-frame matrix and
  projected with optimized matrix multiplication. Duplicate-frame expansion,
  model construction, HRF convolution, filtering, and fitting otherwise retain
  their established behavior.

This is the only standard dense-double path; there is no legacy arithmetic
mode to select. In two development comparisons using 100 deterministic
left-hemisphere voxels each, nonlinear fitting was 3.37x and 4.00x faster. The
newer comparison had identical valid-fit masks, a median parameter and r2
difference of zero, a maximum r2 difference of `6.47e-5`, and no classification
changes at r2 cutoffs from 0.05 through 0.70. Its largest center/width parameter
difference was 5.7% of one stimulus pixel and its worst prediction correlation
was `0.99999127`.

## Progress and quality control

Aggregate progress is reported from the client using completed-voxel counts,
elapsed time, and an estimated remaining time. Parallel workers send compact
completion messages instead of streaming one console line per voxel.

Diagnostic voxels are selected reproducibly across the requested voxel list.
Their summaries use the same formatting as `verbose` output. At completion,
pRFv3 also warns when sampled inputs or fits are unexpectedly identical, when
no sampled fits are valid, or when all finite voxel fits in an analysis are
exactly identical.

Monitoring failures are warning-only and do not invalidate completed fits.

## Stimulus orientation

`pRFGetStimImageFromStimfile` now normalizes a newly reconstructed MGL raster
into the canonical pRF world-coordinate orientation before returning or
saving it. The main pRF analysis consumes that returned image directly.
Therefore:

- If the participant saw the intended world-coordinate stimulus, leave both
  flip checkboxes unchecked.
- Check a flip only to correct an actual horizontal or vertical reversal in
  the displayed experiment.
- A manually generated and saved `pRFStimImage` has the same default
  orientation as one generated during the pRF analysis.

Important migration note: pRFv2 saved its embedded `pRFStimImage` before the
old implicit orientation correction, and the saved structure has no format
version marker. Before using a pRFv2-saved embedded image with pRFv3, reconstruct
and resave it once with `recomputeStimImage` and `saveStimImage`, or use
`saveStimImageOnly`.

## Averaged scans and stimulus reuse

For a non-concatenated scan with multiple linked stimfiles, pRFv3 uses the
generic `averageTSeries` metadata field `params.scanList` to establish averaged
scan provenance. It does not look for a particular subject, session, group
name, filename, or scan number.

It then compares every constituent against the first:

- Existing embedded stimuli require exact equality of `im`, `x`, `y`, and
  `t`.
- Forced reconstructions compare the supported task program, stimulus
  geometry and state, derived task parameters and randomized variables,
  trial-volume sequence, `volTrigRatio`, and screen geometry.

When these values match, the stimulus is constructed once. If saving was
requested, the same canonical `pRFStimImage` is written to every constituent
stimfile. Any missing, unsupported, or unequal information takes the legacy
fallback that constructs every stimulus independently and compares the
resulting images.

Nonzero `shiftList` and `reverseList`, averages of concatenations, and nested
averages remain unsupported. In forced-reconstruction comparison, exact
trigger-derived `stim.t` vectors are not compared; equivalence is based on the
stimulus/task settings and frame content.

Supported reconstruction task files are matched case-insensitively:

- `mglRetinotopy.m`
- `gruRetinotopy.m`
- `offsetRetinotopy.m`
- `mglDoubleBars.m`
- `mglMetalRetinotopy.m`

Stimulus types 1 through 5 are accepted.

## Saving and continuing analyses

When an analysis filename already exists, pRFv3 uses the Command Window rather
than an overwrite dialog. Depending on the mrTools overwrite policy, the user
can merge, overwrite, save a separate copy, or cancel. Temporary policy changes
are kept in memory and restored after the save.

Continuing an analysis uses the name of the analysis actually selected in the
GUI. A stale `params.saveName` stored inside a renamed analysis can no longer
redirect the merge into an unrelated analysis named `pRF`.

## Plotting and voxel interrogation

- ROI pRF outlines are batched into a small number of patch objects and use
  96-point circles. Higher-r2 outlines remain on top.
- RF-center and eccentricity/size point loops are replaced with scatter
  objects.
- Coverage and center plots use equal physical x/y scaling.
- ROI titles display `pRF coverages`, the plotted voxel count, and the r2
  cutoff; underscores are shown literally.
- A voxel that is absent from the saved analysis is fitted only after an
  explicit Yes/No confirmation.
- Shift-click retains the original behavior of suppressing the aggregate ROI
  plot while displaying the single-voxel result.
- A local compatibility guard prevents MATLAB R2025b's legacy `getpts`
  callback from treating a modifier-only key event as a scalar character.

## Models and solvers

Available models remain:

- `gaussian`
- `gaussian-hdr`
- `gaussian-css`
- `gaussian-diffs`
- `gaussian-divs`
- `gaussian-DoG-CSS` (experimental)
- `divNorm` (experimental)

Available solvers are `nelder-mead`, `levenberg-marquardt`,
`nelder-mead-bnd`, and the optional experimental `fmincon` SQP path.

`gaussian-DoG-CSS` and `divNorm` require a bounded solver and validated signs,
positive widths, and positive nonlinear constants. Their neural nonlinearities
are now applied before HRF convolution. These corrections intentionally make
their results different from the corresponding pRFv2 implementations.

## Removed legacy components

pRFv3 does not include the obsolete controller/split execution stack, machine
profiles, passwordless-SSH test helper, or Sherlock support. It also removes
the unused `fminsearchcon` helper, model template, split-merge helper, and the
unimplemented cross-validation GUI option.

## Developer notes

The main performance components are:

- `helpers/pRFFitPrepared.m`: direct prepared-context dispatch.
- `models/pRFPrepareStimulusProjection.m`: exact duplicate-frame/run metadata.
- `models/pRFMakeGaussian.m`: separable Cartesian Gaussian construction with
  a general-grid fallback.
- `models/pRFProjectStimulus.m`: shared projection, expansion, and bounded
  worker-local result cache.
- `helpers/pRFResolveNumWorkers.m`: noninteractive worker resolution.
- `helpers/pRFProgressReporter.m` and `helpers/pRFFormatVoxelFit.m`: client-side progress and
  sampled QC.
- `helpers/pRFSaveAnalysis.m`: command-line conflict handling.
- `helpers/pRFv3ActivatePath.m` and `helpers/pRFv3Startup.m`: version-specific
  path activation.

Projection, prefit, and objective caches are bounded and use exact identity or
exact value keys. Cache and progress telemetry remain available for debugging,
although the main analysis does not save that instrumentation in its output.
