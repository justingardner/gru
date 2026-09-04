% pRF.m
%
%        $Id:$ 
%      usage: pRF(v,params,varargin)
%         by: justin gardner
%       date: 11/20/11
%    purpose: compute pRF analysis on MLR data
%
%             if you just want a default parameter structure you
%             can do:
% 
%             v = newView;
%             [v params] = pRF(v,[],'justGetParams=1','defaultParams=1','scanList=1')
%
%             Note that justGetParams, defaultParams and scanList are independent parameters, so
%             if you want, say to bring up the GUI to set the params, but not run the analysis, you
%             can do:
%             [v params] = pRF(v,[],'justGetParams=1');
%             The exception is the one-shot saveStimImageOnly action: when selected, stimulus
%             images are reconstructed and saved before the returned action flags are cleared.
%
function [v d] = pRF(v,params,varargin)

% check arguments
if nargin < 1
  help pRF
  return
end

% A pRF analysis can be called directly without first installing plugins in
% an mrLoadRet window. Bootstrap from this entry point, then let the shared
% activator discover all helpers from the pRFv3 tree that currently exists.
pluginRoot = fileparts(mfilename('fullpath'));
addpath(genpath(pluginRoot),'-begin');
pRFv3ActivatePath(pluginRoot);

d = [];

% a version number in case we make major changes
pRFVersion = 1;

% params defaults to empty
if nargin < 2,params =[];end

% other arguments
justGetParams=[];defaultParams=[];scanList=[];
groupNum=[];
getArgs(varargin,{'justGetParams=0','defaultParams=0','scanList=[]','groupNum=[]'});

% first get parameters
if isempty(params)
  % get group
  if isempty(groupNum),groupNum = viewGet(v,'curGroup');end
  % put up the gui
  params = pRFGUI('v',v,'groupNum',groupNum,'defaultParams',defaultParams,'scanList',scanList);
end

% Ordinarily justGetParams returns immediately. The save-stimulus-only flag
% is a one-shot action, however: execute it first and return a cleared flag so
% the resulting params can subsequently be passed back to pRF for fitting.
saveStimImageOnlyRequested = false;
if isstruct(params) && isfield(params,'pRFFit') && ...
    isstruct(params.pRFFit) && isscalar(params.pRFFit) && ...
    isfield(params.pRFFit,'saveStimImageOnly')
  stimulusOnlyValue = params.pRFFit.saveStimImageOnly;
  saveStimImageOnlyRequested = ...
    (isnumeric(stimulusOnlyValue) || islogical(stimulusOnlyValue)) && ...
    isreal(stimulusOnlyValue) && isscalar(stimulusOnlyValue) && ...
    isfinite(stimulusOnlyValue) && (stimulusOnlyValue ~= 0);
end
if justGetParams && ~saveStimImageOnlyRequested
  d = params;
  return
end

% Reconcile params with current status of group and ensure that it has
% the required fields. 
params = defaultReconcileParams([],params);

% Abort if params empty
if isempty(params),return,end

% check the params
params = checkPRFparams(params);

% set the group
v = viewSet(v,'curGroup',params.groupName);

% Optionally stop after rebuilding and saving the stimulus images. This must
% run before pool resolution, voxel restriction, overlay allocation, or
% time-series loading so the stimulus-only operation cannot start a fit.
if params.pRFFit.saveStimImageOnly
  global gpRFFitStimImage
  stimParams = params.pRFFit;
  stimParams.saveStimImage = true;
  stimParams.recomputeStimImage = true;
  stimParams.numWorkers = 1;
  fprintf('(pRF) Reconstructing and saving stimulus images for %i selected scan(s).\n', ...
    numel(params.scanNum));
  for scanNum = params.scanNum
    fprintf('(pRF) Reconstructing stimulus for %s:%i.\n',params.groupName,scanNum);
    stim = pRFFit(v,scanNum,[],[],[], ...
      'justGetStimImage=1','fitTypeParams',stimParams);
    if isempty(stim)
      gpRFFitStimImage = [];
      warning('pRF:StimulusOnlyFailed', ...
        ['Could not reconstruct and save the stimulus for %s:%i. ' ...
         'No pRF analysis was created; stimulus files processed before this ' ...
         'failure may already have been updated.'],params.groupName,scanNum);
      return
    end
    clear stim
  end
  % The saved files are now authoritative; do not retain the final movie in
  % the process-global in-memory cache after this explicitly save-only run.
  gpRFFitStimImage = [];
  fprintf('(pRF) Finished saving stimulus images. No pRF analysis was created.\n');
  % All three controls describe stimulus-generation actions. Clear them in
  % returned params so a one-shot justGetParams call can feed directly into
  % a normal fit without reconstructing or saving the stimulus a second time.
  params.pRFFit.saveStimImageOnly = false;
  params.pRFFit.saveStimImage = false;
  params.pRFFit.recomputeStimImage = false;
  if justGetParams
    d = params;
  end
  return
end

% create the parameters for the r2 overlay
dateString = datestr(now);
r2.name = 'r2';
r2.groupName = params.groupName;
r2.function = 'pRF';
r2.reconcileFunction = 'defaultReconcileParams';
r2.data = cell(1,viewGet(v,'nScans'));
r2.date = dateString;
r2.params = cell(1,viewGet(v,'nScans'));
r2.range = [0 1];
r2.clip = [0 1];
% colormap is made with a little bit less on the dark end
r2.colormap = hot(312);
r2.colormap = r2.colormap(end-255:end,:);
r2.alpha = 1;
r2.colormapType = 'setRangeToMax';
r2.interrogator = 'pRFPlot';
r2.mergeFunction = 'pRFMergeParams';

% create the parameters for the polarAngle overlay
polarAngle = r2;
polarAngle.name = 'polarAngle';
polarAngle.range = [-pi pi];
polarAngle.clip = [-pi pi];
polarAngle.colormapType = 'normal';
polarAngle.colormap = hsv(256);

% create the parameters for the eccentricity overlay
eccentricity = r2;
eccentricity.name = 'eccentricity';
eccentricity.range = [0 15];
eccentricity.clip = [0 inf];
eccentricity.colormapType = 'normal';
eccentricity.colormap = copper(256);

% create the parameters for the rfHalfWidth overlay
rfHalfWidth = r2;
rfHalfWidth.name = 'rfHalfWidth';
rfHalfWidth.range = [0 15];
rfHalfWidth.clip = [0 inf];
rfHalfWidth.colormapType = 'normal';
rfHalfWidth.colormap = pink(256);

% Resolve the GUI/requested count once. Store the effective value so an
% unavailable pool falls back to serial for prefit construction as well as
% voxel fitting, without a second startup attempt or worker-count popup.
nProcessors = pRFResolveNumWorkers(params.pRFFit.numWorkers);
params.pRFFit.numWorkers = nProcessors;

% code snippet for clearing precomputed prefit
%global gpRFFitStimImage;gpRFFitStimImage = [];

dispHeader
disp(sprintf('(pRF) Running on scans %s:%s (restrict %s)',params.groupName,num2str(params.scanNum,'%i '),params.restrict ));

for scanNum = params.scanNum
  % see how long it took
  scanStartTic = tic;

  % get voxels that we are restricted to
  [x y z] = getVoxelRestriction(v,params,scanNum);
  if isempty(x)
    disp(sprintf('(pRF) No voxels to analyze with current restriction'));
    return
  end

  % get total number of voxels
  n = length(x);

  % get scan dims
  scanDims = viewGet(v,'scanDims',scanNum);

  % init overlays
  r2.data{scanNum} = nan(scanDims);
  polarAngle.data{scanNum} = nan(scanDims);
  eccentricity.data{scanNum} = nan(scanDims);
  rfHalfWidth.data{scanNum} = nan(scanDims);

  % save pRF parameters
  pRFAnal.d{scanNum}.ver = pRFVersion;
  pRFAnal.d{scanNum}.linearCoords = [];
  pRFAnal.d{scanNum}.params = [];

  % get some information from pRFFit that will be used again in
  % the fits, including concatInfo, stim, prefit, etc.
  fit = pRFFit(v,scanNum,[],[],[],'fitTypeParams',params.pRFFit, ...
    'returnPrefit',true,'returnPreparedContext',true);
  if isempty(fit),return,end
  stim = fit.stim;
  % Projection metadata can contain shared handle objects intended only for
  % the live worker context. Keep saved analyses compatible and compact.
  pRFAnal.d{scanNum}.stim = cellArray(stripStimulusProjectionMetadata(stim));
  pRFAnal.d{scanNum}.stimX = fit.stimX;
  pRFAnal.d{scanNum}.stimY = fit.stimY;
  pRFAnal.d{scanNum}.stimT = fit.stimT;
  concatInfo = fit.concatInfo;
  pRFAnal.d{scanNum}.concatInfo = fit.concatInfo;
  prefit = fit.prefit;
  paramsInfo = fit.paramsInfo;
  pRFAnal.d{scanNum}.paramsInfo = paramsInfo;
  preparedContext = fit.preparedContext;

  % grab all these fields and stick them onto a structure called paramsInfo
  % preallocate some space
  rawParams = nan(fit.nParams,n);
  r = nan(n,fit.concatInfo.n);
  thisr2 = nan(1,n);
  thisPolarAngle = nan(1,n);
  thisEccentricity = nan(1,n);
  thisRfHalfWidth = nan(1,n);
  fitWasReturned = false(1,n);

  % get some info about the scan to pass in (which prevents
  % pRFFit from calling viewGet - which is problematic for distributed computing
  framePeriod = viewGet(v,'framePeriod');
  junkFrames = viewGet(v,'junkFrames',scanNum);

  % compute pRF for each voxel in the restriction
  if params.pRFFit.prefitOnly,algorithm='prefit-only';else,algorithm=params.pRFFit.algorithm;end

  % disp info about fitting
  dispHeader;
  disp(sprintf('(pRF) Scan %s:%i (restrict %s) running on %i processor(s)',params.groupName,scanNum,params.restrict,nProcessors));
  disp(sprintf('(pRF) Computing %s fits using %s for %i voxels',params.pRFFit.rfType,algorithm,n));
  dispHeader;

  % Load enough voxels to keep every worker busy across relatively few
  % parfor launches, while bounding the estimated time-series working set.
  blockSize = chooseVoxelBlockSize(n,fit.concatInfo,nProcessors,params.pRFFit);
  if isfield(params.pRFFit,'verbose') && params.pRFFit.verbose
    disp(sprintf('(pRF) Using adaptive blocks of up to %i voxels',blockSize));
  end

  % Report true completed-voxel counts on the client without forwarding a
  % line of worker output for every fit. Diagnostic voxels are selected with
  % a private random stream, so monitoring cannot change fit RNG state.
  progressReporter = pRFProgressReporter(n,params.pRFFit.progressEvery, ...
    params.pRFFit.diagnosticVoxelCount,params.pRFFit.diagnosticSeed,scanNum);
  diagnosticVoxelMask = false(1,n);
  diagnosticVoxelMask(progressReporter.DiagnosticIndices) = true;
  liveVoxelProgressEnabled = isfinite(progressReporter.ProgressEvery) && ...
    progressReporter.ProgressEvery < n;

  % A Constant transfers the large immutable stimulus/prefit context once
  % per worker and keeps it resident across every voxel block in this scan.
  preparedContextConstant = [];
  if nProcessors > 1 && exist('parallel.pool.Constant','class') == 8
    try
      preparedContextConstant = parallel.pool.Constant( ...
        @() initializePreparedWorkerContext(preparedContext), ...
        @releasePreparedWorkerContext);
    catch contextException
      warning('pRF:PreparedContextConstant', ...
        'Could not keep the prepared context resident on workers: %s', ...
        contextException.message);
    end
  end
  % DataQueue callbacks execute on the client and are the MathWorks-supported
  % way to monitor PARFOR. Cleanup always detaches this listener before
  % releasing the worker-resident context, including during Ctrl+C unwinding.
  progressQueue = [];
  progressListener = [];
  progressQueueEnabled = false;
  if nProcessors > 1 && (liveVoxelProgressEnabled || any(diagnosticVoxelMask)) && ...
      exist('parallel.pool.DataQueue','class') == 8
    try
      progressQueue = parallel.pool.DataQueue;
      progressListener = afterEach(progressQueue,@(payload) progressReporter.record(payload));
      progressQueueEnabled = true;
    catch progressException
      warning('pRF:ProgressQueueUnavailable', ...
        ['Could not start live voxel progress reporting: %s. ' ...
         'Progress will update after each voxel block.'],progressException.message);
    end
  end
  scanResourceCleanup = onCleanup(@() cleanupScanResources( ...
    progressListener,preparedContextConstant));

  fitStartTic = tic;
  % break into blocks of voxels to go easy on memory
  % if blockSize = n then this just does on block at a time.
  for blockStart = 1:blockSize:n

    % display information about what we are doing
    % get blockEnd
    blockEnd = min(blockStart + blockSize-1,n);
    thisBlockSize = blockEnd-blockStart+1;
    
    % load ROI
    loadROI = makeEmptyROI(v,'scanNum',scanNum,'groupNum',params.groupName);
    loadROI.coords(1,1:thisBlockSize) = x(blockStart:blockEnd);
    loadROI.coords(2,1:thisBlockSize) = y(blockStart:blockEnd);
    loadROI.coords(3,1:thisBlockSize) = z(blockStart:blockEnd);
    % load all time series for block, we do this to pass into pRFFit. Generally
    % the purpose here is that if we run on distributed computing, we
    % can't load each voxel's time series one at a time. If this is
    % too large for memory then you can comment this out and not
    % pass it into pRFFit and pRFFit will load the tSeries itself
    loadROI = loadROITSeries(v,loadROI,scanNum,params.groupName,'keepNAN',true);
    % reorder x,y,z coordinates since they can get scrambled in loadROITSeries
    x(blockStart:blockEnd) = loadROI.scanCoords(1,1:thisBlockSize);
    y(blockStart:blockEnd) = loadROI.scanCoords(2,1:thisBlockSize);
    z(blockStart:blockEnd) = loadROI.scanCoords(3,1:thisBlockSize);
    % Broadcast only the numeric time-series block to workers; sending the
    % complete ROI structure also transfers coordinates and other metadata
    % that the fit loop never reads.
    blockTSeries = loadROI.tSeries;
    % keep the linear coords
    pRFAnal.d{scanNum}.linearCoords = [pRFAnal.d{scanNum}.linearCoords sub2ind(scanDims,x(blockStart:blockEnd),y(blockStart:blockEnd),z(blockStart:blockEnd))];

    % now loop over each voxel. Both branches call the same prepared fitter;
    % the Constant branch only changes how its immutable context is shipped.
    if nProcessors <= 1
      % A plain PARFOR can auto-start MATLAB's default pool even when MLR's
      % parallel preference is disabled. Use a true serial loop so
      % numWorkers=1 and failed-pool fallbacks remain noninteractive.
      for i = blockStart:blockEnd
        voxelTSeries = blockTSeries(i-blockStart+1,:)';
        voxelFit = pRFFitPrepared('fit',preparedContext, ...
          voxelTSeries,x(i),y(i),z(i),i,n);
        fitWasReturned(i) = ~isempty(voxelFit);
        if ~isempty(voxelFit)
          thisr2(i) = voxelFit.r2;
          thisPolarAngle(i) = voxelFit.polarAngle;
          thisEccentricity(i) = voxelFit.eccentricity;
          thisRfHalfWidth(i) = voxelFit.std;
          rawParams(:,i) = voxelFit.params(:);
          r(i,:) = voxelFit.r;
        end
        if liveVoxelProgressEnabled || diagnosticVoxelMask(i)
          progressReporter.record(makeVoxelProgressPayload( ...
            i,n,x(i),y(i),z(i),voxelTSeries,voxelFit,diagnosticVoxelMask(i)));
        end
      end
    elseif isempty(preparedContextConstant)
      parfor i = blockStart:blockEnd
        voxelTSeries = blockTSeries(i-blockStart+1,:)';
        voxelFit = pRFFitPrepared('fit',preparedContext, ...
          voxelTSeries,x(i),y(i),z(i),i,n);
        fitWasReturned(i) = ~isempty(voxelFit);
        if ~isempty(voxelFit)
          thisr2(i) = voxelFit.r2;
          thisPolarAngle(i) = voxelFit.polarAngle;
          thisEccentricity(i) = voxelFit.eccentricity;
          thisRfHalfWidth(i) = voxelFit.std;
          rawParams(:,i) = voxelFit.params(:);
          r(i,:) = voxelFit.r;
        end
        if progressQueueEnabled && (liveVoxelProgressEnabled || diagnosticVoxelMask(i))
          sendVoxelProgress(progressQueue, ...
            i,n,x(i),y(i),z(i),voxelTSeries,voxelFit,diagnosticVoxelMask(i));
        end
      end
    else
      parfor i = blockStart:blockEnd
        workerContext = preparedContextConstant.Value;
        voxelTSeries = blockTSeries(i-blockStart+1,:)';
        voxelFit = pRFFitPrepared('fit',workerContext, ...
          voxelTSeries,x(i),y(i),z(i),i,n);
        fitWasReturned(i) = ~isempty(voxelFit);
        if ~isempty(voxelFit)
          thisr2(i) = voxelFit.r2;
          thisPolarAngle(i) = voxelFit.polarAngle;
          thisEccentricity(i) = voxelFit.eccentricity;
          thisRfHalfWidth(i) = voxelFit.std;
          rawParams(:,i) = voxelFit.params(:);
          r(i,:) = voxelFit.r;
        end
        if progressQueueEnabled && (liveVoxelProgressEnabled || diagnosticVoxelMask(i))
          sendVoxelProgress(progressQueue, ...
            i,n,x(i),y(i),z(i),voxelTSeries,voxelFit,diagnosticVoxelMask(i));
        end
      end
    end

    % A completed PARFOR block is authoritative. Reconcile any callbacks
    % still queued on the client and guarantee that every QC sample in this
    % block is printed even if live queue setup was unavailable.
    progressReporter.reconcile(blockStart:blockEnd);
    blockDiagnosticIndices = progressReporter.DiagnosticIndices( ...
      progressReporter.DiagnosticIndices >= blockStart & ...
      progressReporter.DiagnosticIndices <= blockEnd);
    for diagnosticIndex = blockDiagnosticIndices
      blockIndex = diagnosticIndex-blockStart+1;
      storedFit = makeStoredDiagnosticFit(diagnosticIndex,fitWasReturned, ...
        thisr2,thisPolarAngle,thisEccentricity,thisRfHalfWidth,rawParams,r, ...
        params.pRFFit.rfType,params.pRFFit.prefitOnly);
      progressReporter.record(makeVoxelProgressPayload( ...
        diagnosticIndex,n,x(diagnosticIndex),y(diagnosticIndex),z(diagnosticIndex), ...
        blockTSeries(blockIndex,:)',storedFit,true));
    end
  end


  % Stop client callbacks before releasing worker-held context. This ordering
  % avoids progress activity during parallel-resource teardown.
  deleteProgressListener(progressListener);
  progressReporter.finish;
  warnIfAllVoxelFitsIdentical(rawParams,thisr2,fitWasReturned);

  % Release worker-held handle objects and their bounded projection caches
  % before advancing to another scan.
  deletePreparedContextConstant(preparedContextConstant);
  clear scanResourceCleanup

  % Set each final overlay once. This loop intentionally retains the original
  % assignment order, including its behavior if coordinates are duplicated.
  for i = 1:n
    r2.data{scanNum}(x(i),y(i),z(i)) = thisr2(i);
    polarAngle.data{scanNum}(x(i),y(i),z(i)) = thisPolarAngle(i);
    eccentricity.data{scanNum}(x(i),y(i),z(i)) = thisEccentricity(i);
    rfHalfWidth.data{scanNum}(x(i),y(i),z(i)) = thisRfHalfWidth(i);
  end
  % display time update
  dispHeader;
  disp(sprintf('(pRF) Fitting %i voxels took %s.',n,mlrDispElapsedTime(toc(fitStartTic))));
  dispHeader;

  pRFAnal.d{scanNum}.params = rawParams;
  pRFAnal.d{scanNum}.r = r;

  % Obtain one canonical model from the first successful finite voxel fit.
  % Keep the candidate index local to this scan and use the parameters from
  % that same voxel; coordinates and parameters must never get out of sync.
  canonicalCandidates = find(fitWasReturned & isfinite(thisr2) & ...
    all(isfinite(rawParams),1));
  canonicalFit = [];
  for canonicalIndex = canonicalCandidates
    fprintf('\n(pRF) Attempting to obtain canonical from voxel [%d %d %d].', ...
      x(canonicalIndex),y(canonicalIndex),z(canonicalIndex));
    try
      candidateFit = pRFFit(v,scanNum,x(canonicalIndex),y(canonicalIndex), ...
        z(canonicalIndex),'params',rawParams(:,canonicalIndex),'stim',stim, ...
        'concatInfo',concatInfo,'prefit',prefit, ...
        'fitTypeParams',params.pRFFit,'dispIndex',canonicalIndex, ...
        'dispN',n,'framePeriod',framePeriod,'junkFrames',junkFrames, ...
        'paramsInfo',paramsInfo,'getModelResponse',1);
      if ~isempty(candidateFit) && isfield(candidateFit,'canonical') && ...
          ~isempty(candidateFit.canonical)
        canonicalFit = candidateFit;
        fprintf('\n(pRF) Using valid voxel: [%d %d %d].\n', ...
          x(canonicalIndex),y(canonicalIndex),z(canonicalIndex));
        break
      end
    catch canonicalException
      warning('pRF:CanonicalVoxelFailed', ...
        'Could not obtain a canonical model from voxel [%d %d %d]: %s', ...
        x(canonicalIndex),y(canonicalIndex),z(canonicalIndex), ...
        canonicalException.message);
    end
  end
  if isempty(canonicalFit)
    error('pRF:NoCanonicalVoxel', ...
      'Could not obtain a canonical model from any successful finite fit in scan %d.', ...
      scanNum);
  end
  pRFAnal.d{scanNum}.canonicalModel = canonicalFit.canonical;

  iScan = find(params.scanNum == scanNum);
  thisParams.scanNum = params.scanNum(iScan);
  r2.params{scanNum} = thisParams;
  polarAngle.params{scanNum} = thisParams;
  eccentricity.params{scanNum} = thisParams;
  rfHalfWidth.params{scanNum} = thisParams;
  % display how long it took
  disp(sprintf('(pRF) Fitting for %s:%i took in total: %s',params.groupName,scanNum,mlrDispElapsedTime(toc(scanStartTic))));
end

% install analysis
pRFAnal.name = params.saveName;
pRFAnal.type = 'pRFAnal';
pRFAnal.groupName = params.groupName;
pRFAnal.function = 'pRF';
pRFAnal.reconcileFunction = 'defaultReconcileParams';
pRFAnal.mergeFunction = 'pRFMergeParams';
pRFAnal.guiFunction = 'pRFGUI';
pRFAnal.params = params;
pRFAnal.overlays = [r2 polarAngle eccentricity rfHalfWidth];
pRFAnal.curOverlay = 1;
pRFAnal.date = dateString;
v = viewSet(v,'newAnalysis',pRFAnal);

% Save without opening mrTools' overwrite/rename dialogs. Continuing an
% analysis retains its established automatic-merge behavior.
forceMerge = isfield(params,'mergeAnalysis') && params.mergeAnalysis;
[v,savedAnalysisName] = pRFSaveAnalysis(v,pRFAnal.name,forceMerge);
pRFAnal.name = savedAnalysisName;

if ~isempty(viewGet(v,'fignum'))
  refreshMLRDisplay(viewGet(v,'viewNum'));
end

%set(viewGet(v,'figNum'),'Pointer','arrow');drawnow

% for output
if nargout > 1
  for i = 1:length(pRFAnal.d)
    if ~isempty(pRFAnal.d{i})
      pRFAnal.d{i}.r2 = r2.data{i};
    end
  end
  % make d strucutre
  if length(pRFAnal.d) == 1
    d = pRFAnal.d{1};
  else
    d = pRFAnal.d;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getVoxelRestriction    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x y z] = getVoxelRestriction(v,params,scanNum)

x = [];y = [];z = [];

if strncmp(params.restrict,'Base: ',6)
  % get the base name
  baseName = params.restrict(7:end);
  baseNums = [];
  if strcmp(baseName,'ALL')
    for iBase = 1:viewGet(v,'numBase')
      % if the base is a surface or flat then add to the list
      if any(viewGet(v,'baseType',iBase) == [1 2])
	baseNums(end+1) = iBase;
      end
    end
  else
    baseNums = viewGet(v,'baseNum',baseName);
  end
  % cycle through all bases that we are going to run on
  scanCoords = [];
  for iBase = 1:length(baseNums)
    % get the baseNum
    baseNum = baseNums(iBase);
    if isempty(baseNum)
      disp(sprintf('(pRF) Could not find base to restrict to: %s',params.restrict));
      continue
    end
    % get the base
    base = viewGet(v,'base',baseNum);
    if isempty(base)
      disp(sprintf('(pRF) Could not find base to restrict to: %s',params.restrict));
      return;
    end
    % if flat or surface
    if any(base.type == [1 2])
      % get base coordinates from the coordMap
      for corticalDepth = 0:0.1:1
	if base.type == 1
	  % flat map
	  baseCoords = (base.coordMap.innerCoords + corticalDepth * (base.coordMap.outerCoords-base.coordMap.innerCoords));
	  baseCoords = reshape(baseCoords,prod(size(base.data)),3)';
	else
	  % surface
	  baseCoords = (base.coordMap.innerVtcs + corticalDepth * (base.coordMap.outerVtcs-base.coordMap.innerVtcs))';
	end
	% convert to 4xn array
	baseCoords(4,:) = 1;
	% and convert to scan coordinates
	base2scan = viewGet(v,'base2scan',scanNum,params.groupName,baseNum);
	scanCoords = [scanCoords round(base2scan*baseCoords)];
      end
    end
  end
  % check against scandims
  scanDims = viewGet(v,'scanDims',scanNum,params.groupName);
  scanCoords = mrSub2ind(scanDims,scanCoords(1,:),scanCoords(2,:),scanCoords(3,:));
  % remove duplicates and nans
  scanCoords = scanCoords(~isnan(scanCoords));
  scanCoords = unique(scanCoords);
  % convert back to x,y,z coordinates
  [x y z] = ind2sub(scanDims,scanCoords);
elseif strncmp(params.restrict,'ROI: ',5)
  % get the roi name
  roiName = params.restrict(6:end);
  scanCoords = getROICoordinates(v,roiName,scanNum,params.groupName,'straightXform=1');
  if isempty(scanCoords),return,end
  x = scanCoords(1,:);y = scanCoords(2,:);z = scanCoords(3,:);
elseif strncmp(params.restrict,'None',4)
  scanDims = viewGet(v,'scanDims',scanNum,params.groupName);
  [x y z]  = ndgrid(1:scanDims(1),1:scanDims(2),1:scanDims(3));
  x = x(:);y = y(:);z = z(:);
else
  error('pRF:UnknownVoxelRestriction', ...
    'Unknown voxel restriction "%s".',params.restrict);
end

%check if we have already computed Voxels
if isfield(params,'computedVoxels') && (length(params.computedVoxels)>=scanNum) && ~isempty(params.computedVoxels{scanNum})
  % get scan dims
  scanDims = viewGet(v,'scanDims',scanNum,params.groupName);
  % convert x, y, z to linear coords
  linearCoords = sub2ind(scanDims,x,y,z);
  % get new ones
  newLinearCoords = setdiff(linearCoords,params.computedVoxels{scanNum});
  if length(newLinearCoords) ~= length(linearCoords)
    % show what we are doing
    disp(sprintf('(pRF) Dropping %i voxels that have been already computed',length(linearCoords)-length(newLinearCoords)));
    % convert back to x, y, z
    [x y z] = ind2sub(scanDims,newLinearCoords);
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%
%    checkPRFparams    %
%%%%%%%%%%%%%%%%%%%%%%%%
function params = checkPRFparams(params)


% check the pRFFit params
checkFields = {{'stimImageDiffTolerance',5}, ...
               {'saveStimImageOnly',false}, ...
               {'progressEvery',0}, ...
               {'diagnosticVoxelCount',10}, ...
               {'diagnosticSeed',5489}};
for iFit = 1:length(params.pRFFit)

  % set defaults
  for iField = 1:length(checkFields)
    if ~isfield(params.pRFFit(iFit),checkFields{iField}{1})
      params.pRFFit(iFit).(checkFields{iField}{1}) = checkFields{iField}{2};
    end
  end

  % Query the local capacity only for an older parameter structure that
  % genuinely needs the new default. An explicit serial request therefore
  % never touches any Parallel Computing Toolbox API before it is resolved.
  if ~isfield(params.pRFFit(iFit),'numWorkers') || ...
      isempty(params.pRFFit(iFit).numWorkers)
    params.pRFFit(iFit).numWorkers = pRFResolveNumWorkers;
  end

  numWorkers = params.pRFFit(iFit).numWorkers;
  if ~(isnumeric(numWorkers) && isreal(numWorkers) && isscalar(numWorkers) && ...
      isfinite(numWorkers) && numWorkers >= 1 && numWorkers == fix(numWorkers))
    error('pRF:InvalidNumWorkers','pRFFit.numWorkers must be a positive integer.');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    chooseVoxelBlockSize         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function blockSize = chooseVoxelBlockSize(nVoxels,concatInfo,nProcessors,fitParams)

if isfield(fitParams,'voxelBlockSize') && ~isempty(fitParams.voxelBlockSize) && ...
    isnumeric(fitParams.voxelBlockSize) && isscalar(fitParams.voxelBlockSize) && ...
    isfinite(fitParams.voxelBlockSize) && fitParams.voxelBlockSize >= 1
  blockSize = min(nVoxels,max(1,round(fitParams.voxelBlockSize)));
  return
end

nProcessors = max(1,nProcessors);
nTimePoints = max(1,concatInfo.runTransition(end,2));
targetVoxels = max(240,512*nProcessors);
% Allow roughly 512 MiB for the loaded data and temporary per-voxel arrays.
% Six double vectors per voxel is deliberately conservative.
bytesPerVoxel = max(1,6*8*nTimePoints);
% Do not impose the historical 240-voxel floor after applying the memory
% limit: for unusually long scans that floor could exceed the budget it was
% intended to protect. Ordinary scans still use at least 240 voxels through
% targetVoxels above; only the memory-bound case is allowed to go lower.
memoryLimitedVoxels = max(1,floor((512*2^20)/bytesPerVoxel));
blockSize = min(nVoxels,min(targetVoxels,memoryLimitedVoxels));
blockSize = max(1,blockSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stripStimulusProjectionMetadata %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function strippedStim = stripStimulusProjectionMetadata(stim)

strippedStim = stim;
if iscell(strippedStim)
  for iStim = 1:numel(strippedStim)
    if isstruct(strippedStim{iStim}) && isfield(strippedStim{iStim},'projection')
      strippedStim{iStim} = rmfield(strippedStim{iStim},'projection');
    end
  end
elseif isstruct(strippedStim) && isfield(strippedStim,'projection')
  strippedStim = rmfield(strippedStim,'projection');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    releasePreparedWorkerContext %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function releasePreparedWorkerContext(~)

if exist('pRFProjectStimulus','file') == 2
  pRFProjectStimulus('resetCache');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initializePreparedWorkerContext  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function context = initializePreparedWorkerContext(context)

% Drop banks retained while the one-time prefit was constructed; the
% serialized resident context may have a different worker-local handle
% identity even though it represents the same immutable stimulus.
if exist('pRFProjectStimulus','file') == 2
  pRFProjectStimulus('resetCache');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deletePreparedContextConstant   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deletePreparedContextConstant(contextConstant)

if ~isempty(contextConstant)
  try
    if isvalid(contextConstant)
      delete(contextConstant);
    end
  catch
    % Cleanup must not mask a completed fit or the original exception.
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    deleteProgressListener   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deleteProgressListener(progressListener)

if ~isempty(progressListener)
  try
    if isvalid(progressListener)
      delete(progressListener);
    end
  catch
    % Monitoring cleanup must never mask a fit result or original exception.
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      cleanupScanResources      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleanupScanResources(progressListener,preparedContextConstant)

% Keep callback shutdown ahead of worker-context destruction on errors and
% Ctrl+C as well as on the normal completion path.
deleteProgressListener(progressListener);
deletePreparedContextConstant(preparedContextConstant);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    makeVoxelProgressPayload    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function payload = makeVoxelProgressPayload(index,totalVoxels,x,y,z,tSeries,fit,isDiagnostic)

% A scalar is substantially cheaper to serialize than a structure for the
% ordinary completion path, and is also the fail-safe diagnostic fallback.
payload = index;
if ~isDiagnostic,return,end
try
  payload = pRFProgressReporter.makeDiagnosticRecord( ...
    index,totalVoxels,[x y z],tSeries,fit);
catch
  % Monitoring must never interrupt or invalidate a completed voxel fit.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       sendVoxelProgress        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sendVoxelProgress(progressQueue,index,totalVoxels,x,y,z,tSeries,fit,isDiagnostic)

try
  send(progressQueue,makeVoxelProgressPayload( ...
    index,totalVoxels,x,y,z,tSeries,fit,isDiagnostic));
catch
  % The client reconciles the completed block and reconstructs diagnostic
  % records from stored fit arrays if queue delivery becomes unavailable.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    makeStoredDiagnosticFit     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = makeStoredDiagnosticFit(index,fitWasReturned,r2,polarAngle,eccentricity,rfHalfWidth,rawParams,runR,rfType,prefitOnly)

if ~fitWasReturned(index)
  fit = [];
  return
end
fit.r2 = r2(index);
fit.polarAngle = polarAngle(index);
fit.eccentricity = eccentricity(index);
fit.std = rfHalfWidth(index);
fit.params = rawParams(:,index);
fit.r = runR(index,:);
fit.rfType = rfType;
if prefitOnly
  % This marker is present on live prefit-only results and tells the shared
  % formatter to use pRFFit's prefit-only prefix and base output line.
  fit.bestFitVoxel = NaN;
elseif strcmp(rfType,'gaussian-css') && size(rawParams,1) >= 4
  fit.css = rawParams(4,index);
elseif strcmp(rfType,'divNorm') && size(rawParams,1) >= 8
  fit.stdSurround = rawParams(4,index);
  fit.gainCenter = rawParams(5,index);
  fit.gainSurround = rawParams(6,index);
  fit.actConst = rawParams(7,index);
  fit.divConst = rawParams(8,index);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% warnIfAllVoxelFitsIdentical       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function warnIfAllVoxelFitsIdentical(rawParams,r2,fitWasReturned)

validFits = find(fitWasReturned & isfinite(r2) & all(isfinite(rawParams),1));
if isempty(validFits)
  warning('pRF:NoValidVoxelFits','No voxel returned a finite pRF fit.');
  return
end
if numel(validFits) < 3,return,end
referenceParams = rawParams(:,validFits(1));
sameParams = all(all(rawParams(:,validFits) == referenceParams,1));
sameR2 = all(r2(validFits) == r2(validFits(1)));
if sameParams && sameR2
  warning('pRF:AllVoxelFitsIdentical', ...
    ['All %i finite voxel fits returned exactly identical parameters and r2. ' ...
     'Check ROI time-series loading and voxel indexing before trusting this analysis.'], ...
    numel(validFits));
end
