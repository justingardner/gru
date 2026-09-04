classdef pRFProgressReporter < handle
% pRFProgressReporter
%
% Client-side progress and diagnostic reporting for pRF voxel fits. Worker
% iterations send compact completion messages through a DataQueue; only
% rate-limited legacy-style progress and a small, reproducible sample of
% verbose-format fit summaries are printed.

properties (SetAccess=private)
  TotalVoxels
  ProgressEvery
  DiagnosticIndices
  CompletedCount = 0
  ProgressReportCount = 0
end

properties (Access=private)
  Completed
  NextProgress
  StartTic
  LastPrintedCount = 0
  DiagnosticPrinted
  DiagnosticRecords
  Finished = false
  MonitoringWarningIssued = false
end

methods
  function obj = pRFProgressReporter(totalVoxels,progressEvery,diagnosticVoxelCount,diagnosticSeed,scanNum)

    if nargin < 1 || isempty(totalVoxels),totalVoxels = 0;end
    if nargin < 2 || isempty(progressEvery),progressEvery = 0;end
    if nargin < 3 || isempty(diagnosticVoxelCount),diagnosticVoxelCount = 10;end
    if nargin < 4 || isempty(diagnosticSeed),diagnosticSeed = 5489;end
    if nargin < 5 || isempty(scanNum),scanNum = 1;end

    obj.TotalVoxels = max(0,round(double(totalVoxels)));
    if isnumeric(progressEvery) && isscalar(progressEvery) && isfinite(progressEvery) && progressEvery == 0
      obj.ProgressEvery = max(500,ceil(obj.TotalVoxels/100));
    elseif isnumeric(progressEvery) && isscalar(progressEvery) && isfinite(progressEvery) && progressEvery > 0
      obj.ProgressEvery = max(1,round(double(progressEvery)));
    else
      % A negative or nonfinite value explicitly disables progress lines.
      obj.ProgressEvery = inf;
    end

    obj.DiagnosticIndices = pRFProgressReporter.selectDiagnosticIndices( ...
      obj.TotalVoxels,diagnosticVoxelCount,diagnosticSeed,scanNum);
    obj.Completed = false(1,obj.TotalVoxels);
    obj.NextProgress = obj.ProgressEvery;
    obj.DiagnosticPrinted = false(size(obj.DiagnosticIndices));
    obj.DiagnosticRecords = cell(size(obj.DiagnosticIndices));
    obj.StartTic = tic;

  end

  function record(obj,payload)
    % Called on the client for each completed voxel. Monitoring failures are
    % warning-only so they can never abort or change a fit.
    if obj.Finished,return,end
    try
      if isstruct(payload)
        voxelIndex = payload.index;
      else
        voxelIndex = payload;
      end
      wasNew = obj.markCompleted(voxelIndex);
      if isstruct(payload)
        obj.reportDiagnostic(payload);
      end
      if wasNew
        obj.maybePrintProgress;
      end
    catch reportException
      if ~obj.MonitoringWarningIssued
        obj.MonitoringWarningIssued = true;
        warning('pRF:ProgressReporterFailure', ...
          'Progress reporting failed but fitting will continue: %s',reportException.message);
      end
    end
  end

  function reconcile(obj,voxelIndices)
    % PARFOR completion guarantees that all indices in a finished block ran.
    % Reconcile protects against delayed/dropped display callbacks and makes
    % duplicate queue messages harmless.
    if obj.Finished || isempty(voxelIndices),return,end
    voxelIndices = unique(round(double(voxelIndices(:)')));
    voxelIndices = voxelIndices(voxelIndices >= 1 & voxelIndices <= obj.TotalVoxels);
    newIndices = voxelIndices(~obj.Completed(voxelIndices));
    if ~isempty(newIndices)
      obj.Completed(newIndices) = true;
      obj.CompletedCount = obj.CompletedCount + numel(newIndices);
      obj.maybePrintProgress;
    end
  end

  function finish(obj)
    if obj.Finished,return,end
    obj.reconcile(1:obj.TotalVoxels);
    if isfinite(obj.ProgressEvery) && obj.TotalVoxels > 0 && ...
        obj.LastPrintedCount < obj.TotalVoxels
      obj.printProgress(true);
    end
    try
      obj.auditDiagnostics;
    catch reportException
      if ~obj.MonitoringWarningIssued
        obj.MonitoringWarningIssued = true;
        warning('pRF:ProgressReporterFailure', ...
          'Progress reporting failed but fitting will continue: %s',reportException.message);
      end
    end
    obj.Finished = true;
  end

  function summary = getSummary(obj)
    [nDiagnostics,nValidFits,nUniqueInputs,nUniqueFits] = obj.diagnosticDiversity;
    summary.totalVoxels = obj.TotalVoxels;
    summary.progressEvery = obj.ProgressEvery;
    summary.completedCount = obj.CompletedCount;
    summary.progressReportCount = obj.ProgressReportCount;
    summary.diagnosticIndices = obj.DiagnosticIndices;
    summary.diagnosticsReceived = nDiagnostics;
    summary.validDiagnosticFits = nValidFits;
    summary.uniqueDiagnosticInputs = nUniqueInputs;
    summary.uniqueDiagnosticFits = nUniqueFits;
  end
end

methods (Access=private)
  function wasNew = markCompleted(obj,voxelIndex)
    wasNew = false;
    if ~(isnumeric(voxelIndex) && isscalar(voxelIndex) && isfinite(voxelIndex))
      return
    end
    voxelIndex = round(double(voxelIndex));
    if voxelIndex < 1 || voxelIndex > obj.TotalVoxels || obj.Completed(voxelIndex)
      return
    end
    obj.Completed(voxelIndex) = true;
    obj.CompletedCount = obj.CompletedCount + 1;
    wasNew = true;
  end

  function maybePrintProgress(obj)
    if ~isfinite(obj.ProgressEvery) || obj.CompletedCount < obj.NextProgress
      return
    end
    obj.printProgress(obj.CompletedCount >= obj.TotalVoxels);
    obj.NextProgress = (floor(obj.CompletedCount/obj.ProgressEvery)+1)*obj.ProgressEvery;
  end

  function printProgress(obj,~)
    elapsedSeconds = toc(obj.StartTic);
    if obj.CompletedCount > 0
      remainingSeconds = elapsedSeconds*(obj.TotalVoxels/obj.CompletedCount-1);
    else
      remainingSeconds = 0;
    end
    percentDone = 100*obj.CompletedCount/max(1,obj.TotalVoxels);
    progressMessage = sprintf( ...
      '(pRF) %0.1f%% done in %s (Estimated time remaining: %s)', ...
      percentDone,pRFProgressReporter.formatElapsedTime(elapsedSeconds), ...
      pRFProgressReporter.formatElapsedTime(remainingSeconds));
    pRFProgressReporter.displayHeader(progressMessage);
    obj.LastPrintedCount = obj.CompletedCount;
    obj.ProgressReportCount = obj.ProgressReportCount + 1;
  end

  function reportDiagnostic(obj,record)
    samplePosition = find(obj.DiagnosticIndices == round(double(record.index)),1);
    if isempty(samplePosition) || obj.DiagnosticPrinted(samplePosition)
      return
    end
    requiredFields = {'coords','inputFiniteCount','inputSampleCount','inputMean', ...
      'inputStd','inputMin','inputMax','inputTSeries','fitReturned','fitValid', ...
      'r2','polarAngleDegrees','eccentricity','rfHalfWidth','params','r', ...
      'fitDisplay'};
    if ~all(isfield(record,requiredFields)) || numel(record.coords) < 3
      error('pRF:InvalidDiagnosticRecord','Incomplete diagnostic voxel record.');
    end

    disp(record.fitDisplay);
    obj.DiagnosticPrinted(samplePosition) = true;
    obj.DiagnosticRecords{samplePosition} = record;
  end

  function auditDiagnostics(obj)
    [nDiagnostics,nValidFits,nUniqueInputs,nUniqueFits] = obj.diagnosticDiversity;
    if isempty(obj.DiagnosticIndices),return,end
    fprintf('(pRF QC) Received %i/%i samples: %i distinct inputs; %i distinct valid fits.\n', ...
      nDiagnostics,numel(obj.DiagnosticIndices),nUniqueInputs,nUniqueFits);

    if nDiagnostics < numel(obj.DiagnosticIndices)
      warning('pRF:MissingDiagnosticSamples', ...
        ['Only %i of %i selected diagnostic voxels could be reported. ' ...
         'The fits themselves are unaffected.'],nDiagnostics,numel(obj.DiagnosticIndices));
    end
    if nDiagnostics >= 3 && nUniqueInputs == 1
      warning('pRF:IdenticalDiagnosticInputs', ...
        ['All %i sampled voxels had bit-identical input time series. ' ...
         'Check ROI loading, coordinates, and scan selection.'],nDiagnostics);
    end
    if nValidFits >= 3 && nUniqueFits == 1
      warning('pRF:IdenticalDiagnosticFits', ...
        ['All %i valid sampled voxels returned identical parameters and r2. ' ...
         'Check the input data and voxel indexing.'],nValidFits);
    elseif nDiagnostics >= 3 && nValidFits == 0
      warning('pRF:NoValidDiagnosticFits', ...
        'None of the %i sampled voxels returned a valid fit.',nDiagnostics);
    end
  end

  function [nDiagnostics,nValidFits,nUniqueInputs,nUniqueFits] = diagnosticDiversity(obj)
    records = obj.DiagnosticRecords(~cellfun(@isempty,obj.DiagnosticRecords));
    nDiagnostics = numel(records);
    nUniqueInputs = pRFProgressReporter.countUniqueRecordField(records,'inputTSeries');
    validRecords = records(cellfun(@(record) record.fitValid,records));
    nValidFits = numel(validRecords);
    fitSignatures = cell(size(validRecords));
    for iRecord = 1:numel(validRecords)
      fitSignatures{iRecord} = [validRecords{iRecord}.params(:)' validRecords{iRecord}.r2];
    end
    nUniqueFits = pRFProgressReporter.countUniqueValues(fitSignatures);
  end
end

methods (Static)
  function indices = selectDiagnosticIndices(totalVoxels,sampleCount,seed,scanNum)
    % Select one random voxel from each equal-sized index stratum. A private
    % stream preserves the caller's global random-number state.
    if nargin < 2 || isempty(sampleCount),sampleCount = 10;end
    if nargin < 3 || isempty(seed),seed = 5489;end
    if nargin < 4 || isempty(scanNum),scanNum = 1;end
    totalVoxels = max(0,round(double(totalVoxels)));
    if ~(isnumeric(sampleCount) && isscalar(sampleCount) && isfinite(sampleCount))
      sampleCount = 0;
    end
    sampleCount = min(totalVoxels,max(0,round(double(sampleCount))));
    if totalVoxels == 0 || sampleCount == 0
      indices = zeros(1,0);
      return
    end
    if ~(isnumeric(seed) && isscalar(seed) && isfinite(seed)),seed = 5489;end
    if ~(isnumeric(scanNum) && isscalar(scanNum) && isfinite(scanNum)),scanNum = 1;end
    streamSeed = mod(floor(abs(double(seed))) + 104729*floor(abs(double(scanNum))),2^32);
    localStream = RandStream('mt19937ar','Seed',streamSeed);
    edges = floor((0:sampleCount)*totalVoxels/sampleCount);
    indices = zeros(1,sampleCount);
    for iSample = 1:sampleCount
      firstIndex = edges(iSample)+1;
      lastIndex = edges(iSample+1);
      indices(iSample) = firstIndex + floor(rand(localStream,1)*(lastIndex-firstIndex+1));
    end
  end

  function record = makeDiagnosticRecord(index,totalVoxels,coords,tSeries,fit)
    % Constructed only for the small diagnostic sample, so retaining the
    % complete input time series is inexpensive and permits an exact check.
    record.kind = 'pRFDiagnosticVoxel';
    record.index = round(double(index));
    record.totalVoxels = round(double(totalVoxels));
    record.coords = round(double(coords(:)'));
    record.inputTSeries = double(tSeries(:)');
    finiteInput = record.inputTSeries(isfinite(record.inputTSeries));
    record.inputSampleCount = numel(record.inputTSeries);
    record.inputFiniteCount = numel(finiteInput);
    if isempty(finiteInput)
      record.inputMean = NaN;
      record.inputStd = NaN;
      record.inputMin = NaN;
      record.inputMax = NaN;
    else
      record.inputMean = mean(finiteInput);
      record.inputStd = std(finiteInput);
      record.inputMin = min(finiteInput);
      record.inputMax = max(finiteInput);
    end

    record.fitReturned = ~isempty(fit);
    record.fitValid = false;
    record.r2 = NaN;
    record.polarAngleDegrees = NaN;
    record.eccentricity = NaN;
    record.rfHalfWidth = NaN;
    record.params = [];
    record.r = [];
    if record.fitReturned
      record.r2 = pRFProgressReporter.fitField(fit,'r2',NaN);
      record.polarAngleDegrees = 180/pi*pRFProgressReporter.fitField(fit,'polarAngle',NaN);
      record.eccentricity = pRFProgressReporter.fitField(fit,'eccentricity',NaN);
      record.rfHalfWidth = pRFProgressReporter.fitField(fit,'std',NaN);
      record.params = pRFProgressReporter.fitField(fit,'params',[]);
      record.params = double(record.params(:)');
      record.r = pRFProgressReporter.fitField(fit,'r',[]);
      record.r = double(record.r(:)');
      record.fitValid = isnumeric(record.r2) && isscalar(record.r2) && ...
        isfinite(record.r2) && ~isempty(record.params) && ...
        all(isfinite(record.params));
    end

    prefitOnly = record.fitReturned && isstruct(fit) && isfield(fit,'bestFitVoxel');
    prefitOnlyString = '';
    if prefitOnly,prefitOnlyString = ' (prefit only)';end
    displayPrefix = sprintf(sprintf('Voxel %%%i.f/%%i%%s: ', ...
      length(sprintf('%i',record.totalVoxels))),record.index, ...
      record.totalVoxels,prefitOnlyString);
    record.fitDisplay = pRFFormatVoxelFit(fit,record.coords,displayPrefix, ...
      pRFProgressReporter.fitField(fit,'rfType',''),prefitOnly);
  end
end

methods (Static,Access=private)
  function value = fitField(fit,fieldName,defaultValue)
    if isstruct(fit) && isfield(fit,fieldName)
      value = fit.(fieldName);
    else
      value = defaultValue;
    end
  end

  function nUnique = countUniqueRecordField(records,fieldName)
    values = cell(size(records));
    for iRecord = 1:numel(records)
      values{iRecord} = records{iRecord}.(fieldName);
    end
    nUnique = pRFProgressReporter.countUniqueValues(values);
  end

  function nUnique = countUniqueValues(values)
    uniqueValues = cell(1,0);
    for iValue = 1:numel(values)
      isNew = true;
      for iUnique = 1:numel(uniqueValues)
        if isequaln(values{iValue},uniqueValues{iUnique})
          isNew = false;
          break
        end
      end
      if isNew
        uniqueValues{end+1} = values{iValue}; %#ok<AGROW>
      end
    end
    nUnique = numel(uniqueValues);
  end

  function durationString = formatElapsedTime(seconds)
    seconds = max(0,double(seconds));
    if exist('mlrDispElapsedTime','file') == 2
      durationString = mlrDispElapsedTime(seconds);
      return
    end
    hours = floor(seconds/3600);
    minutes = floor((seconds-hours*3600)/60);
    wholeSeconds = floor(seconds-hours*3600-minutes*60);
    milliseconds = floor((seconds-hours*3600-minutes*60-wholeSeconds)*1000);
    durationString = '';
    if hours > 0,durationString = sprintf('%i hours ',hours);end
    if minutes > 0,durationString = sprintf('%s%i min ',durationString,minutes);end
    if wholeSeconds > 0,durationString = sprintf('%s%i sec ',durationString,wholeSeconds);end
    durationString = sprintf('%s%i ms',durationString,milliseconds);
  end

  function displayHeader(message)
    if exist('dispHeader','file') == 2
      dispHeader(message);
    elseif length(message)+2 >= 60
      disp(repmat('=',1,60));
      disp(message);
      disp(repmat('=',1,60));
    else
      fillerLength = (60-(length(message)+2))/2;
      disp(sprintf('%s %s %s',repmat('=',1,floor(fillerLength)), ...
        message,repmat('=',1,ceil(fillerLength))));
    end
  end
end

end
