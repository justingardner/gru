% pRFFit (v3)
%
%      usage: pRFFit(v,scanNum,x,y,s,<dispFit=0>)
%         by: justin gardner
%       date: 11/14/2011
%    purpose: interrogator that fits pRF model to selected voxel
%
function fit = pRFFit(varargin)

fit = [];
preparedCall = false;
preparedContext = [];

% Private entry points used by pRFFitPrepared. The public calling convention
% below is left intact, while scan-invariant parsing and setup can be done
% once before entering a voxel loop.
if ~isempty(varargin) && (ischar(varargin{1}) || (isstring(varargin{1}) && isscalar(varargin{1})))
  preparedCommand = char(varargin{1});
  if strcmp(preparedCommand,'__prepareContext__')
    fit = pRFFit(varargin{2},varargin{3},[],[],[],varargin{4:end}, ...
      'returnPrefit',1,'returnPreparedContext',1);
    if ~isempty(fit) && isfield(fit,'preparedContext')
      fit = fit.preparedContext;
    else
      fit = [];
    end
    return
  elseif strcmp(preparedCommand,'__fitPrepared__')
    preparedContext = varargin{2};
    validatePreparedContext(preparedContext);
    tSeries = varargin{3};
    x = varargin{4};
    y = varargin{5};
    z = varargin{6};
    if length(varargin) >= 7, dispIndex = varargin{7};else,dispIndex = [];end
    if length(varargin) >= 8, dispN = varargin{8};else,dispN = [];end
    v = [];
    scanNum = preparedContext.scanNum;
    fitParams = preparedContext.fitParams;
    fitParams.returnPrefit = false;
    fitParams.returnPreparedContext = false;
    fitParams.getModelResponse = false;
    fitParams.dispstr = makeVoxelDisplayString(fitParams,dispIndex,dispN);
    preparedCall = true;
  end
end

% parse input arguments - note that this is set
% up so that it can also be called as an interrogator
if ~preparedCall
  [v scanNum x y z fitParams tSeries] = parseArgs(varargin);
  if isempty(v),return,end
end
 
% get concat info
if ~preparedCall && (~isfield(fitParams,'concatInfo') || isempty(fitParams.concatInfo))
  fitParams.concatInfo = viewGet(v,'concatInfo',scanNum);
end

% if there is no concatInfo, then make one that will
% treat the scan as a single scan
if ~preparedCall && isempty(fitParams.concatInfo)
  nFrames = viewGet(v,'nFrames',scanNum);
  fitParams.concatInfo.isConcat = false;
  fitParams.concatInfo.n = 1;
  fitParams.concatInfo.whichScan = ones(1,nFrames);
  fitParams.concatInfo.whichVolume = 1:nFrames;
  fitParams.concatInfo.runTransition = [1 nFrames];
  fitParams.concatInfo.junkFrames = viewGet(v,'junkFrames');
  fitParams.concatInfo.totalJunkedFrames = viewGet(v,'totalJunkedFrames',scanNum);
  if numel(fitParams.concatInfo.totalJunkedFrames) > 1
    % first check for consistency in totalJunkedFrames
    if length(unique(fitParams.concatInfo.totalJunkedFrames)) > 1
      disp(sprintf('(pRFFit) totalJunkedFrames are different for different members of component scans - could be an average in which different scans with different number of junked frames were removed. This could cause a problem in computing what the stimulus was for the average. The total junked frames count was: %s, but we will use %i as the actual value for computing the stimulus',num2str(fitParams.concatInfo.totalJunkedFrames),floor(median(fitParams.concatInfo.totalJunkedFrames))));
    end
    fitParams.concatInfo.totalJunkedFrames = floor(median(fitParams.concatInfo.totalJunkedFrames));
  end
elseif ~preparedCall
  fitParams.concatInfo.isConcat = true;
  if ~isfield(fitParams.concatInfo,'totalJunkedFrames')
    fitParams.concatInfo.totalJunkedFrames = viewGet(v,'totalJunkedFrames',scanNum);
  end
end

% get the stimulus movie if it wasn't passed in
if ~preparedCall && (~isfield(fitParams,'stim') || isempty(fitParams.stim))
  fitParams.stim = getStim(v,scanNum,fitParams);
end
if isempty(fitParams.stim),return,end

% if we are being called to just return the stim image 
% then return it here
if fitParams.justGetStimImage
  fit = fitParams.stim;
  return
end
% get the tSeries
if ~isempty(x)
  % if tSeries was not passed in then load it
  if isempty(tSeries)
    % load using loadTSeries
    tSeries = squeeze(loadTSeries(v,scanNum,z,[],x,y));
  end

  % convert to percent tSeries. Note that we  detrend here which is not necessary for concats,
  % but useful for raw/motionCorrected time series. Also, it is very important that
  % the tSeries is properly mean subtracted
  if ~isfield(fitParams.concatInfo,'hipassfilter')
    tSeries = percentTSeries(tSeries,'detrend','Linear', 'spatialNormalization','Divide by mean','subtractMean', 'Yes', 'temporalNormalization', 'No');
  end
  
  % if there are any nans in the tSeries then don't fit
  if any(isnan(tSeries))
    if fitParams.verbose
      if preparedCall
        groupName = preparedContext.groupName;
      else
        groupName = viewGet(v,'groupName');
      end
      disp(sprintf('(pRFFit) Nan found in tSeries for voxel [%i %i %i] in scan %s:%i. Abandoning fit',x,y,z,groupName,scanNum));
    end
    fit=[];return
  end
else
  tSeries = [];
end
% handle junk frames (i.e. ones that have not already been junked)
if ~isempty(fitParams.junkFrames) && ~isequal(fitParams.junkFrames,0)
  % drop junk frames
  disp(sprintf('(pRFFit) Dropping %i junk frames',fitParams.junkFrames));
  tSeries = tSeries(fitParams.junkFrames+1:end);
  if ~isfield(fitParams.concatInfo,'totalJunkedFramesIncludesJunked');
    fitParams.concatInfo.totalJunkedFrames = fitParams.concatInfo.totalJunkedFrames+fitParams.junkFrames;
    fitParams.concatInfo.totalJunkedFramesIncludesJunked = 1;
  end
end

% set up the scan-invariant fit routine params only once. Prepared calls
% receive the exact post-setup structure from pRFFitPrepared.
if ~preparedCall
  fitParams = setFitParams(fitParams);
end

% just return model response for already calculated params
if fitParams.getModelResponse
  % get model fit
  [residual fit.modelResponse fit.rfModel] = getModelResidual(fitParams.params,tSeries,fitParams);
  % get the canonical
  fit.p = getFitParams(fitParams.params,fitParams);
  fit.canonical = getCanonicalHRF(fit.p.canonical,fitParams.framePeriod);
  % return tSeries
  fit.tSeries = tSeries;
  fit.preparedContextUsed = preparedCall;
  fit.objectiveCacheHits = 0;
  fit.objectiveCacheMisses = 0;
  return;
end

% return some fields
fit.stim = fitParams.stim;
fit.stimX = fitParams.stimX;
fit.stimY = fitParams.stimY;
fit.stimT = fitParams.stimT;
fit.concatInfo = fitParams.concatInfo;
fit.nParams = fitParams.nParams;
paramsInfoFields = {'minParams','maxParams','initParams','paramNames','paramDescriptions', ...
  'paramIncDec','paramMin','paramMax', ...
  'stimExtents','stimWidth','stimHeight','stimDeltaT','fixedHRF','optimParams'};
for iField = 1:length(paramsInfoFields)
  if isfield(fitParams,paramsInfoFields{iField})
    fit.paramsInfo.(paramsInfoFields{iField}) = fitParams.(paramsInfoFields{iField});
  end
end

% Warn when the saved stimulus and data lengths differ. The model response
% is deliberately resampled below, so this is diagnostic rather than fatal.
for iScan = 1:fit.concatInfo.n
  sLength = fit.concatInfo.runTransition(iScan,2) - fit.concatInfo.runTransition(iScan,1) + 1;
  if sLength ~= size(fitParams.stim{iScan}.im,3)
    mrWarnDlg(sprintf('(pRFFit) Data length of %i for scan %i (concatNum:%i) does not match stimfile length %i',sLength,scanNum,iScan,size(fitParams.stim{iScan}.im,3)));
  end
end

% do prefit. This computes (or is passed in precomputed) model responses
% for a variety of parameters and calculates the correlation between
% the models and the time series. The one that has the best correlation
% is then used as the initial parameters for the nonlinear fit. This
% helps prevent getting stuck in local minima
if ~isfield(fitParams,'prefitCacheHit'),fitParams.prefitCacheHit = false;end
if ~isfield(fitParams,'prefitCacheKey'),fitParams.prefitCacheKey = '';end
fitParams.prefitBankReused = false;
if isfield(fitParams,'prefit') && ~isempty(fitParams.prefit)
  params = fitParams.initParams;
  % calculate model if not already calculated
  if ~isfield(fitParams.prefit,'modelResponse')
    % first convert the x/y and width parameters into sizes
    % on the actual screen
    fitParams.prefit.x = fitParams.prefit.x*fitParams.stimWidth;
    fitParams.prefit.y = fitParams.prefit.y*fitParams.stimHeight;
    fitParams.prefit.rfHalfWidth = fitParams.prefit.rfHalfWidth*max(fitParams.stimWidth,fitParams.stimHeight);
    % The bank depends only on exact stimulus/model/temporal values, never on
    % the voxel time series. Keep a small process-local cache so repeated
    % analyses of an identical scan can reuse it safely.
    fitParams.prefitCacheKey = makePrefitCacheKey(fitParams,params(4:end));
    [cachedModelResponse fitParams.prefitCacheHit] = ...
      prefitBankCache('get',fitParams.prefitCacheKey);
    if fitParams.prefitCacheHit
      fitParams.prefit.modelResponse = cachedModelResponse;
    else
      % Use the same explicit setting as voxel fitting. In particular,
      % numWorkers=1 takes the ordinary FOR branch without inspecting or
      % starting a parallel pool.
      if fitParams.prefit.n > 1
        nProcessors = pRFResolveNumWorkers(fitParams.numWorkers);
        fitParams.numWorkers = nProcessors;
      else
        nProcessors = 1;
      end
      if fitParams.verbose==1
        disppercent(-inf,sprintf('(pRFFit) Computing %i prefit model responses using %i processors',fitParams.prefit.n,nProcessors));
      end
      % init modelResponse
      allModelResponse = nan(fitParams.prefit.n,fitParams.concatInfo.runTransition(end,end));
      % compute all the model response, using parfor loop
      if nProcessors>1
          parfor i = 1:fitParams.prefit.n
            % fit the model with these parameters
            [residual modelResponse rfModel] = getModelResidual([fitParams.prefit.x(i) fitParams.prefit.y(i) fitParams.prefit.rfHalfWidth(i) params(4:end)],tSeries,fitParams,1);
            % normalize to 0 mean unit length
            allModelResponse(i,:) = (modelResponse-mean(modelResponse))./sqrt(sum(modelResponse.^2))';
            if fitParams.verbose
          disp(sprintf('(pRFFit) Computing prefit model response %i/%i: Center [%6.2f,%6.2f] rfHalfWidth=%5.2f',i,fitParams.prefit.n,fitParams.prefit.x(i),fitParams.prefit.y(i),fitParams.prefit.rfHalfWidth(i)));
            end
          end
      else
          for i = 1:fitParams.prefit.n
            % fit the model with these parameters
            [residual modelResponse rfModel] = getModelResidual([fitParams.prefit.x(i) fitParams.prefit.y(i) fitParams.prefit.rfHalfWidth(i) params(4:end)],tSeries,fitParams,1);
            % normalize to 0 mean unit length
            allModelResponse(i,:) = (modelResponse-mean(modelResponse))./sqrt(sum(modelResponse.^2))';
            if fitParams.verbose
          disp(sprintf('(pRFFit) Computing prefit model response %i/%i: Center [%6.2f,%6.2f] rfHalfWidth=%5.2f',i,fitParams.prefit.n,fitParams.prefit.x(i),fitParams.prefit.y(i),fitParams.prefit.rfHalfWidth(i)));
            end
            if fitParams.verbose,disppercent(i/fitParams.prefit.n);end
          end
      end
      if fitParams.verbose,disppercent(inf);end
      fitParams.prefit.modelResponse = allModelResponse;
      prefitBankCache('put',fitParams.prefitCacheKey,allModelResponse);
      clear allModelResponse;
    end
  else
    fitParams.prefitBankReused = true;
  end
  % save in global, so that when called as an interrogator
  % we don't have to keep computing fitParams
  if ~preparedCall
    global gpRFFitTypeParams
    gpRFFitTypeParams.prefit = fitParams.prefit;
  end
  % return some computed fields
  fit.prefit = fitParams.prefit;
  fit.prefitCacheHit = fitParams.prefitCacheHit;
  fit.prefitBankReused = fitParams.prefitBankReused;
  if fitParams.returnPrefit
    if fitParams.returnPreparedContext
      fit.preparedContext = makePreparedContext(v,scanNum,fitParams);
    end
    return
  end
  % normalize tSeries to 0 mean unit length
  tSeriesNorm = (tSeries-mean(tSeries))/sqrt(sum(tSeries.^2));
  % calculate r for all modelResponse by taking inner product
  r = fitParams.prefit.modelResponse*tSeriesNorm;
  % get best r2 for all the models
  [maxr bestModel] = max(r);
  fitParams.initParams(1) = fitParams.prefit.x(bestModel);
  fitParams.initParams(2) = fitParams.prefit.y(bestModel);
  fitParams.initParams(3) = fitParams.prefit.rfHalfWidth(bestModel);
  if fitParams.prefitOnly
    % return if we are just doing a prefit
    fit = getFitParams(fitParams.initParams,fitParams);
    fit.modelResponse = fitParams.prefit.modelResponse;
    fit.tSeries = tSeries;
    fit.rfType = fitParams.rfType;
    fit.params = fitParams.initParams;
    fit.r2 = maxr^2;
    fit.r = maxr;
    fit.bestFitVoxel = bestModel;
    [fit.polarAngle fit.eccentricity] = cart2pol(fit.x,fit.y);
    % return canonical
    fit.canonicalModel = getCanonicalHRF(fit.canonical,fitParams.framePeriod);
    fit.preparedContextUsed = preparedCall;
    fit.prefitCacheHit = fitParams.prefitCacheHit;
    fit.objectiveCacheHits = 0;
    fit.objectiveCacheMisses = 0;
    % display
    if fitParams.verbose
      disp(pRFFormatVoxelFit(fit,[x y z],fitParams.dispstr,fitParams.rfType,true));
    end
    return
  end
end

% now do nonlinear fit
validateOptimizationParameters(fitParams);
objectiveCacheId = objectiveMemo('create');
objectiveCacheCleanup = onCleanup(@() objectiveMemo('destroy',objectiveCacheId));
if strcmp(lower(fitParams.algorithm),'levenberg-marquardt')
  objective = @(candidateParams) objectiveMemo('evaluate',objectiveCacheId, ...
    candidateParams,tSeries,fitParams,false);
  [params resnorm residual exitflag output lambda jacobian] = lsqnonlin(objective,fitParams.initParams,fitParams.minParams,fitParams.maxParams,fitParams.optimParams);
elseif strcmp(lower(fitParams.algorithm),'nelder-mead')
  objectiveTSeries = (tSeries-mean(tSeries))/var(tSeries.^2);
  objective = @(candidateParams) objectiveMemo('evaluate',objectiveCacheId, ...
    candidateParams,objectiveTSeries,fitParams,false);
  [params fval exitflag output] = fminsearch(objective,fitParams.initParams,fitParams.optimParams);
elseif strcmp(lower(fitParams.algorithm), 'nelder-mead-bnd')
  objectiveTSeries = (tSeries-mean(tSeries))/var(tSeries.^2);
  objective = @(candidateParams) objectiveMemo('evaluate',objectiveCacheId, ...
    candidateParams,objectiveTSeries,fitParams,false);
  [params fval exitflag output] = fminsearchbnd(objective, fitParams.initParams, fitParams.minParams, fitParams.maxParams, fitParams.optimParams);
elseif strcmp(lower(fitParams.algorithm), 'fmincon')
  % fmincon requires a feasible initial point. This affects only the
  % explicitly selected experimental solver; the legacy solvers retain their
  % original initialization and objective paths.
  fminconMinParams = fitParams.minParams;
  degenerateZeroBounds = (fminconMinParams == 0) & (fitParams.maxParams > 0);
  fminconMinParams(degenerateZeroBounds) = sqrt(eps);
  fminconInit = min(max(fitParams.initParams,fminconMinParams),fitParams.maxParams);
  % Use the same time-series scale for optimization and final reporting.
  % The historical Nelder-Mead normalization is retained in its own branch,
  % but it can make degenerate multi-component models appear finite during
  % optimization and NaN when evaluated on the original data scale.
  fminconTSeries = tSeries;
  [fminconInit fminconInitObjective] = makeFiniteFminconInitial( ...
    fminconInit,fminconTSeries,fitParams,fminconMinParams,fitParams.maxParams);
  objective = @(candidateParams) objectiveMemo('evaluate',objectiveCacheId, ...
    candidateParams,fminconTSeries,fitParams,true);
  [params fval exitflag output] = fmincon(objective,fminconInit,[],[],[],[],fminconMinParams,fitParams.maxParams,[],fitParams.optimParams);
  finalFminconObjective = getModelResidual(params,fminconTSeries,fitParams);
  if ~isscalar(finalFminconObjective) || ~isfinite(finalFminconObjective)
    % A singular boundary can occasionally be accepted by SQP after its
    % finite penalty is applied. Return the best finite starting point rather
    % than emitting a NaN fit.
    params = fminconInit;
    fval = fminconInitObjective;
    exitflag = 0;
    output.usedFiniteInitialFallback = true;
  else
    output.usedFiniteInitialFallback = false;
  end
else
  error('pRFFit:UnknownAlgorithm','Unknown optimization algorithm: %s',fitParams.algorithm);
end
objectiveCacheStats = objectiveMemo('stats',objectiveCacheId);

% set output arguments
fit = getFitParams(params,fitParams);
fit.rfType = fitParams.rfType;
fit.params = params;
fit.exitflag = exitflag;
fit.optimOutput = output;
if isstruct(output) && isfield(output,'funcCount')
  fit.funcCount = output.funcCount;
else
  fit.funcCount = [];
end
fit.objectiveCacheHits = objectiveCacheStats.hits;
fit.objectiveCacheMisses = objectiveCacheStats.misses;
fit.preparedContextUsed = preparedCall;
fit.prefitCacheHit = fitParams.prefitCacheHit;
fit.prefitBankReused = fitParams.prefitBankReused;

% get the canonical
fit.canonicalModel = getCanonicalHRF(fit.canonical,fitParams.framePeriod);

% compute r^2
[residual modelResponse rfModel fit.r] = getModelResidual(params,tSeries,fitParams);
if strcmp(lower(fitParams.algorithm),'levenberg-marquardt')
  fit.r2 = 1-sum((residual-mean(residual)).^2)/sum((tSeries-mean(tSeries)).^2);
elseif strcmp(lower(fitParams.algorithm),'nelder-mead')
  fit.r2 = residual^2;
elseif strcmp(lower(fitParams.algorithm), 'nelder-mead-bnd')
  fit.r2 = residual^2;
elseif strcmp(lower(fitParams.algorithm), 'fmincon')
  fit.r2 = residual^2;
end

% compute polar coordinates
[fit.polarAngle fit.eccentricity] = cart2pol(fit.x,fit.y);

% display
if fitParams.verbose
  disp(pRFFormatVoxelFit(fit,[x y z],fitParams.dispstr,fitParams.rfType,false));
end
%%%%%%%%%%%%%%%%%%%%%%
%    setFitParams    %
%%%%%%%%%%%%%%%%%%%%%%
function fitParams = setFitParams(fitParams);

% A prepared scan context is immutable for the lifetime of the voxel loop.
% Build its exact duplicate-frame/run metadata once, before the prefit bank
% and nonlinear fits use the stimulus. Ordinary interrogator calls retain
% their original, unprepared stimulus representation.
if isfield(fitParams,'returnPreparedContext') && fitParams.returnPreparedContext && ...
    exist('pRFPrepareStimulusProjection','file') == 2
  fitParams.stim = pRFPrepareStimulusProjection(fitParams.stim);
end

% set rfType
if ~isfield(fitParams,'rfType') || isempty(fitParams.rfType)
  fitParams.rfType = 'gaussian';
end

% The selected model and constraint UI both need to know the algorithm.
% Older saved parameter structures may not contain it.
if ~isfield(fitParams,'algorithm') || isempty(fitParams.algorithm)
  fitParams.algorithm = 'nelder-mead-bnd';
  disp('(pRFFit) No algorithm provided. Using Default: Nelder-Mead-Bnd');
end

% get stimulus x,y and t
fitParams.stimX = fitParams.stim{1}.x;
fitParams.stimY = fitParams.stim{1}.y;
fitParams.stimT = fitParams.stim{1}.t;
fitParams.gaussianGrid = prepareGaussianGrid(fitParams.stimX,fitParams.stimY);

% set stimulus extents (or reuse the scan-invariant values passed in paramsInfo)
if ~isfield(fitParams,'stimExtents') || isempty(fitParams.stimExtents)
  fitParams.stimExtents(1) = min(fitParams.stimX(:));
  fitParams.stimExtents(3) = max(fitParams.stimX(:));
  fitParams.stimExtents(2) = min(fitParams.stimY(:));
  fitParams.stimExtents(4) = max(fitParams.stimY(:));
end
if ~isfield(fitParams,'stimWidth') || isempty(fitParams.stimWidth)
  fitParams.stimWidth = fitParams.stimExtents(3)-fitParams.stimExtents(1);
end
if ~isfield(fitParams,'stimHeight') || isempty(fitParams.stimHeight)
  fitParams.stimHeight = fitParams.stimExtents(4)-fitParams.stimExtents(2);
end

if ~isfield(fitParams,'initParams')
  % check the rfType to get the correct min/max arrays

  %%%%%%%%%%%%%%%%%%%
  fitParams = prfModel('setParams', fitParams);
  %%%%%%%%%%%%%%%%%%%
  validateModelParameterMetadata(fitParams);
  
  % round constraints
  fitParams.minParams = round(fitParams.minParams*10)/10;
  fitParams.maxParams = round(fitParams.maxParams*10)/10;

  % handle constraints here
  % Check if fit algorithm is one that allows constraints
  algorithmsWithConstraints = {'levenberg-marquardt', 'nelder-mead-bnd', 'fmincon'};
  if any(strcmp(fitParams.algorithm,algorithmsWithConstraints))
    % if constraints allowed then allow user to adjust them here (if they set defaultConstraints)
    if isfield(fitParams,'defaultConstraints') && ~fitParams.defaultConstraints
      % create a dialog to allow user to set constraints
      paramsInfo = {};
      for iParam = 1:length(fitParams.paramNames)
	paramsInfo{end+1} = {sprintf('min%s',fitParams.paramNames{iParam}) fitParams.minParams(iParam) sprintf('Minimum for parameter %s (%s)',fitParams.paramNames{iParam},fitParams.paramDescriptions{iParam}) sprintf('incdec=[%f %f]',-fitParams.paramIncDec(iParam),fitParams.paramIncDec(iParam)) sprintf('minmax=[%f %f]',fitParams.paramMin(iParam),fitParams.paramMax(iParam))};
	paramsInfo{end+1} = {sprintf('max%s',fitParams.paramNames{iParam}) fitParams.maxParams(iParam) sprintf('Maximum for parameter %s (%s)',fitParams.paramNames{iParam},fitParams.paramDescriptions{iParam})  sprintf('incdec=[%f %f]',-fitParams.paramIncDec(iParam),fitParams.paramIncDec(iParam)) sprintf('minmax=[%f %f]',fitParams.paramMin(iParam),fitParams.paramMax(iParam))};
      end
      params = mrParamsDialog(paramsInfo,'Set parameter constraints');
      % if params is not empty then set them
      if isempty(params)
	disp(sprintf('(pRFFit) Using default constraints'));
      else
	% get the parameter constraints back from the dialog entries
	for iParam = 1:length(fitParams.paramNames)
	  fitParams.minParams(iParam) = params.(sprintf('min%s',fitParams.paramNames{iParam}));
	  fitParams.maxParams(iParam) = params.(sprintf('max%s',fitParams.paramNames{iParam}));
	end
      end
    end
    % Now display parameter constraints when detailed logging is requested.
    if fitParams.verbose
      for iParam = 1:length(fitParams.paramNames)
        disp(sprintf('(pRFFit) Parameter %s [min:%f max:%f] (%i:%s)',fitParams.paramNames{iParam},fitParams.minParams(iParam),fitParams.maxParams(iParam),iParam,fitParams.paramDescriptions{iParam}));
      end
    end
  else
    % no constraints allowed
    if fitParams.verbose==1
      disp(sprintf('(pRFFit) !!! Fit constraints ignored for algorithm: %s (if you want to constrain the fits, then use: %s) !!!',fitParams.algorithm,cell2mat(algorithmsWithConstraints)));
    end
  end
end

fitParams.nParams = length(fitParams.initParams);
validateFitParameterVectors(fitParams);
validateExperimentalModelConfiguration(fitParams);

% optimization parameters
if strcmpi(fitParams.algorithm,'fmincon') && exist('fmincon','file') ~= 2
  error('pRFFit:FminconUnavailable','The optional fmincon algorithm requires MATLAB Optimization Toolbox.');
end
if ~isfield(fitParams,'optimParams') || isempty(fitParams.optimParams)
  if strcmpi(fitParams.algorithm,'fmincon')
    fitParams.optimParams = optimoptions('fmincon', ...
      'Algorithm','sqp', ...
      'MaxIterations',inf, ...
      'Display',fitParams.optimDisplay, ...
      'FiniteDifferenceStepSize',1e-2);
  else
    fitParams.optimParams = optimset('MaxIter',inf,'Display',fitParams.optimDisplay);
  end
end

% Cache temporal quantities that are fixed for every objective evaluation.
if ~isfield(fitParams,'stimDeltaT') || isempty(fitParams.stimDeltaT)
  fitParams.stimDeltaT = median(diff(fitParams.stimT));
end
fixedHRFTypes = {'gaussian','gaussian-css','gaussian-diffs','gaussian-divs','gaussian-DoG-CSS','divNorm'};
if any(strcmp(fitParams.rfType,fixedHRFTypes)) && (~isfield(fitParams,'fixedHRF') || isempty(fitParams.fixedHRF))
  fixedP = getFitParams(fitParams.initParams,fitParams);
  fitParams.fixedHRF = getCanonicalHRF(fixedP.canonical,fitParams.stimDeltaT);
end

% parameters for converting the stimulus
params = {'xFlipStimulus','yFlipStimulus','timeShiftStimulus'};
for i = 1:length(params)
  if ~isfield(fitParams,params{i}) || isempty(fitParams.(params{i}))
    fitParams.(params{i}) = 0;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    validateModelParameterMetadata    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function validateModelParameterMetadata(fitParams)

requiredFields = {'paramNames','paramDescriptions','paramIncDec', ...
  'paramMin','paramMax','minParams','maxParams','initParams'};
missingFields = requiredFields(~isfield(fitParams,requiredFields));
if ~isempty(missingFields)
  error('pRFFit:InvalidModelMetadata', ...
    'Model %s did not define required parameter metadata: %s.', ...
    fitParams.rfType,strjoin(missingFields,', '));
end

nParams = numel(fitParams.initParams);
for iField = 1:numel(requiredFields)
  fieldName = requiredFields{iField};
  if numel(fitParams.(fieldName)) ~= nParams
    error('pRFFit:InvalidModelMetadata', ...
      'Model %s field %s has %i entries; expected %i.', ...
      fitParams.rfType,fieldName,numel(fitParams.(fieldName)),nParams);
  end
end
if ~iscellstr(fitParams.paramNames) || ~iscellstr(fitParams.paramDescriptions)
  error('pRFFit:InvalidModelMetadata', ...
    'Model %s parameter names and descriptions must be cell arrays of character vectors.', ...
    fitParams.rfType);
end
numericFields = {'paramIncDec','paramMin','paramMax','minParams','maxParams','initParams'};
for iField = 1:numel(numericFields)
  value = fitParams.(numericFields{iField});
  if ~isnumeric(value) || ~isreal(value) || ~isvector(value) || any(isnan(value(:)))
    error('pRFFit:InvalidModelMetadata', ...
      'Model %s field %s must be a real numeric vector without NaNs.', ...
      fitParams.rfType,numericFields{iField});
  end
end
if any(fitParams.paramMin(:) > fitParams.paramMax(:))
  error('pRFFit:InvalidModelMetadata', ...
    'Model %s has a parameter UI minimum greater than its maximum.',fitParams.rfType);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    validateFitParameterVectors    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function validateFitParameterVectors(fitParams)

requiredFields = {'initParams','minParams','maxParams'};
missingFields = requiredFields(~isfield(fitParams,requiredFields));
if ~isempty(missingFields)
  error('pRFFit:InvalidParameterBounds', ...
    'Missing fit parameter fields: %s.',strjoin(missingFields,', '));
end
nParams = numel(fitParams.initParams);
for iField = 1:numel(requiredFields)
  fieldName = requiredFields{iField};
  value = fitParams.(fieldName);
  if ~isnumeric(value) || ~isreal(value) || ~isvector(value) || ...
      numel(value) ~= nParams || any(isnan(value(:)))
    error('pRFFit:InvalidParameterBounds', ...
      'Fit field %s must be a real %i-element numeric vector without NaNs.', ...
      fieldName,nParams);
  end
end
if any(fitParams.minParams(:) > fitParams.maxParams(:))
  error('pRFFit:InvalidParameterBounds', ...
    'At least one fit parameter minimum is greater than its maximum.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    validateOptimizationParameters  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function validateOptimizationParameters(fitParams)

validateFitParameterVectors(fitParams);
boundedAlgorithms = {'levenberg-marquardt','nelder-mead-bnd','fmincon'};
if any(strcmpi(fitParams.algorithm,boundedAlgorithms)) && ...
    (any(fitParams.initParams(:) < fitParams.minParams(:)) || ...
     any(fitParams.initParams(:) > fitParams.maxParams(:)))
  error('pRFFit:InitialParametersOutsideBounds', ...
    'Initial parameters for %s must lie within the selected model bounds.', ...
    fitParams.algorithm);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    validateExperimentalModelConfiguration    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function validateExperimentalModelConfiguration(fitParams)

experimentalModels = {'gaussian-DoG-CSS','divNorm'};
if ~any(strcmp(fitParams.rfType,experimentalModels))
  return
end

% Response-only calls are used to plot saved analyses. Their actual model
% parameters are validated by the model itself, so historical bound metadata
% does not prevent a valid saved fit from being displayed.
if (isfield(fitParams,'getModelResponse') && fitParams.getModelResponse) || ...
    (isfield(fitParams,'prefitOnly') && fitParams.prefitOnly)
  return
end

if strcmpi(fitParams.algorithm,'nelder-mead')
  error('pRFFit:ExperimentalModelRequiresBoundedAlgorithm', ...
    ['Model %s has essential sign and positivity constraints. Select ' ...
     'nelder-mead-bnd, levenberg-marquardt, or fmincon instead of unconstrained nelder-mead.'], ...
    fitParams.rfType);
end

switch fitParams.rfType
  case 'gaussian-DoG-CSS'
    safeBounds = (fitParams.minParams(3) > 0) && ...
      (fitParams.minParams(4) > 0) && (fitParams.minParams(5) >= 0) && ...
      (fitParams.maxParams(6) <= 0) && (fitParams.minParams(7) > 0);
  case 'divNorm'
    safeBounds = (fitParams.minParams(3) > 0) && ...
      (fitParams.minParams(4) > 0) && all(fitParams.minParams(5:7) >= 0) && ...
      (fitParams.minParams(8) > 0);
end
if ~safeBounds
  error('pRFFit:UnsafeExperimentalModelBounds', ...
    ['Model %s requires positive widths, its documented amplitude signs, ' ...
     'and a positive nonlinear exponent or divisive constant.'],fitParams.rfType);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getModelResidual   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [residual modelResponse rfModel r] = getModelResidual(params,tSeries,fitParams,justGetModel)

if nargin < 4, justGetModel = 0;end

% get the model response
% convert parameter array into a parameter strucutre
p = getFitParams(params,fitParams);

% compute an RF
rfModel = getRFModel(p,fitParams);

% Keep per-run responses separate until the loop is complete. This avoids
% repeatedly reallocating the growing arrays while retaining their order.
modelResponseParts = cell(fitParams.concatInfo.n,1);
residualParts = cell(fitParams.concatInfo.n,1);
modelResponse = [];residual = [];

% The HRF is identical for every run in a concat. For all models except the
% HDR-fitting model it is also identical for every objective evaluation.
if isfield(fitParams,'fixedHRF') && ~isempty(fitParams.fixedHRF)
  hrf = fitParams.fixedHRF;
else
  hrf = getCanonicalHRF(p.canonical,fitParams.stimDeltaT);
end

% create the model for each concat
for i = 1:fitParams.concatInfo.n
  %%%%%%%%%%%%%%%%%%%
  % Get model response, which involves convolving model with stimulus, and with HRF and dropping junk frames.
  thisModelResponse = prfModel('getModelResponse', fitParams, rfModel, hrf, p, i);
  %%%%%%%%%%%%%%%%%%%

  % check if we need to resample, by seeing if the model response has the
  % same number of samples as the concatInfo expects for the scan (which should
  % occur if the stim images are sampled at the volume rate
  nVols = (fitParams.concatInfo.runTransition(i,2)-fitParams.concatInfo.runTransition(i,1)+1);
  if length(thisModelResponse) ~= nVols
    % resample 
    resampledTimeSeries = resample(timeseries(thisModelResponse,fitParams.stimT(1:length(thisModelResponse))),(0:(nVols-1))*fitParams.framePeriod);
    % and pull out of timeseries structure
    thisModelResponse = squeeze(resampledTimeSeries.Data);
  end
  
  % apply concat filtering
  if isfield(fitParams,'applyFiltering') && fitParams.applyFiltering
    thisModelResponse = applyConcatFiltering(thisModelResponse,fitParams.concatInfo,i);
  else
    % with no filtering, just remove mean
    thisModelResponse = thisModelResponse - mean(thisModelResponse);
  end
  
  if ~justGetModel
    % compute correlation of this portion of the model response with time series
    thisTSeries = tSeries(fitParams.concatInfo.runTransition(i,1):fitParams.concatInfo.runTransition(i,2));
    thisTSeries = thisTSeries - mean(thisTSeries);
    
    % check here for length
    if length(thisTSeries) ~= length(thisModelResponse)
      error('pRFFit:ModelDataLengthMismatch', ...
        ['Voxel tSeries length of %i does not match model length of %i. ' ...
         'This can happen if the tSense factor or junk-frame setting is incorrect.'], ...
        length(thisTSeries),length(thisModelResponse));
    end
  
    % Per-run correlations are only part of the four-output reporting call;
    % neither optimizer uses them in its objective.
    if nargout >= 4
      r(i) = corr(thisTSeries(:),thisModelResponse(:));
    end

    if fitParams.betaEachScan
      % scale and offset the model to best match the tSeries
      [thisModelResponse thisResidual] = scaleAndOffset(thisModelResponse',thisTSeries(:));
    else
      thisResidual = [];
    end
  else
    thisResidual = [];
  end
  
  % save each run as a column, preserving the original concatenation order
  modelResponseParts{i} = thisModelResponse(:);
  residualParts{i} = thisResidual(:);
end

modelResponse = vertcat(modelResponseParts{:});
residual = vertcat(residualParts{:});

% return model only
if justGetModel,return,end

% scale the whole time series
if ~fitParams.betaEachScan
  [modelResponse residual] = scaleAndOffset(modelResponse,tSeries(:));
end


% display the fit
if fitParams.dispFit
  dispModelFit(params,fitParams,modelResponse,tSeries,rfModel);
end

% for nelder-mead just compute correlation and return 1-4
if strcmp(lower(fitParams.algorithm),'nelder-mead') || strcmp(lower(fitParams.algorithm), 'nelder-mead-bnd') || strcmp(lower(fitParams.algorithm), 'fmincon')
  residual = -corr(modelResponse,tSeries);
%  disp(sprintf('(pRFFit:getModelResidual) r: %f',residual));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     objectiveMemo        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = objectiveMemo(action,varargin)

% Each fit has an isolated, bounded cache. Parameter keys use MATLAB's
% hexadecimal floating-point representation, so only bit-identical vectors
% are reused. The time series and fit settings are fixed by the session.
persistent sessions nextSessionId
if isempty(sessions)
  sessions = containers.Map('KeyType','char','ValueType','any');
  nextSessionId = uint64(0);
end

switch action
  case 'create'
    nextSessionId = nextSessionId+1;
    sessionId = sprintf('fit-%016x',nextSessionId);
    session.values = containers.Map('KeyType','char','ValueType','any');
    session.order = {};
    session.hits = 0;
    session.misses = 0;
    session.maxEntries = 256;
    sessions(sessionId) = session;
    varargout{1} = sessionId;

  case 'evaluate'
    sessionId = varargin{1};
    params = varargin{2};
    tSeries = varargin{3};
    fitParams = varargin{4};
    finitePenalty = varargin{5};
    session = sessions(sessionId);
    paramsKey = exactParameterKey(params);
    if isKey(session.values,paramsKey)
      residual = session.values(paramsKey);
      session.hits = session.hits+1;
    else
      residual = getModelResidual(params,tSeries,fitParams);
      if finitePenalty && (~isscalar(residual) || ~isfinite(residual))
        residual = 1;
      end
      session.values(paramsKey) = residual;
      session.order{end+1} = paramsKey;
      session.misses = session.misses+1;
      if numel(session.order) > session.maxEntries
        oldestKey = session.order{1};
        session.order(1) = [];
        if isKey(session.values,oldestKey)
          remove(session.values,oldestKey);
        end
      end
    end
    sessions(sessionId) = session;
    varargout{1} = residual;

  case 'stats'
    session = sessions(varargin{1});
    stats.hits = session.hits;
    stats.misses = session.misses;
    stats.entries = session.values.Count;
    varargout{1} = stats;

  case 'destroy'
    sessionId = varargin{1};
    if isKey(sessions,sessionId)
      remove(sessions,sessionId);
    end

  otherwise
    error('pRFFit:InvalidObjectiveCacheAction', ...
      'Unknown objective-cache action: %s',action);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   exactParameterKey        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = exactParameterKey(params)

paramHex = num2hex(params(:));
key = sprintf('%s:%s:%s',class(params),mat2str(size(params)),reshape(paramHex.',1,[]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  makeFiniteFminconInitial    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [finiteParams finiteObjective] = makeFiniteFminconInitial(params,tSeries,fitParams,minParams,maxParams)

finiteParams = params;
finiteObjective = getModelResidual(params,tSeries,fitParams);
if isscalar(finiteObjective) && isfinite(finiteObjective)
  return
end

finiteObjective = inf;
for iParam = 1:numel(params)
  candidate = params;
  if isfinite(minParams(iParam)) && (params(iParam) > minParams(iParam))
    candidate(iParam) = minParams(iParam) + 0.1*(params(iParam)-minParams(iParam));
  elseif isfinite(maxParams(iParam)) && (params(iParam) < maxParams(iParam))
    candidate(iParam) = maxParams(iParam) - 0.1*(maxParams(iParam)-params(iParam));
  else
    continue
  end

  candidateObjective = getModelResidual(candidate,tSeries,fitParams);
  if isscalar(candidateObjective) && isfinite(candidateObjective) && ...
      (candidateObjective < finiteObjective)
    finiteParams = candidate;
    finiteObjective = candidateObjective;
  end
end

if ~isfinite(finiteObjective)
  error('pRFFit:FminconNoFiniteStart', ...
    'The optional fmincon solver could not find a finite starting model.');
end


%%%%%%%%%%%%%%%%%%%%%%
%    dispModelFit    %
%%%%%%%%%%%%%%%%%%%%%%
function dispModelFit(params,fitParams,modelResponse,tSeries,rfModel)

mlrSmartfig('pRFFit_getModelResidual','reuse');
clf
subplot(4,4,[1:3 5:7 9:11 13:15]);
%plot(fitParams.stimT(fitParams.junkFrames+1:end),tSeries,'k-');
plot(tSeries,'k-');
hold on
%plot(fitParams.stimT(fitParams.junkFrames+1:end),modelResponse,'r-');
plot(modelResponse,'r-');
xlabel('Time (sec)');
ylabel('BOLD (% sig change)');
p = getFitParams(params,fitParams);
titleStr = sprintf('x: %s y: %s rfHalfWidth: %s',mlrnum2str(p.x),mlrnum2str(p.y),mlrnum2str(p.std));
titleStr = sprintf('%s\n(timelag: %s tau: %s exponent: %s)',titleStr,mlrnum2str(p.canonical.timelag),mlrnum2str(p.canonical.tau),mlrnum2str(p.canonical.exponent));
if p.canonical.diffOfGamma
  titleStr = sprintf('%s - %s x (timelag2: %s tau2: %s exponent2: %s)',titleStr,mlrnum2str(p.canonical.amplitudeRatio),mlrnum2str(p.canonical.timelag2),mlrnum2str(p.canonical.tau2),mlrnum2str(p.canonical.exponent2));
end
title(titleStr);
axis tight

subplot(4,4,[8 12 16]);
imagesc(fitParams.stimX(:,1),fitParams.stimY(1,:),flipud(rfModel'));
axis image;
hold on
hline(0);vline(0);

subplot(4,4,4);cla
p = getFitParams(params,fitParams);
canonical = getCanonicalHRF(p.canonical,fitParams.framePeriod);
plot(canonical.time,canonical.hrf,'k-')
if exist('myaxis') == 2,myaxis;end

%%%%%%%%%%%%%%%%%%%%%%%%
%    scaleAndOffset    %
%%%%%%%%%%%%%%%%%%%%%%%%
function [modelResponse residual] = scaleAndOffset(modelResponse,tSeries)

designMatrix = modelResponse;
designMatrix(:,2) = 1;

% get beta weight for the modelResponse
if ~any(isnan(modelResponse)) && ~any(isinf(modelResponse))
  beta = pinv(designMatrix)*tSeries;
  beta(1) = max(beta(1),0);
  modelResponse = designMatrix*beta;
  residual = tSeries-modelResponse;
else
  residual = tSeries;
end

%%%%%%%%%%%%%%%%%%%%%%
%%   getFitParams   %%
%%%%%%%%%%%%%%%%%%%%%%
function p = getFitParams(params,fitParams)

p.rfType = fitParams.rfType;
p = prfModel('getFitParams', fitParams, params);

%%%%%%%%%%%%%%%%%%%%%
%%   getGammaHRF   %%
%%%%%%%%%%%%%%%%%%%%%
function fun = getGammaHRF(time,p)

fun = thisGamma(time,1,p.timelag,p.offset,p.tau,p.exponent)/100;
% add second gamma if this is a difference of gammas fit
if p.diffOfGamma
  fun = fun - thisGamma(time,p.amplitudeRatio,p.timelag2,p.offset2,p.tau2,p.exponent2)/100;
end

%%%%%%%%%%%%%%%%%%%
%%   thisGamma   %%
%%%%%%%%%%%%%%%%%%%
function gammafun = thisGamma(time,amplitude,timelag,offset,tau,exponent)

exponent = round(exponent);
% gamma function
gammafun = (((time-timelag)/tau).^(exponent-1).*exp(-(time-timelag)/tau))./(tau*factorial(exponent-1));

% negative values of time are set to zero,
% so that the function always starts at zero
gammafun(find((time-timelag) < 0)) = 0;

% normalize the amplitude
if (max(gammafun)-min(gammafun))~=0
  gammafun = (gammafun-min(gammafun)) ./ (max(gammafun)-min(gammafun));
end
gammafun = (amplitude*gammafun+offset);


%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getCanonicalHRF   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function hrf = getCanonicalHRF(params,sampleRate)

hrf.time = 0:sampleRate:params.lengthInSeconds;

% default to using getGammaHRF
if ~isfield(params,'function')
  params.function = 'getGammaHRF';
end

switch params.function
    case {'rmHrfTwoGammas'}
    % use vistasoft parameterization
    hrf.hrf = rmHrfTwogammas(hrf.time,params.params);
  case {'getGammaHRF'}
    hrf.hrf = getGammaHRF(hrf.time,params);
    % normalize to amplitude of 1
    hrf.hrf = hrf.hrf / max(hrf.hrf);
  otherwise
    error('pRFFit:UnknownCanonicalFunction', ...
      'Unknown canonical function: %s',params.function);
end

%%%%%%%%%%%%%%%%%%%%%%%%
%    rmHrfTwoGammas    %
%%%%%%%%%%%%%%%%%%%%%%%%
function [h]=rmHrfTwogammas(t,params)
% Copied from Vistasoft
%
% Create an HRF based on the SPM two gamma model. 
% 
% [h]=rmHrfTwogammas(t,params)
% 
% t: a range of latencies
% params(1): peak gamma 1
% params(2): fwhm gamma 1
% params(3): peak gamma 2
% params(4): fwhm gamma 2
% params(5): dip
% Final hrf is:   gamma1/max(gamma1)-dip*gamma2/max(gamma2)
% from Glover, NeuroImage, 9:416-429
%
%
% Example 1 Use defaults.
%   tr = 1.5;
%   t  = tr * (0:50);
%   h = rmHrfTwogammas(t);
%   figure; plot(t, h); 
%   xlabel('time (seconds)'); 
%   ylabel('response');
%
% Example 2 Put in your own parameters
%   tr = 1.5;
%   t  = tr * (0:50);
%   params = [3 5 15 20 0.2];
%   h1 = rmHrfTwogammas(t);
%   h2 = rmHrfTwogammas(t, params);
%   figure; plot(t, h1, 'r', t, h2, 'g'); 
%   xlabel('time (seconds)'); 
%   ylabel('response');
%   legend('Default', 'Subject specfic');

% 2.11.2011 JW: Cropped out of rfConvolveTC and made into a separate function 

% If no HRF parameters input, use defauls
if nargin < 2 || isempty(params), 
    params = [5.4 5.2 10.8 7.35 0.35];
end;

% params
peak1 = params(1);
fwhm1 = params(2);
peak2 = params(3);
fwhm2 = params(4);
dip   = params(5);

% sanity check
if(peak1 == 0 || fwhm1 ==0),
    fprintf('[%s]: zero params',mfilename);
    params, %#ok<NOPRT>
    return;
end;

% Taylor:
alpha1=peak1^2/fwhm1^2*8*log(2);
beta1=fwhm1^2/peak1/8/log(2);
gamma1=(t/peak1).^alpha1.*exp(-(t-peak1)./beta1);

if peak2>0 && fwhm2>0
    alpha2=peak2^2/fwhm2^2*8*log(2);
    beta2=fwhm2^2/peak2/8/log(2);
    gamma2=(t/peak2).^alpha2.*exp(-(t-peak2)./beta2);
else
    gamma2=min(abs(t-peak2))==abs(t-peak2);
end
h = gamma1-dip*gamma2;
%h = h./sum(h);

return;

%%%%%%%%%%%%%%%%%%%%
%%   getRFModel   %%
%%%%%%%%%%%%%%%%%%%%
function rfModel = getRFModel(params,fitParams)

% now generate the rfModel
if any(strcmp(fitParams.rfType,{'gaussian','gaussian-hdr', 'gaussian-css', 'gaussian-diffs', 'gaussian-divs', 'gaussian-DoG-CSS', 'divNorm'}))
  rfModel = makeRFGaussian(params,fitParams);
else
  error('pRFFit:UnknownRFType','Unknown rfType: %s',fitParams.rfType);
end

%%%%%%%%%%%%%%%%%%%
%%  prfModel     %%
%%%%%%%%%%%%%%%%%%%
function output = prfModel(varargin)

fitParams = varargin{2};
switch (fitParams.rfType)
  case 'gaussian'
    output = pRF_gaussian(varargin{:});
  case 'gaussian-hdr'
    output = pRF_gaussianhdr(varargin{:});
  case 'gaussian-diffs'
    output = pRF_diffGaussian(varargin{:});
  case 'gaussian-divs'
    output = pRF_divGaussian(varargin{:});
  case 'gaussian-css'
    output = pRF_css(varargin{:});
  case 'gaussian-DoG-CSS'
    output = pRF_DoG_CSS(varargin{:});
  case 'divNorm'
    output = pRF_divNorm(varargin{:});
  otherwise
    error('pRFFit:UnknownRFType','Unknown rfType: %s',fitParams.rfType);
end


%%%%%%%%%%%%%%%%%%%%%%%%
%%   makeRFGaussian   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function rfModel = makeRFGaussian(params,fitParams)

% compute rf
rfModel = pRFMakeGaussian(fitParams,params.x,params.y,params.std);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  prepareGaussianGrid    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grid = prepareGaussianGrid(stimX,stimY)

grid = struct('layout','','firstAxis',[],'secondAxis',[]);
isSupported = isa(stimX,'double') && isa(stimY,'double') && ...
  isreal(stimX) && isreal(stimY) && ~issparse(stimX) && ...
  ~issparse(stimY) && ismatrix(stimX) && ismatrix(stimY) && ...
  isequal(size(stimX),size(stimY)) && ~isempty(stimX);
if ~isSupported
  return
end

% pRFGetStimImageFromStimfile creates an NDGRID layout, but accept MATLAB's
% other common rectilinear convention as well. Comparisons are exact so an
% arbitrary or warped grid always falls back to the original expression.
isNdgrid = all(all(diff(stimX,1,2) == 0)) && ...
  all(all(diff(stimY,1,1) == 0));
if isNdgrid
  grid.layout = 'ndgrid';
  grid.firstAxis = stimX(:,1);
  grid.secondAxis = stimY(1,:);
  return
end

isMeshgrid = all(all(diff(stimX,1,1) == 0)) && ...
  all(all(diff(stimY,1,2) == 0));
if isMeshgrid
  grid.layout = 'meshgrid';
  grid.firstAxis = stimY(:,1);
  grid.secondAxis = stimX(1,:);
end

%%%%%%%%%%%%%%%%%%%
%    parseArgs    %
%%%%%%%%%%%%%%%%%%%
function [v scanNum x y s fitParams tSeries] = parseArgs(args);

v = [];scanNum=[];x=[];y=[];s=[];fitParams=[];tSeries = [];

% check for calling convention from interrogator
if (length(args) >= 7) && isnumeric(args{6})
  v = args{1};
  %overlayNum = args{2};
  scanNum = args{3};
  x = args{4};
  y = args{5};
  s = args{6};
  %roi = args{7};
  fitParams.dispFit = true;
  fitParams.optimDisplay = 'final';
  fitParams.algorithm = 'nelder-mead';
  fitParams.getModelResponse = false;
  fitParams.prefit = [];
  fitParams.xFlipStimulus = 0;
  fitParams.yFlipStimulus = 0;
  fitParams.timeShiftStimulus = 0;
  fitParams.betaEachScan = false;
  fitParams.justGetStimImage = false;
  fitParams.returnPrefit = false;
  fitParams.verbose = 1;
  fitParams.timelag = 1;
  fitParams.tau = 0.6;
  fitParams.exponent = 6;
  fitTypeParams = [];
  getArgs({args{8:end}},{'fitTypeParams=[]'});
  if isempty(fitTypeParams)
    % no fit type params, check if we have them set in
    % the global (this is useful so that when called as an 
    % interrogator we don't have to keep setting them
    global gpRFFitTypeParams
    % if user is holding shift, then reget parameters
    if ~isempty(gcf) && any(strcmp(get(gcf,'CurrentModifier'),'shift'))
      gpRFFitTypeParams = [];
    end
    % get the parameters from the user interface if not already set
    if isempty(gpRFFitTypeParams)
      fitTypeParams = pRFGUI('pRFFitParamsOnly=1','v',v);
      if isempty(fitTypeParams) 
	v = [];
	return
      end
      gpRFFitTypeParams = fitTypeParams;
    else
      % otherwise grab old ones
      disp(sprintf('(pRFFit) Using already set parameters to compute pRFFit. If you want to use different parameters, hold shift down as you click the next voxel'));
      fitTypeParams = gpRFFitTypeParams;
    end
  end
  if ~isempty(fitTypeParams)
    % if fitTypeParams is passed in (usually from pRF / pRFGUI) then
    % grab parameters off that structure
    fitTypeParamsFields = fieldnames(fitTypeParams);
    for i = 1:length(fitTypeParamsFields)
      fitParams.(fitTypeParamsFields{i}) = fitTypeParams.(fitTypeParamsFields{i});
    end
  end
  
% normal calling convention
elseif length(args) >= 5
  v = args{1};
  scanNum = args{2};
  x = args{3};
  y = args{4};
  s = args{5};
  % parse any more arguments
  dispFit=[];stim = [];getModelResponse = [];params = [];concatInfo = [];prefit = [];
  xFlipStimulus=[];yFlipStimulus=[];timeShiftStimulus=[];rfType=[];betaEachScan=[];fitTypeParams = [];
  dispIndex = [];dispN = [];returnPrefit = [];returnPreparedContext=[];tSeries=[];quickPrefit=[];junkFrames=[];
  verbose = [];justGetStimImage = [];framePeriod = [];
  getArgs({args{6:end}},{'dispFit=0','stim=[]','getModelResponse=0','params=[]','concatInfo=[]','prefit=[]','xFlipStimulus=0','yFlipStimulus=0','timeShiftStimulus=0','rfType=gaussian','betaEachScan=0','fitTypeParams=[]','justGetStimImage=[]','verbose=1','dispIndex=[]','dispN=[]','returnPrefit=0','returnPreparedContext=0','quickPrefit=0','tSeries=[]','junkFrames=[]','framePeriod=[]','paramsInfo=[]'});
  % default to display fit
  fitParams.dispFit = dispFit;
  fitParams.stim = stim;
  fitParams.optimDisplay = 'off';
  fitParams.getModelResponse = getModelResponse;
  fitParams.params = params;
  fitParams.concatInfo = concatInfo;
  fitParams.prefit = prefit;
  fitParams.xFlipStimulus = xFlipStimulus;
  fitParams.yFlipStimulus = yFlipStimulus;
  fitParams.timeShiftStimulus = timeShiftStimulus;
  fitParams.rfType = rfType;
  fitParams.betaEachScan = betaEachScan;
  fitParams.justGetStimImage = justGetStimImage;
  fitParams.verbose = verbose;
  fitParams.returnPrefit = returnPrefit;
  fitParams.returnPreparedContext = returnPreparedContext;
  fitParams.junkFrames = junkFrames;
  fitParams.framePeriod = framePeriod;
  % now read in all the fields in the paramsInfo
  if ~isempty(paramsInfo)
    paramsInfoFields = fieldnames(paramsInfo);
    for iField = 1:length(paramsInfoFields)
      fitParams.(paramsInfoFields{iField}) = paramsInfo.(paramsInfoFields{iField});
    end
  end
  if ~isempty(fitTypeParams)
    % if fitTypeParams is passed in (usually from pRF / pRFGUI) then
    % grab parameters off that structure
    fitTypeParamsFields = fieldnames(fitTypeParams);
    for i = 1:length(fitTypeParamsFields)
      fitParams.(fitTypeParamsFields{i}) = fitTypeParams.(fitTypeParamsFields{i});
    end
  end
  if ~isempty(dispIndex) && ~isempty(dispN)
    % create a display string. Note that we use sprintf twice here so that
    % we can create a string with the proper amount of space padding the index
    % so that each row always displays as the same length string
    prefitOnlyStr = '';
    if isfield(fitParams,'prefitOnly') && fitParams.prefitOnly
      prefitOnlyStr = ' (prefit only)';
    end
    fitParams.dispstr = sprintf(sprintf('Voxel %%%i.f/%%i%%s: ',length(sprintf('%i',dispN))),dispIndex,dispN,prefitOnlyStr);
    end
  if getModelResponse && isempty(params)
    disp(sprintf('(pRFFit) Must pass in params when using getModelResponse'));
    fitParams.getModelResponse = false;
  end
else
  help pRFFit;
end

% some default parameters
if ~isfield(fitParams,'prefitOnly') || isempty(fitParams.prefitOnly)
  fitParams.prefitOnly = false;
end
if ~isfield(fitParams,'dispstr')
  fitParams.dispstr = '';
end
if ~isfield(fitParams,'quickPrefit') || isempty(fitParams.quickPrefit)
  fitParams.quickPrefit = false;
end
if ~isfield(fitParams,'verbose') || isempty(fitParams.verbose)
  fitParams.verbose = true;
end
if ~isfield(fitParams,'numWorkers') || isempty(fitParams.numWorkers)
  % Direct/interrogator calls and older saved fit structures use the same
  % all-local-workers default as the main pRF GUI.
  fitParams.numWorkers = pRFResolveNumWorkers;
end
if ~(isnumeric(fitParams.numWorkers) && isreal(fitParams.numWorkers) && ...
    isscalar(fitParams.numWorkers) && isfinite(fitParams.numWorkers) && ...
    fitParams.numWorkers >= 1 && fitParams.numWorkers == fix(fitParams.numWorkers))
  error('pRFFit:InvalidNumWorkers','numWorkers must be a positive integer.');
end
if ~isfield(fitParams,'returnPreparedContext') || isempty(fitParams.returnPreparedContext)
  fitParams.returnPreparedContext = false;
end

% get some info about the scanNum
if ~isfield(fitParams,'framePeriod') || isempty(fitParams.framePeriod)
  fitParams.framePeriod = viewGet(v,'framePeriod');
end
if ~isfield(fitParams,'junkFrames') || isempty(fitParams.junkFrames)
  fitParams.junkFrames = viewGet(v,'junkFrames',scanNum);
end


if isempty(fitParams.prefit) || (fitParams.prefit.quickPrefit ~= fitParams.quickPrefit)
  % set the values over which to first prefit
  % the best of these parameters will then be used 
  % to init the non-linear optimization. Note that the
  % values here are expressed as a factor of the screen
  % dimensions (1 being the width/height of the screen)
  % Later when the prefit is calculated, they will be multiplied
  % by the screenWidth and screenHeight
  if fitParams.quickPrefit
    if fitParams.verbose,disp(sprintf('(pRFFit) Doing quick prefit'));end
    % make sure here that x and y points go through 0 symmetrically
    [prefitx prefity prefitrfHalfWidth] = ndgrid(-0.375:0.125:0.375,-0.375:0.125:0.375,[0.025 0.05 0.15 0.4]);
  else
    [prefitx prefity prefitrfHalfWidth] = ndgrid(-0.4:0.025:0.4,-0.4:0.025:0.4,[0.0125 0.025 0.05 0.1 0.25 0.5 0.75]);
  end
  fitParams.prefit.quickPrefit = fitParams.quickPrefit;
  fitParams.prefit.n = length(prefitx(:));
  fitParams.prefit.x = prefitx(:);
  fitParams.prefit.y = prefity(:);
  fitParams.prefit.rfHalfWidth = prefitrfHalfWidth(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    makePreparedContext       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function context = makePreparedContext(v,scanNum,fitParams)

context.kind = 'pRFPreparedContext';
context.version = 1;
context.scanNum = scanNum;
context.groupName = '';
try
  context.groupName = viewGet(v,'groupName');
catch
  % groupName is used only in a verbose NaN diagnostic.
end
fitParams.returnPrefit = false;
fitParams.returnPreparedContext = false;
fitParams.getModelResponse = false;
context.fitParams = fitParams;

% Do not mirror stimulus, prefit, or cache metadata at the top level:
% scan-invariant state remains available through context.fitParams.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    validatePreparedContext   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function validatePreparedContext(context)

isValid = isstruct(context) && isfield(context,'kind') && ...
  strcmp(context.kind,'pRFPreparedContext') && isfield(context,'version') && ...
  isequal(context.version,1) && isfield(context,'fitParams') && ...
  isfield(context,'scanNum');
if ~isValid
  error('pRFFit:InvalidPreparedContext', ...
    'The supplied value is not a compatible pRFFit prepared context.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    makeVoxelDisplayString    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayString = makeVoxelDisplayString(fitParams,dispIndex,dispN)

displayString = '';
if ~isempty(dispIndex) && ~isempty(dispN)
  prefitOnlyStr = '';
  if isfield(fitParams,'prefitOnly') && fitParams.prefitOnly
    prefitOnlyStr = ' (prefit only)';
  end
  displayString = sprintf(sprintf('Voxel %%%i.f/%%i%%s: ', ...
    length(sprintf('%i',dispN))),dispIndex,dispN,prefitOnlyStr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    makePrefitCacheKey      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = makePrefitCacheKey(fitParams,trailingParams)

% Include the original stimulus values but intentionally ignore derived
% acceleration metadata (unique-frame maps, prepared projections, etc.).
% Those fields change how a response is calculated, not what it is.
stimIdentity = cell(size(fitParams.stim));
stimFields = {'im','x','y','t'};
for iStim = 1:numel(fitParams.stim)
  stimIdentity{iStim} = struct;
  for iField = 1:numel(stimFields)
    fieldName = stimFields{iField};
    if isfield(fitParams.stim{iStim},fieldName)
      stimIdentity{iStim}.(fieldName) = fitParams.stim{iStim}.(fieldName);
    end
  end
end

identity.cacheFormatVersion = 2;
identity.stim = stimIdentity;
identity.stimX = fitParams.stimX;
identity.stimY = fitParams.stimY;
identity.stimT = fitParams.stimT;
identity.concatInfo = fitParams.concatInfo;
identity.prefit.n = fitParams.prefit.n;
identity.prefit.x = fitParams.prefit.x;
identity.prefit.y = fitParams.prefit.y;
identity.prefit.rfHalfWidth = fitParams.prefit.rfHalfWidth;
identity.trailingParams = trailingParams;

% Fields read by the current model-response implementations and temporal
% processing. Presence is encoded as well as value because several paths
% distinguish a missing field from an empty one.
modelFields = {'rfType','framePeriod','stimDeltaT','fixedHRF','applyFiltering', ...
  'timelag','tau','exponent','diffOfGamma','amplitudeRatio','timelag2', ...
  'tau2','exponent2','canonicalFunction','canonicalParams'};
for iField = 1:numel(modelFields)
  fieldName = modelFields{iField};
  identity.modelSettings.(fieldName).present = isfield(fitParams,fieldName);
  if isfield(fitParams,fieldName)
    identity.modelSettings.(fieldName).value = fitParams.(fieldName);
  else
    identity.modelSettings.(fieldName).value = [];
  end
end
key = exactValueHash(identity);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    exactValueHash          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = exactValueHash(value)

bytes = getByteStreamFromArray(value);
try
  digestEngine = java.security.MessageDigest.getInstance('SHA-256');
  digestEngine.update(bytes);
  digest = typecast(digestEngine.digest(),'uint8');
  key = lower(reshape(dec2hex(digest,2).',1,[]));
catch
  % This exact fallback is large, but it is used only in no-JVM MATLAB.
  key = ['raw-' lower(reshape(dec2hex(bytes,2).',1,[]))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    prefitBankCache         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = prefitBankCache(action,key,value)

persistent cacheKeys cacheValues cacheBytes
maxEntries = 4;
maxBytes = 512*2^20;
if isempty(cacheKeys)
  cacheKeys = {};
  cacheValues = {};
  cacheBytes = [];
end

switch action
  case 'get'
    cacheIndex = find(strcmp(cacheKeys,key),1,'last');
    cacheHit = ~isempty(cacheIndex);
    if cacheHit
      cachedValue = cacheValues{cacheIndex};
      % Move the entry to the back to maintain a bounded LRU cache.
      cacheKeys = [cacheKeys(1:cacheIndex-1) cacheKeys(cacheIndex+1:end) cacheKeys(cacheIndex)];
      cacheValues = [cacheValues(1:cacheIndex-1) cacheValues(cacheIndex+1:end) cacheValues(cacheIndex)];
      cacheBytes = [cacheBytes(1:cacheIndex-1) cacheBytes(cacheIndex+1:end) cacheBytes(cacheIndex)];
    else
      cachedValue = [];
    end
    varargout{1} = cachedValue;
    varargout{2} = cacheHit;

  case 'put'
    existingIndex = find(strcmp(cacheKeys,key),1,'last');
    if ~isempty(existingIndex)
      cacheKeys(existingIndex) = [];
      cacheValues(existingIndex) = [];
      cacheBytes(existingIndex) = [];
    end
    valueInfo = whos('value');
    valueBytes = double(valueInfo.bytes);
    if valueBytes > maxBytes
      return
    end
    cacheKeys{end+1} = key;
    cacheValues{end+1} = value;
    cacheBytes(end+1) = valueBytes;
    while numel(cacheKeys) > maxEntries || sum(cacheBytes) > maxBytes
      cacheKeys(1) = [];
      cacheValues(1) = [];
      cacheBytes(1) = [];
    end

  otherwise
    error('pRFFit:InvalidPrefitCacheAction','Unknown prefit cache action: %s',action);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    checkStimForAverages    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stim ignoreMismatchStimfiles] = checkStimForAverages(v,scanNum,groupNum,stim,concatInfo,stimImageDiffTolerance)

ignoreMismatchStimfiles = false;  

% this function will check for some bad casses (like concat of concats etc)
% it will also check that all the component scans of an average have the
% same stim image and warn if they do not. It will then replace the stim cell
% array for the average with a single stim file, so that processing
% can continue as normal for pRFFit

% if not a cell, then ok, return
if ~iscell(stim),return,end

% first check for bad shiftList or refverseLIst
p = viewGet(v,'params',scanNum,groupNum);
if isfield(p,'params') && isfield(p.params,'shiftList') && any(p.params.shiftList~=0)
  error('pRFFit:UnsupportedAverageShift', ...
    'Component scan %s:%i has a nonzero shiftList (%s); shifted averages are not supported.', ...
    viewGet(v,'groupName',groupNum),scanNum,mlrnum2str(p.params.shiftList));
end
if isfield(p,'params') && isfield(p.params,'reverseList') && any(p.params.reverseList~=0)
  error('pRFFit:UnsupportedAverageReversal', ...
    'Component scan %s:%i has a nonzero reverseList (%s); time-reversed averages are not supported.', ...
    viewGet(v,'groupName',groupNum),scanNum,mlrnum2str(p.params.reverseList));
end

% if is a cell, check to see if this is a concat or not
if ~isempty(concatInfo) && (concatInfo.isConcat)
  % this is a concat, so check each one of the elements
  [originalScanNum originalGroupNum] = viewGet(v,'originalScanNum',scanNum,groupNum);
  for i = 1:length(stim)
    % get concatInfo for original scan
    concatInfo = viewGet(v,'concatInfo',originalScanNum(i),originalGroupNum(i));
    if ~isempty(concatInfo)
      error('pRFFit:UnsupportedNestedConcatenation', ...
        'Concatenations containing concatenated component scans are not supported.');
    end
    % check this next scan
    [stim{i} ignoreMismatchStimfiles] = checkStimForAverages(v,originalScanNum(i),originalGroupNum(i),stim{i},concatInfo,stimImageDiffTolerance);
    % if user has accepted all then set stimImageDiffTOlerance to infinity
    if isinf(ignoreMismatchStimfiles),stimImageDiffTolerance = inf;end
    if isempty(stim{i}),stim = [];return,end
  end
else
  % this for orignals
  [originalScanNum originalGroupNum] = viewGet(v,'originalScanNum',scanNum,groupNum);
  % if it is an original than check each element
  if ~isempty(originalScanNum)
    % check that this is not an average of a concat
    for i = 1:length(stim)
      % get concatInfo for original scan
      concatInfo = viewGet(v,'concatInfo',originalScanNum(i),originalGroupNum(i));
      if ~isempty(concatInfo)
	error('pRFFit:UnsupportedAverageOfConcatenation', ...
	  'Averages containing concatenated component scans are not supported.');
      end
      % see if it is an average of an average
      originalOfOriginalScanNum = viewGet(v,'originalScanNum',originalScanNum(i),originalGroupNum(i));
      if length(originalOfOriginalScanNum) > 1
	error('pRFFit:UnsupportedNestedAverage', ...
	  'Averages containing another averaged scan are not supported.');
      end
    end
    % ok, not an average of a concatenation/average so check all the stim files 
    % and warn if there are any inconsistencies
    for i = 1:length(stim)
      if ~isequalwithequalnans(stim{1}.im,stim{i}.im)    
	dispHeader
	disp(sprintf('(pRFFit:checkStimForAverages) !!! Average for %s:%i component scan %i does not match stimulus for other scans. If you wish to continue then this will use the stimfile associated with the first scan in the average !!!',viewGet(v,'groupName',groupNum),scanNum,originalScanNum(i)));
	% display which volumes are different
	diffVols = [];
	for iVol = 1:size(stim{1}.im,3)
	  if ~isequalwithequalnans(stim{1}.im(:,:,iVol),stim{i}.im(:,:,iVol))
	    diffVols(end+1) = iVol;
	  end
	end
	disp(sprintf('(pRFFit) Stimulus files are different at %i of %i vols (%0.1f%%): %s',length(diffVols),size(stim{1}.im,3),100*length(diffVols)/size(stim{1}.im,3),num2str(diffVols)));
	if 100*(length(diffVols)/size(stim{1}.im,3)) < stimImageDiffTolerance
	  disp(sprintf('(pRFFit) This could be from minor timing inconsistencies, so ignoring it. Set stimImageDiffTolerance lower if you want to stop the code when this happens'));
	else
	  % ask user if they want to continue (only if there is a difference of more than 10 vols	  
	  ignoreMismatchStimfiles = askuser('Do you wish to continue',1);
	  if ~ignoreMismatchStimfiles
	    stim = [];
	    return;
	  end
	end
	dispHeader
      end
    end
    % if we passed the above, this is an average of identical
    % scans, so just keep the first stim image since they are all the same
    stim = stim{1};
  end
end

%%%%%%%%%%%%%%%%%
%    getStim    %
%%%%%%%%%%%%%%%%%
function stim = getStim(v,scanNum,fitParams)

% get stimfile
stimfile = viewGet(v,'stimfile',scanNum);
% get volume to trigger ratio
volTrigRatio = viewGet(v,'auxParam','volTrigRatio',scanNum);
% Identify the linked stimulus source as well as the scan number. A view can
% be relinked without changing its group/scan indices.
[stimfileIdentity hasReliableStimfileIdentity] = getStimfileCacheIdentity(stimfile);
% check if global matches
groupNum = viewGet(v,'curGroup');
global gpRFFitStimImage
stimulusCacheVersion = 3;
requiredCacheFields = {'version','scanNum','groupNum','stimfileIdentity', ...
  'volTrigRatio','xFlip','yFlip','timeShift','stim'};
cacheMatches = hasReliableStimfileIdentity && ...
  isstruct(gpRFFitStimImage) && isscalar(gpRFFitStimImage) && ...
  all(isfield(gpRFFitStimImage,requiredCacheFields)) && ...
  isequal(gpRFFitStimImage.version,stimulusCacheVersion) && ...
  isequal(gpRFFitStimImage.scanNum,scanNum) && ...
  isequal(gpRFFitStimImage.groupNum,groupNum) && ...
  isequaln(gpRFFitStimImage.stimfileIdentity,stimfileIdentity) && ...
  isequaln(gpRFFitStimImage.volTrigRatio,volTrigRatio) && ...
  isequal(gpRFFitStimImage.xFlip,fitParams.xFlipStimulus) && ...
  isequal(gpRFFitStimImage.yFlip,fitParams.yFlipStimulus) && ...
  isequal(gpRFFitStimImage.timeShift,fitParams.timeShiftStimulus);
forceStimulusHelper = ...
  (isfield(fitParams,'recomputeStimImage') && fitParams.recomputeStimImage) || ...
  (isfield(fitParams,'saveStimImage') && fitParams.saveStimImage);
if forceStimulusHelper || ~cacheMatches
  disp(sprintf('(pRFFit) Computing stim image'));

  % averageTSeries records its component scans in scanList. If the linked
  % stimfiles describe the same stimulus, reconstruct/load only the first
  % one and let the existing average checker handle the result as before.
  reuseAverageStimulus = false;
  if iscell(stimfile) && numel(stimfile) > 1 && ...
      (~isfield(fitParams.concatInfo,'isConcat') || ~fitParams.concatInfo.isConcat)
    averageInfo = viewGet(v,'params',scanNum,groupNum);
    isAverage = isstruct(averageInfo) && isfield(averageInfo,'params') && ...
      isstruct(averageInfo.params) && isfield(averageInfo.params,'scanList') && ...
      numel(averageInfo.params.scanList) == numel(stimfile);
    if isAverage
      reuseAverageStimulus = pRFGetStimImageFromStimfile(stimfile, ...
        'volTrigRatio',volTrigRatio, ...
        'recomputeStimImage',fitParams.recomputeStimImage, ...
        'compareStimfiles',1);
    end
  end

  if reuseAverageStimulus
    if fitParams.verbose
      disp(sprintf('(pRFFit) Matching stimulus parameters across %i averaged scans; constructing stimulus once',numel(stimfile)));
    end
    if iscell(volTrigRatio)
      firstVolTrigRatio = volTrigRatio{1};
    elseif numel(volTrigRatio) == numel(stimfile)
      firstVolTrigRatio = volTrigRatio(1);
    else
      firstVolTrigRatio = volTrigRatio;
    end
    % The helper returns the final stimulus orientation requested by the GUI.
    % pRFFit must not apply another spatial flip to this returned image. When
    % saving, the helper writes its canonical pre-transform image to every
    % verified matching component stimfile.
    saveStimImageTargets = [];
    if fitParams.saveStimImage
      saveStimImageTargets = stimfile;
    end
    firstStim = pRFGetStimImageFromStimfile(stimfile{1}, ...
      'volTrigRatio',firstVolTrigRatio, ...
      'xFlip',fitParams.xFlipStimulus, ...
      'yFlip',fitParams.yFlipStimulus, ...
      'timeShift',fitParams.timeShiftStimulus, ...
      'verbose',fitParams.verbose, ...
      'saveStimImage',fitParams.saveStimImage, ...
      'recomputeStimImage',fitParams.recomputeStimImage, ...
      'saveStimImageTargets',saveStimImageTargets);
    if isempty(firstStim)
      stim = [];
    else
      % Preserve the legacy cell shape through checkStimForAverages. MATLAB
      % shares the underlying arrays until modified, so this does not make N
      % physical copies of the movie.
      stim = repmat({firstStim},size(stimfile));
    end
  else
    % Exact legacy fallback: construct every supplied component.
    % As above, use the helper output directly without a post-return flip.
    stim = pRFGetStimImageFromStimfile(stimfile,'volTrigRatio',volTrigRatio,'xFlip',fitParams.xFlipStimulus,'yFlip',fitParams.yFlipStimulus,'timeShift',fitParams.timeShiftStimulus,'verbose',fitParams.verbose,'saveStimImage',fitParams.saveStimImage,'recomputeStimImage',fitParams.recomputeStimImage);
  end
  % check for averages
  stim = checkStimForAverages(v,scanNum,viewGet(v,'curGroup'),stim,fitParams.concatInfo,fitParams.stimImageDiffTolerance);
  if isempty(stim),return,end
  % make into cell array
  stim = cellArray(stim);
  % save stim image in global
  gpRFFitStimImage.version = stimulusCacheVersion;
  gpRFFitStimImage.scanNum = scanNum;
  gpRFFitStimImage.groupNum = groupNum;
  gpRFFitStimImage.stimfileIdentity = stimfileIdentity;
  gpRFFitStimImage.volTrigRatio = volTrigRatio;
  gpRFFitStimImage.xFlip = fitParams.xFlipStimulus;
  gpRFFitStimImage.yFlip = fitParams.yFlipStimulus;
  gpRFFitStimImage.timeShift = fitParams.timeShiftStimulus;
  gpRFFitStimImage.stim = stim;
else
  % otherwise load from global
  disp(sprintf('(pRFFit) Using precomputed stim image'));
  stim = gpRFFitStimImage.stim;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getStimfileCacheIdentity     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [identity isReliable] = getStimfileCacheIdentity(stimfile)

identity = [];
isReliable = false;

if iscell(stimfile)
  identity = cell(size(stimfile));
  isReliable = true;
  for iStimfile = 1:numel(stimfile)
    [identity{iStimfile} thisIsReliable] = ...
      getStimfileCacheIdentity(stimfile{iStimfile});
    isReliable = isReliable && thisIsReliable;
  end
  if ~isReliable,identity = [];end
  return
end

if isstring(stimfile) && isscalar(stimfile)
  stimfile = char(stimfile);
end
if ischar(stimfile)
  % Match getStimfile's public filename normalization before consulting file
  % metadata, so extensionless links and explicit .mat links are identical.
  stimfile = setext(stimfile,'mat');
  identity.kind = 'file';
  identity.filename = stimfile;
  fileInfo = dir(stimfile);
  if isscalar(fileInfo) && ~fileInfo.isdir
    identity.filename = fullfile(fileInfo.folder,fileInfo.name);
    identity.bytes = fileInfo.bytes;
    identity.datenum = fileInfo.datenum;
  else
    identity.bytes = [];
    identity.datenum = [];
  end
  isReliable = true;
  return
end

if isstruct(stimfile) && isscalar(stimfile)
  if isfield(stimfile,'filename')
    [identity isReliable] = getStimfileCacheIdentity(stimfile.filename);
  elseif isfield(stimfile,'myscreen') && isstruct(stimfile.myscreen) && ...
      isscalar(stimfile.myscreen) && isfield(stimfile.myscreen,'stimfile')
    [identity isReliable] = getStimfileCacheIdentity(stimfile.myscreen.stimfile);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    applyConcatFiltering    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tSeries = applyConcatFiltering(tSeries,concatInfo,runnum)

% apply the same filter as original data
% check for what filtering was done
tSeries = tSeries(:);

% apply detrending (either if concatInfo does not say what it did or if
% the filterType field has detrend in it)
if ~isfield(concatInfo,'filterType') || ~isempty(findstr('detrend',lower(concatInfo.filterType)))
  tSeries = eventRelatedDetrend(tSeries);
end

% apply hipass filter
if isfield(concatInfo,'hipassfilter') && ~isempty(concatInfo.hipassfilter{runnum})
  % check for length match
  if ~isequal(length(tSeries),length(concatInfo.hipassfilter{runnum}))
    disp(sprintf('(pRFFit:applyConcatFiltering) Mismatch dimensions of tSeries (length: %i) and concat filter (length: %i)',length(tSeries),length(concatInfo.hipassfilter{runnum})));
  else
    tSeries = real(ifft(fft(tSeries) .* repmat(concatInfo.hipassfilter{runnum}', 1, size(tSeries,2)) ));
  end
end

% project out the mean vector
if isfield(concatInfo,'projection') && ~isempty(concatInfo.projection{runnum})
  projectionWeight = concatInfo.projection{runnum}.sourceMeanVector * tSeries;
  tSeries = tSeries - concatInfo.projection{runnum}.sourceMeanVector'*projectionWeight;
end

% now remove mean
tSeries = tSeries-repmat(mean(tSeries,1),size(tSeries,1),1);

% make back into the right dimensions
tSeries = tSeries(:)';
