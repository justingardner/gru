% fitTimecourse.m
%
%        $Id: fitTimecourse.m,v 1.29 2009/02/19 21:18:50 justin Exp $ 
%      usage: d = fitTimecourse(timecourse,stimvol,framePeriod,<varargin>)
%         by: justin gardner
%       date: 07/11/08
%    purpose: fit the passed in timecourses with gamma / or difference
%             of gamma.
%
%             The timecourses argument is a 1xn matrix where n is the #
%             of volumes in the timecourse. If you are using a concatenated
%             time course, you should pass in the concatInfo as an optional
%             argument, so that the proper filtering of the models and
%             respecting of timecourse borders is applied (i.e. add after
%             the first three arguments:
%
%             'concatInfo',d.concatInfo
%  
%             stimvol is a cell array where each cell contains an array of
%             volumes where each event happened.
% 
%             framePeriod is the time between time points in the time series.
%
%             Optional variables:
%
%             'displayFit=1': Will display a figure showing the fit.
%             'checkFit=1': Will allow you to interactively step through
%                each fit to check how good they are.
%             'returnAllFields=1': returns intermediate data fields
%
%              Calculating std over events. You can also calculate the std
%              of response amplitudes if you are using fitType=glm and fitType=nonlin
%              What this does is that it uses the canonicalResponse (either what
%              is passed in or the one computed from all responses -- see fitType=glm)
%              and makes a glm with one column for every single trial and computes
%              the distribution of amplitude of response for each trial type.
%              It then computes the mean and std of this distribution and returns
%              them as d.amplitudeMEAN and d.amplitudeSTD. If you want to compute
%              the std of a group of trials (which may have different mean responses),
%              you can set the variable 'stdGroups' to be a cell array of arrays where
%              each array contains the list of stimulus type numbers that you want to
%              pool together (for instance if you want the std of all trials for stimulus
%              types 3 and 4, you would set 'stdGroups',{[3 4]}. This will then pool
%              the amplitudes of these trial types together after subtracting off the
%              mean response and then computed mean and std (mean should be near 0) and
%              return that in the field d.groupSTD
%
%             'fitType=deconv': Will do a deconvolution on the data
%                and get the amplitude dependent on what amplitudeType
%                is set to:
%                'amplitudeType=fit1': Single gamma function in which the 
%                amplitude is allowed to vary from stimulus type to stimulus
%                type, but all other parameters (tau, timelag, exponent, offset)
%                are the same for all stimulus types.
%                'amplitudeType=fit1each': Single gamma function in which all
%                parameters are allowed to vary across stimulus types.
%                'amplitudeType=fit2': Like fit1 but with a difference of gammas
%                'amplitudeType=fit2each': like fit2each w/difference of gammas
%                Note that the deconv method is the only one that fits to 
%                the data after filtering. All of the rest of the models
%                return amplitude measures that are amplitude of response before
%                filtering. Thus, when you compare deconv vs. glm or nonlin
%                amplitudes, the deconv amplitudes tend to be lower.
%                'amplitudeType=none': Just does deconv and does not fit amplitude
%
%             'fitType=glm': Uses GLM model fit to get amplitude. It
%                does this by computing the average response to all 
%                stimulus types and then allowing that to scale to
%                get the amplitudes. amplitudeType can be:
%                'amplitudeType=area': The mean response is scaled
%                so that the area under the curve is 1% signal change
%                so amplitude values can be interpreted as %signal change
%                'amplitudeType=max': Mean response is scaled so that
%                the maximum response is 1%. The amplitude values can
%                be interpreted as peak response.
%                default is area. You can also pass in your own 'canonicalResponse'
%                which will be used rather than computing the average response
%                to all stimulus types.
%
%             'fitType=nonlin': Does a nonlinear fit to the whole timecourse
%                at each stimulus occurrence a gamma or difference of gamma
%                function is placed, and the parameters are optimized to 
%                minimize the squared residual using lsqnonlin. You can
%                choose which function to fit with using "amplitudeType"
%                (see fitType=deconv, above)
% 
%              For fit1 and fit2, you can set constraints on the parameters by passing
%              in the following:
%
%              minAmplitude=-inf (min amplitude of first gamma)
%              maxAmplitude=inf  (max amplitude of first gamma)
%              minTimelag=0 (min time lag of first gamma)
%              maxTimelag=5 (max time lag of first gamma)
%              minTau=0 (min tau for first gamma)
%              maxTau=inf (max tau for first gamma)
%              minAmplitudeRatio=0 (min amplitude ratio in fit2. Note that the amplitude
%                of the second gamma is always set to -amplitudeRatio times the amplitude
%                of the first gamma to insure that the gammas have opposite sign (usually
%                to fit the positive lobe followed by the longer undershoot).
%              maxAmplitudeRatio=1 (maximum amplitude ratio)
%              minTimelag2=0 ( min time lag for second gamma)
%              maxTimelag2=inf (max time lag for second gamma)
%              minTau2=0 (min tau for second gamma)
%              maxTau2=inf (max tau for second gamma)
% 
function d = fitTimecourse(timecourse,stimvol,framePeriod,varargin)

% display an already processed structure
if (nargin >= 1) && isstruct(timecourse)
  d = [];
  % grab arguments
  args = {};
  if nargin >= 4,args = varargin;end
  if nargin >= 3,args = {framePeriod args{:}};end
  if nargin >= 2,args = {stimvol args{:}};end
  % and display the fit
  dispComputedFit(timecourse,args);
  return
end

% check arguments
if (nargin < 2)
  d = [];
  help fitTimecourse
  return
end

global verbose;

% see if we are passed in argument list (as in a recursive call)
getArgs(varargin,{'args=[]'},'suppressUnknownArgMessage=1');
if isempty(args) args = varargin; end

% get arguments + set defaults
[argNames argValues args] = getArgs(args,{'hdrlen=20','concatInfo=[]','displayFit','fitType=glm','amplitudeType=default','returnAllFields=0','option=none','checkFit=0','canonicalResponse=[]','testRun=0','stdGroups=[]','verbose=1','displayNums=[]'});

% For testing. Call as fitTimecourse([],[],[],'testRun=1) or fitTimecourse([],[],[],'testRun=2');
if testRun
  d = testFitTimecourse(testRun);
  return
end

% setup a d structure to pass around
d = constructD(timecourse,stimvol,framePeriod,hdrlen,concatInfo,option,fitType,verbose);

% check timecourse size
if (ndims(timecourse) == 2) && ~any(size(timecourse)==1)
  % this is a kxn matrix, so need to run for each one.
  nVoxels = size(timecourse,1);
  % if we are just doing deconv with no fit, then pass this on to getr2timecourse
  if strcmp(fitType,'deconv') && strcmp(amplitudeType,'none')
    d.deconv = getr2timecourse(timecourse,d.nhdr,d.hdrlenTR,d.scm,d.framePeriod,d.verbose);
    return
  % some types are implemented for multiple voxel fits (like glmFit) let those through
  % for others, compute voxel by voxel
  elseif ~strcmp(fitType,'glm')
    % otherwise go voxel by voxel
    disppercent(-inf,sprintf('(fitTimecourse) Doing fit for %i voxels',nVoxels));
    for iVoxel = 1:nVoxels
      retd(iVoxel) = fitTimecourse(timecourse(iVoxel,:),stimvol,framePeriod,'args',varargin);
      disppercent(iVoxel/nVoxels);
    end
    disppercent(inf);
    d = retd;
    return
  end
end

% do the fits 
switch lower(fitType)
  case {'deconv','deconvolution'}
    d = deconvFit(d,amplitudeType,displayFit,checkFit,args);
  case {'glm'}
    d = glmFit(d,amplitudeType,option,displayFit,checkFit,canonicalResponse,stdGroups,verbose,displayNums);
  case {'timecourse','nonlin','nonlinear'}
    d = nonlinFit(d,amplitudeType,displayFit,checkFit,option,args);
  otherwise
    disp(sprintf('(fitTimecourse) Unknown type: %s',fitType));
end

% remove unnecessary fields
if ~returnAllFields
  d = removeFields(d,{'dim','stimvol','framePeriod','eachSCM','scm','hdrlen','hdrlenTR','nhdr','timecourse','concatInfo','eTimecourse','applyFiltering','deconvModel','zeroMean','fixedParams','useGLMforSTE'});
  if isfield(d,'deconv')
    d.deconv = removeFields(d.deconv,{'scm','volumes'});
  end
end
d.fitType = fitType;
if ~isfield(d,'amplitudeType'),d.amplitudeType = amplitudeType;end
d = orderfields(d);

%%%%%%%%%%%%%%%%
%%   glmFit   %%
%%%%%%%%%%%%%%%%
function d = glmFit(d,amplitudeType,option,displayFit,checkFit,canonicalResponse,stdGroups,verbose,displayNums)

if nargin < 8,verbose=1;end
% choose amplitude type
if ~any(strcmp(lower(amplitudeType),{'max','area','default'}))
  disp(sprintf('(fitTimecourse:glmFit) Unknown amplitdueType = %s',amplitudeType));
  return
end
% set default method
if strcmp(lower(amplitudeType),'default'),amplitudeType = 'max';end
d.amplitudeType = amplitudeType;

% get all stimvolumes
allStimvol = cell2mat(d.stimvol);

% now get the estimated mean response, i.e. the average one to
% all stimulus types, start by making a new stimulusConvolutionMatrix
% with all the stim volumes. We will call this the canonicalResponse
if ieNotDefined('canonicalResponse')
  d2 = d;d2.stimvol = [];
  d2.stimvol{1} = allStimvol;
  d2.dim(4) = d2.dim(2);d2.dim(1:3) = 1;
  d2 = makescm(d2,d.hdrlenTR,1);

  % then compute the estimated response
  canonicalResponse = getr2timecourse(d.timecourse,1,d.hdrlenTR,d2.scm,d.framePeriod,verbose);
  canonicalResponse = (canonicalResponse.ehdr/100)';

  % Normalize the mean response. This is so that the amplitude measures
  % come out meaning something. We take the area under the curve to be equal
  % to 1% signal change. This means that the amplitude measure is in %signal change/sec
  if any(strcmp(lower(amplitudeType),{'area'}))
    canonicalResponse = canonicalResponse/(sum(canonicalResponse)*length(canonicalResponse)*d.framePeriod);
    % if amplitudeType is set to max, then make the max of response equal to 1. This
    % will make the amplitude measure look like the peak response.
  elseif any(strcmp(lower(amplitudeType),{'max'}))
    canonicalResponse = canonicalResponse/max(canonicalResponse);
  end
else
  % make sure it is a column vector
  if size(canonicalResponse,1) == 1
    canonicalResponse = canonicalResponse';
  end
end

% check canonicalResponse
if size(canonicalResponse,1) ~= d.hdrlenTR
  disp(sprintf('(fitTimecourse:glmFit) Canonical response must be %i volumes long',d.hdrlenTR));
  return
end
d.canonicalResponse = canonicalResponse;

% now multiply the original scm by this mean response to construct
% the GLM
for i = 1:d.nhdr
  % convolve nResponse with each stimulus times to form a column
  if size(canonicalResponse,2) >= i
    % in this case we have a separate canonical for each stimulus type
    glm(:,i) = d.eachSCM{i}*canonicalResponse(:,i);
  else
    % just one canonical
    glm(:,i) = d.eachSCM{i}*canonicalResponse;
  end
  % now filter the column appropriately
  glm(:,i) = applyEventRelatedFiltering(glm(:,i)',d.concatInfo)';
end

% subtract off mean
glm = glm-repmat(mean(glm),size(glm,1),1);

% compute GLM
glmResponse = getr2timecourse(d.timecourse,1,d.nhdr,glm,d.framePeriod,verbose);

% get amplitudes
if d.dim(1) == 1
  d.amplitude = glmResponse.ehdr;
  d.amplitudeSTE = glmResponse.ehdrste;
else
  d.amplitude = glmResponse.ehdr';
  d.amplitudeSTE = glmResponse.ehdrste';
end
d.covar = glmResponse.covar;
d.r2 = glmResponse.r2;
d.glm = glm;
d.displayNums = displayNums;

% next we can compute what the deconvolved response of the glm fit look like
% and compare them to the actual deconvolution
if displayFit 
  dispGLMFit(d)
end

% now try to compute the standard deviation of response. 
if any(strcmp(lower(option),{'std'})) || ~isempty(stdGroups)
  if d.dim(1) == 1,disp(sprintf('(fitTimecourse:glmFit) std not implemented yet for multiple voxel fit'));keyboard;end
  disppercent(-inf,sprintf('(fitTimecourse) Computing std'));
  if isfield(d,'concatInfo'),n = d.concatInfo.n;else,n = 1;end
  allAmplitudes = [];
  for stimTypeNum = 1:d.nhdr
    allAmplitudes{stimTypeNum} = [];
  end
  % We are going to compute this run-by-run since it is more efficient
  for runNum = 1:n
    % first get the data for this run
    if isfield(d,'concatInfo')
      [thisTSeries thisStimvol thisConcatInfo] = extractRuns(runNum,runNum,d.timecourse,d.stimvol,d.concatInfo);
    else
      thisTSeries = d.timecourse;thisStimvol = d.stimvol;thisConcatInfo = d.concatInfo;
    end
    % now compute the glm design matrix for this run.
    glm = [];whichStimTypeNum = [];
    for stimTypeNum = 1:length(thisStimvol)
      % get the canonical response for this stimType
      if size(canonicalResponse,2) >= stimTypeNum
	thisCanonicalResponse = canonicalResponse(:,stimTypeNum);
      else
	thisCanonicalResponse = canonicalResponse;
      end
      % now make the columns of the glm for this stimType (one for 
      % each stimulus occurrence).
      for stimNum = 1:length(thisStimvol{stimTypeNum})
	% compute location where response goes
	responsePos = thisStimvol{stimTypeNum}(stimNum);
	responsePos = responsePos:min(responsePos+length(thisCanonicalResponse)-1,length(thisTSeries));
	% create this columns response by putting canonicalResponse in the correct position
	thisColumn = zeros(length(thisTSeries),1);
	thisColumn(responsePos) = thisCanonicalResponse(1:length(responsePos));
	% stick it into the glm
	glm(:,end+1) = thisColumn;
	% apply filtering
	glm(:,end) = applyEventRelatedFiltering(glm(:,end)',thisConcatInfo)';
	% remember which stimTypeNum this is form
	whichStimTypeNum(end+1) = stimTypeNum;
      end
    end

    if isempty(glm)
      disp(sprintf('(fitTimecourse) Run %i has no events, ignoring this run',runNum));
    else
      % remove mean
      glm = glm-repmat(mean(glm),size(glm,1),1);

      % compute GLM
      glmResponse = getr2timecourse(thisTSeries,1,size(glm,2),glm,d.framePeriod,verbose);
  
      % sort out amplitudes
      for stimNum = 1:length(whichStimTypeNum)
	allAmplitudes{whichStimTypeNum(stimNum)} = [allAmplitudes{whichStimTypeNum(stimNum)} glmResponse.ehdr(stimNum)];
      end
    end
    
    disppercent(runNum/n);
  end
  
  % get mean/std of computed amplitudes
  for stimTypeNum = 1:d.nhdr
    d.amplitudeSTD(stimTypeNum) = std(allAmplitudes{stimTypeNum});
    d.amplitudeMEAN(stimTypeNum) = mean(allAmplitudes{stimTypeNum});
  end

  d.allAmplitudes = allAmplitudes;
  
  % loop through STDgroups - this is a cell array that contains which stimulus types
  % belong to what groups. It is used for computing std across a set of responses by
  % subtracting the mean response for each response type and then putting all the
  % mean subtracted responses together into a big array and computing the standard
  % deviation over that
  for iSTDGroup = 1:length(stdGroups)
    thisGroup = stdGroups{iSTDGroup};
    thisGroupAmplitude = [];
    for iType = 1:length(thisGroup)
      % check response number
      if (thisGroup(iType) < 0) || (thisGroup(iType) > d.nhdr)
	disp(sprintf('(fitTimecourse:glmFit) STD Group %i has a response number %i which is out of range [1 %i]',iSTDGroup,thisGroup(iType),d.nhdr));
      else
	% if it is a valid response number, then subtract mean from the responses for this type
	thisGroupAmplitude = [thisGroupAmplitude (allAmplitudes{thisGroup(iType)}-d.amplitudeMEAN(thisGroup(iType)))];
      end
    end
    
    % when we are done, compute stats
    d.groupSTD.stdGroups{iSTDGroup} = stdGroups{iSTDGroup};
    d.groupSTD.amplitude{iSTDGroup} = thisGroupAmplitude;
    d.groupSTD.amplitudeMEAN(iSTDGroup) = mean(thisGroupAmplitude); % this should be 0
    d.groupSTD.amplitudeSTD(iSTDGroup) = std(thisGroupAmplitude);
  end
  disppercent(inf);
end

% checkthe fit
if checkFit
  % paramsDialog setup
  paramsInfo = {};
  paramsInfo{end+1} = {'responseNum',1,'incdec=[-1 1]',sprintf('minmax=[1 %i]',d.nhdr),'callback',@fitTimecoursePlotCallback,'callbackArg',d};
  paramsInfo{end+1} = {'amplitude',num2cell(d.amplitude),'group=responseNum','type=numeric','editable=0'};
  mrParamsDialog(paramsInfo,'Check fit');
  close(mlrSmartfig('fitTimecourseCheckFit','reuse'));
end
  
%%%%%%%%%%%%%%%%%%%
%%   deconvFit   %%
%%%%%%%%%%%%%%%%%%%
function d = deconvFit(d,amplitudeType,dispfit,checkFit,args)

% check for valid amplitudeType
if ~any(strcmp(lower(amplitudeType),{'fit1','fit1each','fit2','fit2each','default','none'}))
  disp(sprintf('(fitTimecourse:deconvFit) Unknown amplitdueType = %s',amplitudeType));
  return
end
if strcmp(lower(amplitudeType),'default'),amplitudeType = 'fit1';end

% do a deconvolution
if d.verbose,disp(sprintf('(fitTimecourse) Computing event-related analysis on time course'));end
d.deconv = getr2timecourse(d.timecourse-mean(d.timecourse),d.nhdr,d.hdrlenTR,d.scm,d.framePeriod,d.verbose);

if ~strcmp(amplitudeType,'none')

  % and fit the output of the deconvolution
  disppercent(-inf,sprintf('(fitTimecourse) Fitting gamma functions to results of event-related analysis'));

  % make a structure to do the fits in the same way that the nonlin fit works.
  % first make a timecourse that is the concatenation of the ehdrs from the deconv analysis
  fitd.timecourse = d.deconv.ehdr'/100;
  fitd.timecourse = fitd.timecourse(:)';
  fitd.dim = size(fitd.timecourse);
  fitd.nhdr = d.deconv.nhdr;
  fitd.hdrlenTR = d.hdrlenTR;
  fitd.framePeriod = d.framePeriod;
  % make the convolution matrices
  for stimNum = 1:fitd.nhdr
    fitd.stimvol{stimNum} = (stimNum-1)*fitd.hdrlenTR+1;
    stimVol = zeros(1,fitd.dim(2));
    stimVol(fitd.stimvol{stimNum}) = 1;
    fitd.eachSCM{stimNum} = stimconv(stimVol,fitd.hdrlenTR);
  end
  fitd = makescm(fitd,fitd.hdrlenTR);
  
  fitd.hdrlen = d.hdrlen;
  fitd.deconv = d.deconv;
  fitd.concatInfo = [];
  fitd.applyFiltering = 0;
  fitd.zeroMean = 1;
  fitd.deconvModel = 0;
  fitd.useGLMforSTE = 0;

  % now get the fit using nonlinFit routine
  fitd = nonlinFit(fitd,amplitudeType,dispfit,checkFit,[],args);

  % get fits
  d.p = fitd.p;
  d.covar = fitd.covar;
  d.ehdr = fitd.ehdr;
  d.time = fitd.time;
  d.amplitude = fitd.amplitude;
  d.amplitudeSTE = fitd.amplitudeSTE;
  disppercent(inf);
end

%%%%%%%%%%%%%%%%%%%
%%   nonlinFit   %%
%%%%%%%%%%%%%%%%%%%
function d = nonlinFit(d,amplitudeType,displayFit,checkFit,option,args)

% default for option
if ieNotDefined('option'),option='none';end

% declare fit constraints for single (or first gamma in fit2)
minAmplitude=[];
maxAmplitude=[];
minTimelag=[];
maxTimelag=[];
minTau=[];
maxTau=[];

% declare fit constraints for second gamma (only used in fit2)
minAmplitudeRatio=[];
maxAmplitudeRatio=[];
minTimelag2=0;
maxTimelag2=inf;
minTau2=0;
maxTau2=inf;

% check for fit constraints in set options
getArgs(args,{'minAmplitude=0','maxAmplitude=inf','minTimelag=0','maxTimelag=5','minTau=0','maxTau=inf','minAmplitudeRatio=0','maxAmplitudeRatio=1','minTimelag2=5','maxTimelag2=inf','minTau2=0','maxTau2=inf'});

% choose amplitude type
if ~any(strcmp(lower(amplitudeType),{'fit1','fit1each','fit1fixed','fit2','fit2each','default'}))
  disp(sprintf('(fitTimecourse:nonlinFit) Unknown amplitudeType = %s',amplitudeType));
  return
end
% set default method
if strcmp(lower(amplitudeType),'default'),amplitudeType = 'fit2';end
d.amplitudeType = amplitudeType;

% reset the counter
global gComputeResidualCount;gComputeResidualCount = 0;

% initialize parameters for single gamma
amplitude = 1;
timelag = 1;
tau = 0.5;
d.fixedParams.exponent = 6;

% initialize parameters for second gamma, if needed
if strncmp(lower(amplitudeType),'fit2',length('fit2'))
  % the amplitude ratio is the ratio between the two
  % amplitudes, the parameters are set up so that the
  % second gamma is always of opposite sign to the
  % first gamma and the ratio of its amplitude to
  % the first gamma is between 0 and 1
  amplitudeRatio = 0.5;
  timelag2 = 4;
  tau2 = 0.5;
  d.fixedParams.exponent2 = 6;
  % use glm to get standard errors
  if ~isfield(d,'useGLMforSTE'),d.useGLMforSTE = 1;end
end
if ~isfield(d,'useGLMforSTE'),d.useGLMforSTE = 0;end

% now, depending on the fit type we set up the parameters differently
% fit1 is a single gamma fit in which only the amplitude is allowed
% to scale between stimulus types.
if strcmp(lower(amplitudeType),'fit1')
  initParams = [repmat(amplitude,1,d.nhdr) timelag tau];
  minfit = [repmat(minAmplitude,1,d.nhdr) minTimelag minTau];
  maxfit = [repmat(maxAmplitude,1,d.nhdr) maxTimelag maxTau];
% fit1each is a single gamma fit in which everything is allowed
% to scale between stimulus types.
elseif strcmp(lower(amplitudeType),'fit1each')
  initParams = repmat([amplitude;timelag;tau],1,d.nhdr);
  minfit = repmat([minAmplitude;minTimelag;minTau],1,d.nhdr);
  maxfit = repmat([maxAmplitude;maxTimelag;maxTau],1,d.nhdr);
% fit2 is a difference of gamma fit in which only the amplitudes can scale
% between stimulus types
elseif strcmp(lower(amplitudeType),'fit2')
  initParams = [repmat(amplitude,1,d.nhdr) amplitudeRatio timelag timelag2 tau tau2];
  minfit = [repmat(minAmplitude,1,d.nhdr) minAmplitudeRatio minTimelag minTimelag2 minTau minTau2];
  maxfit = [repmat(maxAmplitude,1,d.nhdr) maxAmplitudeRatio maxTimelag maxTimelag2 maxTau maxTau2];
% fit2each is a difference of gamma fit where the parameters are all allowed
% to change for each stimulus type.
elseif strcmp(lower(amplitudeType),'fit2each')
  initParams = repmat([amplitude amplitudeRatio timelag timelag2 tau tau2],1,d.nhdr);
  minfit = repmat([minAmplitude minAmplitudeRatio minTimelag minTimelag2 minTau minTau2],1,d.nhdr);
  maxfit = repmat([maxAmplitude maxAmplitudeRatio maxTimelag maxTimelag2 maxTau maxTau2],1,d.nhdr);
end

% maximum number of interations allowed
if ieNotDefined('maxiter');maxiter = inf;end

% compute deconvolution. This is for making figures
if ~isfield(d,'deconv') && displayFit
  if d.zeroMean
    d.deconv = getr2timecourse(d.timecourse+1,d.nhdr,d.hdrlenTR,d.scm,d.framePeriod,d.verbose);
  else
    d.deconv = getr2timecourse(d.timecourse,d.nhdr,d.hdrlenTR,d.scm,d.framePeriod,d.verbose);  
  end
end

% optimization parameters
optimParams = optimset('MaxIter',maxiter,'Display','off');

% make sure that the timecourse is double, since lsqnonlin chokes if it is not.
d.timecourse = double(d.timecourse);

% use lsqnonlin to find function minimum
[fitParams resnorm residual exitflag output lambda jacobian] = lsqnonlin(@computeResidual,initParams,minfit,maxfit,optimParams,d,amplitudeType,displayFit);

% extract the parameters into a data structure
d.p = getFitParams(fitParams,amplitudeType,d.nhdr,d.fixedParams);

% compute the standard errors of the fits. Note that this definition
% is only good if the measurement errors are normally distributed.
% it also assumes that we have a significant goodness of fit.
% That is, chi-squared is the sum-squared difference between model
% and data divided by the noise std of each measurement. Since
% we don't know this std, we just set it to 1 (i.e. equal across
% all measurements and of value 1) to do the fitting. But this means
% that we can't put the chi-squared value into the chi-squared distribution
% to get the correct p-value for the goodness-of-fit. Instead, we compute
% the noise estimate, by taking the sum-squared residual and dividing it
% by the number of time points minus the number of paramters and use that instead. 
d.covar = inv(jacobian'*jacobian);
noiseVariance = (residual*residual')/(d.dim(end)-length(initParams));
d.pSte = getFitParams(diag(sqrt(noiseVariance * d.covar))',amplitudeType,d.nhdr,d.fixedParams);

% get the last eTimecourse
gComputeResidualCount = 0;
[residual d.eTimecourse] = computeResidual(fitParams,d,amplitudeType,displayFit);
d.r2 = 1-var(residual)/var(d.timecourse);

% now get the fits
if d.deconvModel
  % we deconvolve the model timecourse to get the fit, this is so that it matches
  % what we get in the deconvolution.
  d.fit = getr2timecourse(d.eTimecourse-mean(d.eTimecourse)+1,d.nhdr,d.hdrlenTR,d.scm,d.framePeriod,d.verbose);
  d.ehdr = d.fit.ehdr;
  d.time = d.fit.time;
else
  % just directly get smoother version of the best fit
  d.time = 0:d.framePeriod/20:(d.hdrlen+d.framePeriod);
  for i = 1:d.nhdr
    d.ehdr(i,:) = getThisGamma(d.time,d.p,i)*100;
  end
end

% get canonical resonses (i.e. the found fits to the responses. THis
% is used for amplitude estimation with GLM or for std estimation)
if d.useGLMforSTE || any(strcmp(lower(option),{'std'}))
  pCanonical = d.p;
  for i = 1:d.nhdr
    canonical(:,i) = getThisGamma(d.time,pCanonical,1)'*100;
    % make sure it is right side up
    canonical(:,i) = canonical(:,i)*sign(d.p.amplitude(i));
    %scale it
    canonical(:,i) = canonical(:,i)/max(canonical(:,i));
  end
end

% next we extract the amplitdue values. This can be as easy as just
% using the amplitude parameter from the single gamma fit, but
% if this is a diffOfGammaFit then we get the amplitude by using
% the glm fit with the found canonical response to get the amplitudes
%if d.p.diffOfGammaFit
if d.useGLMforSTE
  % now to compute the amplitudes, get the response we have just computed
  % and pump them through the glm code to get amplitudes and amplitude
  % estimates
  disp(sprintf('(fitTimecourse:nonlinFit) Using glm to get standard errors'));
  glmAmp = glmFit(d,'default',option,0,0,canonical,[],1,[]);
  d.amplitude = glmAmp.amplitude;
  d.amplitudeSTE = glmAmp.amplitudeSTE;
  if isfield(glmAmp,'amplitudeSTD')
    d.amplitudeMEAN = glmAmp.amplitudeMEAN;
    d.amplitudeSTD = glmAmp.amplitudeSTD;
  end
else
  d.amplitude = d.p.amplitude;
  d.amplitudeSTE = d.pSte.amplitude;
end

% display the fit if called for
if displayFit
  titlestr = sprintf('A=%s lag=%s offset=%s tau=%s exponent=%s',num2str(d.p.amplitude),num2str(d.p.timelag),num2str(d.p.offset),num2str(d.p.tau),num2str(d.p.exponent));
  plotFit(titlestr,d.deconv,d.amplitude,d.ehdr,d.time);
end

% if we are to compute the std, and haven't done it already (through the
% GLM estimation above. Do it here.
if any(strcmp(lower(option),{'std'})) && ~isfield(d,'amplitudeSTD')
  disp(sprintf('(fitTimecourse:nonlinFit) Using glm to get std'));
  glmAmp = glmFit(d,'default',option,0,0,canonical);
  d.amplitudeMEAN = glmAmp.amplitudeMEAN;
  d.amplitudeSTD = glmAmp.amplitudeSTD;
end

% display a dialog to allow user to check the fit
if checkFit
  % paramsDialog setup
  paramsInfo = {};
  paramsInfo{end+1} = {'responseNum',1,'incdec=[-1 1]',sprintf('minmax=[1 %i]',d.nhdr),'callback',@fitTimecoursePlotCallback,'callbackArg',d};
  paramsInfo{end+1} = {'amplitude',num2cell(d.p.amplitude),'group=responseNum','type=numeric','editable=0'};
  paramsInfo{end+1} = {'timelag',num2cell(d.p.timelag),'group=responseNum','type=numeric','editable=0'};
  paramsInfo{end+1} = {'offset',num2cell(d.p.offset),'group=responseNum','type=numeric','editable=0'};
  paramsInfo{end+1} = {'tau',num2cell(d.p.tau),'group=responseNum','type=numeric','editable=0'};
  paramsInfo{end+1} = {'exponent',num2cell(d.p.exponent),'group=responseNum','type=numeric','editable=0'};
  if d.p.diffOfGammaFit
    paramsInfo{end+1} = {'amplitude2',num2cell(d.p.amplitude2),'group=responseNum','type=numeric','editable=0'};
    paramsInfo{end+1} = {'timelag2',num2cell(d.p.timelag2),'group=responseNum','type=numeric','editable=0'};
    paramsInfo{end+1} = {'offset2',num2cell(d.p.offset2),'group=responseNum','type=numeric','editable=0'};
    paramsInfo{end+1} = {'tau2',num2cell(d.p.tau2),'group=responseNum','type=numeric','editable=0'};
    paramsInfo{end+1} = {'exponent2',num2cell(d.p.exponent2),'group=responseNum','type=numeric','editable=0'};
  end
  mrParamsDialog(paramsInfo,'Check fit');
  close(mlrSmartfig('fitTimecourseCheckFit','reuse'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%   computeResidual   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function [residual,eTimecourse] = computeResidual(params,d,amplitudeType,displayFit)

global verbose;
% This function takes the timecourse, the stimulus convolution matrix, the
% amplitude, timelag, offset, tau and exponent and convolves a gamma function
% with these parameters with the stimulus timings and computes the difference
% between that ideal time course and the real timecourse. 

% keep track of how often this function is called
global gComputeResidualCount;
gComputeResidualCount = gComputeResidualCount+1;

% get the fit params into a structure
p = getFitParams(params,amplitudeType,d.nhdr,d.fixedParams);

% compute the estimate given these parameters
[eTimecourse ehdr] = computeEstimatedTimecourse(d,p);

% and the residual
residual = d.timecourse-eTimecourse;

% display something
if (mod(gComputeResidualCount,length(params)) == 1) && displayFit
  mlrSmartfig('fitTimecourse','reuse');clf
  % get the fit
  if d.deconvModel
    % normally, when filtering is applied, we grab the fit by deconvolving the
    % eTimecourse, this is done because it will make the fit match the deconvolution
    % properly.
    fit = getr2timecourse(eTimecourse-mean(eTimecourse)+1,d.nhdr,d.hdrlenTR,d.scm,d.framePeriod,verbose);
  else
    % if we don't need to do this, then just get the fits directly from eTimecourse
    fit.ehdr = 100*reshape(eTimecourse,d.hdrlenTR,d.nhdr)';
  end
  titlestr = sprintf('%i: A=%s lag=%s offset=%s tau=%s exponent=%s\n(residual=%0.4f)',gComputeResidualCount,num2str(p.amplitude),num2str(p.timelag),num2str(p.offset),num2str(p.tau),num2str(p.exponent),sqrt(sum(residual.^2)));
  plotFit(titlestr,d.deconv,p.amplitude,fit.ehdr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   computeEstimatedTimecourse   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eTimecourse ehdr] = computeEstimatedTimecourse(d,p)

% now compute the estimated response functions (i.e. the gamma functions we
% are using to fit the timecourse)
time = d.framePeriod/2:d.framePeriod:d.hdrlen+d.framePeriod/2;
eTimecourse = zeros(1,d.dim(end));
for i = 1:d.nhdr
  ehdr = getThisGamma(time,p,i);
  % and convolve those against the stimvols for each stimulus type
  eTimecourse = eTimecourse + (d.eachSCM{i}*ehdr')';
end

% now apply the detrending / hipassfilter
if d.applyFiltering
  eTimecourse = applyEventRelatedFiltering(eTimecourse,d.concatInfo);
end

%%%%%%%%%%%%%%%%%%%%%%
%%   getThisGamma   %%
%%%%%%%%%%%%%%%%%%%%%%
function fun = getThisGamma(time,p,i)

fun = thisGamma(time,p.amplitude(i),p.timelag(i),p.offset(i),p.tau(i),p.exponent(i))/100;
% add second gamma if this is a difference of gammas fit
if p.diffOfGammaFit
  fun = fun + thisGamma(time,p.amplitude2(i),p.timelag2(i),p.offset2(i),p.tau2(i),p.exponent2(i))/100;
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

%%%%%%%%%%%%%%%%%%%%%%
%%   getFitParams   %%
%%%%%%%%%%%%%%%%%%%%%%
function p = getFitParams(params,amplitudeType,nhdr,fixedParams)

% extract the fit params into a more understandable data structure

% get the parameters, this will be different for different amplitudeTypes 
% (see above function nonlinFit for explanation of different amplitudeTypes
if strcmp(lower(amplitudeType),'fit1')
  p.amplitude = params(1:nhdr);
  p.timelag = params(nhdr+1);
  p.offset = 0;
  p.tau = params(nhdr+2);
  p.exponent = fixedParams.exponent;
  p.diffOfGammaFit = 0;
elseif strcmp(lower(amplitudeType),'fit1each')
  params = reshape(params,3,nhdr);
  p.amplitude = params(1,:);
  p.timelag = params(2,:);
  p.offset = 0;
  p.tau = params(3,:);
  p.exponent = fixedParams.exponent;
  p.diffOfGammaFit = 0;
elseif strcmp(lower(amplitudeType),'fit2')
  p.amplitude = params(1:nhdr);
  amplitudeRatio = params(nhdr+1);
  p.timelag = params(nhdr+2);
  p.timelag2 = params(nhdr+3);
  p.offset = 0;
  p.offset2 = 0;
  p.tau = params(nhdr+4);
  p.tau2 = params(nhdr+5);
  p.exponent = fixedParams.exponent;
  p.exponent2 = fixedParams.exponent2;
  p.diffOfGammaFit = 1;
elseif strcmp(lower(amplitudeType),'fit2each')
  params = reshape(params,6,nhdr);
  p.amplitude = params(1,:);
  amplitudeRatio = params(2,:);
  p.timelag = params(3,:);
  p.timelag2 = params(4,:);
  p.offset = 0;
  p.offset2 = 0;
  p.tau = params(5,:);
  p.tau2 = params(6,:);
  p.exponent = fixedParams.exponent;
  p.exponent2 = fixedParams.exponent2;
  p.diffOfGammaFit = 1;
end

% make sure all parameters have settings for each stimulus type
if length(p.amplitude) < nhdr,p.amplitude = repmat(p.amplitude,1,nhdr);end
if length(p.timelag) < nhdr,p.timelag = repmat(p.timelag,1,nhdr);end
if length(p.offset) < nhdr,p.offset = repmat(p.offset,1,nhdr);end
if length(p.tau) < nhdr,p.tau = repmat(p.tau,1,nhdr);end
if length(p.exponent) < nhdr,p.exponent = repmat(p.exponent,1,nhdr);end

% do the same for the other gamma (for difference of gammas)
if p.diffOfGammaFit
  if length(amplitudeRatio) < nhdr,amplitudeRatio = repmat(amplitudeRatio,1,nhdr);end
  % to get the amplitude of the 2nd gamma, we multiply by the first
  % ampitude by the amplitudeRatio (the neg sign is to insure that
  % the second amplitude is always in the opposite direction)
  p.amplitude2 = -p.amplitude.*amplitudeRatio;
  if length(p.timelag2) < nhdr,p.timelag2 = repmat(p.timelag2,1,nhdr);end
  if length(p.offset2) < nhdr,p.offset2 = repmat(p.offset2,1,nhdr);end
  if length(p.tau2) < nhdr,p.tau2 = repmat(p.tau2,1,nhdr);end
  if length(p.exponent2) < nhdr,p.exponent2 = repmat(p.exponent2,1,nhdr);end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   fitTimecoursePlotCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fitTimecoursePlotCallback(d,params)

% callback from checkFit parameters dialog
mlrSmartfig('fitTimecourseCheckFit','reuse');
plotFit('',d.deconv,d.amplitude,d.ehdr,d.time,params.responseNum);

%%%%%%%%%%%%%%%%%
%%   plotFit   %%
%%%%%%%%%%%%%%%%%
function plotFit(titlestr,deconv,amplitude,fitehdr,fittime,whichOne,voxNum)

if ieNotDefined('fittime'),fittime=deconv.time;end
if ieNotDefined('whichOne'),whichOne = 1:deconv.nhdr;end
if ieNotDefined('voxNum'),voxNum = 1;end
if ndims(deconv.ehdr) == 3
  titlestr = sprintf('%s (Voxel %i)',titlestr,voxNum);
end
cla;
legendStr = {};legendSymbol = {}; r2 = [];
for i = whichOne
  color = getSmoothColor(find(i==whichOne),length(whichOne),'hsv');
  if ndims(deconv.ehdr) == 3
    % many voxels, choose the right voxel to plot
    dispData = squeeze(deconv.ehdr(voxNum,i,:));
    dispDataSTE = squeeze(deconv.ehdrste(voxNum,i,:));
    dispFitData = squeeze(fitehdr(voxNum,i,:));
    dispAmp = amplitude(voxNum,i);
  else
    % single voxel
    dispData = deconv.ehdr(i,:);
    dispDataSTE = deconv.ehdrste(i,:);
    dispFitData = fitehdr(i,:);
    dispAmp = amplitude(i);
  end
  % plot
  myerrorbar(deconv.time,dispData,'yError',dispDataSTE,'Symbol','o','MarkerFaceColor',color,'MarkerEdgeColor','w');
  plot(fittime,dispFitData,'-','Color',color);
  % legend
  legendStr{end+1} = sprintf('%f',dispAmp);
  legendSymbol{end+1} = {'ko' color};
  % compute the fit r2
  fitMatch = interp1(fittime,dispFitData,deconv.time);
  r2(end+1) = 1 - var(dispData(:)-fitMatch(:))/var(dispData(:));
end
mylegend(legendStr,legendSymbol);
xlabel('time (sec)');
ylabel('response magnitude');
title(sprintf('%s (r2=%s)',titlestr,num2str(r2,'%0.2f ')));
drawnow

%%%%%%%%%%%%%%%%%%%%
%%   constructD   %%
%%%%%%%%%%%%%%%%%%%%
function d = constructD(timecourse,stimvol,framePeriod,hdrlen,concatInfo,option,fitType,verbose)

% make sure timecourse looks ok.
if (size(timecourse,2) == 1) && (size(timecourse,1) > 1)
  mrWarnDlg(sprintf('(fitTimecourse) Timecourse found to be %ix1 instead of 1x%i. Taking transpose',size(timecourse,1),size(timecourse,1)));
  timecourse = timecourse';
end

% first get stimulus convolution matrix for each condition separately
d.dim(4) = size(timecourse,2);
% note that this hdrlen is different from
% below, becuase it is in # of TR's not
% in seconds.
hdrlenTR = floor(hdrlen/framePeriod)+1;
d.hdrlen = hdrlenTR;
d.concatInfo = concatInfo;
d.nFrames = d.dim(4);

if verbose
  disppercent(-inf,'(fitTimecourse) Constructing stimulus convolution matrices for each condition');
end
for stimNum = 1:length(stimvol)
  d.stimvol{1} = stimvol{stimNum};
  d = makescm(d,[],0);
  eachSCM{stimNum} = d.scm;
  if verbose,disppercent(stimNum/length(stimvol));end
end
if verbose,disppercent(inf);end

% make an scm for the whole design
if verbose
  disppercent(-inf,'(fitTimecourse) Constructing stimulus convolution matrix for whole design');
end
d.stimvol = stimvol;
d = makescm(d,[],1);
allSCM = d.scm;
if verbose,disppercent(inf);end

% make an scm for each stimulus 
% now make a d structure to pass around
clear d;
d.dim = size(timecourse);
d.stimvol = stimvol;
d.framePeriod = framePeriod;
d.concatInfo = concatInfo;
d.eachSCM = eachSCM;
d.scm = allSCM;
d.hdrlen = hdrlen;
d.hdrlenTR = hdrlenTR;
d.nhdr = length(stimvol);
d.timecourse = timecourse;
d.applyFiltering = 1;
d.zeroMean = 0;
d.deconvModel = 1;
d.nFrames = max(d.dim);
d.verbose = verbose;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   testFitTimecourse   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = testFitTimecourse(whichType)

if ieNotDefined('whichType'),whichType = 1;end

if whichType==1
  % make random stimulus volumes
  nStim = 50;
  nStimTypes = 6;
  isi = 25;
  thisStimvol = 1;
  disppercent(-inf,'Making stimulus volumes');
  for i = 1:nStim
    for j = 1:nStimTypes
      d.stimvol{j}(i) = thisStimvol;
      thisStimvol = thisStimvol + isi(ceil(rand*length(isi)));
    end
    disppercent(i/nStim);
  end
  d.dim(4) = thisStimvol+isi(1);
  d.tr = 1;
  d = makescm(d,[],1);
  hdrs = [];
  for i = 1:nStimTypes
    gammaFun = thisGamma(1:d.hdrlen,0.5+i/nStimTypes,0.5,0,0.75,6)/100;
    hdrs = [hdrs gammaFun];
  end
  d.timecourse = (d.scm*hdrs')';
%  d.timecourse = d.timecourse + 0.05*rand(size(d.timecourse))/100;
  d.timecourse = applyEventRelatedFiltering(d.timecourse');

%  d.timecourse = eventRelatedDetrend(d.timecourse')';
%  d.timecourse = d.timecourse - mean(d.timecourse)+1;
  disppercent(inf);
  d = fitTimecourse(d.timecourse,d.stimvol,d.tr,'fitType=glm','amplitudeType=max');
else

  % open a view to test data
  cd ~/data/nyu/jg060918_mr4;
  mrQuit;
  v = newView;
  v = viewSet(v,'curGroup','Concatenation');
  v = viewSet(v,'curScan',1);

  % load the roi tSeries
  roi = loadROITSeries(v,'r_mt',1,3);
  tSeries = mean(roi.tSeries);

  % get concatInfo
  concatInfo = viewGet(v,'concatInfo',1,3);

  % get the stimvols
  scanParams = getEventRelatedVarname(v);
  d = loadScan(v,1,3,0);
  d = getStimvol(d,scanParams{1});
  stimvol = d.stimvol;
  framePeriod = viewGet(v,'framePeriod');
  d = fitTimecourse(tSeries,stimvol,framePeriod,'concatInfo',concatInfo,'returnAllFields','fitType=glm');
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%    dispComputedFIr    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function dispComputedFit(d,args)

getArgs(args,{'subplot1=[1 2 1]','subplot2=[1 2 2]','titleStr=[]'});

if ~isfield(d,'dim')
  disp(sprintf('(fitTimecourse:dispGLMFit) This structure does not have all necessary fields for display. Try running fitTimecourse with the parameter returnAllFields set to true'));
  return
end

% display
if strcmp(d.fitType,'glm')
  dispGLMFit(d,subplot1,subplot2,titleStr);
else
  disp(sprintf('(fitTimecourse:dispComputedFit) Display for fit type %s not implemented yet',d.fitType));
end


%%%%%%%%%%%%%%%%%%%%
%    dispGLMFIt    %
%%%%%%%%%%%%%%%%%%%%
function dispGLMFit(d,subplot1,subplot2,titleStr)

% default subplots
if nargin == 1
  subplot1 = [1 2 1];
  subplot2 = [1 2 2];
  titleStr = [];
end

% now compute an estimated timecourse from the glm fit
glmTimecourse = d.glm*(d.amplitude'/100)+1;

% and deconvolve that
glmDeconv = getr2timecourse(glmTimecourse',d.nhdr,d.hdrlenTR,d.scm,d.framePeriod,d.verbose);

% keep fits
d.ehdr = glmDeconv.ehdr;
d.time = glmDeconv.time;

% get deconv of data
if ~isfield(d,'deconv')
  d.deconv = getr2timecourse(d.timecourse,d.nhdr,d.hdrlenTR,d.scm,d.framePeriod,d.verbose);
end

% display the fit
if nargin == 1
  mlrSmartfig('fitTimecourse','reuse');
  clf;
end
subplot(subplot1(1),subplot1(2),subplot1(3));
[maxVal maxIndex] = max(d.deconv.r2);
if isempty(titleStr),titleStr = 'Deconvolution and GLM fit';end
plotFit(sprintf('%s (r2=%0.2f)',titleStr,maxVal),d.deconv,d.amplitude,glmDeconv.ehdr,glmDeconv.time,d.displayNums,maxIndex);
subplot(subplot2(1),subplot2(2),subplot2(3));
plot(glmDeconv.time,d.canonicalResponse,'k.-');
xlabel('Time (sec)');
ylabel('Response Amplitude');
if ~isfield(d,'canonicalFit') || isempty(d.canonicalFit)
  title('Canonical response');
else
  p = d.canonicalFit.p;
  title(sprintf('Canonical response: A: %0.2f timeLag: %0.2f offset: %0.2f tau: %0.2f exponent: %i',p.amplitude,p.timelag,p.offset,p.tau,p.exponent));
end
drawnow;
%  makeEqualYaxis(1,2);

