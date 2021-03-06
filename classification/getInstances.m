% getInstances.m
%
%        $Id: getInstances.m,v 1.5 2008/10/23 01:52:48 justin Exp $
%      usage: rois = getInstances(v,rois,stimvol,<startLag>,<blockLen>,<type>)
%         by: justin gardner
%       date: 09/29/08
%    purpose: compute instances.
%
%             The default is to take mean after stimulus occurrences, with options:
%               startLag: number of volumes after stimulus occurrence to take mean from (set to 0 for no lag).
%                         Default is to take the stimvol that occurs 3 seconds after stimulus onset.
%               blockLen: number of volumes to compute mean over. Default is to take the median difference
%                         in volumes between stimulus occurences (for a block design as this is intended
%                         for this gives the block length - median is taken to make sure to ignore
%                         any few blocks that might not be a full length - e.g. if the scan ended before
%                         the end of the block).
%               groupTrials: Groups trials into k trials each. That is to compute each instance it
%                            takes the mean across k randomly chosen trials (no replacement). This
%                            might improve snr a bit by averaging noise out from individual trials.
%               minResponseLen: The minimum response length needed to use a trial. Defaults to half
%                                the blockLen
%
%             You can also use a "glm" (or "decon", the same thing) method. This first collapses
%             all of the stimulus volumes into one and computs an average hemodynamic
%             response. It then uses that as the "canonical response" and then computes
%             the amplitude of that response that best accounts for each response. Note
%             that the canonical is computed on a voxel-by-voxel basis. You can use this
%             method with the input argument:
%                type='glm' (or type='decon')
%             Note that there are some other options. The default is to fit the canonical form
%             with a difference of gamma (so that we have a smooth curve). This can be changed:
%                'canonicalType=allfit2'
%             see the function getCanonical for options. You can all set an r2cutoff for
%             using voxels in the mean to compute the canonical response:
%                'r2cutoff=0.1';
%
%             All of the methods above have the following options:
%               n: number of voxels to use (voxels are used in order specified
%                  by the roi field sortindex -- see getSortIndex). Defaults to 100.
%               fieldName: Name of field in which instances are stored (default=instance)
%               verbose: set to 1 for standard messages, set to 2 to display messages about
%                        which trials are being dropped and other more detailed info
%    Note on cross validation:
%         When using the glm/deconv method to obtain instances you would want to first use the above
%         method to obtain train instances, and then pass the train data as an input to obtain
%         test instances. This method will use the canonical response estimated using the training
%         data when computing the test instances. For this usage, do the following:
%         testRois = getInstances(v, testRois, testStimvols, trainRois, <other options>);


function rois = getInstances(v,rois,stimvol,varargin)

% check arguments
if any(nargin == [0])
  help getInstances
  return
end


if isempty(varargin) || ~iscell(varargin{1})
  % get type of instance getting
  % make sure that rois is a cell array of rois and check their validity
  rois = cellArray(rois);
  type = [];
  [argNames argValues args] = getArgs(varargin,{'type=mean'});
  [tf rois] =  validateROIs(rois,stimvol,args);
  if ~tf,return,end
  % and get the instances accordingly
  switch lower(type)
    case {'mean'}
      rois = getInstancesMean(rois,stimvol,args);
    case {'deconv','glm'} %TSL:why have two possible names to mean the same thing?
      rois = getInstancesDeconv(v,rois,stimvol,args);
    otherwise
      disp(sprintf('(getInstance) Unknown type: %s',type));
  end
elseif isroi(varargin{1}{1})
  buildROIs=varargin{1};
  rois=getMatchedInstances(v,rois,stimvol,buildROIs,varargin{2:end});
else
  disp('Input to getInstance not correct, see help.');
  keyboard;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getMatchedInstances   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get test instances based on parameters of build instances.
function rois=getMatchedInstances(v,rois,stimvol,buildROIs,varargin)
getArgs(varargin,{'fieldName=instance'});

for i=1:length(buildROIs)
  if ~strcmp(rois{i}.name, buildROIs{i}.name) || length(rois)~=length(buildROIs)
    disp('(getInstances:getMatchedInstances) Training and test ROIs have different names or number, data need to be properly preprocessed, see splitRuns');
    keyboard;
  end
  
  % get the arguments from the buildROI  and copy necessary fields from buildROI
  thisROI=cellArray(rois{i});
  instanceArgs = buildROIs{i}.(fieldName).info.args;
  [tf thisROI] = validateROIs(thisROI,stimvol,instanceArgs);
  
  thisROI{1}.sortindex=buildROIs{i}.sortindex;

  type=buildROIs{i}.(fieldName).info.type;
  switch lower(type)
    case {'mean'}
      thisROI = getInstancesMean(thisROI,stimvol,instanceArgs);
    case {'deconv','glm'}
      %copy canonical response (and fit) from build to test ROI
      thisROI{1}.canonicalResponse=buildROIs{i}.canonicalResponse; %TSL:no need to pass this in as arg
      if isfield(buildROIs{i},'canonicalFit')
        thisROI{1}.canonicalFit=buildROIs{i}.canonicalFit;
      end
      instanceArgs = {instanceArgs{:} 'setSortIndex=0'};  % setSortIndex has to be 0 (because we don't want alter the sortIndex copied from the build ROI 
      thisROI = getInstancesDeconv(v,thisROI,stimvol,instanceArgs); 
    otherwise
      disp(sprintf('(getInstance) Unknown type: %s',type));
  end
  rois{i}=thisROI{1}; %put the cell array as a struct into rois
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getInstancesDeconv    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rois = getInstancesDeconv(v,rois,stimvol,args)  %TSL: v is not necessary?

% get some default arguments
getArgs(args,{'n=100','fieldName=instance','verbose=1','hdrlen=19','canonicalType=allfit2','r2cutoff=[]','pval=[]','displayFit=0','makeAverage=0','displayNums=[]','saveFitStructures=1','setSortIndex=0'});

% make a stimvol array with one field for each stimulus type
if ~makeAverage  
  stimvolAll = num2cell(cell2mat(stimvol));
else
  stimvolAll = stimvol;
end

% go through each ROI
% canonicalFit = [];
for iROI = 1:length(rois)
  % if canonical is to be calculated across whole roi
  if ~isfield(rois{iROI},'canonicalResponse') || isempty(rois{iROI}.canonicalResponse)
    if any(strcmp(canonicalType,{'all','allfit2','allfit1','none'}))
%       [rois{iROI}.canonicalResponse r2 rois{iROI}.canonicalFit] = getCanonical(v,rois{iROI},stimvol,canonicalType,'hdrlen',hdrlen,'r2cutoff',r2cutoff,'normType=max','pval',pval);
      rois{iROI}.canonicalResponse = getCanonical(v,rois{iROI},stimvol,canonicalType,'hdrlen',hdrlen,'r2cutoff',r2cutoff,'normType=max','pval',pval,'n',n);
    else
      rois{iROI}.canonicalResponse = [];
    end
  else
    disp(sprintf('(getInstances:getInstancesDeconv) Using precomputed canonicalResponse'));
  end
  
  % reject and don't create instances if canonicalResponse was not fit
  if isempty(rois{iROI}.canonicalResponse)
    disp(sprintf('(getInstances) Canonical response was not well fit. Rejecting for ROI %s and moving on',rois{iROI}.name));
    continue;
  end
  
  % set the sort index if desired, use r2 value by default
  if setSortIndex
    rois{iROI}.sortindexType = 'r2';
    [rois{iROI}.r2sorted rois{iROI}.sortindex] = sort(r2(:),1,'descend');
    rois{iROI}.r2 = r2;
  end
  % set up fields
  rois{iROI}.(fieldName).info.type = 'deconv';
  rois{iROI}.(fieldName).info.args = args;
  disppercent(-inf,sprintf('(getInstances) Computing amplitudes for ROI %s',rois{iROI}.name));
  
  % fit glm with the canonical we just computed
  d = fitTimecourse(rois{iROI}.tSeries,stimvolAll,rois{iROI}.framePeriod,'concatInfo',rois{iROI}.concatInfo,'fitType=glm','displayFit',displayFit,'verbose=0',sprintf('hdrlen=%s',mynum2str(hdrlen)),'canonicalResponse',rois{iROI}.canonicalResponse,'displayNums',displayNums,'returnAllFields',saveFitStructures);
  
  % now sort the magnitudes back
  sortindex = rois{iROI}.sortindex(1:min(n,rois{iROI}.n));
  for iType = 1:length(stimvol)
    if makeAverage
      % in this case each instance is made as average across all
      rois{iROI}.(fieldName).instances{iType}(1,:) = d.amplitude(sortindex,iType);
      rois{iROI}.(fieldName).instancesSTE{iType}(1,:) = d.amplitudeSTE(sortindex,iType);
    else
      % normally this is run to make one instance per each showing of the stimulus
      for iInstance = 1:length(stimvol{iType}) %TSL:this loop through instances can be avoided using [temp thisIndex]=ismember(stimvol{iType},cell2mat(stimvolAll); thisAmp=d.amplitude(:,thisIndex);
        % get the stimvol
        thisIndex = find(stimvol{iType}(iInstance) == cell2mat(stimvolAll));
        thisAmplitude = d.amplitude(:,thisIndex);
        thisAmplitudeSTE = d.amplitudeSTE(:,thisIndex);
        % set the field
        rois{iROI}.(fieldName).instances{iType}(iInstance,:) = thisAmplitude(sortindex);
        rois{iROI}.(fieldName).instancesSTE{iType}(iInstance,:) = thisAmplitudeSTE(sortindex);
      end
    end
  end
  % keep some other info
  rois{iROI}.(fieldName).r2 = d.r2(sortindex);
  rois{iROI}.(fieldName).instanceVol = stimvol;
  if saveFitStructures
    % copy fit structure
    rois{iROI}.(fieldName).fit = d;
    if isfield(rois{iROI},'canonicalFit');
      % copy canonicalFit into fit so we can display parameters in fitTimecourse
      rois{iROI}.(fieldName).fit.canonicalFit = rois{iROI}.canonicalFit;  %TSL: why duplicate this?
    end
  end
  
  disppercent(inf);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getInstancesMean    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function rois = getInstancesMean(rois,stimvol,args) 

% get some default arguments
startLag=[];blockLen=[];n=[];minResponseLen=[];fieldName=[];verbose=[];groupTrials = [];
getArgs(args,{'startLag=[]','blockLen=[]','n=100','minResponseLen=[]','fieldName=instance','verbose=1','groupTrials=1'});

% some variables
numTypes = length(stimvol);

% choose startLag if not specified
if isempty(startLag)
  % make the start volume begin 3 seconds after stimulation onset
  for iROI = 1:length(rois)
    startLag(iROI) = round(3./rois{iROI}.framePeriod);
  end
elseif length(startLag) == 1
  % set the start volume the same for each roi
  startLag = repmat(startLag,1,length(rois));
end

% choose block length if necessary 
if isempty(blockLen)
  % get the median time difference between stimulation intervals
  allStimInt=diff(sort(cell2mat(stimvol)));
  blockLen = median(allStimInt);
  if sum(allStimInt>blockLen)/length(allStimInt) > 0.1
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('(getInstances:getInstancesMean) More than 10% of all inter-stim intervals are greater than the automatically determined block length, which might be incorrect. You can specify blockLen explicitly in your call to getInstances.');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
  end
end

if isempty(minResponseLen)
  minResponseLen = blockLen/2;
end

% now compute instances using mean
for iROI = 1:length(rois)
  % set up fields
  rois{iROI}.(fieldName).info.type = 'mean';
  rois{iROI}.(fieldName).info.startLag = startLag(iROI);
  rois{iROI}.(fieldName).info.blockLen = blockLen;
  rois{iROI}.(fieldName).info.minResponseLen = minResponseLen;
  rois{iROI}.(fieldName).info.groupTrials = groupTrials;
  rois{iROI}.(fieldName).instanceVol = {};
  rois{iROI}.(fieldName).info.args = args;

  % remove all trials that go over concatenation boundary
  thisStimvol = removeBoundaryTrials(stimvol,startLag(iROI),blockLen,rois{iROI}.concatInfo,rois{iROI}.nFrames,minResponseLen,verbose);
  
  % group the stimvols into groupSize trials, if asked for
  % this could potentially increase SNR by averaging over
  % multiple trials at once
  [thisStimvol numInstances] = groupStimvol(thisStimvol,groupTrials);
  
  % show percent done
  framePeriod = rois{iROI}.framePeriod;
  disppercent(-1/numTypes,sprintf('(getInstances) Creating %0.1f instances for %s with %i voxels by taking mean from vol %i:%i (%0.2fs-%0.2fs)',mean(numInstances),rois{iROI}.name,n,startLag(iROI),startLag(iROI)+blockLen-1,framePeriod*startLag(iROI),framePeriod*(startLag(iROI)+blockLen-1)));
  
  totalVoxCount=min(n,rois{iROI}.n);
  % cycle over number of stimulus types
  for iType = 1:numTypes
    % and number of voxels
    %     for iVox = 1:min(n,rois{iROI}.n)
    % get the index of the current voxel (this depends on the sort order passed in)
    % usually this should be in r2 order from a localizer
    voxIndex = rois{iROI}.sortindex(1:totalVoxCount);
    % cycle over instances
    instanceCount = 1;
    for iInstance = 1:length(thisStimvol{iType})
      % get the start volume(s)
      thisInstanceStimvol = thisStimvol{iType}{iInstance};
      % now for each start volume (note that this is usually 1, but if
      % you gave groupTrials set to something higher, than each instances
      % will be computed from a group of trials, so there may be multiple
      % start trials), get all the volumes starting from startLag through
      % the end of the blockLen
      thisVols = [];
      for iStartVol = 1:length(thisInstanceStimvol)
        % get the start and end volume
        startVol = thisInstanceStimvol(iStartVol)+startLag(iROI);
        endVol = startVol+blockLen-1;
        % now adjust start/end to allow for blocks that are cutoff by the end of the scan
        [startVol endVol] = checkBlockAgainstConcatInfo(startVol,endVol,rois{iROI}.concatInfo,rois{iROI}.nFrames,verbose);
        thisVols = [thisVols startVol:endVol];
      end
      if isempty(startVol) || isempty(endVol) || (((endVol-startVol)+1)<minResponseLen)
        if verbose >= 1
          disp(sprintf('(getInstances) Dropping trial (type %i, repeat %i) with startVol=%i',iType,iInstance,startVol));
        end
      else
        % get the mean
        rois{iROI}.(fieldName).instances{iType}(instanceCount,1:totalVoxCount) = mean(rois{iROI}.tSeries(voxIndex,thisVols(:)),2)';
        % keep the start volume(s)
        if groupTrials == 1
          rois{iROI}.(fieldName).instanceVol{iType}(instanceCount) = thisInstanceStimvol;
        else
          keyboard; %This part is not implemented yet.
          rois{iROI}.(fieldName).instanceVol{iType}{instanceCount} = thisInstanceStimvol;
        end
        % update instance counter
        instanceCount = instanceCount+1;
      end
    end
    %       disppercent((iType-1)/numTypes,iVox/min(n,rois{iROI}.n));
    %     end
  end
  disppercent(inf);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   checkBlockAgainstConcatInfo   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [startVol endVol] = checkBlockAgainstConcatInfo(startVol,endVol,concatInfo,nFrames,verbose)

if startVol < 1
  if verbose >=1
    disp(sprintf('(getInstances) Stimulus volume %i happens before beginning of scan',startVol));
  end
  startVol = [];endVol = [];
  return
end

% if no concatInfo, then just check against number of frames
if isempty(concatInfo)
  if startVol > nFrames
    if verbose > 1
      disp(sprintf('(getInstances) Stimulus volume %i happens after end of scan',startVol));
    end
    startVol = [];endVol = [];
  end
  endVol = min(endVol,nFrames);
  return
end

% check which scan this volume is from
if startVol <= length(concatInfo.whichScan)
  whichScan = concatInfo.whichScan(startVol);
else
  if verbose > 1
    disp(sprintf('(getInstances) Stimulus volume %i happens after end of scan',startVol));
  end
  startVol = [];endVol = [];
  return
end

% now get the endVolume
endVol = min(endVol,concatInfo.runTransition(whichScan,2));

%%%%%%%%%%%%%%%%%%%%%%
%    validateROIs    %
%%%%%%%%%%%%%%%%%%%%%%
function [tf rois] = validateROIs(rois,stimvol,args)

tf = 1;
% get some default arguments
n=[];verbose=[];fieldName=[];
getArgs(args,{'n=100','fieldName=instance','verbose=1'},'suppressUnknownArgMessage=1');

% go through ROIS
for iROI = 1:length(rois)
  % check sort index
  if ~isfield(rois{iROI},'sortindex') || isempty(rois{iROI}.sortindex)
    if ~isequal(n,inf)
      disp(sprintf('(getInstances) *** Asked for %i voxels, but no sortindex was provided for %s. Using order found in ROI ***',n,rois{iROI}.name));
    end
    rois{iROI}.sortindex = 1:min(rois{iROI}.n,n); 
    rois{iROI}.sortindexType = 'default';
  else
    if ~isfield(rois{iROI},'sortindexType') || isempty(rois{iROI}.sortindexType)
      rois{iROI}.sortindexType = 'user';
    end
  end
  
  % get some info about scan this roi was loaded from
  if isempty(rois{iROI}.concatInfo)
    % not a problem, but just let the user know.
    disp(sprintf('(getInstances) *** tSeries from %s has no concatenation information ***',rois{iROI}.name));
  end
  
  % check for loaded time series
  if ~isfield(rois{iROI},'tSeries') || (size(rois{iROI}.tSeries,1) ~= rois{iROI}.n)
    disp(sprintf('(getInstances) *** tSeries for ROI %s has not been loaded. Use loadROITSeries? ***',rois{iROI}.name));
    tf = 0;
  end
  
  % set up fields
  rois{iROI}.(fieldName).instances = [];
  rois{iROI}.(fieldName).info = [];
%   rois{iROI}.(fieldName).stimvol = stimvol;
end

%%%%%%%%%%%%%%%%%%%%%%
%    groupStimvol    %
%%%%%%%%%%%%%%%%%%%%%%
function [groupedStimvol numInstances] = groupStimvol(stimvol,groupSize,randomSort)

% if groupSize is 1 then just make the stimvols into a cell array
% so we can perform the same processing later on
if groupSize == 1
  for iStimType = 1:length(stimvol)
    groupedStimvol{iStimType} = num2cell(stimvol{iStimType});
    numInstances(iStimType) = length(groupedStimvol{iStimType});
  end
  return
end

% randomly sort or not
if ieNotDefined('randomSort'), randomSort = true;end

% go through each stim type
for iStimType = 1:length(stimvol)
  % get the stimvols for this conditon
  thisStimvol = stimvol{iStimType};
  % randomly sort if asked for
  if randomSort
    thisStimvol = thisStimvol(randperm(length(thisStimvol)));
  end
  iStimGroup = 1;
  
  % and until we run out keep grabbing sets of trials
  while length(thisStimvol) > groupSize
    groupedStimvol{iStimType}{iStimGroup} = thisStimvol(1:groupSize);
    thisStimvol = thisStimvol(groupSize+1:end);
    iStimGroup = iStimGroup+1;
  end
  % compute how many instances we have
  numInstances(iStimType) = length(groupedStimvol{iStimType});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    removeBoundaryTrials    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimvolWithoutBoundaryTrials = removeBoundaryTrials(stimvol,startLag,blockLen,concatInfo,nFrames,minResponseLen,verbose)

for iStimType = 1:length(stimvol)
  % initialize array
  stimvolWithoutBoundaryTrials{iStimType} = [];
  for iStimVol = 1:length(stimvol{iStimType})
    % get the start and end volume
    startVol = stimvol{iStimType}(iStimVol)+startLag;
    endVol = startVol+blockLen-1;
    % check to make sure we have not gone over a concatenation boundary
    [startVol endVol] = checkBlockAgainstConcatInfo(startVol,endVol,concatInfo,nFrames,verbose);
    % if we have, give warning
    if isempty(startVol) || isempty(endVol) || (((endVol-startVol)+1)<minResponseLen)
      % give warning
      if verbose > 1,disp(sprintf('(getInstances) Dropping trial (type %i, repeat %i) with startVol=%i',iStimType,iStimVol,stimvol{iStimType}(iStimVol)+startLag)); end
    else
      % this stimvol passed test, so keep it
      stimvolWithoutBoundaryTrials{iStimType}(end+1) = stimvol{iStimType}(iStimVol);
    end
  end
end

