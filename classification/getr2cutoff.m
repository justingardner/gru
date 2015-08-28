% getr2cutoff.m
%
%        $Id: getr2cutoff.m,v 1.4 2008/08/08 23:08:56 justin Exp $ 
%      usage: cutoff = getr2cutoff(tSeries,stimvol,framePeriod,hdrlen,<'concatInfo',concatInfo>,<'n=1000'>,<returnFullDist=0>,<computeOverAllVoxels=0','<pval=[]>')
%         by: justin gardner
%       date: 08/08/08
%    purpose: permutation test to get r2 cutoff. tSeries is a kxn matrix with k timeseries
%             of length n, stimvol is a cell array of stimvols for each condition.
%             framePeriod is the time between each volume acquisition. hdrlen is the
%             length of each estimated hemodynamic response in *seconds*. Randomization
%             of stimvols is done by resampling the distribution of ISIs from the
%             original stimvol -- thus insuring that the randomized stimvols have
%             the same overall timing (i.e. ISI distibution). Cutoff is a structure
%             that has fields for the p=0.05 and p=0.01 cutoffs (p05 and p01). These
%             fields are *arrays* of lenght k (i.e. each cutoff is computed separately
%             for each timeseries. Set computeOverAllVoxels to true if you want
%             a single distribution over all voxels. Set pval to a value and
%             the field cutoff.pval will contain the appropriate cutoff
%
function cutoff = getr2cutoff(tSeries,stimvol,framePeriod,hdrlen,varargin)

% check arguments
if nargin<3
  help getr2cutoff
  return
end

% get arguments
getArgs(varargin,{'concatInfo=[]','n=1000','returnFullDist=0','verbose=1','computeOverAllVoxels=0','pval=[]'});

% get some fields
nhdr = length(stimvol);
hdrlenTR = round(hdrlen/framePeriod);

% get the number of runs
if isempty(concatInfo)
  numRuns = 1;
else
  numRuns = concatInfo.n;
end

% first collect information about stimvols for each run
for runNum = 1:numRuns
  % get the relevant stimvols + concatInfo
  if ~isempty(concatInfo)
    [thisTSeries thisStimvol thisConcatInfo] = extractRuns(runNum,runNum,tSeries,stimvol,concatInfo);
  else
    thisStimvol = stimvol;
    thisConcatInfo = concatInfo;
    thisTSeries = tSeries;
  end
  % get some info about the original stimvols & tSeries
  firstStimvol(runNum) = first(sort(cell2mat(thisStimvol)));
  numStimvols(runNum) = length(cell2mat(thisStimvol));
  tSeriesLen(runNum) = length(thisTSeries);
  % get the distribution of ISIs
  ISIdist{runNum} = diff(sort(cell2mat(thisStimvol)));
  ISIdistLen(runNum) = length(ISIdist{runNum});
end

% get stimvol lengths
for i = 1:nhdr
  stimvolLens(i) = length(stimvol{i});
end

r2 = zeros(size(tSeries,1),n);
if verbose
  disppercent(-inf,'(getr2cutoff) Computing distribution of r2 expected by chance');
end
for repeatNum = 1:n
  % create a randomized set of stimvols
  allRandStimvol = [];thisStartVol = 1;
  for runNum = 1:numRuns
    % grab a new randomized ISI from the distribution
    randISI = ISIdist{runNum}(ceil(rand(1,numStimvols(runNum))*ISIdistLen(runNum)));
    % get randomized stimvols
    thisRandStimvol = [];
    thisRandStimvol(1) = max(1,round(firstStimvol(runNum)-(randISI(1)/2)+(rand*randISI(1))));
    thisRandStimvol = [thisRandStimvol thisRandStimvol+cumsum(randISI(2:end))];
    thisRandStimvol = thisRandStimvol(thisRandStimvol<=tSeriesLen(runNum));
    % concatenate with rest of stimvols
    allRandStimvol = [allRandStimvol thisRandStimvol+thisStartVol-1];
    thisStartVol = thisStartVol+tSeriesLen(runNum);
  end

  % randomize order of allRandStimvol
  allRandStimvol = allRandStimvol(randperm(length(allRandStimvol)));

  % and sort these into stimvols of randomized lengths
  randStimvol = {};currentRandStimvolPos = 1;
  for i = 1:nhdr
    % get a random number of events for this stimulus type
    randStimvolLen = stimvolLens(ceil(rand*length(stimvol)));
    % figure out which volumes to pull out of the random stimulus volumes
    if i ~= nhdr
      endRandStimvolPos = currentRandStimvolPos+randStimvolLen-1;
    else
      % make sure that we include all volumes for the last stimvol
      endRandStimvolPos = length(allRandStimvol);
    end
    % and get the stimvols
    randStimvol{i} = allRandStimvol(min(end,currentRandStimvolPos):min(end,endRandStimvolPos));
    currentRandStimvolPos = currentRandStimvolPos+randStimvolLen;
  end

  % create scm
  d.dim = size(tSeries);
  d.stimvol = randStimvol;
  d.concatInfo = concatInfo;
  d = makescm(d,hdrlen,1);
  % compute deconvolution
  deconv = getr2timecourse(tSeries,nhdr,hdrlen,d.scm,framePeriod,verbose>1);
  r2(:,repeatNum) = deconv.r2;
  if verbose,disppercent(repeatNum/n);end
end
if verbose,disppercent(inf);end

% sort the randomized r2 distribution
r2 = sort(r2')';

% compute over all voxels
if computeOverAllVoxels
  r2 = sort(r2(:))';
  n = size(r2,2);
end

% return full distribution if called for
if returnFullDist
  cutoff.r2 = r2;
end

% get one-sided p-cutoff values
cutoff.p05 = r2(:,round(n*.95));
cutoff.p01 = r2(:,round(n*.99));

% compute pval
if ~isempty(pval)
  cutoff.pval = r2(:,round(n*(1-pval)));
end


