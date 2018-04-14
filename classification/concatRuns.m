% concatRuns.m
%
%        $Id:$ 
%      usage: [tSeries stimvol concatInfo] = concatRuns(tSeries,stimvol,concatInfo)
%         by: justin gardner
%       date: 10/12/10
%    purpose: This function concats together tSeries/stimvol/concatInfo from different
%             runs (actually usually different sessions) together so that they can
%             be passed into fitTimecourse.
% 
%             The input arguments are all cell arrays of the relevant variables that
%             you want to concat.
%
%             For example, to concat together two time series (tSeries1 and tSeries2) and
%             their stimvols and concatInfo to make a single tSeries, stimvol and concatInfo:
%
%             [tSeries stimvol concatInfo] = concatRuns({tSeries1 tSeries2},{stimvol1 stimvol2},{concatInfo1 concatInfo2});
%
%             Each tSeries (e.g. tSeries1, tSeries2) should be a kxn row array (where k is
%             the number of voxels and n is the number of volumes).
%
%             See also extractRuns
%
function [tSeriesOut stimvolOut concatInfoOut] = concatRuns(tSeries,stimvol,concatInfo,varargin)

% check arguments
tSeriesOut = [];stimvolOut = [];concatInfoOut = [];
if nargin < 3
  help concatRuns
  return
end

% parse arguments
verbose = [];
getArgs(varargin,{'verbose=1'});

% make sure everything is a cell array
tSeries = cellArray(tSeries);
stimvol = cellArray(stimvol);
concatInfo = cellArray(concatInfo);

% check that the lengths all match
n = length(tSeries);
if (length(stimvol) ~= n)
  disp(sprintf('(concatRuns) Length of stimvol (%i) does not much number of tSeries (%i)',length(stimvol),n));
  return
end
if (length(concatInfo) ~= n)
  disp(sprintf('(concatRuns) Length of concatInfo (%i) does not much number of tSeries (%i)',length(concatInfo),n));
  return
end

% now check length of tSeries against concatInfo to make sure everything matches
for i = 1:n
  thisN = size(tSeries{i},2);
  concatInfoThisN = length(concatInfo{i}.whichVolume);
  if thisN ~= concatInfoThisN
    disp(sprintf('(concatRuns) tSeries %i has %i volumes, but concatInfo says it should have %i volumes',i,thisN,concatInfoThisN));
    return
  end
end

% concat the tSeries together
tSeriesOut = cell2mat(tSeries);
if verbose,disp(sprintf('(concatRuns) New tSeries is %ix%i',size(tSeriesOut,1),size(tSeriesOut,2)));end

% concat the stimvols together
stimvolOut = stimvol{1};
nStimvolTypes = length(stimvolOut);
nVolumes = size(tSeries{1},2);
for i = 2:n
  % make sure that we have a matching number of stimvol types
  if length(stimvol{i}) ~= nStimvolTypes
    disp(sprintf('(concatRuns) Number of stimvol types is %i in run %i which does not match number inf run 1 (%i)',length(stimvol{i}),i,nStimvolTypes));
    return
  end
  % now concat them together, adding the number of volumes from the last run
  for j = 1:length(stimvol{i})
    stimvolOut{j} = [stimvolOut{j} stimvol{i}{j}+nVolumes];
  end
  % update number of volumes
  nVolumes = nVolumes + size(tSeries{i},2);
end


% concat the concatInfo together
concatInfoOut = concatInfo{1};
% check for junk frames
if any(concatInfo{1}.junkFrames~=0)
  disp(sprintf('(concatRuns) *** Found non-zero junk frames in runs %s of session %i. If you got stimvols with getStimvol this should be ok, but in general you do not need to set junkFrames for eventRelated runs ***',mynum2str(find(concatInfo{1}.junkFrames),'sigfigs=0'),1));
end
nVolumes = size(tSeries{1},2);
for i = 2:n
  % concat the whichScan field (adding the number of total scans we have so far
  concatInfoOut.whichScan = [concatInfoOut.whichScan concatInfo{i}.whichScan+concatInfoOut.n];
  % concat the volume number together
  concatInfoOut.whichVolume = [concatInfoOut.whichVolume concatInfo{i}.whichVolume];
  % nifti headers, filenames and path
  concatInfoOut.nifti = [concatInfoOut.nifti concatInfo{i}.nifti];
  concatInfoOut.filename = {concatInfoOut.filename{:} concatInfo{i}.filename{:}};
  concatInfoOut.path = {concatInfoOut.path{:} concatInfo{i}.path{:}};
  % junk frames
  concatInfoOut.junkFrames = [concatInfoOut.junkFrames concatInfo{i}.junkFrames];
  if any(concatInfo{i}.junkFrames~=0)
    disp(sprintf('(concatRuns) *** Found non-zero junk frames in runs %s of session %i. If you got stimvols with getStimvol this should be ok, but in general you do not need to set junkFrames for eventRelated runs ***',mynum2str(find(concatInfo{i}.junkFrames),'sigfigs=0'),i));
  end
  % hipassfilter
  concatInfoOut.hipassfilter = {concatInfoOut.hipassfilter{:} concatInfo{i}.hipassfilter{:}};
  % now update runTransitions
  concatInfoOut.runTransition = [concatInfoOut.runTransition ; concatInfo{i}.runTransition+nVolumes];
  % update number of scans
  concatInfoOut.n = concatInfoOut.n+concatInfo{i}.n;
  % update number of volumes
  nVolumes = nVolumes + size(tSeries{i},2);
end


return



% test code
% test on single time series
testOutput = fitTimecourse(j1.mean_roiTSeries, j1.stimvols, 1.2, 'concatInfo',j1.concatInfo, 'amplitudeType',j.amplitudeType, 'fitType',j.fitType, 'hdrlen',j.hdrlen);
% concat the same time series twice and test
[tSeries stimvol concatInfo] = concatRuns({j1.mean_roiTSeries j1.mean_roiTSeries},{j1.stimvols j1.stimvols},{j1.concatInfo j1.concatInfo});
testOutput = fitTimecourse(tSeries, stimvol, 1.2, 'concatInfo',concatInfo, 'amplitudeType',j.amplitudeType, 'fitType',j.fitType, 'hdrlen',j.hdrlen);
% concat two different time series
[tSeries stimvol concatInfo] = concatRuns({j1.mean_roiTSeries j2.mean_roiTSeries},{j1.stimvols j2.stimvols},{j1.concatInfo j2.concatInfo});
testOutput = fitTimecourse(tSeries, stimvol, 1.2, 'concatInfo',concatInfo, 'amplitudeType',j.amplitudeType, 'fitType',j.fitType, 'hdrlen',j.hdrlen);








