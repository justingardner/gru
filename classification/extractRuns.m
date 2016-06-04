% extractRuns.m
%
%        $Id: extractRuns.m,v 1.1 2008/08/08 22:47:18 justin Exp $ 
%      usage: [tSeries stimvol concatInfo] = extractRuns(startRun,endRun,tSeries,stimvol,concatInfo)
%         by: justin gardner
%       date: 08/08/08
%    purpose: function that extracts startRun:endRun of the runs in the tSeries
%             stimvol concatInfo.
%
%             see also concatRuns
%
function [tSeries stimvol concatInfo] = extractRuns(startRun,endRun,tSeries,stimvol,concatInfo)

% get the volumes associated with those runs
volNums = concatInfo.runTransition(startRun,1):concatInfo.runTransition(endRun,2);
minVol = min(volNums);maxVol = max(volNums);
nVols = length(volNums);

% get time series
if ~isempty(tSeries)
  tSeries = tSeries(volNums);
end

% extract the concatInfo pieces associated with these runs
concatInfo.n = endRun-startRun+1;

fieldNames = {'whichScan','whichVolume'};
for i = 1:length(fieldNames)
  concatInfo.(fieldNames{i}) = concatInfo.(fieldNames{i})(volNums);
end
concatInfo.whichScan = concatInfo.whichScan-startRun+1;
fieldNames = {'nifti','junkFrames'};
for i = 1:length(fieldNames)
  concatInfo.(fieldNames{i}) = concatInfo.(fieldNames{i})(startRun:endRun);
end
fieldNames = {'filename','path','hipassfilter'};
for i = 1:length(fieldNames)
  for j = 1:concatInfo.n
    tmp{j} = concatInfo.(fieldNames{i}){startRun+(j-1)};
  end
  concatInfo.(fieldNames{i}) = tmp;
end
concatInfo.runTransition = concatInfo.runTransition(startRun:endRun,:)-minVol+1;

% get the correct stimvols
for i = 1:length(stimvol)
  stimvol{i} = stimvol{i}((stimvol{i}>=minVol)&(stimvol{i}<=maxVol))-minVol+1;
end

