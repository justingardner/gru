% pmMLR.m
%
%      usage: results = pmMLR(homedir, stimfile, datafile, stimradius)
%         by: justin gardner
%       date: 01/25/20
%    purpose: Wrapper for Gari's validation call
%
function results = pmMLR(homedir, stimfile, datafile, stimradius)

% set up return value, default to failure
results.status = -1;
results.errorStr = '';

% check arguments
if ~any(nargin == [4]) || (nargout ~= 1)
  help pmMLR
  return
end

% name for temporary MLR directory to process data
dataDir = fullfile(homedir,'pmMLR');

% now make that as a temporary directory
makeEmptyMLRDir(dataDir,'description=Temporary directory for pmMLR to process pRF data','defaultParams=1');

% get a view from this directory
cd(dataDir);
v = newView;

% save file as a Raw/TSeries
v = saveNewTSeries(v,stimfile);

% check to see if we got a scan
if viewGet(v,'nScans') ~= 1
  results.status = -1;
  results.errorString = sprintf('(pmMLR) Could not get scan %s',stimfile);
  return
end

% setup the pRF analysis
[v params] = pRF(v,[],'justGetParams=1','defaultParams=1');

% set the stimulus file

keyboard


