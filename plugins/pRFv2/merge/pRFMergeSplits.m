% pRFMergeSplits.m
%
%     Usage: pRFMergeSplits(varargin)
%      Date: 06/29/2017
%        by: Akshay Jagadeesh
%
function [rawParams, r, overlays] = pRFMergeSplits(analysisName)

disp('***** Merging Split Analyses *****');

%Load master struct
m = load(sprintf('Splits/%s_master.mat', analysisName));
suid = m.suid;
sherlockSessionPath = m.sherlockSessionPath;

%% First, pull the analyses into the local directory.
splits = dir('Splits/');
analyses = dir(sprintf('Splits/Analysis/%s_split*_Anal.mat', analysisName));

% Check to make sure we have all the analyses
if length(analyses) == m.numSplits
  disp(sprintf('All splits found. Merging %d splits.', m.numSplits));
else
  disp(sprintf('WARNING: All splits are not in the local analyses folder. Only found %d out of %d splits.', length(analyses), m.numSplits));
  keyboard
end

% Load the overlays from the master struct
r2 = m.overlays.r2;
polarAngle = m.overlays.polarAngle;
eccentricity = m.overlays.eccentricity;
rfHalfWidth = m.overlays.rfHalfWidth;
pRFAnal = m.pRFAnal;
x = m.x; y = m.x; z = m.z;
scanNum = m.scanNum;
fit = m.fit;
params = m.params;
%v = m.v;
v = getMLRView;
scanDims = viewGet(v, 'scanDims');

% replace 3 with fit.nParams
rawParams = nan(fit.nParams, length(x));
r = nan(length(x), fit.concatInfo.n);

pRFAnal.d{scanNum}.linearCoords2 = [];

for ai = 1:length(analyses)
  l1 = load(sprintf('Splits/Analysis/%s_split%d_Anal.mat', analysisName, ai));
  nVox = length(l1.splits.r2);
  if ai == 1
    nLastVox = nVox;
  end
  startIndex = 1+nLastVox*(ai-1);
  %startIndex = ceil(length(x) / length(analyses))*(ai-1);
  rawParams(:,startIndex:(startIndex+nVox-1)) = l1.splits.params;
  
  %Get scan coords
  x = l1.splits.scanCoords(1,:); y = l1.splits.scanCoords(2,:); z = l1.splits.scanCoords(3,:);

  pRFAnal.d{scanNum}.linearCoords2 = [pRFAnal.d{scanNum}.linearCoords2 sub2ind(scanDims,x,y,z)];
 
  % Set overlays
  for vi = 1:nVox
    iMaster = startIndex-1+vi;
    r2.data{scanNum}(x(vi), y(vi), z(vi)) = l1.splits.r2(vi);
    polarAngle.data{scanNum}(x(vi), y(vi), z(vi)) = l1.splits.polarAngle(vi);
    eccentricity.data{scanNum}(x(vi), y(vi), z(vi)) = l1.splits.eccentricity(vi);
    rfHalfWidth.data{scanNum}(x(vi), y(vi), z(vi)) = l1.splits.rfHalfWidth(vi);
    r(iMaster, :) = l1.splits.r(vi, :);
  end
  nLastVox=nVox;

end

pRFAnal.d{scanNum}.params = rawParams;
pRFAnal.d{scanNum}.r = r;
iScan = find(params.scanNum == scanNum);
thisParams.scanNum = params.scanNum(iScan);
r2.params{scanNum} = thisParams;
polarAngle.params{scanNum} = thisParams;
eccentricity.params{scanNum} = thisParams;
rfHalfWidth.params{scanNum} = thisParams;
overlays = [r2 polarAngle eccentricity rfHalfWidth];

if ~ieNotDefined('fromController')
  return
end

% install analysis
pRFAnal.name = analysisName;
pRFAnal.type = 'pRFAnal';
pRFAnal.groupName = params.groupName;
pRFAnal.function = 'pRF';
pRFAnal.reconcileFunction = 'defaultReconcileParams';
pRFAnal.mergeFunction = 'pRFMergeParams';
pRFAnal.guiFunction = 'pRFGUI';
pRFAnal.params = params;
pRFAnal.overlays = [r2 polarAngle eccentricity rfHalfWidth];
pRFAnal.curOverlay = 1;
pRFAnal.date = datestr(now);
v = viewSet(v,'newAnalysis',pRFAnal);

if isfield(params, 'mergeAnalysis') && params.mergeAnalysis
  saveMethod = mrGetPref('overwritePolicy');
  mrSetPref('overwritePolicy', 'Merge');
end
saveAnalysis(v,pRFAnal.name);
if isfield(params, 'mergeAnalysis') && params.mergeAnalysis
  mrSetPref('overwritePolicy', saveMethod);
end

if ~isempty(viewGet(v, 'fignum'))
  refreshMLRDisplay(viewGet(v, 'viewNum'));
end
