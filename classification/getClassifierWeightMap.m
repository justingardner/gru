% getClassifierWeightMap.m
%
%        $Id: getClassifierWeightMap.m,v 1.2 2008/10/03 12:00:31 justin Exp $ 
%      usage: [weightMaps weightMapNames] = getClassifierWeightMap(v,rois,stimvol)
%         by: justin gardner
%       date: 10/01/08
%    purpose: returns classifier weight maps that can be installed
%             in MLR with mrDispOverlay. e.g.
%
% mrQuit(0);
% v = newView;
% rois = loadROITSeries(v,{'l_mt','r_mt'});
% stimvol = getStimvol(v,'localdir');
% [weightMaps weightMapNames] = getClassifierWeightMap(v,rois,stimvol);
% v = loadAnat(v,'jg_left_occipital.hdr');
% mrDispOverlay(weightMaps,rois{1}.scanNum,rois{1}.groupNum,v,weightMapNames,'colormapType','normal','cmap',splitcmap,'range=[-1 1]','clip=[-1 1]');
% overlay = refreshMLRDisplay(viewGet(v,'viewNum'));
%
function [weightMaps weightMapNames] = getClassifierWeightMap(v,rois,stimvol,varargin)

% check arguments
if any(nargin == [0])
  help getClassifierWeightMap
  return
end

% get arguments
fieldName = [];startLag = [];blockLen = [];
getArgs(varargin,{'fieldName=classify','startLag=[]','blockLen=[]'});

% init some variables
weightMaps = {};weightMapNames = {};
scanDims = [];

% check scan dims
rois = cellArray(rois);
for iROI = 1:length(rois)
  % get scan dims
  thisScanDims = viewGet(v,'scanDims',rois{iROI}.scanNum,rois{iROI}.groupNum);
  if ~isempty(scanDims) && ~isequal(scanDims,thisScanDims)
    disp(sprintf('(getClassifierWeightMap) Mismatched scan dims for %s',rois{iROI}.name));
    keyboard
  end
  scanDims = thisScanDims;
end

% compute instances
rois = getInstances(v,rois,stimvol,'startLag',startLag,'blockLen',blockLen,'fieldName',fieldName,'n=inf');

% now build the classifiers
rois = buildClassifier(rois,'verbose=1');

% init the weight maps
numClasses = length(stimvol);
numWeightMaps = nchoosek(numClasses,2);
for i = 1:numWeightMaps
  weightMaps{i} = nan(scanDims);
end

% cycle through each ROI
for iROI = 1:length(rois)
  % check to see if we have scan coordinates
  if ~isfield(rois{iROI},'scanCoords')
    rois{iROI}.scanCoords = getROICoordinates(getMLRView,rois{iROI});
  end
  % get the scan linear coords
  if ~isfield(rois{iROI},'scanLinearCoords')
    thisScanDims = viewGet(v,'scanDims',rois{iROI}.scanNum,rois{iROI}.groupNum);
    rois{iROI}.scanLinearCoords = sub2ind(thisScanDims,rois{iROI}.scanCoords(1,:),rois{iROI}.scanCoords(2,:),rois{iROI}.scanCoords(3,:));
  end
  
  weightMapNum = 1;
  % get the sort index that is being used
  if isfield(rois{iROI},'sortindex')
    sortindex = rois{iROI}.sortindex;
  else
    disp(sprintf('(getClassifierWeightMap) No sortindex has been set, ordering in roi order (this should be correct)'));
    sortindex = 1:rois{iROI}.n;
  end
  % and cycle through classes retrieving the weights
  for iClass = 1:numClasses
    for jClass = iClass+1:numClasses
      weightMapNames{weightMapNum} = sprintf('%i vs %i',iClass,jClass);
      % get the weights and enter them into the weight map
      weightMaps{weightMapNum}(rois{iROI}.scanLinearCoords(sortindex)) = rois{iROI}.classify.classifier.svm(iClass,jClass).w;
      weightMapNum = weightMapNum+1;
    end
  end
end



  

