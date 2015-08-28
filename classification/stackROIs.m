% stackROIs.m
%
%        $Id: stackROIs.m,v 1.3 2009/02/17 00:08:30 justin Exp $ 
%      usage: stackedROIs = stackROIs(rois1,rois2,<varargin>)
%         by: justin gardner
%       date: 10/14/08
%    purpose: 
%
function rois = stackROIs(rois1,rois2,varargin)

% check arguments
if nargin < 2
  help stackROIs
  return
end

type = [];fieldName = [];
getArgs(varargin,{,'type=stackTimecourse','fieldName=classify'});
% stack rois
roiNum = 0;rois = {};
disppercent(-inf,'(stackROIs) Stacking ROIs');
for roiNum1 = 1:length(rois1)
  % check for which roi in rois2 matches the roi name of roi in rois1
  roiMatchNum = [];
  for roiNum2 = 1:length(rois2)
    if strcmp(rois1{roiNum1}.name,rois2{roiNum2}.name)
      roiMatchNum = roiNum2;
      break;
    end
  end
  % no match, give error
  if isempty(roiMatchNum)   
    disp(sprintf('(stackROIs) No matching ROI for %s',rois1{roiNum1}.name));
  else
    % stack the two 
    if isequal(type,'stackTimecourse')
      stackedROI = stackROI(rois1{roiNum1},rois2{roiMatchNum});
    elseif isequal(type,'stackInstances')
      stackedROI = stackROIInstances(rois1{roiNum1},rois2{roiMatchNum},fieldName);
    else
      disp(sprintf('(stackROIs) Unknown stack type: %s',type));
      return
    end
    if ~isempty(stackedROI)
      roiNum = roiNum+1;
      rois{roiNum} = stackedROI;
    end
      
  end
  disppercent(roiNum1/length(rois1),sprintf('(stackROIs) Stacking %s',rois1{roiNum1}.name));
end
disppercent(inf,'(stackROIs) Stacking ROIs');

%%%%%%%%%%%%%%%%%%
%%   stackROI   %%
%%%%%%%%%%%%%%%%%%
function roi = stackROIInstances(roi1,roi2,fieldname)

roi = [];

% check for sortindex
if ~isfield(roi1,'sortindex') || ~isfield(roi2,'sortindex')
  disp(sprintf('(stackROIs) No sortindex for %s (need to run getSortIndex)',roi1.name));
  return
end

% check for sortindex
if ~isfield(roi1,fieldname) || ~isfield(roi2,fieldname)
  disp(sprintf('(stackROIs) No field %s for %s (need to run getInstances)',fieldname,roi1.name));
  return
end

% get the stimvols
stimvol1 = roi1.(fieldname).stimvol;
stimvol2 = roi2.(fieldname).stimvol;

if length(stimvol1) ~= length(stimvol2)
  disp(sprintf('(stackROIs) Number of stimulus types in rois (%s:%i vs %s:%i) does not match',roi1.name,length(stimvol1),roi2.name,length(stimvol2)));
  return
end

roi.name = strcmp('%s %s',roi1.name,roi2.name);

for stimNum = 1:length(stimvol1)
  % get the minimum number of repetitions and minimum number of dimensions (i.e. voxels)
  minInstances = min(size(roi1.(fieldname).instances{stimNum},1),size(roi2.(fieldname).instances{stimNum},1));
  minDimensions = min(size(roi1.(fieldname).instances{stimNum},2),size(roi2.(fieldname).instances{stimNum},2));
  % now concat together 
  roi.(fieldname).instances{stimNum} = [roi1.(fieldname).instances{stimNum}(1:minInstances,1:minDimensions) roi2.(fieldname).instances{stimNum}(1:minInstances,1:minDimensions)];
end  
roi.(fieldname).stimvol


%%%%%%%%%%%%%%%%%%
%%   stackROI   %%
%%%%%%%%%%%%%%%%%%
function roi = stackROI(roi1,roi2)

roi = [];

% check for sortindex
if ~isfield(roi1,'sortindex') || ~isfield(roi2,'sortindex')
  disp(sprintf('(stackROIs) No sortindex for %s (need to run getSortIndex)',rois1{roiNum1}.name));
  return
end

% init by making the stacked roi be the first roi
roi = roi1;

minNVols = min(size(roi1.tSeries,2),size(roi2.tSeries,2));

% now add the second rois data
roi.coords = [roi.coords roi2.coords];
roi.scanCoords = [roi.scanCoords roi2.scanCoords];
roi.scanNum = [roi.scanNum roi2.scanNum];
roi.groupNum = [roi.groupNum roi2.groupNum];
roi.tSeries = [roi.tSeries(:,1:minNVols) ; roi2.tSeries(:,1:minNVols)];
roi.scanLinearCoords = [roi.scanLinearCoords roi2.scanLinearCoords];

% for now, just keep one of the rois scan/group
roi.scanNum = roi1.scanNum;
roi.groupNum = roi1.groupNum;

% make the sourceROI array if it doesn't exist already
% this keep tracks of which stackedROI, each voxel comes from
if ~isfield(roi1,'sourceROI')
  roi1.sourceROI(1:roi1.n) = 1;
end
roi2SourceNum = max(roi1.sourceROI)+1;

% check for stacked ROI in roi 2
if isfield(roi2,'sourceROI')
  disp(sprintf('(stackROIs) ROI %s should not be a stacked ROI',roi2.name));
  roi = [];
end

theseSourceROI = [];roi2VoxNum = 0;
roi.sortindex = [];roi.r2 = [];roi.sourceROI = [];
% now we update the sort index. We do this so that
% the voxels from each roi are interleaved with each other
for roi1VoxNum = 1:roi1.n
  % see which source ROI this voxel came from
  theseSourceROI(end+1) = roi1.sourceROI(roi1VoxNum);
  % check to see if we have seen this source ROI
  % on this cycle
  if length(unique(theseSourceROI)) ~= length(theseSourceROI)
    theseSourceROI = theseSourceROI(end);
    % add a voxel from roi2
    if roi2.n > roi2VoxNum
      roi2VoxNum = roi2VoxNum+1;
      roi.sortindex(end+1) = roi2.sortindex(roi2VoxNum)+roi1.n;
      roi.r2(end+1) = roi2.r2(roi2VoxNum);
      roi.sourceROI(end+1) = roi2SourceNum;
    end
  end
  roi.sortindex(end+1) = roi1.sortindex(roi1VoxNum);
  roi.r2(end+1) = roi1.r2(roi1VoxNum);
  roi.sourceROI(end+1) = roi1.sourceROI(roi1VoxNum);
end
% now add all remaining roi2 voxels (if there are any that haven't
% been interleaved above
for roi2VoxNum = roi2VoxNum+1:roi2.n
  roi.sortindex(end+1) = roi2.sortindex(roi2VoxNum)+roi1.n;
  roi.r2(end+1) = roi2.r2(roi2VoxNum);
  roi.sourceROI(end+1) = roi2SourceNum;
end

% set the number of voxels
roi.n = length(roi.sortindex);

% just make sure everything adds up
if length(roi.sortindex) ~= (roi1.n+roi2.n)
  disp(sprintf('(stackROIs) Number of voxels in combined ROI %i does not match individual rois %i %i',roi.n,roi1.n,roi2.n));
end
