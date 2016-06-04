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

