% getSortIndex.m
%
%        $Id: getSortIndex.m,v 1.2 2008/12/12 19:43:22 justin Exp $ 
%      usage: rois = getSortIndex(v,rois,r2)
%             rois = getSortIndex(v,rois,r2Left,r2Right)
%         by: justin gardner
%       date: 09/30/08
%    purpose: adds a sort index to the rois based on the r2 values passed in.
%             if two r2 matrices are passed in then this is a split design
%             with different sort values for ROIs that start with 'l' or start with 'r'
%             the first r2 matrix will be used for ROIs that start with 'l' and
%             the second will be used for ROIs that start with 'r'. The r2 matrices
%             are of the dimensions of the scan
function rois = getSortIndex(v,rois,r2Left,r2Right)

% check arguments
if ~any(nargin == [2 3 4])
  help getSortIndex
  return
end

% no sort index, in this case we will just add af field for sorting in
% order of what is in the roi
if nargin == 2
  r2Left = nan;
  r2Right = nan;
  issplit = 0;
elseif nargin == 3
  % this is not a left/right split design
  r2Right = r2Left;
  issplit = 0;
else
  issplit = 1;
end

rois = cellArray(rois);

for iROI = 1:length(rois)
  % get linear scan coordinates
  scanDims = viewGet(v,'scanDims',rois{iROI}.scanNum,rois{iROI}.groupNum);
  rois{iROI}.scanLinearCoords = sub2ind(scanDims,rois{iROI}.scanCoords(1,:),rois{iROI}.scanCoords(2,:),rois{iROI}.scanCoords(3,:));
  % choose which localizer to use for the sortindex
  if issplit
    if rois{iROI}.name(1) == 'l'
      disp(sprintf('(getSortIndex) Using right localizer for %s sortindex',rois{iROI}.name));
      [rois{iROI}.r2 rois{iROI}.sortindex] = sort(r2Left(rois{iROI}.scanLinearCoords),'descend');
    elseif rois{iROI}.name(1) == 'r'
      disp(sprintf('(getSortIndex) Using left localizer for %s sortindex',rois{iROI}.name));
      [rois{iROI}.r2 rois{iROI}.sortindex] = sort(r2Right(rois{iROI}.scanLinearCoords),'descend');
    else
      disp(sprintf('(getSortIndex) %s is neither left nor right *not* using localizer to sort',rois{iROI}.name));
      rois{iROI}.r2 = nan(1:rois{iROI}.n);
      rois{iROI}.sortindex = 1:rois{iROI}.n;
    end
  else
    % this is not a split left right. 
    if isnan(r2Left)
      % this case is for when we just want in ROI order (i.e. we don't care)
      rois{iROI}.r2 = nan(1,rois{iROI}.n);
      rois{iROI}.sortindex = 1:rois{iROI}.n;
    else
      % this case we are acually pulling off values
      [rois{iROI}.r2 rois{iROI}.sortindex] = sort(r2Left(rois{iROI}.scanLinearCoords),'descend');
    end
  end
end
  

