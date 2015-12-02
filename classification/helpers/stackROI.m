
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
