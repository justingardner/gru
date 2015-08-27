% mlrComparePRF.m
%
%        $Id:$ 
%      usage: mlrComparePRF()
%         by: justin gardner
%       date: 08/26/15
%    purpose: 
%
function retval = mlrComparePRF()

% check arguments
if ~any(nargin == [0])
  help mlrComparePRF
  return
end

%vals1 = getOverlayValues('ROI4','Bars Concatenation',1,'pRF','r2');
vals1 = getOverlayValues('ROI4','MotionComp',3,'pRFtest','r2');
vals2 = getOverlayValues('ROI4','MotionComp',3,'pRFfull','r2');
plot(vals1,vals2,'ko','MarkerFaceColor','k');
xlabel('pRF r2');
ylabel('pRFFull r2');
keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getOVerlayValues    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function vals = getOverlayValues(roiName,groupName,scanNum,analysisName,overlayName)

% init view
v = newView;
% set group and scan
v = viewSet(v,'curGroup',groupName);
v = viewSet(v,'curScan',scanNum);
% load analysis
v = loadAnalysis(v,sprintf('pRFAnal/%s.mat',stripext(analysisName)));
% load ROI
v = loadROI(v,roiName);
% get coordinates of ROI
roi = loadROITSeries(v,roiName,[],[],'loadType=none');
scanDims = viewGet(v,'scanDims',scanNum);
linearCoords = sub2ind(scanDims,roi.scanCoords(1,:),roi.scanCoords(2,:),roi.scanCoords(3,:));
% get overlay
overlayNum = viewGet(v,'overlayNum',overlayName);
overlay = viewGet(v,'overlayData',scanNum,overlayNum);
% get values
vals = overlay(linearCoords);
% delete view
deleteView(v);
