% computeCombinedCorrect.m
%
%        $Id: computeCombinedCorrect.m,v 1.7 2009/01/19 19:11:52 justin Exp $ 
%      usage: [combined rois] = computeCombinedCorrect(v,rois,stimvol1,stimvol2,<startLag>,<blockLen>)
%         by: justin gardner
%       date: 09/30/08
%    purpose: Computes the combined correct across hemispheres for l/r split designs. Display
%             with displayCombinedCorrect.
%
%             use loadROITSeries:
%               rois = loadROITSeries(v,{'l_mt','r_mt','l_v3a','r_v3a','l_v7','r_v7'},1,'Concatenation')
%             then add sortIndex with getSortIndex.
% 
%             Get stimvols using getStimvol
function [combined rois] = computeCombinedCorrect(v,rois,stimvolLeft,stimvolRight,varargin)

% check arguments
if any(nargin == [0 1 2 3])
  help computeCombinedCorrect
  return
end

% get arguments
startLag=[];blockLen=[];type=[];kernelfun=[];kernelargs=[];C=[];n=[];groupTrials=[];
getArgs(varargin,{'startLag=[]','blockLen=[]','type=[]','kernelfun=[]','kernelargs=[]','C=[]','n=100','groupTrials=1'});

% get instances
if iscell(rois) && (length(rois) > 0) && iscell(rois{1})
  for stackNum = 1:length(rois)
    rois = getInstances(v,rois,stimvolLeft,'startLag',startLag,'blockLen',blockLen,'fieldName=classifyLeft','n',n,'groupTrials',groupTrials);
    rois = getInstances(v,rois,stimvolRight,'startLag',startLag,'blockLen',blockLen,'fieldName=classifyRight','n',n,'groupTrials',groupTrials);
  end
else
  rois = getInstances(v,rois,stimvolLeft,'startLag',startLag,'blockLen',blockLen,'fieldName=classifyLeft','n',n,'groupTrials',groupTrials);
  rois = getInstances(v,rois,stimvolRight,'startLag',startLag,'blockLen',blockLen,'fieldName=classifyRight','n',n,'groupTrials',groupTrials);
end

% do classification
rois = leaveOneOut(rois,'hailString=Stimulus left ','fieldName=classifyLeft','type',type,'kernelfun',kernelfun,'kernelargs',kernelargs,'C',C);
rois = leaveOneOut(rois,'hailString=Stimulus right ','fieldName=classifyRight','type',type,'kernelfun',kernelfun,'kernelargs',kernelargs,'C',C);
  
% now combine across left and right
iArea = 0;clear combined;
for iROI = 1:2:length(rois)
  % make sure the pair is corect
  if ~strcmp(rois{iROI}.name(2:end),rois{iROI+1}.name(2:end))
    disp(sprintf('(doit) %s and %s are not left/right matched',rois{iROI}.name,rois{iROI+1}.name));
  else
    % get which ROI is left and which one is right
    if rois{iROI}.name(1) == 'l'
      if rois{iROI+1}.name(1) == 'r'
	iLeftROI = iROI;
	iRightROI = iROI+1;
      else
	disp(sprintf('(computeCombinedCorrect) Expecting l_ and r_ rois, but found %s and %s',rois{iROI}.name,rois{iROI+1}.name));
      end
    elseif rois{iROI}.name(2) == 'r'
      if rois{iROI+1}.name(1) == 'l'
	iLeftROI = iROI+1;
	iRightROI = iROI;
      else
	disp(sprintf('(computeCombinedCorrect) Expecting l_ and r_ rois, but found %s and %s',rois{iROI}.name,rois{iROI+1}.name));
      end
    else
      disp(sprintf('(computeCombinedCorrect) Expecting l_ and r_ rois, but found %s and %s',rois{iROI}.name,rois{iROI+1}.name));
    end
      
    iArea = iArea+1;
    % remove the leading l or r from the name
    if any(strcmp(rois{iROI}.name(2),{'_',' '}))
      combined.name{iArea} = rois{iROI}.name(3:end);
    else
      combined.name{iArea} = rois{iROI}.name(2:end);
    end
    combined.contra(iArea) = (rois{iRightROI}.classifyLeft.leaveOneOut.correct + rois{iLeftROI}.classifyRight.leaveOneOut.correct)/2;
    combined.contraSTE(iArea) = sqrt((rois{iRightROI}.classifyLeft.leaveOneOut.correctSTE^2 + rois{iLeftROI}.classifyRight.leaveOneOut.correctSTE^2)/2);
    combined.ipsi(iArea) = (rois{iRightROI}.classifyRight.leaveOneOut.correct + rois{iLeftROI}.classifyLeft.leaveOneOut.correct)/2;
    combined.ipsiSTE(iArea) = sqrt((rois{iRightROI}.classifyRight.leaveOneOut.correctSTE^2 + rois{iLeftROI}.classifyLeft.leaveOneOut.correctSTE^2)/2);
  end
end
for i = 1:length(stimvolLeft)
  % get number of trials
  combined.nLeft(i) = length(stimvolLeft{i});
  % get number of instnaces (this can be different if you have grouped trials to compute instances)
  combined.nLeftInstances(i) = size(rois{1}.classifyLeft.instances{i},1);
end
for i = 1:length(stimvolRight)
  % get number of trials
  combined.nRight(i) = length(stimvolRight{i});
  % get number of instnaces (this can be different if you have grouped trials to compute instances)
  combined.nRightInstances(i) = size(rois{1}.classifyRight.instances{i},1);
end

combined.groupTrials = groupTrials;
combined.blockLen = blockLen;
combined.startLag = startLag;
  
