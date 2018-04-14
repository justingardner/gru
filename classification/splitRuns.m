%This function separates a concatenated time series data into training and test partition,
%used for cross validation in MVPA. Test runs are specified in testRunLabels,
%and non-test data are defined as all remaining data other than test data.
%You can specify arbitrary runs to be extracted, e.g., testRunLabels could be [1 3 8].
%There are two usages of this function. If passed with two variables,
%assume the first one is a data structure (s) used in cross validation. This
%structure should contain ... (see leaveOneRunOut for specification). The
%second variable is testRunLabels. Used this way, the function will call
%itself recursively so that it modifies the roi field of s. i.e. 
%Usage 1: [sTrain sTest]=splitRuns(s,testRunLabels);
%In the second usage, user can pass in timeseries, testRunLabels,stimvol,
%and concatInfo, and it'll return corresponding pieces for training and
%test data, i.e.,
%Usage 2: [trainTSeries trainStimvols trainConcatInfo testTSeries testStimvols testConcatInfo]=splitRuns(tSeries,testRunLabels,stimvol,concatInfo);
%by Taosheng Liu 7/14/2014

% function [trainTSeries trainStimvols trainConcatInfo testTSeries testStimvols testConcatInfo]=splitRuns(dat,testRunLabels,stimvol,concatInfo)
function [arg1 arg2 arg3 arg4 arg5 arg6]=splitRuns(dat,testRunLabels,stimvol,concatInfo)

if nargin==2 && nargout==2 %if user passed two args, assume the first one is an cell array of ROI structs
  sTrain=dat; %initialize train and test data struct, copy all the content from the original
  sTest=dat;
  allROIs=dat.roi;  %extract the ROI data, as we only need to modify the ROI struct
  for iROI=1:length(allROIs)
    trainROI=allROIs{iROI};
    testROI=allROIs{iROI};
    [trainTSeries trainStimvols trainConcatInfo testTSeries testStimvols testConcatInfo]=splitRuns(allROIs{iROI}.tSeries,testRunLabels,dat.stimvol,allROIs{iROI}.concatInfo);
    trainROI.tSeries=trainTSeries;
    trainROI.concatInfo=trainConcatInfo;
    trainROI.nFrames=size(trainTSeries,2);
    
    testROI.tSeries=testTSeries;
    testROI.concatInfo=testConcatInfo;
    testROI.nFrames=size(testTSeries,2);
    
    sTrain.roi{iROI}=trainROI;
    sTest.roi{iROI}=testROI;
    sTrain.stimvol=trainStimvols; %TSL:note this is repeatedly assigned across ROIs, but should be the same value
    sTest.stimvol=testStimvols;
  end
  arg1=sTrain;
  arg2=sTest;
elseif nargin==4 && nargout==6
  tSeries=dat;
  
  allRuns=1:concatInfo.n; %get all of the run labels
  if ~all(ismember(testRunLabels, allRuns))
    error('some test runs are not present in the original data');
    return;
  end
  
  if length(testRunLabels)~=length(unique(testRunLabels))
    warning('Some test runs have the same run labels, ignored.');
  end
  testRunLabels=unique(testRunLabels); %this sort the run labels and return the unique ones
  
  %extract test runs one run at a time, and concatenate them
  tempTSeries={}; tempVols={}; tempConcatInfo={};
  %extraction
  for i=1:length(testRunLabels)
    [tempTSeries{i},tempVols{i},tempConcatInfo{i}]=extractRuns(testRunLabels(i),testRunLabels(i),tSeries,stimvol,concatInfo);
  end
  %concatenation
  [testTSeries,testStimvols,testConcatInfo]=concatRuns(tempTSeries,tempVols,tempConcatInfo);
  
  %now extract training runs one run at a time, and concatenate them
  trainingRunLabels=setdiff(allRuns, testRunLabels);
  tempTSeries={}; tempVols={}; tempConcatInfo={};
  %extraction
  for i=1:length(trainingRunLabels)
    [tempTSeries{i},tempVols{i},tempConcatInfo{i}]=extractRuns(trainingRunLabels(i),trainingRunLabels(i),tSeries,stimvol,concatInfo);
  end
  %concatenation
  [trainTSeries,trainStimvols,trainConcatInfo]=concatRuns(tempTSeries,tempVols,tempConcatInfo);
  
  arg1=trainTSeries;
  arg2=trainStimvols;
  arg3=trainConcatInfo;
  arg4=testTSeries;
  arg5=testStimvols;
  arg6=testConcatInfo;
  
else
  error('Wrong usage of splitRuns, see help');
end
