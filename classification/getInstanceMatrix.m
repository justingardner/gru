% getInstanceMatrix.m
%
%        $Id:$ 
%      usage: rois = getInstanceMatrix(v,varname,rois,...)
%         by: justin gardner
%       date: 04/10/14
%    purpose: get an instance matrix and paramsVals array. For example,
%             if you run with the variable 'myRandomDir' then it will
%             put in the classify field of the passed in roi cell array
%             a field instanceMatrix which is just a big matrix of instances x voxels for each roi
%             and a paramVals array which is what the value of myRandomDir was for each instance
% 
%             v = newView;
%             v = viewSet(v,'curGroup','Concatenation');
%             rois = loadROITSeries(v,'lMT');
%             rois = getInstanceMatrix(v,'myRandomDir',rois,'taskNum=2');
%
%             Note that this function is just a wrapper around getStimvol and getInstances.
%             You can pass taskNum,phaseNum and segmentNum as optional arguments and those
%             will be passed on to getStimvol.
%             Any other argument will get passed to getInstances, so for example, if you want
%             to set startLag for getInstances, then you would do:
%
%             rois = getInstanceMatrix(v,'myRandomDir',rois,'taskNum=2','startLag=3');
%
%
function rois = getInstanceMatrix(v,varname,rois,varargin)

% check arguments
if nargin < 3
  help getInstanceMatrix
  return
end

% first get arguments for getStimvol
[argNames argValues args] = getArgs(varargin,{'taskNum=[]','phaseNum=[]','segmentNum=[]'});

% get the instances
[stimvol stimNames] = getStimvol(v,varname,'taskNum',taskNum,'phaseNum',phaseNum,'segmentNum',segmentNum);
if isempty(stimvol),return,end

% get the stimulus values out of the stimNames string
for iStim = 1:length(stimNames)
  % get the string from stimNames - this is usally used as a legend title,
  % but we are going to extract values from it.
  stimNameStr = stimNames{iStim};
  thisStimVal = [];
  % go through and check for varname. This is because there may be multiple
  % values listed in the string like: 'val1=3 val2=4'
  while ~isempty(stimNameStr)
    [stimName stimNameStr] = strtok(stimNameStr,'=');
    if strcmp(stimName,varname)
      thisStimVal = mrStr2num(stimNameStr(2:end));
    end
  end
  % see if we got something
  if isempty(thisStimVal);
    disp(sprintf('(getInstanceMatrix) !!! Could not find name of variable in stimNames returned by getStimvol. Setting value of %s to nan !!!',varname));
    stimVal(iStim) = nan;
  else
    stimVal(iStim) = thisStimVal;
  end
end

% now get instances
rois = getInstances(v,rois,stimvol,args{:});

% now sort into one big matrix
for iROI = 1:length(rois)
  % initalize matrices
  rois{iROI}.classify.instanceMatrix = [];rois{iROI}.classify.paramVals = [];
  % set the paramName field
  rois{iROI}.classify.paramName = varname;
  for iInstance = 1:length(rois{iROI}.classify.instances)
    % concatenate the instance matrix
    rois{iROI}.classify.instanceMatrix = [rois{iROI}.classify.instanceMatrix ; rois{iROI}.classify.instances{iInstance}];
    % concatenate the param vals
    rois{iROI}.classify.paramVals = [rois{iROI}.classify.paramVals ; repmat(stimVal(iInstance),size(rois{iROI}.classify.instances{iInstance},1),1)];
  end
end
