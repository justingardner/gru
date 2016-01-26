% preprocessInstances.m
%
%        $Id:$ 
%      usage: instances = preprocessInstances(instances)
%         by: justin gardner
%       date: 10/25/13
%    purpose: Runs various pre-processing steps on instances. You can:
%             remove the mean across all instances (regardless of type):
%               instances = preprocessInstances(instances,'demean=1');
%             z-score the instances across all instances
%               instances = preprocessInstances(instances,'zecore=1');
%             do pca and keep only components that account for a certain amount of variance
%             across instances. For example to keep components accounting for 60% of variance
%               [instances pSettings] = preprocessInstances(instances,'pca=0.6');
%             note that the variable pSettings will contain the amount of variance accounted 
%             for by each component pSettings.v and the pc components pSettings.pc. The value
%             pSettings.nPC tells you how many PC components were kept (this calls getInstancesPCA)
%             You can also keep a fixed number of pc's. for exmaple to keep 5:
%               [instances pSettings] = preprocessInstances(instances,'pca=5');
%
%             This function is usually called by buildClassifier and classifyInstance. The settings
%             used by buildClassifier are remembered and applied to classifyInstance. So for example:
%             In build the call could look like this:
%               [instances pSettings] = preprocessInstances(instances,'pca=5');
%             and in test:
%               [instances] = preprocessInstances(instances,pSettings);
%             Note that this will apply the PC xform found in the build to the test (i.e. it doesn't
%             compute a new PC basis for the test instances).
%             
%             
function [instances pSettings] = preprocessInstances(instances,varargin)


% second argument can be a pSettings with settings already set
if (length(varargin)>=1) && isstruct(varargin{1})
  pSettings = varargin{1};
  if length(varargin) > 1
    varargin = varargin{2:end};
  else
    varargin = {};
  end
else
  pSettings = [];
end

% now parse args (allow passing of an args cell array with arguments - from buildClassifier)
[~,~,args2] = getArgs(varargin,{'args=[]','suppressUnknwonArgMessage=1'});
if ~isempty(args),varargin={args{:} args2{:}};end
zscore=[];pca=[];demean=[];
getArgs(varargin,{'demean=[]','zscore=0','pca=0','args=[]'});

% if demean is default then
if isempty(demean)
  % set it to true if we are zscoring otherwise false
  if zscore,demean = true;else,demean = false;end
end

% set pSettings variable
if ~isfield(pSettings,'demean'),pSettings.demean = demean;end
if ~isfield(pSettings,'zscore'),pSettings.zscore = zscore;end
if ~isfield(pSettings,'pca'),pSettings.pca = pca;end

% for calls like classifyInstance that send a single instance rather than a cell
if ~iscell(instances)
  instances = cellArray(instances);
  returnAsSingleInstance = true;
else
  returnAsSingleInstance = false;
end

% demean and zscore processing
if pSettings.demean || pSettings.zscore
  % make data matrix which should be one row for each instance, one column for each voxel
  % so d = n x k where n is number of instances and k is number of voxels
  d = [];nClasses = length(instances);
  for iInstance = 1:nClasses
    d(end+1:end+size(instances{iInstance},1),:) = instances{iInstance};
    % remember how many instances there originaly was
    nInstances(iInstance) = size(instances{iInstance},1);
  end
  instances = {};

  % demean
  if pSettings.demean
    % remove mean across instances (i.e. each voxels response will be 0 across all instances)
    meand = mean(d);
    d = bsxfun(@minus,d,meand);
  end

  % z-score
  if pSettings.zscore
    % remove std across instances (i.e. each voxels response will have std=1 across all instances)
    stdd = std(d);
    d = bsxfun(@rdivide,d,stdd);
    disp('(ppInstances) zscore');
  end

  % now sort back into instances
  thisRow = 1;
  for iInstance = 1:nClasses
    thisEndRow = thisRow + nInstances(iInstance) - 1;
    instances{iInstance} = d(thisRow:thisEndRow,:);
    thisRow = thisEndRow+1;
  end
end

% pca basis for pSettings
if ~isempty(pSettings.pca) && ~isequal(pSettings.pca,0)
  if ~isfield(pSettings,'pc')
    % do pca decomposition
    [instances pSettings.pcaInfo] = getInstancesPCA(instances,'pcaComponents',pSettings.pca);
  else
    % PC already computed, so use that
    instances = getInstancesPCA(instances,pSettings.pcaInfo);
  end
end

% some calls expect a single instance to be returned (classifyInstance)
if returnAsSingleInstance
  instances = instances{:};
end

