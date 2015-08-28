% buildClassifier.m
%
%        $Id: buildClassifier.m,v 1.3 2008/10/04 12:00:19 justin Exp $ 
%      usage: classifier = buildClassifier(instances,varargin)
%         by: justin gardner
%       date: 10/01/08
%    purpose: Creates a classifier for later classification.
%             instances can be those returned from getInstances
%             it is a cell array for each category there is 
%             a cell with an nxk matrix in which n is the number
%             of repeats and k is the number of dimensions (e.g. voxels).
%
%             type: Can be any one of 'fisher', 'svm' or 'mahalanobis'
%             kernelfun,C: For svm classification, see getsvm
%             
%
function classifier = buildClassifier(instances,varargin)

classifier = [];
% check arguments
if any(nargin == [0])
  help buildClassifier
  return
end

% parse input arguments
type = [];kernelargs = [];C=[];verbose = [];zscore=[];pca=[];demean=[];
[~,~,preprocessArgs] = getArgs(varargin,{'type=fisher','kernelfun=[]','kernelargs=[]','C=[]','fieldName=classify','verbose=0','projectionLine=[]','biasPoint=[]','w=[]','saveInstances=1'});

% see if we are passed in a cell array of rois. If so, then call buildClassifier
% sequentially on each roi and put the output into the field specified by classField
if isfield(instances{1},fieldName) && isfield(instances{1},'name')
  for iROI = 1:length(instances)
    if ~isfield(instances{iROI}.(fieldName),'instances')
      disp(sprintf('(buildClassifier) No instances found in %s for %s',fieldName,instances{iROI}.name));
    else
      % put the output into the roi with the field specified by classField
      disppercent(-inf,sprintf('(buildClassifier) Building %s classifier for ROI %s',type,instances{iROI}.name));
      instances{iROI}.(fieldName).classifier = buildClassifier(instances{iROI}.(fieldName).instances,'type',type,'kernelfun',kernelfun,'kernelargs',kernelargs,'C',C,'saveInstances',saveInstances);
      disppercent(inf);
    end
  end
  classifier = instances;
  return
end

% preprocess instances
[instances classifier] = preprocessInstances(instances,'args',preprocessArgs);

% build classifier
if any(strcmp(type,{'fisher','svm'}))
  classifier = buildBinaryLinearClassifier(classifier,instances,type,kernelfun,kernelargs,C);
elseif strcmp(type,'mahalanobis')
  classifier = buildMahalanobisClassifier(classifier,instances);
elseif strcmp(lower(type),'userdefinedlinear')
  classifier = buildUserDefinedLinearClassifier(classifier,instances,w,projectionLine,biasPoint);
else
  disp(sprintf('(buildClassifier) Unknown classifier type %s',type));
end

% save the instnaces if called for
if saveInstances
  classifier.instances = instances;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    buildUserDefineLinearClassifier    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function classifier = buildUserDefinedLinearClassifier(classifier,instances,w,projectionLine,biasPoint)

% set classifier type
classifier.type = 'userDefinedLinear';

% number of classes
classifier.numClasses = length(instances);
if (classifier.numClasses ~= 2)
  disp(sprintf('(buildClassifier) User defined linear classifier only implemented for two classes'));
  classifier = [];
  return
end

% number of dimensions
nDims1 = size(instances{1},2);
nDims2 = size(instances{2},2);
% make sure they match
if (nDims1 ~= nDims2)
  disp(sprintf('(buildClassifier) Instances must all have the same number of dimensions (i.e. 2nd dim)'));
  classifier = [];
  return
end
nDims = nDims1;

if ~isempty(projectionLine) && isempty(w)
  % projectionLine should be a 2xn where each row is a point in the n dimensional input
  % space that defines the projection line. So first check that the passed in dimensions
  % are correct
  if size(projectionLine,2) ~= nDims
    disp(sprintf('(buildClassifier) Number of dimensions in projectionLine should be %i. That is, it should be an array of two points in the instance space',nDims));
    classifier = [];
    return
  end

  % keep projection line
  classifier.projectionLine = projectionLine;

  % create the projection line based on the points
  classifier.w = (projectionLine(2,:)-projectionLine(1,:));
  
% user passed in weights
elseif ~isempty(w) && isempty(projectionLine)
  classifier.w = w(:)';
  % check weights are of right length
  if length(classifier.w) ~= nDims
    disp(sprintf('(buildClassifier) Number of dimensions in weight vector w should be %i.',nDims));
    classifier = [];
    return
  end
    
else
  disp(sprintf('(buildClassifier) One of w or projectionLine needs to be passed in for userDefinedLinear classifier. If w is passed in then that is used as the weight vector. If projectionLine is passed in then that is used to define the weights'));
  classifier = [];
  return
end
  
    

% get bias point on projection line if it has not been set
if isempty(biasPoint)
  classifier.biasPoint = 0;
  [class classRaw projection1] = classifyInstance(classifier,instances{1});
  [class classRaw projection2] = classifyInstance(classifier,instances{2});
  % now put the bias point in betweeen the two distributions
  % assuming equal standard deviation
  biasPoint = (mean(projection1)+(mean(projection2)-mean(projection1))/2);
end
classifier.biasPoint = -biasPoint;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   buildMahalanobisClassifier   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function classifier = buildMahalanobisClassifier(classifier,instances)

% set classifier type
classifier.type = 'mahalanobis';

% number of classes
classifier.numClasses = length(instances);

% create mean and covariance matrix for each class
for iClass = 1:classifier.numClasses
  % compute the mean
  classifier.mean(iClass,:) = mean(instances{iClass});
  % compute the covarinace matrix
  classifier.inverseCovar{iClass} = pinv(cov(instances{iClass}));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   buildBinaryLinearClassifier   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function classifier = buildBinaryLinearClassifier(classifier,instances,type,kernelfun,kernelargs,C)

% if fisher then set kernelfun to fisher
if strcmp(type,'fisher'),kernelfun='fisher';end 

% if svm set defaults here
if strcmp(type,'svm')
  % default to linear
  if isempty(kernelfun),kernelfun = 'linear';end
  % default ot C value of 1;
  if isempty(C),C=1;end
end

% set up variable
classifier.type = type;
classifier.params.kernelfun = kernelfun;
classifier.params.kernelargs = kernelargs;
classifier.params.C = C;

% number of classes
classifier.numClasses = length(instances);

% build classifiers on each pair of instances
for iClass = 1:(classifier.numClasses-1)
  for jClass = iClass+1:classifier.numClasses
    classifier.svm(iClass,jClass) = getsvm(instances{iClass},instances{jClass},kernelfun,kernelargs,C);
  end
end

