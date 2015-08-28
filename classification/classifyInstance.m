% classifyInstance.m
%
%        $Id: classifyInstance.m,v 1.3 2009/01/27 13:01:08 justin Exp $ 
%      usage: [class rawClassValue rawClassMatrix] = classifyInstance(classifier,instance)
%         by: justin gardner
%       date: 10/01/08
%    purpose: classify a novel instance. Classifier is the 
%             structure returned by buildClassifier.
%             instance is a 1xk array which has k dimensions (e.g. num voxels)
% 
%             class returns what class the instance is from
%             rawClassValue returns the value of the classifier for that best class (for linear
%                classifiers is the projection on to the line perpendicular to the decision boundary)
%             rawClassMatrix returns a matrix of all classification values which is kxk where k is
%                the number of classes
% 
%
function [class classifierOut rawClassificationOutput] = classifyInstance(classifier,instance)

% check arguments
if any(nargin == [0])
  help classifyInstance
  return
end

% see if this is a single instance, or if it is an array of instances
if iscell(instance)
  % preprocess the instances as called for by the classifier variable
  instance = preprocessInstances(instance,classifier);

  % this is a cell array. Do each cell 
  for i = 1:length(instance)
    [class{i} classifierOut{i} rawClassificationOutput{i}] = classifyInstance(classifier,instance{i});
  end
  % create confusion matrix
  numClasses = length(class);
  confusionMatrix = zeros(numClasses,numClasses);
  for iClass = 1:numClasses
    numInstances = length(class{iClass});
    for iInstance = 1:numInstances
      confusionMatrix(iClass,class{iClass}(iInstance)) = confusionMatrix(iClass,class{iClass}(iInstance)) + 1;
    end
  end
  return
elseif isnumeric(instance)
  % if we have an array of instances, classify each one and return
  if size(instance,1) > 1
    for i = 1:size(instance,1)
      [class(i) classifierOut(i) rawClassificationOutput(i,:,:)] = classifyInstance(classifier,instance(i,:));
    end
    return
  end
end

% choose which classifier to use
if strcmp(classifier.type,'mahalanobis')
  [class classifierOut rawClassificationOutput] = classifyInstanceMahalanobis(classifier,instance);
elseif strcmp(classifier.type,'userDefinedLinear')
  [class classifierOut rawClassificationOutput] = classifyInstanceUserDefinedLinear(classifier,instance);
else
  [class classifierOut rawClassificationOutput] = classifyInstanceLinear(classifier,instance);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   classifyInstanceMahalanobis   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [class classifierOut rawClassificationOutput] = classifyInstanceMahalanobis(classifier,instance)

% calculate mahalanobis distance to each cluster
for iClass = 1:classifier.numClasses
  distanceToMean = instance-classifier.mean(iClass,:);
  mahalanobisDistance(iClass) = sqrt(distanceToMean*classifier.inverseCovar{iClass}*distanceToMean');
end
[minDistance class] = min(mahalanobisDistance);

% FIX: is this really what we want to return
classifierOut = minDistance;

rawClassificationOutput = mahalanobisDistance;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   classifyInstanceLinear   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [class classifierOut rawClassificationOutput] = classifyInstanceLinear(classifier,instance)

% init the classifier output matrix. This will hold the output
% of every binary classification.
binaryClassificationOutput = zeros(classifier.numClasses,classifier.numClasses);

% build classifiers on each pair of instances
for iClass = 1:(classifier.numClasses-1)
  for jClass = iClass+1:classifier.numClasses
    binaryClassificationOutput(iClass,jClass) = getsvm(instance,classifier.svm(iClass,jClass));
  end
end

% make into a (anti) symmetric matrix
binaryClassificationOutput = binaryClassificationOutput - binaryClassificationOutput';

% sum the classification across each pair
for iClass = 1:classifier.numClasses
  fullClassificationOutput(iClass) = sum(binaryClassificationOutput(iClass,:));
end

% and pick the maximum as the output class
[maxValue class] = max(fullClassificationOutput);

classifierOut = maxValue;
rawClassificationOutput = flipud(binaryClassificationOutput);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    classifyInstanceUserDefinedLinear    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [class classifierOut rawClassificationOutput] = classifyInstanceUserDefinedLinear(classifier,instance)

% project
rawClassificationOutput = sum(instance .* classifier.w) + classifier.biasPoint;

% return best class
classifierOut = rawClassificationOutput;
class = 1+(sign(classifierOut)+1)/2;
classifierOut = abs(classifierOut);
