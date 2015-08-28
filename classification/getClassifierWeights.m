% getClassifierWeights.m
%
%        $Id:$ 
%      usage: [weights bias] = getClassifierWeights(c,class1,class2)
%         by: justin gardner
%       date: 07/22/13
%    purpose: gets the weights/bias for the classifier built by builtClassifier
%             used to classify class1 vs class2
%       e.g.: [w bias] = getClassifierWeights(c,1,2);
%             Also the value otherVals is a structure with other classifier specific info
%             [w bias otherVals] = getClassifierWeights(c,1,2);
%
function [w bias otherVals] = getClassifierWeights(c,class1,class2)

% default return values
w = [];bias = [];otherVals = [];

% check arguments
if ~any(nargin == [3])
  help getClassifierWeights
  return
end

% check structure
if ~isstruct(c) || ~isfield(c,'type')
  disp(sprintf('(getClassifierWeights) c must be a structure returned by buildClassifier'));
  return
end

% weights and bais may be tucked away in different fields for different classifiers
if any(strcmp(c.type,{'svm','fisher'}))
  % check class bounds
  if (max(class1,class2) > c.numClasses) || (min(class1,class2) < 1)
    disp(sprintf('(getClassifierWeights) Class %i vs %i out of range [1:%i]',class1,class2,c.numClasses));
    return
  end
  % make sure this is a linear classifier
  if strcmp(c.type,'svm') && ~any(strcmp(c.svm(class1,class2).kernelfun,{'linear','fisher'}))
    disp(sprintf('(getClassifierWeights) Weights for non-linear classifier not returned'));
    return
  end
  if strcmp(c.svm(class1,class2).kernelfun,'fisher')
    % get weights and bias for fisher discriminant
    w = c.svm(class1,class2).w;
    bias = c.svm(class1,class2).bias;
    otherVals = c.svm(class1,class2);
  else
    % get weights and bias for svm classification
    w = c.svm(class1,class2).w;
    bias = c.svm(class1,class2).b;
    otherVals = c.svm(class1,class2);
  end
elseif any(strcmp(c.type,{'userDefinedLinear'}))
  % get weights and bias for user defined linear classification
  w = c.w;
  bias = c.biasPoint;
end

