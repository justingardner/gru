% testClassifier.m
%
%        $Id:$ 
%      usage: [pCorrect cInfo] = testClassifier(testInstances,classifier)
%         by: justin gardner
%       date: 11/07/13
%    purpose: Returns pCorrect over all instnaces plus added information in cInfo.
%       e.g.: testClassifier(testInstances,buildClassifier(buildInstances));
%             
%
function [pCorrect cInfo] = testClassifier(testInstances,classifier)

% check arguments
if ~any(nargin == [2])
  help testClassifier
  return
end

% classify
cInfo.class = classifyInstance(classifier,testInstances);

% compute how many we get correct
for iInstance = 1:length(testInstances)
  cInfo.misClassifiedN(iInstance) = sum(cInfo.class{iInstance}~=iInstance);
  cInfo.totalN(iInstance) = length(cInfo.class{iInstance});
end
% compute pCorrect
pCorrect = (sum(cInfo.totalN)-sum(cInfo.misClassifiedN))/sum(cInfo.totalN);

% get the classifier name
if classifier.pca
  cInfo.classifierName =sprintf('%s (PCA basis: %i components accounting for %0.2f%% of variance)',classifier.type,classifier.pcaInfo.nPC,100*sum(classifier.pcaInfo.v(1:classifier.pcaInfo.nPC)));
else
  cInfo.classifierName = sprintf('%s: (demean=%i, zscore=%i)',classifier.type,classifier.demean,classifier.zscore);
end


