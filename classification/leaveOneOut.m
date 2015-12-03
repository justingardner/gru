% leaveOneOut.m
%
%        $Id: leaveOneOut.m,v 1.4 2008/10/02 12:00:32 justin Exp $
%      usage: leaveOneOut(instances)
%         by: justin gardner
%       date: 09/11/08
%    purpose: Do leaveOneOut cross-validated classification on the passed in instances
%             Instances is a cell array which has a kxn matrix for each class, where
%             k = number of repetitions and n = number of dimensions (e.g. num voxels).
%
%             can also be called on ROIS that have instances computed
%             rois = getInstances(v,rois,stimvol);
%             rois = leaveOneOut(rois);
%
%             type = one of: 'fisher', 'svm', 'mahalonobis'
%
%             parfor added by dan 2015/12/02
%
function retval = leaveOneOut(instances,varargin)

% check arguments
if any(nargin == [0])
    help leaveOneOut
    return
end

% get arguments
type = [];kernelfun = [];kernelargs = [];C=[];fieldName=[];hailString=[];
getArgs(varargin,{'type=fisher','kernelfun=[]','kernelargs=[]','C=[]','fieldName=classify','hailString=[]'});

% see if we are passed in a cell array of rois. If so, then call leaveOneOut
% sequentially on each roi and put the output into the field specified by classField
if isfield(instances{1},fieldName) && isfield(instances{1},'name')
    for iROI = 1:length(instances)
        if ~isfield(instances{iROI}.(fieldName),'instances')
            disp(sprintf('(leaveOneOut) No instances found in %s for %s',fieldName,instances{iROI}.name));
        else
            
            keyboard
            % put the output into the roi with the field specified by classField
            instances{iROI}.(fieldName).leaveOneOut = leaveOneOut(instances{iROI}.(fieldName).instances,'type',type,'kernelfun',kernelfun,'kernelargs',kernelargs,'C',C,sprintf('hailString=%s%s: ',hailString,instances{iROI}.name));
        end
    end
    retval = instances;
    return
end

% number of classes we are going to classify into
numClasses = length(instances);
% number of repeats we have in each class
for iClass = 1:numClasses
    numReps(iClass) = size(instances{iClass},1);
    numDims(iClass) = size(instances{iClass},2);
    disp(sprintf('(leaveOneOut) Class %i has %i instances with %i dimensions',iClass,numReps(iClass),numDims(iClass)));
end

% check for dimensions being bad
if length(unique(numDims)) ~= 1
    disp(sprintf('(leaveOneOut) All instances must have the same number of dimensions'));
    return
end

% save info on how the classification was done
retval.type = type;
retval.classifierParams.kernelfun = kernelfun;
retval.classifierParams.kernelargs = kernelargs;
retval.classifierParams.C = C;

% check for NaN values, and remove
for iClass = 1:numClasses
    naninst = isnan(instances{iClass});
    vox = any(naninst,1);
    if any(vox)
        disp(sprintf('(leaveOneOut) Some voxels include NaNs: Removing these...'));
        instances{iClass} = instances{iClass}(:,~vox);
        numReps(iClass) = size(instances{iClass},1);
        numDims(iClass) = size(instances{iClass},2);
        disp(sprintf('(leaveOneOut) Class %i now has %i instances with %i dimensions',iClass,numReps(iClass),numDims(iClass)));
    end
end

% cycle through class and repetitions, removing a single
% instance, and building the classifier on the remaining
% instances and testing on that single instance.
retval.whichClass = cell(1,numClasses);
disppercent(-1/numClasses,sprintf('(leaveOneOut) %sPerforming leave-one-out cross-validation with classifier %s',hailString,retval.type));
for iClass = 1:numClasses
    % setup variables for parallel loop, it's important to do this now
    % otherwise parfor complains.
    whichClass = zeros(1,numReps(iClass));
    classifierOut = zeros(1,numReps(iClass));
    inst = instances{iClass};
    numRep = numReps(iClass);
    parfor iRep = 1 : numReps(iClass)
        % get the test instance
        testInstance = inst(iRep,:);
        % cerate the training instances, by removing just the testInstance
        trainingInstances = instances;
        trainingInstances{iClass} = instances{iClass}(setdiff(1:numRep,iRep),:);
        % now build the classifier
        thisClassifier = buildClassifier(trainingInstances,sprintf('type=%s',type),'kernelfun',kernelfun,'kernelargs',kernelargs,'C',C);
        % and try to classify the instance
        [whichClass(iRep), classifierOut(iRep)] = classifyInstance(thisClassifier,testInstance);
        % update disppercent
%         disppercent((iClass-1)/numClasses,iRep/numRep);
    end
    disppercent((iClass-1)/numClasses);
    % copy parallelized outputs back into retval
    retval.whichClass{iClass} = whichClass;
    retval.classifierOut{iClass} = classifierOut;
    % compute how many were correct for this class
    correctByClass(iClass) = sum(retval.whichClass{iClass}==iClass);
    % and compute the confusion matrix row for this class
    for jClass = 1:numClasses
        retval.confusionMatrix(iClass,jClass) = sum(retval.whichClass{iClass}==jClass)/numReps(iClass);
    end
end

% now make into percent correct
retval.correct = sum(correctByClass)/sum(numReps);

disppercent(inf,sprintf('(leaveOneOut) %s%s classifier produced %0.2f%% correct and',hailString,retval.type,retval.correct*100));

retval.correctSTE = sqrt(retval.correct*(1-retval.correct)/sum(numReps));

