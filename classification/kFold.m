% kFold.m
%
%        $Id: kFold.m
%      usage: kFold(instances)
%         by: steeve laquitaine
%       date: 12/24/16
%    purpose: Do k fold cross-validated classification on the passed in instances
%             Instances is a cell array which has a kxn matrix for each class, where
%             k = number of repetitions and n = number of dimensions (e.g. num voxels).
%
%             can also be called on ROIS that have instances computed
%             rois = getInstances(v,rois,stimvol);
%             rois = kFold(rois);
%
%             type = one of: 'fisher', 'svm', 'mahalonobis'
%
%
%Additional tags:
%           
%           
%           'permutationUnBal=1': to shuffle the classes and keep them unbalanced [steeve 151202]
%             'permutationBal=1': to shuffle the classes and balance them [steeve 151202]
%             'balancByBootSt=1': to balance dataset by boostrapping [steeve 151203]
%             'balancByRemovI=1': to balance datset by removing instances [steeve 151203]
%
%                  'numFolds=10': number of folds for k-folds (defaults 10) [steeve 151224]
%                  
%
%
%note: use simInstances.m to test

function retval = kFold(instances,varargin)

%check arguments
if any(nargin == [0])
    help kFold
    return
end

%get arguments
type = [];kernelfun = [];kernelargs = [];C=[];fieldName=[];hailString=[];permutation=[];
balancByBootSt=[];balancByRemovI=[];numFolds=[];
getArgs(varargin,{'type=fisher','kernelfun=[]','kernelargs=[]','C=[]',...
    'fieldName=classify','hailString=[]','permutationBal=0',...
    'permutationUnBal=0','balancByBootSt=0','balancByRemovI=0','numFolds=10'});


% see if we are passed in a cell array of rois. If so, then call kFold
% sequentially on each roi and put the output into the field specified by classField
if isfield(instances{1},fieldName) && isfield(instances{1},'name')
    for iROI = 1:length(instances)
        if ~isfield(instances{iROI}.(fieldName),'instances')
            disp(sprintf('(kFold) No instances found in %s for %s',fieldName,instances{iROI}.name));
        else
            %case we want to balance
            %unbalanced dataset
            if balancByBootSt == 1
                fprintf('%s \n','(kFold)','Bootstrapping to balance dataset')
                %get classes
                nClasses = length(instances{1}.classify.instances);
                for ci = 1 : nClasses
                    ni(ci) = size(instances{iROI}.(fieldName).instances{ci},1);
                end
                [nInew,maxci] = max(ni);
                %get classes to bootstrap
                class2Boot = setdiff(1:nClasses,maxci);
                %bootstrap each class to get
                %maxci instances per class
                for ci = 1 : length(class2Boot)
                    tmp = instances{iROI}.(fieldName).instances{class2Boot(ci)};
                    ipos = randi(ni(class2Boot(ci)),nInew,1);
                    instances{iROI}.(fieldName).instances{class2Boot(ci)} = tmp(ipos,:);
                end
            end
            
            %case we want to balance
            %unbalanced dataset
            if balancByRemovI == 1
                fprintf('%s \n','(kFold)','Removing instances to balance dataset')
                %get classes
                nClasses = length(instances{1}.classify.instances);
                for ci = 1 : nClasses
                    ni(ci) = size(instances{iROI}.(fieldName).instances{ci},1);
                end
                nInew = min(ni);
                for ci = 1 : nClasses
                    instances{iROI}.(fieldName).instances{ci} = instances{iROI}.(fieldName).instances{ci}(1:nInew,:);
                end
            end
            
            %case we want to permutate
            %the classes
            if permutationUnBal == 1
                fprintf('%s \n','(kFold)','Suffling instance classes')
                %get classes
                nClasses = length(instances{1}.classify.instances);
                %# of instances per class
                for ci = 1 : nClasses
                    ni(ci) = size(instances{iROI}.(fieldName).instances{ci},1);
                end
                niend = cumsum(ni);
                nist = [1 niend(1:end-1)+1];
                %stack classes instances
                stackedi = cell2mat([instances{iROI}.(fieldName).instances]');
                %shuffle instances position
                shf = randperm(size(stackedi,1));
                stackedi = stackedi(shf,:);
                %feed instances back to a class
                for ci = 1 : nClasses
                    tm{ci} = stackedi(nist(ci):niend(ci),:);
                end
                instances{iROI}.(fieldName).instances = tm;
            end
            
            %case we want to permutate and balance dataset
            %the classes
            if permutationBal == 1
                fprintf('%s \n','(kFold)','Suffling instance classes')
                %get classes
                nClasses = length(instances{1}.classify.instances);
                %stack classes instances
                stackedi = cell2mat([instances{iROI}.(fieldName).instances]');
                %shuffle instances position
                shf = randperm(size(stackedi,1));
                stackedi = stackedi(shf,:);
                %calculate new # of instances per class
                ni = repmat(floor(size(stackedi,1)/nClasses),nClasses,1);
                niend = cumsum(ni);
                nist = [1 niend(1:end-1)+1];
                %feed instances back to a class
                for ci = 1 : nClasses
                    if ci == nClasses
                        tm{ci} = stackedi(nist(ci):end,:);
                    else
                        tm{ci} = stackedi(nist(ci):niend(ci),:);
                    end
                end
                instances{iROI}.(fieldName).instances = tm;
            end
            
            %put the output into the roi with the field specified by classField
            instances{iROI}.(fieldName).kFold = kFold(instances{iROI}.(fieldName).instances,'type',type,'kernelfun',kernelfun,'kernelargs',kernelargs,'C',C,sprintf('hailString=%s%s: ',hailString,instances{iROI}.name));
        end
    end
    retval = instances;
    return
end


%-------------------------- preprocess dataset -------------------
%case we want to balance
%unbalanced dataset
if balancByBootSt == 1
    fprintf('%s \n','(kFold)','Bootstrapping to balance dataset')
    %get classes
    nClasses = length(instances);
    for ci = 1 : nClasses
        ni(ci) = size(instances{ci},1);
    end
    [nInew,maxci] = max(ni);
    %get classes to bootstrap
    class2Boot = setdiff(1:nClasses,maxci);
    %bootstrap each class to get
    %maxci instances per class
    for ci = 1 : length(class2Boot)
        tmp = instances{class2Boot(ci)};
        ipos = randi(ni(class2Boot(ci)),nInew,1);
        instances{class2Boot(ci)} = tmp(ipos,:);
    end
end
%case we want to balance
%unbalanced dataset
if balancByRemovI == 1
    fprintf('%s \n','(kFold)','Removing instances to balance dataset')
    %get classes
    nClasses = length(instances);
    for ci = 1 : nClasses
        ni(ci) = size(instances{ci},1);
    end
    nInew = min(ni);
    for ci = 1 : nClasses
        instances{ci} = instances{ci}(1:nInew,:);
    end
end

%case we want to permute
%the classes
if permutationUnBal == 1
    fprintf('%s \n','(kFold)','Suffling instance classes')
    %get classes
    nClasses = length(instances);
    %# of instances per class
    for ci = 1 : nClasses
        ni(ci) = size(instances{ci},1);
    end
    niend = cumsum(ni);
    nist = [1 niend(1:end-1)+1];
    %stack classes instances
    stackedi = cell2mat(instances');
    %shuffle instances position
    shf = randperm(size(stackedi,1));
    stackedi = stackedi(shf,:);
    %feed instances back to a class
    for ci = 1 : nClasses
        tm{ci} = stackedi(nist(ci):niend(ci),:);
    end
    instances = tm;
end

%case we want to permute and balance classes
if permutationBal == 1
    fprintf('%s \n','(kFold)','Shuffling and balancing instance classes')
    %get classes
    nClasses = length(instances);
    %stack classes instances
    stackedi = cell2mat(instances');
    %shuffle instance position
    shf = randperm(size(stackedi,1));
    stackedi = stackedi(shf,:);
    %calculate new # of instances per class
    ni = repmat(floor(size(stackedi,1)/nClasses),nClasses,1);
    niend = cumsum(ni);
    nist = [1;niend(1:end-1)+1];
    %feed instances back to a class
    for ci = 1 : nClasses
        if ci == nClasses
            tm{ci} = stackedi(nist(ci):end,:);
        else
            tm{ci} = stackedi(nist(ci):niend(ci),:);
        end
    end
    instances = tm;
end
%-------------------------------------------------------------------

% number of classes we are going to classify into
numClasses = length(instances);
% number of repeats we have in each class
for iClass = 1:numClasses
    numReps(iClass) = size(instances{iClass},1);
    numDims(iClass) = size(instances{iClass},2);
    disp(sprintf('(kFold) Class %i has %i instances with %i dimensions',iClass,numReps(iClass),numDims(iClass)));
    numRepByFold(iClass) = floor(numReps(iClass)/numFolds);
    numRepsTrue(iClass) = numRepByFold(iClass) * numFolds;
    instances{iClass} = instances{iClass}(1:numRepsTrue(iClass),:);
    disp(sprintf('(kFold) Class %i has now %i instances with %i dimensions (%i folds of %i instances)',iClass,numRepsTrue(iClass),numDims(iClass),numFolds,numRepByFold(iClass)));
end

% check for dimensions being bad
if length(unique(numDims)) ~= 1
    disp(sprintf('(kFold) All instances must have the same number of dimensions'));
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
        disp(sprintf('(kFold) Some voxels include NaNs: Removing these...'));
        instances{iClass} = instances{iClass}(:,~vox);
        numReps(iClass) = size(instances{iClass},1);
        numDims(iClass) = size(instances{iClass},2);
        disp(sprintf('(kFold) Class %i has %i instances with %i dimensions',iClass,numReps(iClass),numDims(iClass)));
        numRepByFold(iClass) = floor(numReps(iClass)/numFolds);
        numRepsTrue(iClass) = numRepByFold(iClass) * numFolds;
        instances{iClass} = instances{iClass}(1:numRepsTrue(iClass),:);
        disp(sprintf('(kFold) Class %i now has %i instances with %i dimensions (%i folds of %i instances)',iClass,numRepsTrue(iClass),numDims(iClass),numFolds,numRepByFold(iClass)));
    end
end

% cycle through class and repetitions, removing 1 fold of 
% instances, and building the classifier on the remaining
% folds and testing on that removed fold.

for iClass = 1 : numClasses
    % setup variables for parallel loop, it's important to do this now
    % otherwise parfor complains.
    numRep = numRepsTrue(iClass);
    numRepsByFoldThisClass = numRepByFold(iClass);
    inst = instances{iClass};
    %test instances positions within class    
    %last fold might contain more or less data
    testIxSt = 1 : numRepByFold(iClass) : numRep;
    testIxEnd = testIxSt-1;  
    testIxEnd(1) = [];
    if length(testIxEnd(end))<=length(testIxSt)
        testIxEnd(end+1) = numRep;
    end    
    %build on training and classify instances from remaining fold    
    parfor iFold = 1 : numFolds
        %get this fold's instances
        thisFold = testIxSt(iFold):testIxEnd(iFold);
        % get the test fold instances        
        testInstances = inst(thisFold,:);
        % cerate the training instances, by removing just the testInstances
        trainingInstances = instances;
        trainingInstances{iClass} = instances{iClass}(setdiff(1:numRep,thisFold),:);
        % now build the classifier
        thisClassifier = buildClassifier(trainingInstances,sprintf('type=%s',type),'kernelfun',kernelfun,'kernelargs',kernelargs,'C',C);
        % and try to classify these instances
        for iTest = 1 : numRepsByFoldThisClass
            [whichClass{iFold}(iTest), classifierOut{iFold}(iTest)] = classifyInstance(thisClassifier,testInstances(iTest,:));
        end        
        % update disppercent
        % disppercent((iClass-1)/numClasses,iFold/numRep);
    end    
    disppercent((iClass-1)/numClasses);
    % copy parallelized outputs back into retval
    retval.whichClass{iClass} = [whichClass{:}];
    retval.classifierOut{iClass} = [classifierOut{:}];
    % compute how many were correct for this class
    correctByClass(iClass) = sum(retval.whichClass{iClass}==iClass);
    % and compute the confusion matrix row for this class
    for jClass = 1:numClasses
        retval.confusionMatrix(iClass,jClass) = sum(retval.whichClass{iClass}==jClass)/numRepsTrue(iClass);
    end
end

% now make into percent correct
retval.correct = sum(correctByClass)/sum(numRepsTrue);
disppercent(inf,sprintf('(kFold) %s%s classifier produced %0.2f%% correct and',hailString,retval.type,retval.correct*100));
retval.correctSTE = sqrt(retval.correct*(1-retval.correct)/sum(numRepsTrue));
