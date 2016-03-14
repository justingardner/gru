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
%             NaN check added by dan
%               
%
%Additional tags:
%
%           'permutationUnBal=1': to shuffle the classes and keep them unbalanced [steeve 151202]
%             'permutationBal=1': to shuffle the classes and balance them [steeve 151202]
%             'balancByBootSt=1': to balance dataset by boostrapping [steeve 151203]
%             'balancByRemovI=1': to balance datset by removing instances [steeve 151203]

function retval = leaveOneOut(instances,varargin)

%check arguments
if any(nargin == [0])
    help leaveOneOut
    return
end

%get arguments
type = [];kernelfun = [];kernelargs = [];C=[];fieldName=[];hailString=[];permutation=[];
balancByBootSt=[];balancByRemovI=[];fullBoot=[];uselibsvm=[]; permutationUnBal=[]; permutationBal=[];
getArgs(varargin,{'type=fisher','kernelfun=[]','kernelargs=[]','C=[]','fieldName=classify','hailString=[]','permutation=0','balancByBootSt=0','fullBoot=0','balancByRemovI=0','permutationUnBal=0','permutationBal=0','uselibsvm=0'});

% see if we are passed in a cell array of rois. If so, then call leaveOneOut
% sequentially on each roi and put the output into the field specified by classField
if isfield(instances{1},fieldName) && isfield(instances{1},'name')
    for iROI = 1:length(instances)
        if ~isfield(instances{iROI}.(fieldName),'instances')
            disp(sprintf('(leaveOneOut) No instances found in %s for %s',fieldName,instances{iROI}.name));
        else
            %case we want to balance
            %unbalanced dataset
            if balancByBootSt == 1
                fprintf('%s \n','(leaveOneOut)','Bootstrapping to balance dataset')
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
                
                %
                if fullBoot
                    % bootstrap the ENTIRE class that is short instances
                    for ci = 1 : length(class2Boot)
                        tmp = instances{iROI}.(fieldName).instances{class2Boot(ci)};
                        ipos = randi(ni(class2Boot(ci)),nInew,1);
                        instances{iROI}.(fieldName).instances{class2Boot(ci)} = tmp(ipos,:);
                    end
                else
                    % add bootstrapped instances, but use all the instances
                    % that we have already
                    for ci = 1 : length(class2Boot)
                        % get all the instances that we currently have
                        tmp = instances{iROI}.(fieldName).instances{class2Boot(ci)};
                        ipos = [(1:size(tmp,1))' ; randi(ni(class2Boot(ci)),nInew-size(tmp,1),1)];
                        instances{iROI}.(fieldName).instances{class2Boot(ci)} = tmp(ipos,:);
                    end
                    
                end
            end
            
            %case we want to balance
            %unbalanced dataset
            if balancByRemovI == 1
                fprintf('%s \n','(leaveOneOut)','Removing instances to balance dataset')
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
                fprintf('%s \n','(leaveOneOut)','Suffling instance classes')
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
                fprintf('%s \n','(leaveOneOut)','Suffling instance classes')
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
            instances{iROI}.(fieldName).leaveOneOut = leaveOneOut(instances{iROI}.(fieldName).instances,'type',type,'kernelfun',kernelfun,'kernelargs',kernelargs,'C',C,sprintf('hailString=%s%s: ',hailString,instances{iROI}.name));
        end
    end
    retval = instances;
    return
end

%-------------------------- preprocess dataset -------------------
%case we want to balance
%unbalanced dataset
if balancByBootSt == 1
    fprintf('%s \n','(leaveOneOut)','Bootstrapping to balance dataset')
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
    fprintf('%s \n','(leaveOneOut)','Removing instances to balance dataset')
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
    fprintf('%s \n','(leaveOneOut)','Suffling instance classes')
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
    fprintf('%s \n','(leaveOneOut)','Shuffling and balancing instance classes')
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
% parfor_progress(numClasses*numReps);
% disppercent(-1/numClasses,sprintf('(leaveOneOut) %sPerforming leave-one-out cross-validation with classifier %s',hailString,retval.type));
for iClass = 1:numClasses
    % setup variables for parallel loop, it's important to do this now
    % otherwise parfor complains.
    whichClass = zeros(1,numReps(iClass));
    classifierOut = zeros(1,numReps(iClass));
    inst = instances{iClass};
    numRep = numReps(iClass);
    if uselibsvm
        if length(instances)>2
            disp('libsvm=1 flag does not work with non-binary groups');
            return
        end 
    end
    parfor iRep = 1 : numReps(iClass)
        % get the test instance
        testInstance = inst(iRep,:);
        % cerate the training instances, by removing just the testInstance
        trainingInstances = instances;
        trainingInstances{iClass} = instances{iClass}(setdiff(1:numRep,iRep),:);
        if uselibsvm
            % we will ignore the 'type' argument and just use a libsvm svm with
            % default arguments
            % stack trainingInstances and make groups
            trainData = [trainingInstances{1};trainingInstances{2}];
            groupData = [ones(size(trainingInstances{1},1),1);2*ones(size(trainingInstances{2},1),1)];
            csvm = svmtrain(groupData,trainData,'-q');
            [whichClass(iRep)] = svmpredict(1,testInstance,csvm,'-q');
        else
            % now build the classifier
            thisClassifier = buildClassifier(trainingInstances,sprintf('type=%s',type),'kernelfun',kernelfun,'kernelargs',kernelargs,'C',C);
            % and try to classify the instance
            [whichClass(iRep), classifierOut(iRep)] = classifyInstance(thisClassifier,testInstance);
        end
        % update disppercent
%         disppercent(sum(whichClass>0)/numClasses,iRep/numRep);
%         disp(sprintf('rep %i/%i done',iRep,numReps(iClass)));
%         parfor_progress;
    end
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
% parfor_progress(0);

% now make into percent correct
retval.correct = sum(correctByClass)/sum(numReps);

disp(sprintf('(leaveOneOut) %s%s classifier produced %0.2f%% correct and',hailString,retval.type,retval.correct*100));
% disppercent(inf,sprintf('(leaveOneOut) %s%s classifier produced %0.2f%% correct and',hailString,retval.type,retval.correct*100));
retval.correctSTE = sqrt(retval.correct*(1-retval.correct)/sum(numReps));

