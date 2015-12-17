% leaveOneOutNpermut.m
%
%      usage: leaveOneOutOverPermutations(instances)
%         by: steeve laquitaine
%       date: 2015/12/16
%    purpose: leaveOneOut of N permutations of classes to get chance
%             distribution of accuracies
%
%             instances = leaveOneOutNpermut(instances,'nPerm=20');
%               
%
%
%Additional tags: see leaveOneOut

function retval = leaveOneOutNpermut(instances,varargin)

%permute N times
fieldName=[];
getArgs(varargin,{'nPerm','fieldName=classify'});

%check # number of permutations
%is defined
if ieNotDefined('nPerm')
    fprintf('%s \n','(leaveOneOutNpermut) Please define nPerm')
    keyboard
end

%init
correctThisPerm = nan(1,nPerm);
nRois = length(instances);

%check instance fields
if isfield(instances{1},fieldName) && isfield(instances{1},'name')    
    %each roi
    for roi = 1 : length(instances)
        if ~isfield(instances{roi}.(fieldName),'instances')
            disp(sprintf('(leaveOneOut) No instances found in %s for %s',fieldName,instances{iROI}.name));
        else
            %percent correct by permutation
            inst = instances(roi);
            parfor i = 1 : nPerm
                tmp = leaveOneOut(inst,'permutationUnBal=1');%classify
                correctThisPerm(i) = tmp{1}.classify.leaveOneOut.correct;
            end
            %percent correct by roi and permutation
            retval{roi}.corrects = correctThisPerm;
            
            %correct stats
            retval{roi}.correct = mean(retval{roi}.corrects);
            retval{roi}.correctSTE = std(retval{roi}.corrects);
            
            fprintf('%s %.2f %s \n','(leaveOneOutNpermut) fisher classifier produced',retval{roi}.correctSTE,' on average')
            fprintf('%s %.3f \n','(leaveOneOutNpermut) std over permutation :',retval{roi}.correctSTE)
        end
    end
    return
end

%percent correct by permutation
parfor i = 1 : nPerm
    tmp = leaveOneOut(instances,'permutationUnBal=1');%classify
    correctThisPerm(i) = tmp.correct;
end
%percent correct by roi and permutation
retval.corrects = correctThisPerm;

%correct stats
retval.correctSTE = mean(retval.corrects);
retval.correctSTE = std(retval.corrects);

fprintf('%s %i %s \n','(leaveOneOutNpermut) fisher classifier produced',retval.correctSTE ,' on average')
fprintf('%s %i \n','(leaveOneOutNpermut) std over permutation :',retval.correctSTE)





