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
            %95% confidence interval
            retval{roi}.CI95 = retval{roi}.correct + 1.96*retval{roi}.correctSTE/sqrt(nPerm);

            fprintf('%s %.3f %s \n','(leaveOneOutNpermut) fisher classifier produced',retval{roi}.correct,' on average')
            fprintf('%s %.3f \n','(leaveOneOutNpermut) std over permutation :',retval{roi}.correctSTE)
            fprintf('%s %.3f \n','(leaveOneOutNpermut) 95% CI:',retval{roi}.CI95)
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
retval.correct = mean(retval.corrects);
retval.correctSTE = std(retval.corrects);

%95% confidence interval
retval.CI95 = retval.correct + 1.96*retval.correctSTE/sqrt(nPerm);

%output
fprintf('%s %.3f %s \n','(leaveOneOutNpermut) fisher classifier produced',retval.correct,' on average')
fprintf('%s %.3f  \n','(leaveOneOutNpermut) std over permutation :',retval.correctSTE)
fprintf('%s %.3f \n','(leaveOneOutNpermut) 95% CI:',retval.CI95)





