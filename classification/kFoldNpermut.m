% kFoldNpermut.m
%
%      usage: kFoldNpermut(instances)
%         by: steeve laquitaine
%       date: 2015/12/16
%    purpose: kFold of N permutations of classes to get chance
%             distribution of accuracies
%
%             instances = kFoldNpermut(instances,'nPerm=20','numFolds=10');
%               
%
%
%Additional tags: see leaveOneOut

function retval = kFoldNpermut(instances,varargin)

%permute N times
fieldName=[];numFolds=[];
getArgs(varargin,{'nPerm','fieldName=classify','numFolds=10'});

%check # number of permutations
%is defined
if ieNotDefined('nPerm')
    fprintf('%s \n','(kFoldNpermut) Please define nPerm')
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
            disp(sprintf('(kFoldNpermut) No instances found in %s for %s',fieldName,instances{iROI}.name));
        else
            %percent correct by permutation
            inst = instances(roi);                                    
            parfor i = 1 : nPerm
                %classify                
                tmp = kFold(inst,'permutationUnBal=1',['numFolds=' num2str(numFolds)]);
                correctThisPerm(i) = tmp{1}.classify.kFold.correct;
            end                                   
            %percent correct by roi and permutation
            retval{roi}.corrects = correctThisPerm;            
            %correct stats
            retval{roi}.correct = mean(retval{roi}.corrects);
            retval{roi}.correctSTE = std(retval{roi}.corrects);
            %95% confidence interval
            errormargin = 1.96*retval{roi}.correctSTE/sqrt(nPerm);
            retval{roi}.CI95updown = retval{roi}.correct + [-errormargin errormargin];            
            %print
            fprintf('%s %.3f %s \n','(kFoldNpermut) fisher classifier produced',retval{roi}.correct,' on average')
            fprintf('%s %.3f \n','(kFoldNpermut) std over permutation :',retval{roi}.correctSTE)
            fprintf('%s %.3f %.3f %s \n','(kFoldNpermut) 95% CI: [',retval{roi}.CI95updown,']')
        end
    end
    return
end

%percent correct by permutation
parfor i = 1 : nPerm
    %classify
    tmp = kFold(instances,'permutationUnBal=1',['numFolds=' num2str(numFolds)]);
    correctThisPerm(i) = tmp.correct;
end
%percent correct by roi and permutation
retval.corrects = correctThisPerm;

%correct stats
retval.correct = mean(retval.corrects);
retval.correctSTE = std(retval.corrects);

%95% confidence interval
errormargin = 1.96*retval.correctSTE/sqrt(nPerm);
retval.CI95updown = retval.correct + [-errormargin errormargin];

%output
fprintf('%s %.3f %s \n','(kFoldNpermut) fisher classifier produced',retval.correct,' on average')
fprintf('%s %.3f  \n','(kFoldNpermut) std over permutation :',retval.correctSTE)
fprintf('%s %.3f %.3f %s \n','(kFoldNpermut) 95% CI: [',retval.CI95updown,']')





