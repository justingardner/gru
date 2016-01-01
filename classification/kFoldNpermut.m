% kFoldNpermut.m
%
%      usage: kFoldNpermut(instances)
%         by: steeve laquitaine
%       date: 2015/12/16
%    purpose: kFold of N permutations of classes to get chance
%             distribution of accuracies
%
%      usage:
%e.g., 
%
%             instances = kFoldNpermut(instances,'nPerm=20','numFolds=10','type=svm');
%               
%
%
%Additional tags: see leaveOneOut

function retval = kFoldNpermut(instances,varargin)

%get arguments
type = [];kernelfun = [];kernelargs = [];C=[];fieldName=[];
numFolds=[];fieldName=[];numFolds=[];
permutationUnBal=[]; permutationBal=[]; balancByBootSt=[];balancByRemovI=[];
getArgs(varargin,{'type=fisher','kernelfun=[]','kernelargs=[]','C=[]',...
    'fieldName=classify','permutationBal=0',...
    'permutationUnBal=0','balancByBootSt=0','balancByRemovI=0','numFolds=10',...
    'nPerm=1'});

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
                tmp = kFold(inst,['numFolds=' num2str(numFolds)],'type',type,'kernelfun',kernelfun,'kernelargs',kernelargs,'C',C,varargin{:});
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
            fprintf('%s %.3f %s \n','(kFoldNpermut) fisher classifier produced',retval{roi}.correct*100,'% correct on average')
            fprintf('%s %.3f %s \n','(kFoldNpermut) std over permutation :',retval{roi}.correctSTE*100,'%')
            fprintf('%s %.3f %.3f %s \n','(kFoldNpermut) 95% CI: [',retval{roi}.CI95updown*100,']')
        end
    end
    return
end

%percent correct by permutation
parfor i = 1 : nPerm
    %classify
    tmp = kFold(instances,['numFolds=' num2str(numFolds)],'type',type,'kernelfun',kernelfun,'kernelargs',kernelargs,'C',C,varargin{:});    
    correctThisPerm(i) = tmp.correct;
    confusionMatrix{i} = tmp.confusionMatrix;
end
%percent correct by roi and permutation
retval.corrects = correctThisPerm;

%confusion matrices
retval.confusionMatrix = confusionMatrix;

%correct stats
retval.correct = mean(retval.corrects);
retval.correctSTE = std(retval.corrects);

%95% confidence interval
errormargin = 1.96*retval.correctSTE/sqrt(nPerm);
retval.CI95updown = retval.correct + [-errormargin errormargin];

%output
fprintf('%s %.3f %s \n','(kFoldNpermut) fisher classifier produced',retval.correct*100,'% correct on average')
fprintf('%s %.3f %s \n','(kFoldNpermut) std over permutation :',retval.correctSTE*100,'%')
fprintf('%s %.3f %.3f %s \n','(kFoldNpermut) 95% CI: [',retval.CI95updown*100,']')




