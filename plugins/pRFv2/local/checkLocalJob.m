% checkSherlockJob.m
% 
%

function isJobDone = checkLocalJob(split)

global controller

isJobDone = false;

if isfile(fullfile(controller.curPath,'Splits/Analysis',sprintf('%s_split%i_Anal.mat',controller.prfName,split.num)))
    disp('Analysis found. Job has successfully completed running! Shutting down parallel processor...');
    isJobDone = true;
    cancel(split.localProcessor);
end