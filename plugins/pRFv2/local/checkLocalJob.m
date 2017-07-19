% checkSherlockJob.m
% 
%

function isJobDone = checkLocalJob(split)

global controller

isJobDone = false;

if isfile(fullfile(controller.curPath,'Splits/Analysis',sprintf('split%i_Anal.mat',split.num)))
    disp('Analysis found. Job has successfully completed running! Shutting down parallel loop...');
    isJobDone = true;
    cancel(split.localProcessor);
end