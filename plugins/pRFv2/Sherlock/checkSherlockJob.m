% checkSherlockJob.m
% 
%

function isJobDone = checkSherlockJob(split)

global controller

isJobDone = false;

splitName = sprintf('split%i',split.num);

[~,out] = system(sprintf('ssh %s@sherlock.stanford.edu "if [ -f %s/Splits/Analysis/%s_Anal.mat ]; then echo exists; else echo doesNotExist; fi"', controller.suid, controller.sherlockSessionPath, splitName));

if strcmp(deblank(out), 'exists')
    disp('Analysis found. Job has successfully completed running! Pulling files to local machine now...');
    system(sprintf('rsync -q %s@sherlock.stanford.edu:%s/Splits/Analysis/%s_Anal.mat Splits/Analysis/.', controller.suid, controller.sherlockSessionPath, splitName));
    isJobDone = true;
end
