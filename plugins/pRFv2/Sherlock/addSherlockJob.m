% addSherlockJob.m
%
%   Adds jobs to Sherlock queue, given the split file handle.
%
%      by: akshay jagadeesh
%    date: 07/13/2017
function split = addSherlockJob(split)

global controller

disp('Generating batch scripts');
system(sprintf('sh ~/proj/gru/plugins/pRFv2/Sherlock/generateBatchScripts.sh "%s" "%s" "%s" "%s"', controller.prfName, controller.sherlockSessionPath, controller.suid, num2str(split.num)));

disp('Transferring batch scripts to Sherlock and running');
system(sprintf('rsync -q Splits/Scripts/%s.sbatch %s@sherlock.stanford.edu:%s/Splits/Scripts/.', sprintf('%s_%s',controller.prfName,split.splitName), controller.suid, controller.sherlockSessionPath));
system(sprintf('ssh %s@sherlock.stanford.edu "cd %s/Splits/Scripts/; sbatch %s.sbatch"', controller.suid, controller.sherlockSessionPath, sprintf('%s_%s',controller.prfName,split.splitName)));

