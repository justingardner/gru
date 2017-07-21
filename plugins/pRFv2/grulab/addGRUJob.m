% addSherlockJob.m
%
%   Adds jobs to Sherlock queue, given the split file handle.
%
%      by: dan birman
%    date: 07/13/2017
function split = addGRUJob(split)

warning('This is an example function: replace GRU with the machine name and save into the appropriate folders');

machine = 'gru';

global controller

disp(sprintf('Running on machine %s',machine));
script = sprintf('bash&; echo pid; matlab -nodesktop; cd ~/proj/gru; startup; cd %s; pRFRunSplits(''%s'',%i);',controller.curPath,controller.prfName,split.num);
[~,pid] = system(sprintf('ssh gru@%s.stanford.edu "%s"',machine,script));

split.pid = pid;