% addLocalJob.m
%
%   Adds jobs to Sherlock queue, given the split file handle.
%
%      by: dan birman
%    date: 07/18/2017
function split = addLocalJob(split)

global controller

disp('Running local split');

split.localProcessor = parfeval(@pRFRunSplits,0,controller.prfName,split.num);