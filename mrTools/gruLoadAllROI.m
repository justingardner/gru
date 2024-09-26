
% gruLoadAllROI.m
%
%      usage: gruLoadAllROI(d,myROIsPath)
%         by: steeve laquitaine from justin gardner codes
%       date: 09/12/2015
%    purpose: loads all roi in myROIsPath and displays on base anatomy
%
%       usage: 
%
%               v = gruLoadAllROI(v,'~/data/mlrAnatDB/s0025/mlrROIs/',1)


function v = gruLoadAllROI(v,myROIsPath,varargin)

%get all roi in directory
roisAll = dir([myROIsPath '/*.mat']);
roiAll  = {roisAll.name};
varargin{end+1} = myROIsPath;

%load roi
v = loadROI(v,roiAll,varargin{:});

%display in mrLoadRet
viewSet(v,'showrois','all perimeter');
roiNum = viewGet(v,'roiGroup');         
for j = 1 : numel(roiNum)
    v = viewSet(v,'roiColor','white',roiNum(j));
end
refreshMLRDisplay(v);