% checkMatchingVar.m
%
%        $Id: checkMatchingVar.m,v 1.1 2009/02/12 19:18:39 justin Exp $ 
%      usage: checkMatchingVar(v,varname,scanNums,<'groupName=?'>,<'taskNum=?'>,<phaseNum=?'>,<'segmentNum=?'>)
%         by: justin gardner
%       date: 02/10/09
%    purpose: Check that the variable varname had the same values in all scans listed in scanNums
%             This can be used if, for example, you have run your stimulus program
%             with the same random seed and you want to make sure that the event
%             times were the same
%
function retval = checkMatchingVar(v,varname,scanNums,varargin)

retval = true;

% check arguments
if ~any(nargin == [3 4 5 6 7 8])
  help checkMatchingVar
  return
end

% get arguments
validVars = {'groupName',viewGet(v,'curGroup'),'taskNum=[]','phaseNum=[]','segmentNum=[]'};
getArgs(varargin,validVars);

% get group name
groupName = viewGet(v,'groupName',groupName);

% get info for each scan
for iScan = 1:length(scanNums)
  var(iScan) = getVarInfo(v,varname,scanNums(iScan),groupName,taskNum,phaseNum,segmentNum);
end

% make sure all scans have the same number of runs
for iScan = 2:length(scanNums)
  if var(1).nRuns ~= var(iScan).nRuns
    disp(sprintf('(checkMatchingVar) Scan %s:%i has %i runs which DOES NOT MATCH %i runs in scan %i',groupName,scanNums(iScan),var(iScan).nRuns,var(1).nRuns,scanNums(1)));
    retval = false;
  end
  for iRun = 1:min(var(iScan).nRuns,var(1).nRuns)
    if isequal(var(1).settings{iRun},var(iScan).settings{iRun})
      disp(sprintf('(checkMatchingVar) Scan %s:%i, run %i matches %i trials with %s set to %s',groupName,scanNums(iScan),iRun,length(var(iScan).settings{iRun}),varname,num2str(unique(var(iScan).settings{iRun}))));
    else
      disp(sprintf('(checkMatchingVar) Scan %s:%i, run %i DOES NOT MATCH settings of %s in Scan %i',groupName,scanNums(iScan),iScan,iRun,varname,scanNums(1)));
      retval = false;
    end
  end
end

% check stimvols
for iScan = 2:length(scanNums)
  if isequal(var(1).stimvol,var(iScan).stimvol)
    disp(sprintf('(checkMatchingVar) Scan %s:%i matches stimvols from Scan %i',groupName,scanNums(iScan),scanNums(1)));
  else
    disp(sprintf('(checkMatchingVar) Scan %s:%i DOES NOT MATCH stimvols from Scan %i',groupName,scanNums(iScan),scanNums(1)));
    retval = false;
  end
end

% say that they are the same
if retval
  disp(sprintf('(checkMatchingVar) Variable %s completely matches for %s:%s',varname,groupName,num2str(scanNums)));
end


%%%%%%%%%%%%%%%%%%%%
%%   getVarInfo   %%
%%%%%%%%%%%%%%%%%%%%
function var = getVarInfo(v,varname,scanNum,groupName,taskNum,phaseNum,segmentNum)

disppercent(-inf,sprintf('(checkMatchingVar) Getting settings for %s in scan %i',varname,scanNum));

% set the group and scan
v = viewSet(v,'curGroup',groupName);
v = viewSet(v,'curScan',scanNum);

% keep settings
var.scanNum = viewGet(v,'curScan');
var.groupName = viewGet(v,'groupName',viewGet(v,'curGroup'));

% get stimvols
var.stimvol = getStimvol(v,varname,'taskNum',taskNum,'phaseNum',phaseNum,'segmentNum',segmentNum,'verbose',false);

% get stimfiles
var.stimfiles = viewGet(v,'stimfile');
var.nRuns = length(var.stimfiles);

% now get variable setting
for i = 1:length(var.stimfiles)
  % get task parameters
  e = getTaskParameters(var.stimfiles{i}.myscreen,var.stimfiles{i}.task);
  if ~isempty(taskNum),e = e{taskNum};end
  if ~isempty(phaseNum),e = e(phaseNum);end
  % now get the var settins
  var.settings{i} = getVarFromParameters(varname,e);
end

disppercent(inf);


