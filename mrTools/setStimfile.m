% setStimfile.m
%
%      usage: setStimfile()
%         by: justin gardner
%       date: 11/08/06
%    purpose: 
%
function retval = setStimfile()

% check arguments
if ~any(nargin == [0 1])
  help setStimfile
  return
end

% check mrSession
if ~isfile('mrSession.mat') 
  disp(sprintf('(setStimfile) No mrSession.mat'));
  return
end
mrSession = load('mrSession.mat');
if ~isfield(mrSession,'groups') || ~isfield(mrSession,'session')
  disp(sprintf('(setStimfile) Unrecogonized mrSession format'));
  return
end

view = newView;

% get the stim files
[stimFileNums stimFileNames] = getStimFiles;
numStimFiles = length(stimFileNums);

% make sure we have some
if isempty(stimFileNums)
  disp(sprintf('(setStimfile) Could not find any stimfiles in Etc'));
  deleteView(view);
  return
end

% get the initial association
for scanNum = 1:viewGet(view,'nScans',1)
  % get the current associated stimfilename
  stimFileName = viewGet(view,'stimFileName',scanNum,1);
  % find it in current list
  whichStimFile = [];
  if length(stimFileName)==1
    stimFileName = sprintf('%s.mat',stripext(getLastDir(stimFileName{1}))); 
    whichStimFile = find(strcmp(stimFileName,stimFileNames));
  end
  % if it is not empty then that is what we want to use
  if ~isempty(whichStimFile)
    scanStimFileNums(scanNum) = stimFileNums(whichStimFile);
  else
    % otherwise, take the stimfile that matches the scanNum
    if (length(stimFileNums) >= scanNum)
      scanStimFileNums(scanNum) = stimFileNums(scanNum);
    else
      scanStimFileNums(scanNum) = 0;
    end
  end
end

% show info
groupInfo('Raw');
listStimfiles(view,scanStimFileNums);

% ask the user if they want to change
disp('To change the association, make an array that contains the');
disp('stimfile number associated with each scan. Enter a 0 for a');
disp('scan that is not associated with any stimfile.');
r = input('Enter association array or (y to accept current, q to quit): ','s');
while (r ~= 'y')
  if r == 'q',deleteView(view);return,end
  scanStimFileNums = str2num(r);
  listStimfiles(view,scanStimFileNums);
  r = input('Enter association array or (y to accept current, q to quit): ','s');
end

disp(sprintf('(setStimfile) Saving stimfile association'));

% now go through and make association
for scanNum = 1:viewGet(view,'nScans',1)
  if scanNum > length(scanStimFileNums)
    % nothing for this scan, so set to empty
    whichStimFile = [];
  else
    % get which file from the array
    whichStimFile = find(scanStimFileNums(scanNum)==stimFileNums);
  end
  % then set the file appropriately
  if ~isempty(whichStimFile)
    view = viewSet(view,'stimFileName',stimFileNames{whichStimFile},scanNum,1);
  else
    view = viewSet(view,'stimFileName','',scanNum,1);
  end
end
saveSession;
deleteView(view);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display the stimfile info and who they are linked to
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function listStimfiles(view, scanStimFileNums)

% now find out what stim data exists (these are files created by)
% the stimulation program that contain stimulus times.
dirList = dir('Etc');
stimFileNum = [];
printString = '';
firstline = 1;
for i = 1:length(dirList);
  name = dirList(i).name;
  % name should be of form yymmdd_stimnn.mat (nn is sequence number)
  if (regexp(name,'\d\d\d\d\d\d_\w\w\w\w\d\d.mat'))
    % get the stimfile number
    stimFileNum(end+1) = str2num(name(12:13));
    % see if that stimFileNum is taken or not
    thisScanNums = [];
    for j = 1:length(scanStimFileNums)
      if scanStimFileNums(j) == stimFileNum(end);
	thisScanNums(end+1) = j;
      end
    end
    % postpend a line break
    if ~firstline
      printString = sprintf('%s\n',printString);
    else
      firstline = 0;
    end
    % load the stimfile
    load(sprintf('Etc/%s',name));
    % make a string of some info myscreen
    stimfileStr = '';
    if isfield(myscreen,'starttime')
      stimfileStr = sprintf('%sTime: %s ',stimfileStr,myscreen.starttime);
    end
    if isfield(myscreen,'endtime')
      stimfileStr = sprintf('%s(End: %s) ',stimfileStr,myscreen.endtime);
    end
    if isfield(myscreen,'volnum')
      stimfileStr = sprintf('%sNumvols: %i ',stimfileStr,myscreen.volnum);
    end
    % printout with brackets if already used
    % set string to show number, time and volume count
    if ~isempty(thisScanNums)
      thisScanNumsStr = '';thisScanVolsStr = '';
      for k = 1:length(thisScanNums)
	thisScanNumsStr = sprintf('%s%i ',thisScanNumsStr,thisScanNums(k));
	thisScanVolsStr = sprintf('%s%i ',thisScanVolsStr,viewGet(view,'nFrames',thisScanNums(k)));
      end
      printString = sprintf('%s%s->%s%sScan: %s',printString,name,thisScanNumsStr,stimfileStr,thisScanVolsStr);
    else
      printString = sprintf('%s%s %s',printString,name,stimfileStr);
    end
  end
end
disp(printString);


%%%%%%%%%%%%%%%%%%%%
% get stimFileNums
%%%%%%%%%%%%%%%%%%%%
function [stimFileNums stimFileNames] = getStimFiles

dirList = dir('Etc');
stimFileNums = [];stimFileNames = {};
for i = 1:length(dirList);
  name = dirList(i).name;
  % name should be of form yymmdd_stimnn.mat (nn is sequence number)
  if (regexp(name,'\d\d\d\d\d\d_\w\w\w\w\d\d.mat'))
    % get the stimfile number
    stimFileNums(end+1) = str2num(name(12:13));
    stimFileNames{end+1} = name;
  end
end

