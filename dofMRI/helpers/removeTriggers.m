% removeTriggers.m
%
%        $Id:$ 
%      usage: removeTriggers(stimfileName,removeTriggers)
%         by: justin gardner
%       date: 04/09/15
%    purpose: Function to remove triggers from stim file. Based on removeTriggers.
%
%             e.g.  To remove volume acquisition triggers 1 and 4
%             removeTriggers('150403_stim01.mat',[1 4]);
% 
%
%             'verbose=1': Display verbose info
%       
%
function stimfile = removeTriggers(stimfileName,removeTriggers,varargin)

% check arguments
if nargin < 2
  help removeTriggers
  return
end

% get arguments
getArgs(varargin,{'verbose=1'});

% if passed in a stimfile struct then just work on that
% and return it
saveStimfile = true;
if isstruct(stimfileName)
  % set the variable
  stimfile = stimfileName;
  saveStimfile = false;
else
  % check for already run
  stimfile = [];
  stimfileNameOut = sprintf('%sremoveTriggerOriginal.mat',stripext(stimfileName));
  if isfile(stimfileNameOut)
    if askuser(sprintf('(removeTriggers) Found original stimfile %s. Remove triggers from that?',stimfileNameOut));
      stimfile = load(stimfileNameOut);
    end
  end

  % load the stimfile
  if isempty(stimfile)
    if isstr(stimfileName)
      stimfileName = setext(stimfileName,'mat');
      if ~isfile(stimfileName)
	disp(sprintf('(removeTriggers) Could not find stimfile %s',stimfileName));
	return
      end
      stimfile = load(stimfileName);
    else
      disp(sprintf('(removeTriggers) Stimfile input should be the name of a stimfile'));
      return
    end
  end
end

% check for valid stimfile
if ~isfield(stimfile,'myscreen')
  disp(sprintf('(removeTriggers) File is not a stimfile (missing myscreen)'));
  return
end

% keep original
stimfileOriginal = stimfile;

% number of vols and events
nVols = stimfile.myscreen.volnum;
numEvents = stimfile.myscreen.events.n;

% go find the trigger events
allRemove = [];
for i = 1:length(removeTriggers)
  % get trigNum to remove
  trigNum = removeTriggers(i);
  % negative numbers mean from end
  if removeTriggers(i) < 1
    trigNum = nVols+removeTriggers(i);
  end
  % go find that trigger
  thisRemove = first(sort(intersect(find(stimfile.myscreen.events.volnum==(trigNum-1)),find(stimfile.myscreen.events.tracenum==1))));
  if isempty(thisRemove)
    disp(sprintf('(removeTriggers) Could not find trigger event %i',removeTriggers(i)));
  else
    allRemove(end+1) = thisRemove;
  end
end

allRemove = sort(allRemove);
% now remove them, updating all the events afterwords as if the volume never occurred
eventFieldnames = fieldnames(stimfile.myscreen.events);
% first update volume counts by subtracting one after event to be removed
for i = 1:length(allRemove)
  % and update every volnum after the event as if it never happened
  stimfile.myscreen.events.volnum(allRemove(i):numEvents) = stimfile.myscreen.events.volnum(allRemove(i):numEvents)-1;
end
% now remove the events
for i = 1:length(allRemove)
  % decrement n
  stimfile.myscreen.events.n = stimfile.myscreen.events.n-1;
  for iField = 1:length(eventFieldnames);
    % remove the event from every field (except n)
    if ~strcmp(eventFieldnames{iField},'n')
      stimfile.myscreen.events.(eventFieldnames{iField}) = stimfile.myscreen.events.(eventFieldnames{iField})([1:(allRemove(i)-1) (allRemove(i)+1):end]);
    end
  end
  % now decrement allRemove to account for the event that was just removed
  allRemove = allRemove-1;
end

% update volnum
stimfile.myscreen.volnum = stimfile.myscreen.volnum - length(allRemove);

if saveStimfile
  % save the file back
  stimfile.myscreen.removeTriggers = removeTriggers;
  if ~isfile(stimfileNameOut)
    save(stimfileNameOut,'-struct','stimfileOriginal');
  end
  save(stimfileName,'-struct','stimfile');
end



