% unignoreTriggers.m
%
%        $Id:$ 
%      usage: stimfile = unignoreTriggers(stimfile,triggersToUnignore)
%         by: justin gardner
%       date: 05/26/15
%    purpose: Unignores triggers that have been ignored. First argument (stimfile) is a stimfile struct (or filename)
%             2nd argument is the triggers you wish to unignore
%
function stimfile = unignoreTriggers(stimfileName,triggersToUnignore,varargin)

stimfile = [];
% check arguments
if nargin < 2
  help unignoreTriggers
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

% get volume trace
volumeTrace = find(strcmp('volume',stimfile.myscreen.traceNames));
if isempty(volumeTrace)
  disp(sprintf('(unignoreTriggers) Could not find volume trace'));
  keyboard
end

% get ignoredVolumes trace
ignoredTrace = find(strcmp('ignoredVolumes',stimfile.myscreen.traceNames));
if isempty(ignoredTrace)
  disp(sprintf('(unignoreTriggers) Could not find ignoredVolumes trace'));
  keyboard
end

% get events
e = stimfile.myscreen.events;

% get ignored Events
ignoredEvents = find(e.tracenum==ignoredTrace);

% go through and unignore
for iIgnore = 1:length(triggersToUnignore)
  % validate that it exists
  if (triggersToUnignore(iIgnore) <= length(ignoredEvents)) && (triggersToUnignore(iIgnore)>0)
    % get which one to ignore
    unignoreEvent = ignoredEvents(triggersToUnignore(iIgnore));
    % found it, so set it to a volume event
    e.tracenum(unignoreEvent) = volumeTrace;
    % set all subsequent events to have one more volume
    if length(e.volnum) > unignoreEvent
      e.volnum(unignoreEvent+1:end) = e.volnum(unignoreEvent+1:end)+1;
    end
    % decrement the ignoredInitialVols counter 
    stimfile.myscreen.ignoredInitialVols = stimfile.myscreen.ignoredInitialVols-1;
    % increment the volume counter
    stimfile.myscreen.volnum = stimfile.myscreen.volnum+1;
  else
    disp(sprintf('(unignoreTriggers) !!! Could not find ignored volume: %i !!!',unignoreEvent));
  end
end

% save back the events
stimfile.myscreen.events = e;

% done save if passed in filename
if saveStimfile
  % save the file back
  stimfile.myscreen.removeTriggers = removeTriggers;
  if ~isfile(stimfileNameOut)
    save(stimfileNameOut,'-struct','stimfileOriginal');
  end
  save(stimfileName,'-struct','stimfile');
end

