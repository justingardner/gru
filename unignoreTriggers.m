% unignoreTriggers.m
%
%        $Id:$ 
%      usage: stimfile = unignoreTriggers(stimfile,triggersToUnignore)
%         by: justin gardner
%       date: 05/26/15
%    purpose: Unignores triggers that have been ignored. First argument (stimfile) is a stimfile struct (or filename)
%             2nd argument is the triggers you wish to unignore
%
function stimfile = unignoreTriggers(stimfile,triggersToUnignore,varargin)

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

% ok, go to work now. 
disp(sprintf('(unignoreTriggers) Not implemented yet!'));
keyboard

% done save if passed in filename
if saveStimfile
  % save the file back
  stimfile.myscreen.removeTriggers = removeTriggers;
  if ~isfile(stimfileNameOut)
    save(stimfileNameOut,'-struct','stimfileOriginal');
  end
  save(stimfileName,'-struct','stimfile');
end

