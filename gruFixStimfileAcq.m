% gruFixStimfileAcq.m
%
%        $Id:$ 
%      usage: gruFixStimfileAcq(stimfileName,'removeAcq')
%         by: justin gardner
%       date: 04/09/15
%    purpose: Function to remove triggers from stim file. Based on addTriggers
%             (check there for implementation of addingTriggers - that code
%             just adds one or more triggers equally spaced after the 
%             existing ones)..
%
%             Note that you can run this multiple times - since it
%             keeps a copy of the original file and always runs from
%             that original one. You can also revert
%             to the original, by not setting any remove or add acq:
%             gruFixStimfileAcq(stimfilename);
%
%             e.g.  To remove volume acquisition triggers 1 and 4
%             gruFixStimfileAcq('150403_stim01.mat','removeAcq=[1 4]');
%
%             e.g.  To remove volume acquisition triggers the last
%             two acquisition triggers
%             gruFixStimfileAcq('150403_stim01.mat','remove=Acq[-1 -2]');
%
%             e.g. To run on all scans in Raw on a session
%             v = newView;
%             gruFixStimfileAcq(v,'remove=[-1 -2]');
%
%             To run on a subset of scans from a session
%             gruFixStimfileAcq(v,'scanNums=[3 5]','remove=[-1 -2]');
%
%
function gruFixStimfileAcq(stimfileName,varargin)

% check arguments
if nargin < 1
  help gruFixStimfileAcq
  return
end

% get arguments
getArgs(varargin,{'verbose=1','removeAcq=[]','scanNums=[]'});

% check for view argument
if isview(stimfileName)
  % passed in view
  v = stimfileName;
  % get scanNums to run this for
  if isempty(scanNums)
    scanNums = 1:viewGet(v,'nScans');
  end
  % run for each scans stimfiles one-by-one
  for i = 1:length(scanNums)
    stimfileName = viewGet(v,'stimfilename',scanNums(i));
    if isempty(stimfileName)
      disp(sprintf('(gruFixStimfileAcq) No associated stimfile for scan %i',scanNums(i)));
    else
      args = {stimfileName{1} varargin{:}};
      disp(sprintf('(gruFixStimfileAcq) Fixing %s',stimfileName{1}));
      gruFixStimfileAcq(args{:});
    end
  end
  return
end

% check for already run
stimfile = [];
stimfileNameOriginal = sprintf('%s_original.mat',stripext(stimfileName));
if isfile(stimfileNameOriginal)
  stimfile = load(stimfileNameOriginal);
end

% load the stimfile
if isempty(stimfile)
  if isstr(stimfileName)
    stimfileName = setext(stimfileName,'mat');
    if ~isfile(stimfileName)
      disp(sprintf('(gruFixStimfileAcq) Could not find stimfile %s',stimfileName));
      return
    end
    stimfile = load(stimfileName);
  else
    disp(sprintf('(gruFixStimfileAcq) Stimfile input should be the name of a stimfile'));
    return
  end
end

% check for valid stimfile
if ~isfield(stimfile,'myscreen')
  disp(sprintf('(gruFixStimfileAcq) File is not a stimfile (missing myscreen)'));
  return
end

% keep original
stimfileOriginal = stimfile;

% not changed yet
stimfileChanged = false;

% number of vols and events
nVols = stimfile.myscreen.volnum;
numEvents = stimfile.myscreen.events.n;

if ~isempty(removeAcq)
  % go find the trigger events
  allRemove = [];
  for i = 1:length(removeAcq)
    % get trigNum to remove
    trigNum = removeAcq(i);
    % negative numbers mean from end
    if removeAcq(i) < 1
      trigNum = nVols+removeAcq(i);
    end
    % go find that trigger
    thisRemove = first(sort(intersect(find(stimfile.myscreen.events.volnum==(trigNum-1)),find(stimfile.myscreen.events.tracenum==1))));
    if isempty(thisRemove)
      disp(sprintf('(gruFixStimfileAcq) Could not find trigger event %i',removeAcq(i)));
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
    stimfile.myscreen.events.volnum(allRemove(i):numEvents) = max(stimfile.myscreen.events.volnum(allRemove(i):numEvents)-1,0);
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
  
  % mark changed
  stimfileChanged = true;
  
  % mark what was done to the stimfile
  stimfile.myscreen.removeAcq = removeAcq;
end

if stimfileChanged
  % save the file back
  if ~isfile(stimfileNameOriginal)
    save(stimfileNameOriginal,'-struct','stimfileOriginal');
  end
  save(stimfileName,'-struct','stimfile');
else
  save(stimfileName,'-struct','stimfileOriginal');
  if isfile(stimfileNameOriginal)
    delete(stimfileNameOriginal);
  end
end


