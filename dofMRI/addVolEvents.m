% addVolEvent.m
%
%        $Id:$ 
%      usage: stimfile = addVolEvent(stimfile,addVolTimes)
%         by: justin gardner
%       date: 11/07/15
%    purpose: Pulled from dofmricni adds a vol event at specified time to the stimfile
%
function stimfile = addVolEvents(stimfile,addVolTimes)

% shortcut
e = stimfile.myscreen.events;

% cycle over all vol times to be added
for i = 1:length(addVolTimes)
  volTime = addVolTimes(i);
  % get the event number where we want to insert the event
  addNum = last(find(e.time(1:e.n)<volTime))+1;
  lastEntry = addNum-1;
  nextEntry = min(addNum,e.n);
  if isempty(addNum),addNum=1;lastEntry=1;end
  % get time of this event relative to the ones before it
  % and after it (i.e. 0 if same as last event, 1 if same
  % as next event.
  thisFraction = (volTime-e.time(lastEntry))/(e.time(nextEntry)-e.time(lastEntry));
  % add the event
  e.n = e.n+1;
  e.tracenum = [e.tracenum(1:addNum-1) 1 e.tracenum(addNum:end)];
  e.data = [e.data(1:addNum-1) 1 e.data(addNum:end)];
  % get interpolated ticknum
  ticknum = round(e.ticknum(lastEntry)+thisFraction*(e.ticknum(nextEntry)-e.ticknum(lastEntry)));
  % ticknum may be nan if we are placing an event after the last event in the stimfile
  if isnan(ticknum)
    % get ticknum relative to last entry - find difference in time, and then multipy by framesPerSecond
    ticknum = e.ticknum(lastEntry) + round((volTime - e.time(lastEntry))*stimfile.myscreen.framesPerSecond);
  end
  % save the event
  e.ticknum = [e.ticknum(1:addNum-1) ticknum e.ticknum(addNum:end)];
  e.volnum = [e.volnum(1:addNum-1) e.volnum(lastEntry) e.volnum(addNum:end)+1];
  e.time = [e.time(1:addNum-1) volTime e.time(addNum:end)];
  e.force = [e.force(1:addNum-1) 1 e.force(addNum:end)];
end

% fix the volume number of all the volume events (other events are fixed by code above)
% but volume events get updated afterwords - i.e. they are logged starting at 0. 
e.volnum(e.tracenum==1) = 0:(sum(e.tracenum==1)-1);

% put update events back into structure
stimfile.myscreen.events = e;

% record times of events that have been added
if ~isfield(stimfile.myscreen,'addVolTimes')
  stimfile.myscreen.addVolTimes = [];
end
stimfile.myscreen.addVolTimes(end+1:end+length(addVolTimes)) = addVolTimes;

% update number of volumes
stimfile.myscreen.volnum = stimfile.myscreen.volnum+length(addVolTimes);

