 function e = vfAnalysis(stimfileNames,varargin)
 
% default return argument
fit = [];

% default to working on current directory
if nargin < 1, stimfileNames = [];end

% parse arguments
getArgs(varargin,{'dispFit=1','combineData=0'});
% get filenames and path
[e.path stimfileNames] = getStimfileNames(stimfileNames);
if isempty(e.path),return,end

% check for valid files as we go through
% nFiles will contain how many vaild files we have
e.nFiles = 0;
e.visualStaircase = {};
e.auditoryStaircase = {};

% cycle through all files
for iFile = 1:length(stimfileNames)
  % display what is happening
  dispHeader(sprintf('(vfAnalysis) Analyzing file (%i/%i): %s        ',iFile,length(stimfileNames),stimfileNames{iFile}));
 
  % load and parse the stimfile
  d = loadStimfile(fullfile(e.path,stimfileNames{iFile}));

  % valid file, so keep its information
  if ~isempty(d)
    if combineData
      % if combining, then a stimfile that contains the same
      % conditions as a previous one will be combined together
      [isNewFile e] = combineStimfiles(e,d);
    else
      % if not combining then treat every file as new
      isNewFile = 1;
    end
    if isNewFile
      % update count
      e.nFiles = e.nFiles + 1;
      % and list of filenames
      e.filenames{e.nFiles} = stimfileNames{iFile};
      % and the data
      e.d{e.nFiles} = d;
      % see if this is a staircase or psychometric function
    end
  end
end

% if no valid files found return
if e.nFiles == 0
  disp(sprintf('(estimationAnalysis) No files found'));
  return
end

k=2


%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getStimfileNames    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stimfilePath stimfileNames] = getStimfileNames(stimfileNames)

%  check if we are examining a single mat file
if isfile(setext(stimfileNames,'mat'))
  % make sure the extension is .mat
  stimfileNames = setext(stimfileNames,'mat');
  % first check to see if it has a path
  if ~isempty(fileparts(stimfileNames))
    stimfilePath = fileparts(stimfileNames);
    stimfileNames = getLastDir(stimfileNames);
  else
    stimfilePath = pwd;
  end
else
  % not a single file, so go look for the path
  if isempty(stimfileNames)
    % get current directory
    stimfilePath = pwd;
    % check if it is a directory
  elseif isdir(stimfileNames)
    % then use that as the path
    stimfilePath = stimfileNames;
    % see if it is a subject ID
  elseif ((length(stimfileNames)>1) && (isequal(lower(stimfileNames(1)),'s'))) || isnumeric(stimfileNames)
    % turn a numeric into a string SID
    if isnumeric(stimfileNames)
      stimfileNames = sprintf('s%03i',stimfileNames);
    end
    % get the path for this subject
    stimfilePath = fullfile('~/data/estimation',stimfileNames);
  else
    disp(sprintf('(estimationAnalysis) Could not find %s',stimfileNames));
    stimfilePath = '';
    stimfileNames = '';
    return
  end
  % get everything in the directory that is a mat file
  matfiles = dir(fullfile(stimfilePath,'*.mat'));
  stimfileNames = {matfiles(:).name};
end

% make sure we are returning a cell array
stimfileNames = cellArray(stimfileNames);


%%%%%%%%%%%%%%%%%%%%%%%%%
%    loadStimfile    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function d = loadStimfile(stimfileName)

% default to empty
d = [];

% add mat extension
stimfileName = setext(stimfileName,'mat');

% see if it exists
if ~isfile(stimfileName)
  disp(sprintf('(estimationAnalysis:loadStimfile) Could not find stimfile: %s',stimfileName));
  return
end

% load the file
s = load(stimfileName);
if ~isfield(s,'myscreen') || ~isfield(s,'task')
  disp(sprintf('(estimationAnalysis:loadStimfile) No myscreen or task in stimfile: %s',stimfileName));
  return
end

% set task name
taskFilename = s.task{1}{1}.taskFilename;

% parse into parameters
d = getTaskParameters(s.myscreen,s.task);
d = d{1};

% print what we found
disp(sprintf('(vfAnalysis:loadStimfile) Found task: %s (%i trials) SID: %s Date: %s',taskFilename,d.nTrials,s.myscreen.SID,s.myscreen.starttime));

% get the variables
d.myscreen = s.myscreen;
d.task = s.task;
d.stimulus = s.stimulus;

%%%%%% graph %%%%%%

d.task{1}{1}.parameter.calculated.cond = nan
%conditions 1-4 are focused, 5-8 distributed. coordinated Q1 -> clockwise
for i = 1:length(d.task{1}{1}.parameter.calculated.att);
    if d.task{1}{1}.parameter.calculated.att(i) == 1;
        d.task{1}{1}.parameter.calculated.cond(i) = d.task{1}{1}.parameter.calculated.quad(i);
    elseif d.task{1}{1}.parameter.calculated.att(i) == 2;
         d.task{1}{1}.parameter.calculated.cond(i) = d.task{1}{1}.parameter.calculated.quad(i)+4;
    end
end
delta = {[] [] [] [] [] [] [] []};
corr = {[] [] [] [] [] [] [] []};
rt = {[] [] [] [] [] [] [] []};
for k = 2:length(d.task{1}{1}.parameter.calculated.att);
    for i = 1:8;
        if d.task{1}{1}.parameter.calculated.cond(k) == i;
            delta{i} = [delta{i} 1-(d.parameter.contrast/d.task{1}{1}.parameter.calculated.stairDiff(k))/d.parameter.contrast];
            corr{i} = [corr{i} d.randVars.correct(k-1)];
            rt{i} = [rt{i} d.reactionTime(k-1)];
        end
    end
end

react = {rt{1} rt{5} rt{2} rt{6} rt{3} rt{7} rt{4} rt{8}};

%change to delta values
for i = 1:8;
    for k = 1:length(d.stimulus.stair(i).s.strength);
        d.stimulus.stair(i).s.strength(k) = 1-(d.parameter.contrast/d.stimulus.stair(i).s.strength(k))/d.parameter.contrast;
    end
end

for q = 1:4
    figure(q); subplot(2,2,1)
    a = 1:length(d.stimulus.stair(2*q-1).s.strength);
    scatter(a,d.stimulus.stair(2*q-1).s.strength); ylim([-.1,1.1]); hold on; line([0,length(a)],[(1-(d.parameter.contrast/d.stimulus.stair(2*q-1).s.threshold)/d.parameter.contrast),(1-(d.parameter.contrast/d.stimulus.stair(2*q-1).s.threshold)/d.parameter.contrast)],'Color','k');
    set(gca, 'YDir','reverse')
    title(sprintf('Quadrant %u: Attended Threshold',q)), xlabel('Trial number'); ylabel('Relative contrast difference');
    subplot(2,2,3), scatter(d.stimulus.stair(2*q-1).s.strength,react{2*q-1});
    title(sprintf('Quadrant %u: Attended Reaction Time',q)), xlabel('Relative contrast difference'); ylabel('Reaction time (ms)'); xlim([-.05,.55]);
    hold on; c = polyfit(d.stimulus.stair(2*q-1).s.strength,react{2*q-1},1); y_est = polyval(c,d.stimulus.stair(2*q-1).s.strength); plot(d.stimulus.stair(2*q-1).s.strength,y_est);
    hold on; subplot(2,2,2)
    b = 1:length(d.stimulus.stair(2*q).s.strength);
    scatter(b,d.stimulus.stair(2*q).s.strength); ylim([-.1,1.1]); hold on; line([0,length(a)],[(1-(d.parameter.contrast/d.stimulus.stair(2*q).s.threshold)/d.parameter.contrast),(1-(d.parameter.contrast/d.stimulus.stair(2*q).s.threshold)/d.parameter.contrast)],'Color','k');
    set(gca, 'YDir','reverse')
    title(sprintf('Quadrant %u: Distributed Threshold',q)), xlabel('Trial number'); ylabel('Relative contrast difference');
    subplot(2,2,4), scatter(d.stimulus.stair(2*q).s.strength,react{2*q});
    title(sprintf('Quadrant %u: Attended Reaction Time',q)), xlabel('Relative contrast difference'); ylabel('Reaction time (ms)'); xlim([-.05,.55]);
    hold on; c = polyfit(d.stimulus.stair(2*q).s.strength,react{2*q},1); y_est = polyval(c,d.stimulus.stair(2*q).s.strength); plot(d.stimulus.stair(2*q).s.strength,y_est)
end
    k=2