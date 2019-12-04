% alaisBurrAnalysis.m
%
%      usage: alaisBurrAnalysis()
%         by: justin gardner
%       date: 10/11/19
%    purpose: analyze data from Alais & Burr replication (experiment alaisburr.m)
%
function e = alaisBurrAnalysis(stimfileNames,varargin)

% default return argument
fit = [];

% default to working on current directory
if nargin < 1, stimfileNames = [];end

% parse arguments
getArgs(varargin,{'dispFit=1'});

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
  dispHeader(sprintf('(alaisBurrAnalysis) Analyzing file (%i/%i): %s        ',iFile,length(stimfileNames),stimfileNames{iFile}));
  
  % load and parse the stimfile
  d = loadStimfile(fullfile(e.path,stimfileNames{iFile}));

  % valid file, so keep its information
  if ~isempty(d)
    [isNewFile e] = combineStimfiles(e,d);
    if isNewFile
      % update count
      e.nFiles = e.nFiles + 1;
      % and list of filenames
      e.filenames{e.nFiles} = stimfileNames{iFile};
      % and the data
      e.d{e.nFiles} = d;
      % see if this is a staircase or psychometric function
      if ~e.d{e.nFiles}.isStaircase
	% keep which ones are psychometric functions
	e.isPsycho(e.nFiles) = 1;
      else
	% collect staircases
	if e.d{e.nFiles}.stimulus.visual
	  % add to the visual staircases
	  e.visualStaircase{end+1} = e.d{e.nFiles}.stimulus.stair;
	else
	  % add to the auditory staircases
	  e.auditoryStaircase{end+1} = e.d{e.nFiles}.stimulus.stair;
	end
	% not a psychometric function
	e.isPsycho(e.nFiles) = 0;
      end
    end
  end
end

% now cycle through and fit functions
for iFile = 1:e.nFiles
  if ~e.d{iFile}.isStaircase
    % fit psychometric function
    e.d{iFile} = fitPsychometricFunction(e.d{iFile});
  end
end

% convert to struct
e.visualStaircase = cell2mat(e.visualStaircase);
e.auditoryStaircase = cell2mat(e.auditoryStaircase);

% if no valid files found return
if e.nFiles == 0
  disp(sprintf('(alaisBurrAnalysis) No files found'));
  return
end

% display the fits
if dispFit
  % plot each function
  for iFile = find(e.isPsycho)
    %diplay fit
    dispFits(e.d{iFile});
  end
  % display visual staircase
  if ~isempty(e.visualStaircase)
    mlrSmartfig('alaisBurrAnalysis_visualStaircase','reuse');clf;
    doStaircase('threshold',e.visualStaircase,'type=weibull','dispFig=1','titleStr=Visual staircase','useCurrentFig=1');
  end
  % display auditory staircase
  if ~isempty(e.auditoryStaircase)
    mlrSmartfig('alaisBurrAnalysis_auditoryStaircase','reuse');clf;
    doStaircase('threshold',e.auditoryStaircase,'type=weibull','dispFig=1','titleStr=Auditory staircase','useCurrentFig=1');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    combineStimfiles    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [isNewFile e] = combineStimfiles(e,d)

% default to adding to list
isNewFile = true;

% only combine bimodal conditions
if d.stimulus.bimodal
  % look for matching stimfiles
  for iStimfile = 1:length(e.d)
    % found a match
    if isequal(e.d{iStimfile}.experimentName,d.experimentName)
      e.d{iStimfile} = concatStimfile(e.d{iStimfile},d);
      isNewFile = false;
      return
    end
  end
  % no matches
end

%%%%%%%%%%%%%%%%%%%%%%%%
%    concatStimfile    %
%%%%%%%%%%%%%%%%%%%%%%%%
function d = concatStimfile(d1,d2)

% set the number of trials
d.nTrials = d1.nTrials + d2.nTrials;

% concat fields
d.reactionTime = [d1.reactionTime d2.reactionTime];
d.response = [d1.response d2.response];

% copy these fileds
copyFields = {'parameter','randVars'};
for iField = 1:length(copyFields)
  fieldsToConcat = fieldnames(d1.(copyFields{iField}));
  for iConcatField = 1:length(fieldsToConcat)
    d.(copyFields{iField}).(fieldsToConcat{iConcatField}) = [d1.(copyFields{iField}).(fieldsToConcat{iConcatField}) d2.(copyFields{iField}).(fieldsToConcat{iConcatField})];
  end
end

% grab fields
grabFields = {'stimulusType','visualWidth','displacement','experimentName','nCond','condNames','isStaircase','condWidth','condDisplacement'};
for iField = 1:length(grabFields)
  d.(grabFields{iField}) = d1.(grabFields{iField});
end

% concat the condTrialNums being careful to add the correct number of trials
for iCond = 1:d.nCond
  d.condTrialNums{iCond} = [d1.condTrialNums{iCond} (d2.condTrialNums{iCond}+d1.nTrials)];
end

%%%%%%%%%%%%%%%%
% display fits %
%%%%%%%%%%%%%%%%
function dispFits(d)

% open figure
mlrSmartfig(sprintf('alaisBurrAnalysis_psychometricfits_%s',d.experimentName),'reuse');clf;

nPlots = length(d.visualWidth);
if nPlots > 1
  % plot all together
  nPlots = nPlots+1;
  subplot(1,nPlots,1);
  % display the psychometric functions all together
  dispPsychometricFunction(d,1:d.nCond);
  % now do each width separately
  for iWidth = 1:length(d.visualWidth)
    subplot(1,nPlots,iWidth+1);
    dispPsychometricFunction(d,find(d.condWidth == d.visualWidth(iWidth)));
  end
else
  % display the psychometric functions all together
  dispPsychometricFunction(d,1:d.nCond);
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    display a psychometric function    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispPsychometricFunction(d,whichConds)

% set colors to display in
if (length(whichConds) == 1)
  dataColors = {'k'};
  fitColors = {'r'};
else
  for iCond = 1:length(whichConds)
    dataColors{iCond} = getSmoothColor(iCond,length(whichConds),'cool');
    fitColors{iCond} = getSmoothColor(iCond,length(whichConds),'cool');
  end
end

% title string
titleStr = d.experimentName;

for iCond = 1:length(whichConds)
  % plot fit
  plot(d.fit(whichConds(iCond)).fitX,d.fit(whichConds(iCond)).fitY*100,'-','Color',fitColors{iCond});hold on
  
  % plot in percentile
  myerrorbar(d.cond(whichConds(iCond)).uniquePosDiff,100*d.cond(whichConds(iCond)).correctBinned,'yError',100*d.cond(whichConds(iCond)).correctBinnedError,'Symbol','o','MarkerFaceColor',dataColors{iCond});
  xlabel('Position difference (deg)');
  yaxis(0,100);
  ylabel('Percent rightwards choices (100%%)');

  % append fit parameters to title
  titleStr = sprintf('%s\nMean: %0.2f Std: %0.2f lambda: %0.2f',titleStr,d.fit(whichConds(iCond)).mean,d.fit(whichConds(iCond)).std,d.fit(whichConds(iCond)).lambda);
  
end

% display title
title(titleStr);

% display a legend for more than one
if length(whichConds) > 1
  % display legend
  hLegend = mylegend({d.condNames{whichConds}},dataColors);
  set(hLegend,'Location','northwest');
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    fitPsychometricFunction    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = fitPsychometricFunction(d)

% fit for each set of data in the d structure
for iCond = 1:d.nCond
  % get these trial nums
  trialNums = d.condTrialNums{iCond};

  % get the differences
  d.cond(iCond).posDiff = d.parameter.posDiff(trialNums);
  d.cond(iCond).uniquePosDiff = unique(d.cond(iCond).posDiff);

  % compute whether answer are correct
  correct = d.parameter.centerWhich(trialNums) == d.response(trialNums);

  % bin and average
  for iVal = 1:length(d.cond(iCond).uniquePosDiff)
    % get the trials with the setting of posDiff
    whichTrials = find(d.cond(iCond).posDiff == d.cond(iCond).uniquePosDiff(iVal));
    % compute how many trials
    nTrials = length(whichTrials);
    % compute correct
    d.cond(iCond).correctBinned(iVal) = sum(correct(whichTrials))/nTrials;
    % compute ste
    d.cond(iCond).correctBinnedError(iVal) = d.cond(iCond).correctBinned(iVal)*(1-d.cond(iCond).correctBinned(iVal))/sqrt(nTrials);
  end

  % fit a cumulative gaussian to data
  d.fit(iCond) = fitCumulativeGaussian(d.cond(iCond).uniquePosDiff,d.cond(iCond).correctBinned);
end

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
  disp(sprintf('(alaisBurrAnalysis:loadStimfile) Could not find stimfile: %s',stimfileName));
  return
end

% load the file
s = load(stimfileName);
if ~isfield(s,'myscreen') || ~isfield(s,'task')
  disp(sprintf('(alaisBurrAnalysis:loadStimfile) No myscreen or task in stimfile: %s',stimfileName));
  return
end

% now check to see if this has the alaisBurr experiment in it
if ~iscell(s.task) || ~iscell(s.task{1})
  disp(sprintf('(alaiasBurrAnalysis:loadStimfile) Task variable has incorrect phases in stimfile: %s',stimfileName));
  return
end

% check task filename
taskFilename = s.task{1}{1}.taskFilename;
if isempty(strfind(lower(taskFilename),'alaisburr'))
  disp(sprintf('(alaiasBurrAnalysis:loadStimfile) Incorrect task in stimfile: %s',taskFilename));
  return
end

% parse into parameters
d = getTaskParameters(s.myscreen,s.task);
d = d{1};

% print what we found
disp(sprintf('(alaiasBurrAnalysis:loadStimfile) Found task: %s (%i trials) SID: %s Date: %s',taskFilename,d.nTrials,s.myscreen.SID,s.myscreen.starttime));

% get the variables
d.myscreen = s.myscreen;
d.task = s.task;
d.stimulus = s.stimulus;

% get experiment type
d.stimulusType = '';
if d.stimulus.visual
  d.stimulusType = sprintf('%s%s ',d.stimulusType,'Visual');
end
if d.stimulus.auditory
  d.stimulusType = sprintf('%s%s ',d.stimulusType,'Auditory');
end
if d.stimulus.bimodal
  d.stimulusType = sprintf('%s%s ',d.stimulusType,'Bimodal');
end
d.stimulusType = strtrim(d.stimulusType);

% get width
d.visualWidth = d.stimulus.width;
if isfield(s.task{1}{1}.parameter,'displacement')
  d.displacement = s.task{1}{1}.parameter.displacement;
else
  d.displacement = [];
end

% set experiment name and conditions
if d.stimulus.bimodal
  % set experiment name
  d.experimentName = sprintf('%s: [%s] width: %s',d.stimulusType,mlrnum2str(d.displacement),mlrnum2str(d.visualWidth));
  % set number of conditions
  d.nCond = length(d.displacement) * length(d.visualWidth);
  % get trials for each condition
  iCond = 1;
  for iDisplacementCond = 1:length(d.displacement)
    for iWidthCond = 1:length(d.visualWidth)
      % set condition names
      d.condNames{iCond} = sprintf('Displacement: %s width: %s',mlrnum2str(d.displacement(iDisplacementCond)),mlrnum2str(d.visualWidth(iWidthCond)));
      % get trials for this condition
      d.condTrialNums{iCond} = find((d.parameter.displacement == d.displacement(iDisplacementCond)) & (d.parameter.width == d.visualWidth(iWidthCond)));
      % remember the parameters
      d.condWidth(iCond) = d.visualWidth(iWidthCond);
      d.condDisplacement(iCond) = d.displacement(iDisplacementCond);
      % update iCond
      iCond = iCond + 1;
    end
  end
elseif d.stimulus.visual
  % set experiment name
  d.experimentName = sprintf('%s: width: %s',d.stimulusType,mlrnum2str(d.visualWidth));
  % set number of conditions
  d.nCond = length(d.visualWidth);
  for iCond = 1:d.nCond
    d.condNames{iCond} = sprintf('Visual: width: %s',mlrnum2str(d.visualWidth(iCond)));
    % set trial nums to all
    d.condTrialNums{iCond} = find(d.parameter.width == d.visualWidth(iCond));
    % remember the parameters
    d.condWidth(iCond) = d.visualWidth(iCond);
  end
else
  % set experiment name
  d.experimentName = sprintf('%s',d.stimulusType);
  % set number of conditions
  d.nCond = 1;
  % set condition name
  d.condNames{1} = sprintf('Auditory');
  % set trial nums to all
  d.condTrialNums{1} = 1:d.nTrials;
end  

% check if this is a staircase
if isfield(d.stimulus,'useStaircase') && d.stimulus.useStaircase
  d.isStaircase = 1;
else
  d.isStaircase = 0;
end

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
    stimfilePath = fullfile('~/data/alaisburr',stimfileNames);
  else
    disp(sprintf('(alaisBurrAnalysis) Could not find %s',stimfileNames));
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


  
  
  

