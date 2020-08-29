function e = estimationAnalysis(stimfileNames,varargin)

% default return argument
fit = [];

% default to working on current directory
if nargin < 1, stimfileNames = [];end

% parse arguments
getArgs(varargin,{'dispFit=1','combineData=1','numBins=33'});

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
  dispHeader(sprintf('(estimationAnalysis) Analyzing file (%i/%i): %s        ',iFile,length(stimfileNames),stimfileNames{iFile}));
  
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

% if no valid files found return
if e.nFiles == 0
  disp(sprintf('(estimationAnalysis) No files found'));
  return
end

for iFile = 1:length(stimfileNames)

%remove junk first term from response list
e.d{iFile}.task{1}{1}.randVars.calculated.est = e.d{iFile}.task{1}{1}.randVars.calculated.est(2:end)

[resp, dists] = binData(e.d{iFile})

shift = 1/e.d{iFile}.task{1}{1}.parameter.numberOffsets; 

graphDists(resp, dists, e, shift, numBins, iFile)

end
k=2
               
               
%%%%%%%%%%%%%%%%%%%%%%%%%%
%    combineStimfiles    %        
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [isNewFile e] = combineStimfiles(e,d)

% default to adding to list
isNewFile = true;

% first file
if ~isfield(e,'d'),return,end

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
function dispFits(d,iFit)

% open figure
mlrSmartfig(sprintf('%i_%s',iFit,fixBadChars(d.experimentName)),'reuse');clf;

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
  titleStr = sprintf('%s\nMean: %0.2f Std: %0.2f lambda: %0.2f goodness: %g',titleStr,d.fit(whichConds(iCond)).mean,d.fit(whichConds(iCond)).std,d.fit(whichConds(iCond)).lambda,d.fit(whichConds(iCond)).percent);

  
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
    % remember nTrials
    d.cond(iCond).nTrials(iVal) = nTrials;
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
  disp(sprintf('(estimationAnalysis:loadStimfile) Could not find stimfile: %s',stimfileName));
  return
end

% load the file
s = load(stimfileName);
if ~isfield(s,'myscreen') || ~isfield(s,'task')
  disp(sprintf('(estimationAnalysis:loadStimfile) No myscreen or task in stimfile: %s',stimfileName));
  return
end

% now check to see if this has the estimation experiment in it
if ~iscell(s.task) || ~iscell(s.task{1})
  disp(sprintf('(estimationAnalysis:loadStimfile) Task variable has incorrect phases in stimfile: %s',stimfileName));
  return
end

% check task filename
taskFilename = s.task{1}{1}.taskFilename;
if isempty(strfind(lower(taskFilename),'estimation')) & isempty(strfind(lower(taskFilename),'estimation'))
  disp(sprintf('(estimationAnalysis:loadStimfile) Incorrect task in stimfile: %s',taskFilename));
  return
end

% parse into parameters
d = getTaskParameters(s.myscreen,s.task);
d = d{1};

% print what we found
disp(sprintf('(estimationAnalysis:loadStimfile) Found task: %s (%i trials) SID: %s Date: %s',taskFilename,d.nTrials,s.myscreen.SID,s.myscreen.starttime));

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


%% binData %%
%% create a matrix of responses filtered by trial parameters %%
%% matrix is build downwards:  center offset -> AV discrepancy -> width (av1 c1 w1 -> av1 c1 w2)
%% each row is a different condition, 1000 is a blank
function [resp, dists] = binData(d)
if d.stimulusType(1) == 'B'
d.displacement = unique(d.displacement)
dists = length(d.originalTaskParameter.posDiff) * length(d.visualWidth) * length(d.originalTaskParameter.displacement)
resp(1:dists,1:length(d.task{1}{1}.randVars.calculated.est)) = 1000;
for iPos = 1:length(d.originalTaskParameter.posDiff)
   for iDiscrep = 1:length(d.originalTaskParameter.displacement)
       for iWidth = 1:length(d.visualWidth)
           for iTrial = 1:length(d.task{1}{1}.randVars.calculated.est)
               if ((d.parameter.displacement(iTrial) == d.originalTaskParameter.displacement(iDiscrep)) ...
                   && (d.parameter.posDiff(iTrial) == d.originalTaskParameter.posDiff(iPos)) ...
                   && (d.parameter.width(iTrial) == d.originalTaskParameter.width(iWidth)))
                   resp((iPos-1)*length(d.originalTaskParameter.displacement)*length(d.visualWidth) + (iDiscrep-1)*length(d.visualWidth)+iWidth, iTrial) = d.task{1}{1}.randVars.calculated.est(iTrial);
               else end
           end
       end         
   end
end
end
if d.stimulusType ~= 'B'
d.displacement = unique(d.displacement)
dists = length(d.originalTaskParameter.posDiff) * length(d.visualWidth)
resp(1:dists,1:length(d.task{1}{1}.randVars.calculated.est)) = 1000;
for iPos = 1:length(d.originalTaskParameter.posDiff)
       for iWidth = 1:length(d.visualWidth)
           for iTrial = 1:length(d.task{1}{1}.randVars.calculated.est)
               if  (d.parameter.posDiff(iTrial) == d.originalTaskParameter.posDiff(iPos)) ...
                   && (d.parameter.width(iTrial) == d.originalTaskParameter.width(iWidth))
                   resp((iPos-1)*length(d.visualWidth)+iWidth, iTrial) = d.task{1}{1}.randVars.calculated.est(iTrial);
               else end
           end
       end         
end
end


function graphDists(resp, dists, e, shift, numBins, iFile)
%initialize empty arrays for summary statistics
posOffs = []
estAvg = []
estSig = []
%pull and graph estimates at each offset
if e.d{iFile}.stimulusType(1) == 'B' %%bimodal data
for iGraph = 4:(dists-3)
    i = length(e.d{iFile}.originalTaskParameter.displacement)*length(e.d{iFile}.originalTaskParameter.width);
    k = [resp(iGraph-(2*i),1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est))+2*shift resp(iGraph-i,1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est))+shift resp(iGraph,1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est)) resp(iGraph+i,1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est))-shift resp(iGraph+(2*i),1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est))-2*shift];
    j = find(k < 101);
    estimateValues = ones(1,length(j));
    for iValue = 1:length(j);
        estimateValues(iValue) = k(j(iValue));
    end
    
    %titles
    conditions = ones(3,dists)
    %width
    conditions(1,1:dists) = repmat(e.d{iFile}.originalTaskParameter.width, 1, (length(e.d{iFile}.originalTaskParameter.posDiff)*length(e.d{iFile}.originalTaskParameter.displacement)));
    %AV discrepancy
    conditions(2,1:dists) = repmat(repelem(e.d{iFile}.originalTaskParameter.displacement, length(e.d{iFile}.originalTaskParameter.width)), 1, length(e.d{iFile}.originalTaskParameter.posDiff));
    %center offset
    conditions(3,1:dists) = repelem(e.d{iFile}.originalTaskParameter.posDiff, (length(e.d{iFile}.originalTaskParameter.displacement)*length(e.d{iFile}.originalTaskParameter.width)));
    conditions(3,1:dists) = conditions(3,1:dists)*(.5/e.d{iFile}.task{1}{1}.parameter.rightCue)+.50
    %create hist
    
    figure(iGraph-3)
    subplot(2,2,2*iFile-1)
    hist(estimateValues,(0:(1/numBins):1))
    ylim([0 20])
    titleStr = sprintf('Width: %0.2f AV diff: %0.2f Center offset: %0.2f',conditions(1,iGraph),conditions(2,iGraph),conditions(3,iGraph));
    title(titleStr)
    xlabel('Offset estimate')
    ylabel('Number of Judgements')
    hold on
    L1 = scatter(.5+((conditions(3,iGraph)+conditions(2,iGraph))/(2*e.d{iFile}.task{1}{1}.parameter.rightCue)),0,'black','DisplayName','Auditory Cue')
    L2 = scatter(.5+((conditions(3,iGraph)-conditions(2,iGraph))/(2*e.d{iFile}.task{1}{1}.parameter.rightCue)),0,'green','DisplayName','Visual Cue')
    %pdf
    [m,s] = normfit(estimateValues)
    pdf = normpdf((0:.005:1),m,s)
    plot((0:.005:1),pdf)
    L3 = scatter(m,0,'red','DisplayName','Estimate Average')
    legend([L1,L2,L3], 'Auditory location','Visual location', 'Estimate Average')
    %qqplot
    subplot(2,2,2*iFile)
    qqplot(estimateValues);
    [h, p, kstat] = lillietest(estimateValues);
    titleStr = sprintf('Normal Data: %0.2f; P = %0.2f',abs(1-h),p);
    title(titleStr)
    %label axis
    posOffs = [posOffs conditions(3,iGraph)]
    estAvg = [estAvg m]
    estSig = [estSig s]
end
end
if e.d{iFile}.stimulusType(1) ~= 'B' %% unimodal data
    for iGraph = 4:(dists-3)
    i = length(e.d{iFile}.originalTaskParameter.width);
    k = [resp(iGraph-(2*i),1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est))+2*shift resp(iGraph-i,1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est))+shift resp(iGraph,1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est)) resp(iGraph+i,1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est))-shift resp(iGraph+(2*i),1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est))-2*shift];
    j = find(k < 101);
    estimateValues = ones(1,length(j));
    for iValue = 1:length(j);
        estimateValues(iValue) = k(j(iValue));
    end
    
    %titles
    conditions = ones(3,dists)
    %width
    conditions(1,1:dists) = repmat(e.d{iFile}.originalTaskParameter.width, 1, (length(e.d{iFile}.originalTaskParameter.posDiff)));
    %AV discrepancy
    %% conditions(2,1:dists) = repmat(repelem(d.originalTaskParameter.displacement, length(d.originalTaskParameter.width)), 1, length(d.originalTaskParameter.posDiff));
    %center offset
    conditions(3,1:dists) = repelem(e.d{iFile}.originalTaskParameter.posDiff, (length(e.d{iFile}.originalTaskParameter.width)));
    conditions(3,1:dists) = conditions(3,1:dists)*(.5/e.d{iFile}.task{1}{1}.parameter.rightCue)+.5
    
    %create hist
    figure(iGraph-3)
    subplot(2,2,2*iFile-1)
    hist(estimateValues,(0:(1/numBins):1))
    ylim([0 20])
    [m,s] = normfit(estimateValues)
    %title graphs
    if e.d{iFile}.stimulusType(1) == 'V'
    titleStr = sprintf('%s: Width: %0.2f // Center offset: %0.2f // N: %0.2f // Mu: %.02f Sigma: %0.2f',e.d{iFile}.stimulusType,conditions(1,iGraph),conditions(3,iGraph),length(estimateValues),m,s);
    title(titleStr)
    end
    if e.d{iFile}.stimulusType(1) == 'A'
    titleStr = sprintf('%s: Center offset: %0.2f // N: %.02f // Mu: %.02f Sigma: %0.2f',e.d{iFile}.stimulusType,conditions(3,iGraph),length(estimateValues),m,s);
    title(titleStr)
    end
    %label axis
    xlabel('Offset estimate')
    ylabel('Number of Judgements')
    hold on
    L1 = scatter((conditions(3,iGraph)),0,'black','DisplayName','Stimulus Offset')
    %pdf
    [m,s] = normfit(estimateValues)
    pdf = normpdf((0:.005:1),m,s)
    plot((0:.005:1),pdf)
    L2 = scatter(m,0,'red','DisplayName','Estimate Average')
    legend([L1,L2], 'Cue Location', 'Estimate Average')
    %% qq plot
    subplot(2,2,2*iFile)
    qqplot(estimateValues);
    [h, p, kstat] = lillietest(estimateValues);
    titleStr = sprintf('Normal Data: %0.2f; P = %0.2f',abs(1-h),p);
    title(titleStr)
    %grab offest parameters for summary statistics
    posOffs = [posOffs conditions(3,iGraph)]
    estAvg = [estAvg m]
    estSig = [estSig s]
    end
    
 %summary statistic graphs
 figure(100)
 subplot(2,2,2*iFile-1)
 scatter(posOffs,estAvg)
 ylim([0 1])
 titleStr = sprintf('Unimodal %s Estimates',e.d{iFile}.stimulusType);
 title(titleStr)
 xlabel('Stimulus Offset')
 ylabel('Average Response')
 subplot(2,2,2*iFile)
 scatter(posOffs,estSig)
 ylim([.04 .15])
 titleStr = sprintf('Unimodal %s Variation',e.d{iFile}.stimulusType);
 title(titleStr)
 xlabel('Stimulus Offset')
 ylabel('Response Standard Deviation')
end
