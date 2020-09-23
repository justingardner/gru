
 function e = estimationAnalysis(stimfileNames,varargin)

% default return argument
fit = [];

% default to working on current directory
if nargin < 1, stimfileNames = [];end

% parse arguments
getArgs(varargin,{'dispFit=1','combineData=1','numBins=100','Neighbors=2','pNorm=0','dispGraphs=1'});
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

%switch order in which we do analysis (cleaner imo but can change)
if e.nFiles == 3
    for iFile = 1:3
        if e.d{iFile}.stimulusType(1) == 'A'
            a = e.d{iFile}
        elseif e.d{iFile}.stimulusType(1) == 'V'
            v = e.d{iFile}
        elseif e.d{iFile}.stimulusType(1) == 'B'
            b = e.d{iFile}
        end
    end
    e.d{1} = a
    e.d{2} = v
    e.d{3} = b
    e.d{3}.originalTaskParameter.displacement = unique(e.d{3}.originalTaskParameter.displacement)
    numSubs = e.nFiles+length(e.d{3}.originalTaskParameter.displacement)-1
end
if e.nFiles < 3
    numSubs = e.nFiles
end
numSkips = Neighbors

%% bin the data and graph responses in each condition %%
for iFile = 1:e.nFiles
e.d{iFile}.task{1}{1}.randVars.calculated.est = e.d{iFile}.task{1}{1}.randVars.calculated.est(2:end) %trim files (the last response (nan) is put first, but parameter isnt)

[resp, dists] = binData(e, iFile) %organize responses
for val = 1:length(resp);
    if (resp(val) > 1) && (resp(val) < 2); resp(val) = 1; end
    if resp(val) < 0; resp(vap) = 0; end
end
shift = 1/(e.d{iFile}.task{1}{1}.parameter.numberOffsets-1); %graph input

[e] = graphDists(resp, dists, e, shift, numBins, iFile, numSubs, numSkips, Neighbors, pNorm)

end

%% calculate and graph log likelihood fits of models
[loglikes] = modelCompare(e,numSkips)
graphLikelihoods(loglikes,numSkips,e,pNorm)
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
d.task{1}{1}.randVars.calculated.est = [d1.task{1}{1}.randVars.calculated.est d2.task{1}{1}.randVars.calculated.est(2:end)]
d.originalTaskParameter.posDiff =  d1.originalTaskParameter.posDiff
d.originalTaskParameter.displacement = d1.originalTaskParameter.displacement
d.originalTaskParameter.width = d1.originalTaskParameter.width
d.task{1}{1}.parameter.numberOffsets = d1.task{1}{1}.parameter.numberOffsets
d.task{1}{1}.parameter.rightCue = d1.task{1}{1}.parameter.rightCue
d.parameter.posDiff = [d1.parameter.posDiff(1:end-1) d2.parameter.posDiff]
d.parameter.width = [d1.parameter.width(1:end-1) d2.parameter.width]
d.parameter.displacement = [d1.parameter.displacement(1:end-1) d2.parameter.displacement]


% copy these fileds
%copyFields = {'parameter','randVars'};
%for iField = 1:length(copyFields)
%  fieldsToConcat = fieldnames(d1.(copyFields{iField}));
%  for iConcatField = 1:length(fieldsToConcat)
%    d.(copyFields{iField}).(fieldsToConcat{iConcatField}) = [d1.(copyFields{iField}).(fieldsToConcat{iConcatField}) d2.(copyFields{iField}).(fieldsToConcat{iConcatField})];
%  end
%end

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
function [resp, dists] = binData(e, iFile)
if e.d{iFile}.stimulusType(1) == 'B'
dists = length(e.d{iFile}.originalTaskParameter.posDiff) * length(e.d{iFile}.visualWidth) * length(e.d{iFile}.originalTaskParameter.displacement)
resp(1:dists,1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est)) = 1000;
for iPos = 1:length(e.d{iFile}.originalTaskParameter.posDiff)
   for iDiscrep = 1:length(e.d{iFile}.originalTaskParameter.displacement)
       for iWidth = 1:length(e.d{iFile}.visualWidth)
           for iTrial = 1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est)
               if ((e.d{iFile}.parameter.displacement(iTrial) == e.d{iFile}.originalTaskParameter.displacement(iDiscrep)) ...
                   && (e.d{iFile}.parameter.posDiff(iTrial) == e.d{iFile}.originalTaskParameter.posDiff(iPos)) ...
                   && (e.d{iFile}.parameter.width(iTrial) == e.d{iFile}.originalTaskParameter.width(iWidth)))
                   resp((iPos-1)*length(e.d{iFile}.originalTaskParameter.displacement)*length(e.d{iFile}.visualWidth) + (iDiscrep-1)*length(e.d{iFile}.visualWidth)+iWidth, iTrial) = e.d{iFile}.task{1}{1}.randVars.calculated.est(iTrial);
               else end
           end
       end         
   end
end
end
if e.d{iFile}.stimulusType ~= 'B'
e.d{iFile}.displacement = unique(e.d{iFile}.displacement)
dists = length(e.d{iFile}.originalTaskParameter.posDiff) * length(e.d{iFile}.visualWidth)
resp(1:dists,1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est)) = 1000;
for iPos = 1:length(e.d{iFile}.originalTaskParameter.posDiff)
       for iWidth = 1:length(e.d{iFile}.visualWidth)
           for iTrial = 1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est)
               if  (e.d{iFile}.parameter.posDiff(iTrial) == e.d{iFile}.originalTaskParameter.posDiff(iPos)) ...
                   && (e.d{iFile}.parameter.width(iTrial) == e.d{iFile}.originalTaskParameter.width(iWidth))
                   resp((iPos-1)*length(e.d{iFile}.visualWidth)+iWidth, iTrial) = e.d{iFile}.task{1}{1}.randVars.calculated.est(iTrial);
               else end
           end
       end         
end
end


function [e] = graphDists(resp, dists, e, shift, numBins, iFile, numSubs, numSkips, Neighbors, pNorm)
%initialize empty arrays for summary statistics
posOffs = []
estAvg = []
estSig = []

%pull and graph estimates at each offset

%% bimodal data %%
if e.d{iFile}.stimulusType(1) == 'B'
e.d{iFile}.respMatrix = ones(length(e.d{3}.originalTaskParameter.displacement)*(51-2*numSkips),250)+4; %%hardcoded for 51 offsets and max 250 responses
e.d{iFile}.normality = []
for iGraph = numSkips*length(e.d{3}.originalTaskParameter.displacement)+1:(dists-numSkips*length(e.d{3}.originalTaskParameter.displacement))
    i = length(e.d{iFile}.originalTaskParameter.displacement)*length(e.d{iFile}.originalTaskParameter.width);
    k = [];
    for neighbor = Neighbors:(-1):(-Neighbors);
        k = [k resp(iGraph-(neighbor*i),1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est))+neighbor*shift];
    end
    %k = [resp(iGraph-(2*i),1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est))+2*shift resp(iGraph-i,1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est))+shift resp(iGraph,1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est)) resp(iGraph+i,1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est))-shift resp(iGraph+(2*i),1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est))-2*shift];
    j = find(k < 101);
    estimateValues = ones(1,length(j));
    for iValue = 1:length(j);
        estimateValues(iValue) = k(j(iValue));
        e.d{iFile}.respMatrix(iGraph-numSkips*length(e.d{3}.originalTaskParameter.displacement),iValue) = k(j(iValue));
    end
    
    %titles
    conditions = ones(3,dists);
    %width
    conditions(1,1:dists) = repmat(e.d{iFile}.originalTaskParameter.width, 1, (length(e.d{iFile}.originalTaskParameter.posDiff)*length(e.d{iFile}.originalTaskParameter.displacement)));
    %AV discrepancy
    conditions(2,1:dists) = repmat(repelem(e.d{iFile}.originalTaskParameter.displacement, length(e.d{iFile}.originalTaskParameter.width)), 1, length(e.d{iFile}.originalTaskParameter.posDiff));
    conditions(2,1:dists) = conditions(2,1:dists)*(.5/e.d{iFile}.task{1}{1}.parameter.rightCue);
    %center offset
    conditions(3,1:dists) = repelem(e.d{iFile}.originalTaskParameter.posDiff, (length(e.d{iFile}.originalTaskParameter.displacement)*length(e.d{iFile}.originalTaskParameter.width)));
    conditions(3,1:dists) = conditions(3,1:dists)*(.5/e.d{iFile}.task{1}{1}.parameter.rightCue)+.50;
    e.d{iFile}.conditions = conditions;
    %create hist
    
    figure(ceil((iGraph-numSkips*length(e.d{3}.originalTaskParameter.displacement))/length(e.d{3}.originalTaskParameter.displacement)))
    subplot(2,numSubs,2*iFile-1+2*mod(iGraph-1, length(e.d{3}.originalTaskParameter.displacement)))
    hist(estimateValues,(0:(1/numBins):1.02))
    %ylim([0 20])
    xlim([-.1 1.1])
    [m,s] = normfit(estimateValues);
    titleStr = sprintf('Bimodal: Discrepancy: %0.2f Position: %0.2f // N: %0.2f // Mu: %.02f Sigma: %0.2f',conditions(2,iGraph),conditions(3,iGraph),length(estimateValues),m,s);
    title(titleStr)
    xlabel('Offset estimate')
    ylabel('Number of Judgements')
    hold on
    L1 = scatter(conditions(3,iGraph)-conditions(2,iGraph),0,'black','DisplayName','Auditory Cue');
    L2 = scatter(conditions(3,iGraph)+conditions(2,iGraph),0,'green','DisplayName','Visual Cue');
    %pdf
    pdf = normpdf((0:.005:1),m,s)
    %plot((0:.005:1),pdf,'red')
    L3 = scatter(m,0,'red','DisplayName','Estimate Average');
    legend([L1,L2,L3],'Auditory location','Visual location','Estimate Average')
    %qqplot
    subplot(2,numSubs,2*iFile++2*mod(iGraph-1, length(e.d{3}.originalTaskParameter.displacement)))
    qqplot(estimateValues);
    [h, p, kstat] = lillietest(estimateValues);
    e.d{iFile}.normality = [e.d{iFile}.normality p]
    titleStr = sprintf('QQ plot: P = %0.2f',p);
    title(titleStr)
    %label axis
    posOffs = [posOffs conditions(3,iGraph)];
    estAvg = [estAvg m];
    estSig = [estSig s];
end
figure(100)

 for val = 1:2:length(posOffs) %only plot normal data - HARD CODED for 2 discrepancies and 2.4 offset
     if (e.d{1}.normality(ceil(val/2)) < pNorm) || (e.d{2}.normality(ceil(val/2)) < pNorm)
         estAvg(val) = NaN
         estSig(val) = NaN
         posOff(val) = NaN
     end
 end    
  for val = 8:2:(length(posOffs)-6) %only plot normal data
     if (e.d{1}.normality(ceil(val/2)-3) < pNorm) || (e.d{2}.normality(ceil(val/2)+3) < pNorm)
         estAvg(val) = NaN
         estSig(val) = NaN
         posOff(val) = NaN
     end
 end   
 
 subplot(3,2,2*iFile-1)
 for x = 1:length(e.d{3}.originalTaskParameter.displacement)
 scatter(posOffs(x:length(e.d{3}.originalTaskParameter.displacement):end),estAvg(x:length(e.d{3}.originalTaskParameter.displacement):end));
 hold on
 end
 ylim([0 1])
 titleStr = sprintf('Bimodal Estimates',e.d{iFile}.stimulusType);
 title(titleStr)
 xlabel('Stimulus Offset')
 ylabel('Average Response')
 %%%legend
 if length(e.d{3}.originalTaskParameter.displacement) == 2
 leg1 = legend(['Discrepancy: ' num2str(e.d{3}.originalTaskParameter.displacement(1))],['Discrepancy: ' num2str(e.d{3}.originalTaskParameter.displacement(2))])
 end
 if length(e.d{3}.originalTaskParameter.displacement) == 3
 leg1 = legend(['Discrepancy: ' num2str(e.d{3}.originalTaskParameter.displacement(1))],['Discrepancy: ' num2str(e.d{3}.originalTaskParameter.displacement(2))],['Discrepancy: ' num2str(e.d{3}.originalTaskParameter.displacement(2))])  
 end
 subplot(3,2,2*iFile)
 for x = 1:length(e.d{3}.originalTaskParameter.displacement)
 scatter(posOffs(x:length(e.d{3}.originalTaskParameter.displacement):end),estSig(x:length(e.d{3}.originalTaskParameter.displacement):end))
 hold on
 end
 %ylim([.04 .15])
 titleStr = sprintf('Unimodal %s Variation',e.d{iFile}.stimulusType);
 title(titleStr)
 xlabel('Stimulus Offset')
 ylabel('Response Standard Deviation')
 %%%legend
 if length(e.d{3}.originalTaskParameter.displacement) == 2
 legend(['Discrepancy: ' num2str(e.d{3}.originalTaskParameter.displacement(1))],['Discrepancy: ' num2str(e.d{3}.originalTaskParameter.displacement(2))])
 end
 if length(e.d{3}.originalTaskParameter.displacement) == 3
 legend(['Discrepancy: ' num2str(e.d{3}.originalTaskParameter.displacement(1))],['Discrepancy: ' num2str(e.d{3}.originalTaskParameter.displacement(2))],['Discrepancy: ' num2str(e.d{3}.originalTaskParameter.displacement(2))])  
 end

 figure(99)
 subplot(1,3,iFile)
 hold on
 bimodalbias = estAvg-posOffs
 a = scatter(posOffs(1:2:end),bimodalbias(1:2:end)) %hard coded
 ylim([ -.2 .2])
 b = scatter(posOffs(2:2:end),bimodalbias(2:2:end))
 hold on
 plot([0 1],[0 0],'k')
 plot([.5 .5],[-.2 .2],'k')
 xlabel('Stimulus Offset')
 ylabel('Estimate Bias')
 legend([a,b],'No discrepancy','.06 discrepancy')
 titleStr = sprintf('Bimodal Bias',e.d{iFile}.stimulusType);
 title(titleStr)
 
end

%% unimodal data %%
if e.d{iFile}.stimulusType(1) ~= 'B'
e.d{iFile}.respMatrix = ones(51-2*numSkips,250)+4; %%hardcoded for 51 offsets and max 250 responses
e.d{iFile}.normality = []
    for iGraph = numSkips+1:(dists-numSkips)
    i = length(e.d{iFile}.originalTaskParameter.width);
    k = [];
    for neighbor = Neighbors:(-1):(-Neighbors);
        k = [k resp(iGraph-(neighbor*i),1:length(e.d{iFile}.task{1}{1}.randVars.calculated.est))+neighbor*shift];
    end
    j = find(k < 101);
    estimateValues = ones(1,length(j));
    for iValue = 1:length(j);
        estimateValues(iValue) = k(j(iValue));
        e.d{iFile}.respMatrix(iGraph-numSkips,iValue) = k(j(iValue));
    end
    
    %titles
    conditions = ones(3,dists);
    %width
    conditions(1,1:dists) = repmat(e.d{iFile}.originalTaskParameter.width, 1, (length(e.d{iFile}.originalTaskParameter.posDiff)));
    %center offset
    conditions(3,1:dists) = repelem(e.d{iFile}.originalTaskParameter.posDiff, (length(e.d{iFile}.originalTaskParameter.width)));
    conditions(3,1:dists) = conditions(3,1:dists)*(.5/e.d{iFile}.task{1}{1}.parameter.rightCue)+.5;
    e.d{iFile}.conditions = conditions;
    
    %create hist
    figure(iGraph-numSkips)
    subplot(2,numSubs,2*iFile-1)
    hist(estimateValues,(0:(1/numBins):1))
    %ylim([0 20])
    xlim([-.1 1.1])
    [m,s] = normfit(estimateValues);
    %title graphs
    if e.d{iFile}.stimulusType(1) == 'V'
    titleStr = sprintf('%s: Position: %0.2f // N: %0.2f // Mu: %.02f Sigma: %0.2f',e.d{iFile}.stimulusType,conditions(3,iGraph),length(estimateValues),m,s);
    title(titleStr)
    end
    if e.d{iFile}.stimulusType(1) == 'A'
    titleStr = sprintf('%s: Position: %0.2f // N: %.02f // Mu: %.02f Sigma: %0.2f',e.d{iFile}.stimulusType,conditions(3,iGraph),length(estimateValues),m,s);
    title(titleStr)
    end
    %label axis
    xlabel('Offset estimate')
    ylabel('Number of Judgements')
    hold on
    L1 = scatter((conditions(3,iGraph)),0,'black','DisplayName','Stimulus Offset');
    %pdf
    [m,s] = normfit(estimateValues);
    pdf = normpdf((0:.005:1),m,s);
    plot((0:.005:1),pdf)
    L2 = scatter(m,0,'red','DisplayName','Estimate Average');
    legend([L1,L2], 'Cue Location', 'Estimate Average')
    %% qq plot
    subplot(2,numSubs,2*iFile)
    qqplot(estimateValues);
    [h, p, kstat] = lillietest(estimateValues);
    e.d{iFile}.normality = [e.d{iFile}.normality p]
    titleStr = sprintf('QQ plot: P = %0.2f',p);
    title(titleStr)
    %grab offest parameters for summary statistics
    posOffs = [posOffs conditions(3,iGraph)];
    estAvg = [estAvg m];
    estSig = [estSig s];
    end
    
 %summary statistic graphs
 figure(100)
 for val = 1:length(posOffs) %only plot normal data
     if e.d{iFile}.normality(val) < pNorm
         estAvg(val) = NaN
         estSig(val) = NaN
         posOff(val) = NaN
     end
 end       
 subplot(3,2,2*iFile-1)
 scatter(posOffs,estAvg)
 ylim([0 1])
 titleStr = sprintf('Unimodal %s Estimates',e.d{iFile}.stimulusType);
 title(titleStr)
 xlabel('Stimulus Offset')
 ylabel('Average Response')
 subplot(3,2,2*iFile)
 scatter(posOffs,estSig)
 ylim([.04 .15])
 titleStr = sprintf('Unimodal %s Variation',e.d{iFile}.stimulusType);
 title(titleStr)
 xlabel('Stimulus Offset')
 ylabel('Response Standard Deviation')
 
 figure(99)
 hold on
 subplot(1,3,iFile)
 scatter(posOffs,(estAvg-posOffs))
 ylim([ -.2 .2])
 hold on
 plot([0 1],[0 0],'k')
 plot([.5 .5],[-.2 .2],'k')
 xlabel('Stimulus Offset')
 ylabel('Estimate Bias')
 titleStr = sprintf('Unimodal %s Bias',e.d{iFile}.stimulusType);
 title(titleStr)
end

function [loglikes] = modelCompare(e,numSkips)  %%for now, hard coded for the 2 offsets and 51 disps. Need to change when we collect a lot of data, and needs to be the same conditions.
loglikes = [] %%first row OI, 2nd OS, 3rd visual capture, 4th auditory capture
for offset = 7:(102-numSkips*4-6) %hard coded for delta=2.4 with 51 offsets:: 51*2=102 total conditions, take 4 off for each skip (2 deltas each side) and 6 (2.4/50) each side. if you change the discrepnacy, need to change all the 6's. change number of discreps, change 4.
    
    %%%%%%% optimal integration %%%%%%%
    if mod(offset,2) %no discrepancy
        %%unimodal inputs
        audResp = e.d{1}.respMatrix(ceil(offset/2),1:250)
        audResp = audResp(audResp<2)
        [muA,sA] = normfit(audResp)
        visResp = e.d{2}.respMatrix(ceil(offset/2),1:250)
        visResp = visResp(visResp<2)
        [muV,sV] = normfit(visResp)
        %bimodal predictions
        weightV = (sA*sA)/(sA*sA+sV*sV)
        weightA = 1-weightV
        oiMu = weightV*muV + weightA*muA
        oiS = sqrt(sA*sA*sV*sV/(sA*sA+sV*sV))
        oiDist = normpdf((0:.01:1),oiMu,oiS)/100
        %grab bimodal responses
        OIloglike = 0
        biResp = e.d{3}.respMatrix(offset,1:250)
        biResp = biResp(biResp<2)
        %fix some reporting errors where matlab gives values >1 or 0>
        for val = 1:length(biResp)
            if biResp(val) < 0
                biResp(val) = 0
            end
            if biResp(val) > 1
                biResp(val) = 1
            end
        end
        %calculate log likelihood
        for val = 1:length(biResp)
            OIloglike = OIloglike + log(oiDist(round(biResp(val)*100+1)))
        end
        loglikes(1,offset-6) = OIloglike
        figure(ceil(offset/2)); subplot(2,4,5); OIn = plot(0:.01:1,oiDist*100,'cyan'); % graph model prediction at appropriate condition, all hardcoded rn for number of conditions
    end
     if ~mod(offset,2) %discrepancy
        %%unimodal inputs
        audResp = e.d{1}.respMatrix(ceil(offset/2)-3,1:250)
        audResp = audResp(audResp<2)
        [muA,sA] = normfit(audResp)
        visResp = e.d{2}.respMatrix(ceil(offset/2)+3,1:250)
        visResp = visResp(visResp<2)
        [muV,sV] = normfit(visResp)
        %bimodal predictions
        weightV = (sA*sA)/(sA*sA+sV*sV)
        weightA = 1-weightV
        oiMu = weightV*muV + weightA*muA
        oiS = sqrt(sA*sA*sV*sV/(sA*sA+sV*sV))
        oiDist = normpdf((0:.01:1),oiMu,oiS)/100
        %grab bimodal responses
        OIloglike = 0
        biResp = e.d{3}.respMatrix(offset,1:250)
        biResp = biResp(biResp<2)
        %fix some reporting errors where matlab gives values >1 or 0>
        for val = 1:length(biResp)
            if biResp(val) < 0
                biResp(val) = 0
            end
            if biResp(val) > 1
                biResp(val) = 1
            end
        end
        %calc log likelihood
        for val = 1:length(biResp)
            OIloglike = OIloglike + log(oiDist(round(biResp(val)*100+1)))
        end
        loglikes(1,offset-6) = OIloglike
        figure(ceil(offset/2)); subplot(2,4,7); OIy = plot(0:.01:1,oiDist*100,'cyan'); % graph model prediction at appropriate condition
     end
     
     %%%%%%%%%% %optimal switching %%%%%%%%%%%%
    if mod(offset,2) %no discrepancy
        %%unimodal inputs
        audResp = e.d{1}.respMatrix(ceil(offset/2),1:250)
        audResp = audResp(audResp<2)
        [muA,sA] = normfit(audResp)
        Adist = normpdf((0:.01:1),muA,sA)/100
        visResp = e.d{2}.respMatrix(ceil(offset/2),1:250)
        visResp = visResp(visResp<2)
        [muV,sV] = normfit(visResp)
        Vdist = normpdf((0:.01:1),muV,sV)/100
        %bimodal predictions
        weightV = (sA*sA)/(sA*sA+sV*sV)
        weightA = 1-weightV
        osDist = weightV*Vdist+weightA*Adist
        %grab bimodal responses
        OSloglike = 0
        biResp = e.d{3}.respMatrix(offset,1:250)
        biResp = biResp(biResp<2)
         %fix some reporting errors where matlab gives values >1 or 0>
        for val = 1:length(biResp)
            if biResp(val) < 0
                biResp(val) = 0
            end
            if biResp(val) > 1
                biResp(val) = 1
            end
        end
        %calculate log likelihood
        for val = 1:length(biResp)
            OSloglike = OSloglike + log(osDist(round(biResp(val)*100+1)))
        end
        loglikes(2,offset-6) = OSloglike
        figure(ceil(offset/2)); subplot(2,4,5); OSn = plot(0:.01:1,osDist*100,'magenta'); % graph model prediction at appropriate condition
    end
    if ~mod(offset,2) %discrepancy
        %%unimodal inputs
        audResp = e.d{1}.respMatrix(ceil(offset/2)-3,1:250)
        audResp = audResp(audResp<2)
        [muA,sA] = normfit(audResp)
        Adist = normpdf((0:.01:1),muA,sA)/100
        visResp = e.d{2}.respMatrix(ceil(offset/2)+3,1:250)
        visResp = visResp(visResp<2)
        [muV,sV] = normfit(visResp)
        Vdist = normpdf((0:.01:1),muV,sV)/100
        %bimodal predictions
        weightV = (sA*sA)/(sA*sA+sV*sV)
        weightA = 1-weightV
        osDist = weightV*Vdist+weightA*Adist
        %grab bimodal responses
        OSloglike = 0
        biResp = e.d{3}.respMatrix(offset,1:250)
        biResp = biResp(biResp<2)
        %fix some reporting errors where matlab gives values >1 or 0>
        for val = 1:length(biResp)
            if biResp(val) < 0
                biResp(val) = 0
            end
            if biResp(val) > 1
                biResp(val) = 1
            end
        end
        %calculate log likelihood
        for val = 1:length(biResp)
            OSloglike = OSloglike + log(osDist(round(biResp(val)*100+1)))
        end
        loglikes(2,offset-6) = OSloglike
        figure(ceil(offset/2)); subplot(2,4,7); OSy = plot(0:.01:1,osDist*100,'magenta'); % graph model prediction at appropriate condition
    end
    
    
    %%%%%%% visual capture %%%%%%%
    if mod(offset,2) %no discrepancy
        %%unimodal inputs
        audResp = e.d{1}.respMatrix(ceil(offset/2),1:250)
        audResp = audResp(audResp<2)
        [muA,sA] = normfit(audResp)
        visResp = e.d{2}.respMatrix(ceil(offset/2),1:250)
        visResp = visResp(visResp<2)
        [muV,sV] = normfit(visResp)
        %bimodal predictions
        weightV = 1
        weightA = 0
        vcMu = weightV*muV + weightA*muA
        vcS = sV
        vcDist = normpdf((0:.01:1),vcMu,vcS)/100
        %grab bimodal responses
        vcloglike = 0
        biResp = e.d{3}.respMatrix(offset,1:250)
        biResp = biResp(biResp<2)
        %fix some reporting errors where matlab gives values >1 or 0>
        for val = 1:length(biResp)
            if biResp(val) < 0
                biResp(val) = 0
            end
            if biResp(val) > 1
                biResp(val) = 1
            end
        end
        %calculate log likelihood
        for val = 1:length(biResp)
          vcloglike = vcloglike + log(vcDist(round(biResp(val)*100+1)))
        end
        loglikes(3,offset-6) = vcloglike
        figure(ceil(offset/2)); subplot(2,4,5); VCn = plot(0:.01:1,vcDist*100,'green'); % graph model prediction at appropriate condition
    end
     if ~mod(offset,2) %discrepancy
        %%unimodal inputs
        audResp = e.d{1}.respMatrix(ceil(offset/2)-3,1:250)
        audResp = audResp(audResp<2)
        [muA,sA] = normfit(audResp)
        visResp = e.d{2}.respMatrix(ceil(offset/2)+3,1:250)
        visResp = visResp(visResp<2)
        [muV,sV] = normfit(visResp)
        %bimodal predictions
        weightV = 1
        weightA = 0
        vcMu = weightV*muV + weightA*muA
        vcS = sV
        vcDist = normpdf((0:.01:1),vcMu,vcS)/100
        %grab bimodal responses
        vcloglike = 0
        biResp = e.d{3}.respMatrix(offset,1:250)
        biResp = biResp(biResp<2)
        %fix some reporting errors where matlab gives values >1 or 0>
        for val = 1:length(biResp)
            if biResp(val) < 0
                biResp(val) = 0
            end
            if biResp(val) > 1
                biResp(val) = 1
            end
        end
        %calc log likelihood
        for val = 1:length(biResp)
            vcloglike = vcloglike + log(vcDist(round(biResp(val)*100+1)))
        end
        loglikes(3,offset-6) = vcloglike
        figure(ceil(offset/2)); subplot(2,4,7); VCy = plot(0:.01:1,vcDist*100,'green'); % graph model prediction at appropriate condition
     end
  
     %%%%%%%%%% auditory capture %%%%%%%%%%%%
    if mod(offset,2) %no discrepancy
        %%unimodal inputs
        audResp = e.d{1}.respMatrix(ceil(offset/2),1:250)
        audResp = audResp(audResp<2)
        [muA,sA] = normfit(audResp)
        visResp = e.d{2}.respMatrix(ceil(offset/2),1:250)
        visResp = visResp(visResp<2)
        [muV,sV] = normfit(visResp)
        %bimodal predictions
        weightV = 0
        weightA = 1
        acMu = weightV*muV + weightA*muA
        acS = sA
        acDist = normpdf((0:.01:1),acMu,acS)/100
        %grab bimodal responses
        acloglike = 0
        biResp = e.d{3}.respMatrix(offset,1:250)
        biResp = biResp(biResp<2)
        %fix some reporting errors where matlab gives values >1 or 0>
        for val = 1:length(biResp)
            if biResp(val) < 0
                biResp(val) = 0
            end
            if biResp(val) > 1
                biResp(val) = 1
            end
        end
        %calculate log likelihood
        for val = 1:length(biResp)
          acloglike = acloglike + log(acDist(round(biResp(val)*100+1)))
        end
        loglikes(4,offset-6) = acloglike
        figure(ceil(offset/2)); subplot(2,4,5); ACn = plot(0:.01:1,acDist*100,'black'); % graph model prediction at appropriate condition
        leg2 = legend([OIn,OSn,VCn,ACn],'Optimal Integration','Optimal Switching','Visual Capture','Auditory Capture')
    end
     if ~mod(offset,2) %discrepancy
        %%unimodal inputs
        audResp = e.d{1}.respMatrix(ceil(offset/2)-3,1:250)
        audResp = audResp(audResp<2)
        [muA,sA] = normfit(audResp)
        visResp = e.d{2}.respMatrix(ceil(offset/2)+3,1:250)
        visResp = visResp(visResp<2)
        [muV,sV] = normfit(visResp)
        %bimodal predictions
        weightV = 0
        weightA = 1
        acMu = weightV*muV + weightA*muA
        acS = sA
        acDist = normpdf((0:.01:1),acMu,acS)/100
        %grab bimodal responses
        acloglike = 0
        biResp = e.d{3}.respMatrix(offset,1:250)
        biResp = biResp(biResp<2)
        %fix some reporting errors where matlab gives values >1 or 0>
        for val = 1:length(biResp)
            if biResp(val) < 0
                biResp(val) = 0
            end
            if biResp(val) > 1
                biResp(val) = 1
            end
        end
        %calc log likelihood
        for val = 1:length(biResp)
            acloglike = acloglike + log(acDist(round(biResp(val)*100+1)))
        end
        loglikes(4,offset-6) = acloglike
        figure(ceil(offset/2)); subplot(2,4,7); ACy = plot(0:.01:1,acDist*100,'black'); % graph model prediction at appropriate condition
        leg2 = legend([OIy,OSy,VCy,ACy],'Optimal Integration','Optimal Switching','Visual Capture','Auditory Capture')
     end

end


function graphLikelihoods(loglikes,numSkips,e,pNorm)
%%log likelihoods by cue position

%%filter by normality
for model = 1:4 % change to number of models
    for val = 1:2:(102-numSkips*4-6-6)
        if (e.d{1}.normality(ceil(val/2)+3) < pNorm) || (e.d{2}.normality(ceil(val/2)+3) < pNorm)
            loglikes(model,val) = NaN
            conditions(3,ceil(val/2)+numSkips+3) = NaN
        end
    end
end     
for model = 1:4 % change to number of models
    for val = 2:2:(102-numSkips*4-6-6)
        if (e.d{1}.normality(ceil(val/2)) < pNorm) || (e.d{2}.normality(ceil(val/2)+6) < pNorm)
            loglikes(model,val) = NaN
            conditions(3,ceil(val/2)+numSkips+3) = NaN
        end
    end
end     

figure(101)
subplot(2,2,1) %optimal integration
OIno = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),loglikes(1,1:2:(102-numSkips*4-6-6)))
hold on
OIyo = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),loglikes(1,2:2:(102-numSkips*4-6-6)))
xlim([0 1])
legend([OIno,OIyo],'Discrepancy','No Discrepancy')
titleStr = sprintf('Optimal Integration: negative log likelihoods by cue position'); title(titleStr); xlabel('Cue offset'); ylabel('log likelihood');
subplot(2,2,2) %optimal switching
OSno = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),loglikes(2,1:2:(102-numSkips*4-6-6)))
hold on
OSyo = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),loglikes(2,2:2:(102-numSkips*4-6-6)))
xlim([0 1])
legend([OSno,OSyo],'Discrepancy','No Discrepancy')
titleStr = sprintf('Optimal Switching: negative log likelihoods by cue position'); title(titleStr); xlabel('Cue offset'); ylabel('log likelihood');
subplot(2,2,3) %visual capture
vcno = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),loglikes(3,1:2:(102-numSkips*4-6-6)))
hold on
vcyo = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),loglikes(3,2:2:(102-numSkips*4-6-6)))
xlim([0 1])
legend([vcno,vcyo],'Discrepancy','No Discrepancy')
titleStr = sprintf('Visual Capture: negative log likelihoods by cue position'); title(titleStr); xlabel('Cue offset'); ylabel('log likelihood');
subplot(2,2,4) %Auditory capture
acno = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),loglikes(4,1:2:(102-numSkips*4-6-6)))
hold on
acyo = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),loglikes(4,2:2:(102-numSkips*4-6-6)))
xlim([0 1])
legend([acno,acyo],'Discrepancy','No Discrepancy')
titleStr = sprintf('Auditory capture: negative log likelihoods by cue position'); title(titleStr); xlabel('Cue offset'); ylabel('log likelihood');


%%%%%%%%%%% loglikehood ratios between models %%%%%%%%%%%%%%%
figure(102) %% optimal integration
OIvsOS = [] 
for val = 1:length(loglikes(1,1:end))
    OIvsOS = [OIvsOS loglikes(1,val)-loglikes(2,val)]
end
subplot(1,2,1)
OIvOSno = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),OIvsOS(1:2:(102-numSkips*4-6-6))); hold on; xlim([0 1]);
xlabel('Cue offset'); ylabel('Negative log likelihood ratio'); titleStr = sprintf('Optimal Integration Ratios: no AV discrepancy'); legend([OIvOSno],'Optimal Switching'); title(titleStr);
subplot(1,2,2)
OIvOSyo = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),OIvsOS(2:2:(102-numSkips*4-6-6))); hold on; xlim([0 1]);
xlim([0 1]); xlabel('Cue offset'); ylabel('Negative log likelihood ratio'); titleStr = sprintf('Optimal Integration Ratios: AV discrepancy = 6'); legend([OIvOSyo],'Optimal Switching'); title(titleStr)

OIvsVC = []
for val = 1:length(loglikes(1,1:end))
    OIvsVC = [OIvsVC loglikes(1,val)-loglikes(3,val)]
end
subplot(1,2,1)
OIvVCno = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),OIvsVC(1:2:(102-numSkips*4-6-6))); hold on; xlim([0 1]);
legend([OIvVCno],'Visual Capture');
subplot(1,2,2)
OIvVCyo = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),OIvsVC(2:2:(102-numSkips*4-6-6)));
legend([OIvVCyo],'Visual Capture');

OIvsAC = []
for val = 1:length(loglikes(1,1:end))
    OIvsAC = [OIvsAC loglikes(1,val)-loglikes(4,val)]
end
subplot(1,2,1)
OIvACno = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),OIvsAC(1:2:(102-numSkips*4-6-6))); hold on; xlim([0 1]);
legend([OIvOSno,OIvVCno,OIvACno],'Optimal Switching','Visual Capture','Auditory Capture');
subplot(1,2,2)
OIvACyo = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),OIvsAC(2:2:(102-numSkips*4-6-6)));
legend([OIvOSyo,OIvVCyo,OIvACyo],'Optimal Switching','Visual Capture','Auditory Capture');



figure(103) %% optimal switching
OSvsOI = []
for val = 1:length(loglikes(1,1:end))
    OSvsOI = [OSvsOI loglikes(2,val)-loglikes(1,val)]
end
subplot(1,2,1)
OSvOIno = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),OSvsOI(1:2:(102-numSkips*4-6-6))); hold on; xlim([0 1]);
xlabel('Cue offset'); ylabel('Negative log likelihood ratio'); titleStr = sprintf('Optimal Switching Ratios: no AV discrepancy'); legend([OSvOIno],'Optimal Integration'); title(titleStr);
subplot(1,2,2)
OSvOIyo = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),OSvsOI(2:2:(102-numSkips*4-6-6))); hold on; xlim([0 1]);
xlim([0 1]); xlabel('Cue offset'); ylabel('Negative log likelihood ratio'); titleStr = sprintf('Optimal Switching Ratios: AV discrepancy = 6'); legend([OSvOIyo],'Optimal Integration'); title(titleStr)
        
OSvsVC = []
for val = 1:length(loglikes(1,1:end))
    OSvsVC = [OSvsVC loglikes(2,val)-loglikes(3,val)]
end
subplot(1,2,1)
OSvVCno = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),OSvsVC(1:2:(102-numSkips*4-6-6))); hold on; xlim([0 1]);
legend([OSvVCno],'Visual Capture');
subplot(1,2,2)
OSvVCyo = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),OSvsVC(2:2:(102-numSkips*4-6-6)));
legend([OSvVCyo],'Visual Capture');


OSvsAC = []
for val = 1:length(loglikes(1,1:end))
    OSvsAC = [OSvsAC loglikes(2,val)-loglikes(4,val)]
end
subplot(1,2,1)
OSvACno = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),OSvsAC(1:2:(102-numSkips*4-6-6))); hold on; xlim([0 1]);
legend([OSvOIno,OSvVCno,OSvACno],'Optimal Integration','Visual Capture','Auditory Capture');
subplot(1,2,2)
OSvACyo = scatter(e.d{3}.conditions(3,(1+numSkips*2+6):2:(102-numSkips*2-6)),OSvsAC(2:2:(102-numSkips*4-6-6)));
legend([OSvOIyo,OSvVCyo,OSvACyo],'Optimal Integration','Visual Capture','Auditory Capture');        




