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
getArgs(varargin,{'dispFit=1','combineData=1'});
 
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
 
% calculate/formalize stuff we use for analysis later
[aPSE,aSTD,vPSE,vSTD,bPSE,bSTD,delta,widthArray,deltaArray,visID,audID,bID,aErr,vErr,bErr] = barData(e); %easy-access parameters, summary statistics+their error, and d file specification
[audR2Percent,visR2Percent,bR2Percent] = goodnessOfFit(e,audID,visID,bID) %curve confidence measures of fit
%add curve confidence as d file quality
e.d{audID}.fit.percent = audR2Percent
for i = 1:length(e.d{visID}.visualWidth)
e.d{visID}.fit(i).percent = visR2Percent(i)
end
for i = 1:length(e.d{bID}.condWidth)
e.d{bID}.fit(i).percent = bR2Percent(i)
end
 
% display the fits
if dispFit
  % plot each psychometric function
  for iFile = find(e.isPsycho)
    %diplay fit
    dispFits(e.d{iFile},iFile);
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
 
%bar graph of model comparison
for whichWidth = 1:length(vPSE)
    [mle, suboptimal, SCS, scsSubVisW, scsSubAudW] = calcModelThresholds(whichWidth,vSTD,vPSE,aSTD,aPSE,bPSE,delta)
    % bar graph of different model predictions
    graphModelThresholds(aSTD,vSTD,mle,SCS,suboptimal,bSTD,widthArray,whichWidth,visID,audID,bID,e,aErr,vErr,bErr)
end
 
% graph observered PSEs
figure(4+length(vPSE))
nplots = 3
graphPSEs(bPSE,vPSE,deltaArray,widthArray,bErr)
hold on
%grap predicted PSEs
graphPredictedPSE(scsSubVisW,scsSubAudW,aSTD,vSTD,deltaArray,widthArray,vPSE,aPSE,bPSE,delta)
graphNoShiftPredictedPSE(aSTD,vSTD,deltaArray,widthArray,vPSE,aPSE,bPSE,bErr)
 
end
 
 
 
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
titleStr = d.stimulusType;
 
for iCond = 1:length(whichConds)
  % plot fit
  plot(d.fit(whichConds(iCond)).fitX,d.fit(whichConds(iCond)).fitY*100,'-','Color',fitColors{iCond});hold on
  
  % plot in percentile
  myerrorbar(d.cond(whichConds(iCond)).uniquePosDiff,100*d.cond(whichConds(iCond)).correctBinned,'yError',100*d.cond(whichConds(iCond)).correctBinnedError,'Symbol','o','MarkerFaceColor',dataColors{iCond});
  xlabel('Probe Offset (visual degrees)');
  yaxis(0,100);
  ylabel('Frequency Identified Rightwards (%)');
  % append fit parameters to title
  if d.stimulusType(1) == 'V'
      titleStr = sprintf('%s\nWidth: %d; Mean: %0.2f; Std: %0.2f; P: %g',titleStr,d.visualWidth(whichConds(iCond)),d.fit(whichConds(iCond)).mean,d.fit(whichConds(iCond)).std,1-d.fit(whichConds(iCond)).percent);
  end
  if d.stimulusType(1) == 'A'
      titleStr = sprintf('%s\nMean: %0.2f; Std: %0.2f; P: %g',titleStr,d.fit(whichConds(iCond)).mean,d.fit(whichConds(iCond)).std,1-d.fit(whichConds(iCond)).percent);
  end
  if d.stimulusType(1) == 'B'
      titleStr = sprintf('%s\nWidth: %d; Displacement: %d; Mean: %0.2f; Std: %0.2f; P: %g',titleStr,d.visualWidth(whichConds(ceil(iCond/length(d.displacement)))),d.condDisplacement(whichConds(iCond)),d.fit(whichConds(iCond)).mean,d.fit(whichConds(iCond)).std,1-d.fit(whichConds(iCond)).percent);
  end
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
if isempty(strfind(lower(taskFilename),'alaisburr')) & isempty(strfind(lower(taskFilename),'estimation'))
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
  d.experimentName = sprintf('%s: [%s] Width: %s',d.stimulusType,mlrnum2str(d.displacement),mlrnum2str(d.visualWidth));
  % set number of conditions
  d.nCond = length(d.displacement) * length(d.visualWidth);
  % get trials for each condition
  iCond = 1;
  for iDisplacementCond = 1:length(d.displacement)
    for iWidthCond = 1:length(d.visualWidth)
      % set condition names
      d.condNames{iCond} = sprintf('Displacement: %s Width: %s',mlrnum2str(d.displacement(iDisplacementCond)),mlrnum2str(d.visualWidth(iWidthCond)));
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
    d.condNames{iCond} = sprintf('Width: %s',mlrnum2str(d.visualWidth(iCond)));
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
 
 
%%%%%%%%%%%%%%%%%%%%
%      barData     %
%%%%%%%%%%%%%%%%%%%%
function [aPSE,aSTD,vPSE,vSTD,bPSE,bSTD,delta,widthArray,deltaArray,visID,audID,bID,aErr,vErr,bErr] = barData(e)
 
% extract data from every d (type) file
for iFile = 1:e.nFiles
    
    % grab auditory PSE/STD
    if strcmp(e.d{iFile}.stimulusType,'Auditory')
        aPSE = e.d{iFile}.fit.mean;
        aSTD = e.d{iFile}.fit.std;
        aErr = e.d{iFile}.fit.covar;
        % save id of auditory file for later
        audID = iFile
    end
    
    % grab visual PSEs/STDs
    if strcmp(e.d{iFile}.stimulusType,'Visual')
        % returns array of values ordered by ascending visual width
        vPSE = [e.d{iFile}.fit.mean];
        vSTD = [e.d{iFile}.fit.std];
        vErr = [e.d{iFile}.fit.covar];
        % save array of widths for later usage
        widthArray = e.d{iFile}.visualWidth
        % save id of visual file for later
        visID = iFile
    end
    
    % grab bimodal PSEs/STDs
    if strcmp(e.d{iFile}.stimulusType,'Bimodal')
        % returns array ordered by ascending visual width within ascending ordering by delta
        bPSE = [e.d{iFile}.fit.mean];
        bSTD = [e.d{iFile}.fit.std];
        bErr = [e.d{iFile}.fit.covar];
        % grab delta value and array (assumes only 3 discrepansies)
        delta = max(e.d{iFile}.displacement)
        deltaArray = [e.d{iFile}.displacement]
        % save id of bimodal file for later
        bID = iFile
    end
end
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     calcModelThresholds    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mle,suboptimal,SCS,scsSubVisW,scsSubAudW] = calcModelThresholds(whichWidth,vSTD,vPSE,aSTD,aPSE,bPSE,delta)
 
% calculate the MLE threshold 
mle = sqrt((vSTD(whichWidth)*vSTD(whichWidth)*aSTD*aSTD)/(vSTD(whichWidth)*vSTD(whichWidth)+aSTD*aSTD))
 
% calculate weights for SCS/suboptimal models using conflict condition
scsSubVisW = (abs(((aPSE-5)-bPSE(whichWidth))/10)+abs(((aPSE+5)-bPSE(whichWidth+2*length(vPSE)))/10))/2
scsSubAudW = 1 - scsSubVisW
 
% suboptimal threshold
suboptimal = sqrt(scsSubVisW*vSTD(whichWidth)^2+(scsSubAudW*aSTD)^2)
 
% SCS threshold
SCS = sqrt((scsSubVisW*(vSTD(whichWidth)*vSTD(whichWidth)))+(scsSubAudW*(aSTD*aSTD))+scsSubVisW*scsSubAudW*(vPSE(whichWidth)-aPSE)^2)
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   graphModelThresholds   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graph thresholds for unimodals, each bimodal conflict condition, and model-predicted thresholds
function graphModelThresholds(aSTD,vSTD,mle,SCS,suboptimal,bSTD,widthArray,whichWidth,visID,bID,audID,e,aErr,vErr,bErr)
 
% formalize bar graph inputs
avgbSTD = (bSTD(whichWidth)+bSTD(whichWidth+length(vSTD))+bSTD(whichWidth+length(vSTD)+length(vSTD)))/3
graphStats = [vSTD(whichWidth) aSTD avgbSTD mle suboptimal SCS];
graphWidth = widthArray(whichWidth);
 
% initiate figure
figure(whichWidth+3);
nplots = 4
conditionString = 'Visual Width Condition: %g'
conditionTitle = sprintf(conditionString,graphWidth)
title(conditionTitle)
subplot(1,4,1)
 
% call bar graph
bar(graphStats);
hold on
 
% show error bars, taken from covar matrix (bimodal averaged)
avgBerr = sqrt(((bErr(2,whichWidth*2))^2 + (bErr(2,(2*length(vSTD)+whichWidth*2)))^2 + (bErr(2,(3*length(vSTD)+whichWidth*2)))^2)/3)
errorbar(graphStats,[vErr(2,(whichWidth*2)) aErr(2,2) avgBerr 0 0 0],'.')
 
% calculate p values for model predictions
[mleH mleP] = ztest(mle,avgbSTD,avgBerr,'tail','both')
[subH subP] = ztest(suboptimal,avgbSTD,avgBerr,'tail','both')
[scsH scsP] = ztest(SCS,avgbSTD,avgBerr,'tail','both')
format shortE
text(3.55,mle+.2,[num2str(mleP)],'fontsize',7)
text(4.55,suboptimal+.2,[num2str(subP)],'fontsize',7)
text(5.55,SCS+.2,[num2str(scsP)],'fontsize',7)
 
% label
set(gca,'XTickLabel',{'V','A','Bimodal','OI','SI','SCS'})
ylabel('Threshold (Visual Degrees)')
xlabel('Observed Data/Model Predictions')
intString = 'Threshold Comparison, Width: %g'
graphTitle = sprintf(intString,graphWidth)
title(graphTitle)
hold off
 
% display visual curve
subplot(1,4,2)
dispPsychometricFunction(e.d{visID},[whichWidth])
 
% display auditory curve (indexed as bimodal for some reason)
subplot(1,4,3)
dispPsychometricFunction(e.d{bID},[1])
 
% display bimodal curves (indexed as bimodal for some reason)
subplot(1,4,4)
dispPsychometricFunction(e.d{audID},[(whichWidth) (whichWidth+length(vSTD)) (whichWidth+2*length(vSTD))])
 
 
%%%%%%%%%%%%%%%%%%%%%%%%
%     graphPSEs (Observed)   %
%%%%%%%%%%%%%%%%%%%%%%%%
function graphPSEs(bPSE,vPSE,deltaArray,widthArray,bErr)
 
% plot scatterplots and best fit lines
subplot(1,3,1)
y1 = [bPSE(1) bPSE(1+length(vPSE)) bPSE(1+2*length(vPSE))];
s1 = scatter(deltaArray-.5,y1,200,getSmoothColor(2,10,hot),'filled')
hold on
j = errorbar(deltaArray-.5,y1,[bErr(1,1) bErr(1,1+2*length(vPSE)) bErr(1,1+4*length(vPSE))],'.')
j.Color = 'black'
hold on
p = polyval(polyfit(deltaArray,y1,1),deltaArray);
plot(deltaArray-.5,p,'color',getSmoothColor(2,10,hot))
hold on
intString = string('width: %g; slope: %g')
t1 = sprintf(intString,widthArray(1),(p(3)-p(1))/2)
 
% adding lines for each width condition (only works if 5> different width conditions)
if length(vPSE) > 1
    y2 = [bPSE(2) bPSE(2+length(vPSE)) bPSE(2+2*length(vPSE))]
    s2 = scatter(deltaArray-.166,y2,200,getSmoothColor(4,10,hot),'filled')
    hold on
    j = errorbar(deltaArray-.166,y2,[bErr(1,3) bErr(1,3+2*length(vPSE)) bErr(1,3+4*length(vPSE))],'.')
    j.Color = 'black'
    hold on
    p2 = polyval(polyfit(deltaArray,y2,1),deltaArray);
    plot(deltaArray-.166,p2,'color',getSmoothColor(4,10,hot))
    hold on
    intString = string('width: %g; slope: %g')
    t2 = sprintf(intString,widthArray(2),(p2(3)-p2(1))/2)
end
if length(vPSE) > 2
    y3 = [bPSE(3) bPSE(3+length(vPSE)) bPSE(3+2*length(vPSE))]
    s3 = scatter(deltaArray+.166,y3,200,getSmoothColor(6,10,hot),'filled')
    hold on
    j = errorbar(deltaArray+.166,y3,[bErr(1,5) bErr(1,5+2*length(vPSE)) bErr(1,5+4*length(vPSE))],'.')
    j.Color = 'black'
    hold on
    p3 = polyval(polyfit(deltaArray,y3,1),deltaArray);
    plot(deltaArray+.166,p3,'color',getSmoothColor(6,10,hot))
    hold on
    intString = string('width: %g; slope: %g')
    t3 = sprintf(intString,widthArray(3),(p3(3)-p3(1))/2)
end
if length(vPSE) > 3
    y4 = [bPSE(4) bPSE(4+length(vPSE)) bPSE(4+2*length(vPSE))]
    s4 = scatter(deltaArray+.5,y4,200,getSmoothColor(7,10,hot),'filled')
    hold on
    j = errorbar(deltaArray+.5,y4,[bErr(1,7) bErr(1,7+2*length(vPSE)) bErr(1,7+4*length(vPSE))],'.')
    j.Color = 'black'
    hold on
    p4 = polyval(polyfit(deltaArray,y4,1),deltaArray);
    plot(deltaArray+.5,p4,'color',getSmoothColor(7,10,hot))
    intString = string('width: %g; slope: %g')
    t4 = sprintf(intString,widthArray(4),(p4(3)-p4(1))/2)
end
 
% plot expected PSEs (hard-coded for now)
plot([-5 0 5],[5 0 -5],['--','black'])
plot([-5 0 5],[-5 0 5],['--','cyan'])
 
% label and title
title('Observed PSEs by width')
xlabel('Audio-Visual discrepancy (delta)')
ylabel('Observed PSE')
 
% create legend (again, only for 5> conditions) (iterated poorly, will clean this up at some point)
if length(vPSE) == 1
legend([s1],t1)
end
if length(vPSE) == 2
legend([s1 s2],t1,t2)
end
if length(vPSE) == 3
legend([s1 s2 s3],t1,t2,t3)
end
if length(vPSE) == 4
legend([s1 s2 s3 s4],t1,t2,t3,t4)
end
hold off
 
 
%%%%%%%%%%%%%%%%%%%%%%%%
%     graphPredictedPSEs   %
%%%%%%%%%%%%%%%%%%%%%%%%
function graphPredictedPSE(scsSubVisW,scsSubAudW,aSTD,vSTD,deltaArray,widthArray,vPSE,aPSE,bPSE,delta)
subplot(1,3,2)
% plot scatterplots and best fit lines
mleVw = (aSTD^2)/((aSTD^2)+vSTD(1)^2)
mleAw = 1 - mleVw
y11 = [(mleVw*(5+vPSE(1))+mleAw*(-5+aPSE)) (mleVw*vPSE(1)+mleAw*aPSE) (mleVw*(-5+vPSE(1))+mleAw*(5+aPSE))];
s11 = scatter(deltaArray-.5,y11,200,getSmoothColor(2,10,hot),'filled')
hold on
p11 = polyval(polyfit(deltaArray-.5,y11,1),deltaArray);
plot(deltaArray-.5,p11,'color',getSmoothColor(2,10,hot))
hold on
intString = string('width: %g; slope: %g; Vweight: %g')
t11 = sprintf(intString,widthArray(1),(p11(3)-p11(1))/2,mleVw)
 
% adding lines for each width condition (only works if 5> different width conditions)
if length(vPSE) > 1
   mleVw = (aSTD^2)/((aSTD^2)+vSTD(2)^2)
   mleAw = 1 - mleVw
   y22 = [(mleVw*(5+vPSE(2))+mleAw*(-5+aPSE)) (mleVw*vPSE(2)+mleAw*aPSE) (mleVw*(-5+vPSE(2))+mleAw*(5+aPSE))];
   s22 = scatter(deltaArray-.166,y22,200,getSmoothColor(4,10,hot),'filled')
   hold on
   p22 = polyval(polyfit(deltaArray,y22,1),deltaArray);
   plot(deltaArray-.166,p22,'color',getSmoothColor(4,10,hot))
   hold on
   intString = string('width: %g; slope: %g; Vweight: %g')
   t22 = sprintf(intString,widthArray(2),(p22(3)-p22(1))/2,mleVw)
   end
if length(vPSE) > 2
   mleVw = (aSTD^2)/((aSTD^2)+vSTD(3)^2)
   mleAw = 1 - mleVw
   y33 = [(mleVw*(5+vPSE(3))+mleAw*(-5+aPSE)) (mleVw*vPSE(3)+mleAw*aPSE) (mleVw*(-5+vPSE(3))+mleAw*(5+aPSE))];
   s33 = scatter(deltaArray+.166,y33,200,getSmoothColor(6,10,hot),'filled')
   hold on
   p33 = polyval(polyfit(deltaArray,y33,1),deltaArray);
   plot(deltaArray+.166,p33,'color',getSmoothColor(6,10,hot))
   hold on
   intString = string('width: %g; slope: %g; Vweight: %g')
   t33 = sprintf(intString,widthArray(3),(p33(3)-p33(1))/2,mleVw)
   end
if length(vPSE) > 3
   mleVw = (aSTD^2)/((aSTD^2)+vSTD(4)^2)
   mleAw = 1 - mleVw
   y44 = [(mleVw*(5+vPSE(4))+mleAw*(-5+aPSE)) (mleVw*vPSE(4)+mleAw*aPSE) (mleVw*(-5+vPSE(4))+mleAw*(5+aPSE))];
   s44 = scatter(deltaArray+.5,y44,200,getSmoothColor(7,10,hot),'filled')
   hold on
   p44 = polyval(polyfit(deltaArray,y44,1),deltaArray);
   plot(deltaArray+.5,p44,'color',getSmoothColor(7,10,hot))
   hold on
   intString = string('width: %g; slope: %g; Vweight: %g')
   t44 = sprintf(intString,widthArray(4),(p44(3)-p44(1))/2,mleVw)
   end
 
% plot expected PSEs (hard-coded for now)
plot([-5 0 5],[5 0 -5],['--','black'])
plot([-5 0 5],[-5 0 5],['--','cyan'])
 
% label and title
title('MLE expected PSEs')
xlabel('Audio-Visual discrepancy (delta)')
 
% create legend (again, only for 5> conditions) (iterated poorly, will clean this up at some point)
if length(vPSE) == 1
legend([s11],t11)
end
if length(vPSE) == 2
legend([s11 s22],t11,t22)
end
if length(vPSE) == 3
legend([s11 s22 s33],t11,t22,t33)
end
if length(vPSE) == 4
legend([s11 s22 s33 s44],t11,t22,t33,t44)
end
hold off
 
 
 
 
 
% graph scs suboptimal
subplot(1,3,3)
 
scsSubVisW = (abs(((aPSE-5)-bPSE(1))/10)+abs(((aPSE+5)-bPSE(1+2*length(vPSE)))/10))/2
scsSubAudW = 1 - scsSubVisW
y111 = [(scsSubVisW*(5+vPSE(1))+scsSubAudW*(-5+aPSE)) (scsSubVisW*vPSE(1)+scsSubAudW*aPSE) (scsSubVisW*(-5+vPSE(1))+scsSubAudW*(5+aPSE))];
s111 = scatter(deltaArray,y111,200,getSmoothColor(2,10,hot),'filled')
hold on
p111 = polyval(polyfit(deltaArray-.5,y111,1),deltaArray);
plot(deltaArray-.5,p111,'color',getSmoothColor(2,10,hot))
hold on
intString = string('width: %g; slope: %g; Vweight: %g')
t111 = sprintf(intString,widthArray(1),(p111(3)-p111(1))/2,scsSubVisW)
 
% adding lines for each width condition (only works if 5> different width conditions)
if length(vPSE) > 1
   scsSubVisW = (abs(((aPSE-5)-bPSE(2))/10)+abs(((aPSE+5)-bPSE(2+2*length(vPSE)))/10))/2
   scsSubAudW = 1 - scsSubVisW
   y222 = [(scsSubVisW*(5+vPSE(2))+scsSubAudW*(-5+aPSE)) (scsSubVisW*vPSE(2)+scsSubAudW*aPSE) (scsSubVisW*(-5+vPSE(2))+scsSubAudW*(5+aPSE))];
   s222 = scatter(deltaArray-.166,y222,200,getSmoothColor(4,10,hot),'filled')
   hold on
   p222 = polyval(polyfit(deltaArray,y222,1),deltaArray);
   plot(deltaArray-.166,p222,'color',getSmoothColor(4,10,hot))
   hold on
   intString = string('width: %g; slope: %g; Vweight: %g')
   t222 = sprintf(intString,widthArray(2),(p222(3)-p222(1))/2,scsSubVisW)
   end
if length(vPSE) > 2
   scsSubVisW = (abs(((aPSE-5)-bPSE(3))/10)+abs(((aPSE+5)-bPSE(3+2*length(vPSE)))/10))/2
   scsSubAudW = 1 - scsSubVisW
   y333 = [(scsSubVisW*(5+vPSE(3))+scsSubAudW*(-5+aPSE)) (scsSubVisW*vPSE(3)+scsSubAudW*aPSE) (scsSubVisW*(-5+vPSE(3))+scsSubAudW*(5+aPSE))];
   s333 = scatter(deltaArray+.166,y333,200,getSmoothColor(6,10,hot),'filled')
   hold on
   p333 = polyval(polyfit(deltaArray,y333,1),deltaArray);
   plot(deltaArray+.166,p333,'color',getSmoothColor(6,10,hot))
   hold on
   intString = string('width: %g; slope: %g; Vweight: %g')
   t333 = sprintf(intString,widthArray(3),(p333(3)-p333(1))/2,scsSubVisW)
   end
if length(vPSE) > 3
   scsSubVisW = (abs(((aPSE-5)-bPSE(4))/10)+abs(((aPSE+5)-bPSE(4+2*length(vPSE)))/10))/2
   scsSubAudW = 1 - scsSubVisW
   y444 = [(scsSubVisW*(5+vPSE(4))+scsSubAudW*(-5+aPSE)) (scsSubVisW*vPSE(4)+scsSubAudW*aPSE) (scsSubVisW*(-5+vPSE(4))+scsSubAudW*(5+aPSE))];
   s444 = scatter(deltaArray+.5,y444,200,getSmoothColor(7,10,hot),'filled')
   hold on
   p444 = polyval(polyfit(deltaArray,y444,1),deltaArray);
   plot(deltaArray+.5,p444,'color',getSmoothColor(7,10,hot))
   hold on
   intString = string('width: %g; slope: %g; Vweight: %g')
   t444 = sprintf(intString,widthArray(4),(p444(3)-p444(1))/2,scsSubVisW)
   end
 
% plot expected PSEs (hard-coded for now)
plot([-5 0 5],[5 0 -5],['--','black'])
plot([-5 0 5],[-5 0 5],['--','cyan'])
 
% label and title
title('SCS/Suboptimal expected PSEs')
xlabel('Audio-Visual discrepancy (delta)')
 
% create legend (again, only for 5> conditions) (iterated poorly, will clean this up at some point)
if length(vPSE) == 1
legend([s111],t111)
end
if length(vPSE) == 2
legend([s111 s222],t111,t222)
end
if length(vPSE) == 3
legend([s111 s222 s333],t111,t222,t333)
end
if length(vPSE) == 4
legend([s111 s222 s333 s444],t111,t222,t333,t444)
end
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%
%     graphNoShiftPredictedPSEs   %
%%%%%%%%%%%%%%%%%%%%%%%%
function graphNoShiftPredictedPSE(aSTD,vSTD,deltaArray,widthArray,vPSE,aPSE,bPSE,bErr)
figure(5+length(vPSE))
nplots = 2
 
 
%noshiftmle
subplot(1,2,1)
mleVw = (aSTD^2)/((aSTD^2)+vSTD(1)^2)
mleAw = 1 - mleVw
y11 = [(mleVw*(5)+mleAw*(-5)) 0 (mleVw*(-5)+mleAw*(5))];
s11 = scatter(deltaArray-.5,y11,200,getSmoothColor(2,10,hot),'filled')
hold on
p11 = polyval(polyfit(deltaArray,y11,1),deltaArray);
plot(deltaArray-.5,p11,'color',getSmoothColor(2,10,hot))
hold on
intString = string('width: %g; slope: %g; Vweight: %g')
t11 = sprintf(intString,widthArray(1),(p11(3)-p11(1))/2,mleVw)
 
% adding lines for each width condition (only works if 5> different width conditions)
if length(vPSE) > 1
   mleVw = (aSTD^2)/((aSTD^2)+vSTD(2)^2)
   mleAw = 1 - mleVw
   y22 = [(mleVw*(5)+mleAw*(-5)) 0 (mleVw*(-5)+mleAw*(5))];
   s22 = scatter(deltaArray-.166,y22,200,getSmoothColor(4,10,hot),'filled')
   hold on
   p22 = polyval(polyfit(deltaArray,y22,1),deltaArray);
   plot(deltaArray-.166,p22,'color',getSmoothColor(4,10,hot))
   hold on
   intString = string('width: %g; slope: %g; Vweight: %g')
   t22 = sprintf(intString,widthArray(2),(p22(3)-p22(1))/2,mleVw)
   end
if length(vPSE) > 2
   mleVw = (aSTD^2)/((aSTD^2)+vSTD(3)^2)
   mleAw = 1 - mleVw
   y33 = [(mleVw*(5)+mleAw*(-5)) 0 (mleVw*(-5)+mleAw*(5))];
   s33 = scatter(deltaArray+.166,y33,200,getSmoothColor(6,10,hot),'filled')
   hold on
   p33 = polyval(polyfit(deltaArray,y33,1),deltaArray);
   plot(deltaArray+.166,p33,'color',getSmoothColor(6,10,hot))
   hold on
   intString = string('width: %g; slope: %g; Vweight: %g')
   t33 = sprintf(intString,widthArray(3),(p33(3)-p33(1))/2,mleVw)
   end
if length(vPSE) > 3
   mleVw = (aSTD^2)/((aSTD^2)+vSTD(4)^2)
   mleAw = 1 - mleVw
   y44 = [(mleVw*(5)+mleAw*(-5)) 0 (mleVw*(-5)+mleAw*(5))];
   s44 = scatter(deltaArray+.5,y44,200,getSmoothColor(7,10,hot),'filled')
   hold on
   p44 = polyval(polyfit(deltaArray,y44,1),deltaArray);
   plot(deltaArray+.5,p44,'color',getSmoothColor(7,10,hot))
   hold on
   intString = string('width: %g; slope: %g; Vweight: %g')
   t44 = sprintf(intString,widthArray(4),(p44(3)-p44(1))/2,mleVw)
   end
 
% plot expected PSEs (hard-coded for now)
plot([-5 0 5],[5 0 -5],['--','black'])
plot([-5 0 5],[-5 0 5],['--','cyan'])
 
% label and title
title('MLE expected PSEs, no shift')
xlabel('Audio-Visual discrepancy (delta)')
 
% create legend (again, only for 5> conditions) (iterated poorly, will clean this up at some point)
if length(vPSE) == 1
legend([s11],t11)
end
if length(vPSE) == 2
legend([s11 s22],t11,t22)
end
if length(vPSE) == 3
legend([s11 s22 s33],t11,t22,t33)
end
if length(vPSE) == 4
legend([s11 s22 s33 s44],t11,t22,t33,t44)
end
hold off
 
 
 
 
%noshiftscs
subplot(1,2,2)
 
scsSubVisW = (abs(((aPSE-5)-bPSE(1))/10)+abs(((aPSE+5)-bPSE(1+2*length(vPSE)))/10))/2
scsSubAudW = 1 - scsSubVisW
y111 = [(scsSubVisW*(5)+scsSubAudW*(-5)) 0 scsSubVisW*(-5)+scsSubAudW*(5)];
s111 = scatter(deltaArray-.5,y111,200,getSmoothColor(2,10,hot),'filled')
hold on
p111 = polyval(polyfit(deltaArray-.5,y111,1),deltaArray);
plot(deltaArray,p111,'color',getSmoothColor(2,10,hot))
hold on
intString = string('width: %g; slope: %g; Vweight: %g')
t111 = sprintf(intString,widthArray(1),(p111(3)-p111(1))/2,scsSubVisW)
 
% adding lines for each width condition (only works if 5> different width conditions)
if length(vPSE) > 1
   scsSubVisW = (abs(((aPSE-5)-bPSE(2))/10)+abs(((aPSE+5)-bPSE(2+2*length(vPSE)))/10))/2
   scsSubAudW = 1 - scsSubVisW
   y222 = [(scsSubVisW*(5)+scsSubAudW*(-5)) 0 (scsSubVisW*(-5)+scsSubAudW*(5))];
   s222 = scatter(deltaArray-.166,y222,200,getSmoothColor(4,10,hot),'filled')
   hold on
   p222 = polyval(polyfit(deltaArray,y222,1),deltaArray);
   plot(deltaArray-.166,p222,'color',getSmoothColor(4,10,hot))
   hold on
   intString = string('width: %g; slope: %g; Vweight: %g')
   t222 = sprintf(intString,widthArray(2),(p222(3)-p222(1))/2,scsSubVisW)
   end
if length(vPSE) > 2
   scsSubVisW = (abs(((aPSE-5)-bPSE(3))/10)+abs(((aPSE+5)-bPSE(3+2*length(vPSE)))/10))/2
   scsSubAudW = 1 - scsSubVisW
   y333 = [(scsSubVisW*(5)+scsSubAudW*(-5)) 0 (scsSubVisW*(-5)+scsSubAudW*(5))];
   s333 = scatter(deltaArray+.166,y333,200,getSmoothColor(6,10,hot),'filled')
   hold on
   p333 = polyval(polyfit(deltaArray,y333,1),deltaArray);
   plot(deltaArray+.166,p333,'color',getSmoothColor(6,10,hot))
   hold on
   intString = string('width: %g; slope: %g; Vweight: %g')
   t333 = sprintf(intString,widthArray(3),(p333(3)-p333(1))/2,scsSubVisW)
   end
if length(vPSE) > 3
   scsSubVisW = (abs(((aPSE-5)-bPSE(4))/10)+abs(((aPSE+5)-bPSE(4+2*length(vPSE)))/10))/2
   scsSubAudW = 1 - scsSubVisW
   y444 = [(scsSubVisW*(5)+scsSubAudW*(-5)) 0 (scsSubVisW*(-5)+scsSubAudW*(5))];
   s444 = scatter(deltaArray+.5,y444,200,getSmoothColor(7,10,hot),'filled')
   hold on
   p444 = polyval(polyfit(deltaArray,y444,1),deltaArray);
   plot(deltaArray+.5,p444,'color',getSmoothColor(7,10,hot))
   hold on
   intString = string('width: %g; slope: %g; Vweight: %g')
   t444 = sprintf(intString,widthArray(4),(p444(3)-p444(1))/2,scsSubVisW)
   end
 
% plot expected PSEs (hard-coded for now)
plot([-5 0 5],[5 0 -5],['--','black'])
plot([-5 0 5],[-5 0 5],['--','cyan'])
 
% label and title
title('SCS/Suboptimal expected PSEs, no shift')
xlabel('Audio-Visual discrepancy (delta)')
 
% create legend (again, only for 5> conditions) (iterated poorly, will clean this up at some point)
if length(vPSE) == 1
legend([s111],t111)
end
if length(vPSE) == 2
legend([s111 s222],t111,t222)
end
if length(vPSE) == 3
legend([s111 s222 s333],t111,t222,t333)
end
if length(vPSE) == 4
legend([s111 s222 s333 s444],t111,t222,t333,t444)
end
 
%alais burr figure
figure(6+length(vPSE))
nplots = 2
subplot(1,2,1)
y1 = [bPSE(1) bPSE(1+length(vPSE)) bPSE(1+2*length(vPSE))];
s1 = scatter(deltaArray-.5,y1,200,getSmoothColor(2,10,hot),'filled')
hold on
j = errorbar(deltaArray-.5,y1,[bErr(1,1) bErr(1,1+2*length(vPSE)) bErr(1,1+4*length(vPSE))],'.')
j.Color = 'black'
hold on
plot(deltaArray-.5,p11,'color',getSmoothColor(2,10,hot))
hold on
intString = string('width: %g; slope: %g')
t1 = sprintf(intString,widthArray(1),(p11(3)-p11(1))/2)
 
% adding lines for each width condition (only works if 5> different width conditions)
if length(vPSE) > 1
    y2 = [bPSE(2) bPSE(2+length(vPSE)) bPSE(2+2*length(vPSE))]
    s2 = scatter(deltaArray-.166,y2,200,getSmoothColor(4,10,hot),'filled')
    hold on
    j = errorbar(deltaArray-.166,y2,[bErr(1,3) bErr(1,3+2*length(vPSE)) bErr(1,3+4*length(vPSE))],'.')
    j.Color = 'black'
    hold on
    plot(deltaArray-.166,p22,'color',getSmoothColor(4,10,hot))
    hold on
    intString = string('width: %g; slope: %g')
    t2 = sprintf(intString,widthArray(2),(p22(3)-p22(1))/2)
end
if length(vPSE) > 2
    y3 = [bPSE(3) bPSE(3+length(vPSE)) bPSE(3+2*length(vPSE))]
    s3 = scatter(deltaArray+.166,y3,200,getSmoothColor(6,10,hot),'filled')
    hold on
    j = errorbar(deltaArray+.166,y3,[bErr(1,5) bErr(1,5+2*length(vPSE)) bErr(1,5+4*length(vPSE))],'.')
    j.Color = 'black'
    hold on
    plot(deltaArray+.166,p33,'color',getSmoothColor(6,10,hot))
    hold on
    intString = string('width: %g; slope: %g')
    t3 = sprintf(intString,widthArray(3),(p33(3)-p33(1))/2)
end
if length(vPSE) > 3
    y4 = [bPSE(4) bPSE(4+length(vPSE)) bPSE(4+2*length(vPSE))]
    s4 = scatter(deltaArray+.5,y4,200,getSmoothColor(7,10,hot),'filled')
    hold on
    j = errorbar(deltaArray+.5,y4,[bErr(1,7) bErr(1,7+2*length(vPSE)) bErr(1,7+4*length(vPSE))],'.')
    j.Color = 'black'
    hold on
    plot(deltaArray+.5,p44,'color',getSmoothColor(7,10,hot))
    intString = string('width: %g; slope: %g')
    t4 = sprintf(intString,widthArray(4),(p44(3)-p44(1))/2)
end
 
% plot expected PSEs (hard-coded for now)
plot([-5 0 5],[5 0 -5],['--','black'])
plot([-5 0 5],[-5 0 5],['--','cyan'])
 
% label and title
title('Observed PSEs with Optimal Integration Prediction')
xlabel('Audio-Visual discrepancy (delta)')
ylabel('Observed PSE')
 
% create legend (again, only for 5> conditions) (iterated poorly, will clean this up at some point)
if length(vPSE) == 1
legend([s1],t1)
end
if length(vPSE) == 2
legend([s1 s2],t1,t2)
end
if length(vPSE) == 3
legend([s1 s2 s3],t1,t2,t3)
end
if length(vPSE) == 4
legend([s1 s2 s3 s4],t1,t2,t3,t4)
end
hold off
 
 
%alais burr figure but scs predicts
subplot(1,2,2)
y1 = [bPSE(1) bPSE(1+length(vPSE)) bPSE(1+2*length(vPSE))];
s1 = scatter(deltaArray-.5,y1,200,getSmoothColor(2,10,hot),'filled');
hold on
j = errorbar(deltaArray-.5,y1,[bErr(1,1) bErr(1,1+2*length(vPSE)) bErr(1,1+4*length(vPSE))],'.');
j.Color = 'black'
hold on
plot(deltaArray,p111,'color',getSmoothColor(2,10,hot));
hold on
intString = string('width: %g; slope: %g');
t1 = sprintf(intString,widthArray(1),(p111(3)-p111(1))/2);
 
% adding lines for each width condition (only works if 5> different width conditions)
if length(vPSE) > 1
    y2 = [bPSE(2) bPSE(2+length(vPSE)) bPSE(2+2*length(vPSE))];
    s2 = scatter(deltaArray-.166,y2,200,getSmoothColor(4,10,hot),'filled');
    hold on
    j = errorbar(deltaArray-.166,y2,[bErr(1,3) bErr(1,3+2*length(vPSE)) bErr(1,3+4*length(vPSE))],'.');
    j.Color = 'black'
    hold on
    plot(deltaArray-.166,p222,'color',getSmoothColor(4,10,hot));
    hold on
    intString = string('width: %g; slope: %g')
    t2 = sprintf(intString,widthArray(2),(p222(3)-p222(1))/2);
end
if length(vPSE) > 2
    y3 = [bPSE(3) bPSE(3+length(vPSE)) bPSE(3+2*length(vPSE))];
    s3 = scatter(deltaArray+.166,y3,200,getSmoothColor(6,10,hot),'filled');
    hold on
    j = errorbar(deltaArray+.166,y3,[bErr(1,5) bErr(1,5+2*length(vPSE)) bErr(1,5+4*length(vPSE))],'.')
    j.Color = 'black'
    hold on
    plot(deltaArray+.166,p333,'color',getSmoothColor(6,10,hot));
    hold on
    intString = string('width: %g; slope: %g');
    t3 = sprintf(intString,widthArray(3),(p333(3)-p333(1))/2);
end
if length(vPSE) > 3
    y4 = [bPSE(4) bPSE(4+length(vPSE)) bPSE(4+2*length(vPSE))];
    s4 = scatter(deltaArray+.5,y4,200,getSmoothColor(7,10,hot),'filled');
    hold on
    j = errorbar(deltaArray+.5,y4,[bErr(1,7) bErr(1,7+2*length(vPSE)) bErr(1,7+4*length(vPSE))],'.')
    j.Color = 'black'
    hold on
    plot(deltaArray+.5,p444,'color',getSmoothColor(7,10,hot));
    intString = string('width: %g; slope: %g');
    t4 = sprintf(intString,widthArray(4),(p444(3)-p444(1))/2);
end
 
% plot expected PSEs (hard-coded for now)
plot([-5 0 5],[5 0 -5],['--','black'])
plot([-5 0 5],[-5 0 5],['--','cyan'])
 
% label and title
title('Observed PSEs with SCS & Suboptimal Integration Predictions')
xlabel('Audio-Visual discrepancy (delta)')
ylabel('Observed PSE')
 
% create legend (again, only for 5> conditions) (iterated poorly, will clean this up at some point)
if length(vPSE) == 1
legend([s1],t1)
end
if length(vPSE) == 2
legend([s1 s2],t1,t2)
end
if length(vPSE) == 3
legend([s1 s2 s3],t1,t2,t3)
end
if length(vPSE) == 4
legend([s1 s2 s3 s4],t1,t2,t3,t4)
end
hold off
 
 
%%%%%%%%%%%%%%%%%%%%%%%
% Curve goodnessOfFit %
%%%%%%%%%%%%%%%%%%%%%%%
function [audR2Percent,visR2Percent,bR2Percent] = goodnessOfFit(e,audID,visID,bID)
perms = 5
 
%auditory
r2RandAUD = []
r2Aud = e.d{audID}.fit.r2;
for i = 1:perms
    %randomize correct binned
    cbRand = randsample(e.d{audID}.cond.correctBinned, length(e.d{audID}.cond.correctBinned), false);
    %fit the gaussian and grab r2 value
    randFits = fitCumulativeGaussian(e.d{audID}.cond.uniquePosDiff,cbRand);
    r2RandAUD(i) = randFits.r2;
end
figure(13);
histogram(r2RandAUD);
audR2Percent = sum(r2RandAUD < r2Aud)/perms;
str = sprintf('Auditory, r2 = %g, percentile = %g',r2Aud,audR2Percent);
title(str);
 
%visual
r2RandVIS = ones(length(e.d{visID}.visualWidth),perms);
r2VIS = [];
for k = 1:length(e.d{visID}.visualWidth)
    r2VIS(k) = e.d{visID}.fit(k).r2;
    for i = 1:perms
        cbRand = randsample(e.d{visID}.cond(k).correctBinned, length(e.d{visID}.cond(k).correctBinned), false);
        randFits = fitCumulativeGaussian(e.d{visID}.cond(k).uniquePosDiff,cbRand);
        r2RandVIS(k,i) = randFits.r2;
    end
figure(13+k)
histogram(r2RandVIS(k,1:perms));
visR2Percent(k) = sum(r2RandVIS(k,1:perms) < r2VIS(k))/perms;
str = sprintf('Visual width: %g, r2 = %g, percentile = %g',e.d{visID}.visualWidth(k),r2VIS(k),visR2Percent(k));
title(str);
end
 
%bimodal
r2RandB = ones(length(e.d{bID}.condWidth),perms);
r2B = [];
for k = 1:length(e.d{bID}.condWidth);
    r2B(k) = e.d{bID}.fit(k).r2;
    for i = 1:perms;
        cbRand = randsample(e.d{bID}.cond(k).correctBinned, length(e.d{bID}.cond(k).correctBinned), false);
        randFits = fitCumulativeGaussian(e.d{bID}.cond(k).uniquePosDiff, cbRand);
        r2RandB(k,i) = randFits.r2;
    end
figure (17+k)
histogram(r2RandB(k,1:perms));
bR2Percent(k) = sum(r2RandB(k,1:perms) < r2B(k))/perms;
str = sprintf('Bimodal width: %g, discrepancy: %g, r2 = %g, percentile = %g',e.d{bID}.condWidth(k),e.d{bID}.condDisplacement(k),r2B(k),bR2Percent(k))
title(str);
end

