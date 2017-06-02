% motionEnergyModelTest.m
%
%      usage: motionEnergyModelTest()
%         by: justin gardner
%       date: 06/01/17
%    purpose: 
%
function retval = motionEnergyModelTest(varargin)

getArgs(varargin,{'recompute=0','dataDir=~/Google Drive/motionEnergy','coherence=[0:0.1:1]','direction',[0:45:359],'n=100'});

% name of files
filename = sprintf('co%idir%in%i',length(coherence),length(direction),n);
stimulusFileName = fullfile(dataDir,'stimuli',filename);
responseFileName = fullfile(dataDir,'response',filename);

if ~recompute && isfile(setext(stimulusFileName,'mat'))
  % load precomputed stimulus file
  disp(sprintf('(motionEnergyModelTest) Loading stimulus: %s',stimulusFileName));
  load(stimulusFileName);
else
  % compute stimulus file
  [s msc] = motionEnergyModelMakeStimulus('screenName=offscreen','coherence',coherence,'direction',direction,'n',n);
  save(stimulusFileName,'s','msc');
  dispHeader(sprintf('(motionEnergyModelTest) Saving stimulus: %s',stimulusFileName));
end

if ~recompute && isfile(setext(responseFileName,'mat'))
  disp(sprintf('(motionEnergyModelTest) Loading response: %s',responseFileName));
  load(responseFileName);
else
  % compute response
  m = motionEnergyModel(s,'myscreen',msc,'dispFigures=0','removeFilters=1');
  % save
  disp(sprintf('(motionEnergyModelTest) Saving response: %s',responseFileName));
  save(responseFileName,'m');
end

keyboard

