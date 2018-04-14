% motionEnergyModelTest.m
%
%      usage: motionEnergyModelTest()
%         by: justin gardner
%       date: 06/01/17
%    purpose: 
%
function retval = motionEnergyModelTest(varargin)

getArgs(varargin,{'recompute=0','dataDir=~/Desktop/motionEnergy','coherence=[0.06 0.12 0.24]','direction',5,'n=1'});

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

% cycle over direction and coherence, collecting simulations
meanResponse = [];stdResponse = [];
mlrSmartfig('motionEnergyModelTest');
for iDirection = 1:length(direction)
  for iCoherence = 1:length(coherence)
    % reset tresponse
    r = [];iResponse = 1;
    % find all s that match the current direction and coherence
    for iStimulus = 1:length(s)
      if isequal(s{iStimulus}.dir,direction(iDirection)) && isequal(s{iStimulus}.coherence,coherence(iCoherence));
	r(iResponse,:,:,:) = m.r{iStimulus}.meanResponse;
	iResponse = iResponse+1;
      end
    end
    % compute mean and ste
    for iTF = 1:m.nTF
      for iSF = 1:m.nSF
	% compute mean and standard error
	meanResponse = squeeze(mean(r(:,iTF,iSF,:),1));
	stdResponse = squeeze(std(r(:,iTF,iSF,:),1));
	myerrorbar(m.orientationPreference,meanResponse,'yError',stdResponse,'polarPlot=1','Color',getSmoothColor(length(coherence)-iCoherence+1,length(coherence)));
	hold on
      end
    end
  end
end

  
keyboard

