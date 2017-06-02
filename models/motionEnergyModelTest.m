% motionEnergyModelTest.m
%
%      usage: motionEnergyModelTest()
%         by: justin gardner
%       date: 06/01/17
%    purpose: 
%
function retval = motionEnergyModelTest()

% check arguments
if ~any(nargin == [0])
  help motionEnergyModelTest
  return
end


%[s msc] = motionEnergyModelMakeStimulus('screenName=offscreen','coherence=[0:0.1:1]','direction=0:45:359','n=10');
[s msc] = motionEnergyModelMakeStimulus('screenName=offscreen','coherence=1','direction=180','n=1');
save '~/Google Drive/motionEnergy/motionEnergyStimulus' s msc


for iStim = 1:length(s)
  dispHeader(sprintf('(motionEnergyModelTest) Computing model responses for coherence: %0.2f and direction: %0.2f n=%i',s{iStim}.coherence,s{iStim}.dir,s{iStim}.n));
  m(iStim) = motionEnergyModel(s{iStim}.s,'myscreen',msc,'dispFigures=0');
  save '~/Google Drive/motionEnergy/motionEnergyResponse' m s msc
end


keyboard