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

n = 100;
coherence = 0:0.2:1;
direction = [0 90 180 270];

for i = 1:n
  for iCoherence = 1:length(coherence)
    for iDirection = 1:length(direction)
      [s msc] = motionEnergyModelMakeStimulus('screenName','offscreen','coherence',coherence(iCoherence),'direction',direction(iDirection));
      m(iCohernece,iDirection,i) = motionEnergyModel(s,'myscreen',myscreen);
    end
  end
end

keyboard