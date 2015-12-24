% sdtPlot.m
%
%        $Id:$ 
%      usage: sdtPlot()
%         by: justin gardner
%       date: 05/01/15
%    purpose: Make signal detection plots
%
%             Default is to draw signal and noise with means 1 and -1, respectively and standard deviation 1
%             Criterion at 0
%
%             To draw with different means and standard deviations, you can set like the following
%             sdtPlot('signalMean=3','signalSTD=4','noiseMean=0','noiseSTD=2','criterion=1.5');
%
function retval = sdtPlot(varargin)

% parse arguments
getArgs(varargin,{'signalMean=0.5','noiseMean=-0.5','signalSTD=1','noiseSTD=1','criterion=0','lineWidth=2','signalColor','k','noiseColor','r','hitAlpha=0.8','missAlpha=0.2','correctRejectAlpha=0.2','falseAlarmAlpha=0.8','plotROC=[]','pHits=[]','pFalseAlarms=[]','titleStr=[]','p2AFC=[]','dispROC=1','dispSignal=1','dispNoise=1','lowerCriterion=[]','dispCriterion=1','dispTitle=1'});

% if we are given either pHits or pFalseAlarms then check that we have both
if ~isempty(pHits) || ~isempty(pFalseAlarms)
  % compute the needed signalMean and noiseMean from pHits and pFalseAlarms
  [signalMean noiseMean titleStr] = useHitsFalseAlarms(pHits,pFalseAlarms,signalSTD,noiseSTD,criterion);
  if isempty(signalMean),return,end
end

% if we are given 2AFC then compute necessary signalMean and noiseMean to generate that
if ~isempty(p2AFC)
  % compute the needed signalMean and noiseMean from pHits and pFalseAlarms
  [signalMean noiseMean titleStr] = use2AFC(p2AFC,signalSTD,noiseSTD,criterion);
  if isempty(signalMean),return,end
end

% figure out the range of signal strengths over which to compute the distributions
% for signal and noise (this is a bit arbitrary, we just want to be able to see
% the distributions visually in the plots and be able to compute stuff
% out to where the distributions are nearly zero).
signalStrengths = getSignalStrengths(signalMean,signalSTD,noiseMean,noiseSTD);

% compute the signal and noise probability distribution functions
% in this case, gaussian distributions with the appropriate mean
% and standrad deviation
signalPDF = gaussDist(signalStrengths,signalMean,signalSTD);
noisePDF = gaussDist(signalStrengths,noiseMean,noiseSTD);

% compute ROC curve from distributions. That is, compute the hitRate and falseAlarmRate
% as a function of criterion and these will be plotted against each other
for iSignalStrength = 1:length(signalStrengths)
  hitRate(iSignalStrength) = sum(signalPDF(iSignalStrength:end)/sum(signalPDF));
  falseAlarmRate(iSignalStrength) = sum(noisePDF(iSignalStrength:end)/sum(signalPDF));
end

% calculate area under ROC (roughly intergate by computing the rectangles
% under the curve - each one has height hitRate and width is the difference
% between falseAlarmRates on either side. Then sum all that together
areaUnderROC = sum(hitRate(2:end).*-diff(falseAlarmRate));

% plot signal and noise distribution
f = mlrSmartfig('sdtPlot','reuse');clf;set(f,'Render','OpenGL');
subplot(1,1+dispROC,1);hold on
% z is just to make these plot in front of the filled pieces
z = 0.5*ones(1,length(signalStrengths));
if dispSignal,plot3(signalStrengths,signalPDF,z,'-','LineWidth',lineWidth,'Color',signalColor);end
if dispNoise,plot3(signalStrengths,noisePDF,z,'-','LineWidth',lineWidth,'Color',noiseColor);end

% get the index in signalStrengths where the criterion lives (this is just for
% convenience in drawing the filled distributions)
[~,criterionIndex] = min(abs(signalStrengths-criterion));

% draw criterion (set zData to make it plot in front of the distributions)
if (criterion>min(signalStrengths)) && (criterion < max(signalStrengths)) && dispCriterion
  h = vline(criterion,'w-');set(h,'LineWidth',lineWidth);set(h,'ZData',[1 1]);
end

% plot the hit and false alarm rates filled in regions
if dispSignal,fill(signalStrengths([criterionIndex criterionIndex:end end criterionIndex]),[0 signalPDF(criterionIndex:end) 0 0],signalColor,'FaceAlpha',hitAlpha,'EdgeAlpha',0);end
% if lower criterion is set, then we have a special
% case in which we are going to only plot a rectangle
% between lowerCriterion and criterion in the noise 
% distribution
if isempty(lowerCriterion)
  % drawing the usual criterion to infinity
  if dispNoise,fill(signalStrengths([criterionIndex criterionIndex:end end criterionIndex]),[0 noisePDF(criterionIndex:end) 0 0],noiseColor,'FaceAlpha',falseAlarmAlpha,'EdgeAlpha',0);end
else
  % drawing the special case of lowerCriterion to criterion
  [~,lowerCriterionIndex] = min(abs(signalStrengths-lowerCriterion));
  if dispNoise,fill(signalStrengths([lowerCriterionIndex lowerCriterionIndex:criterionIndex criterionIndex lowerCriterionIndex]),[0 noisePDF(lowerCriterionIndex:criterionIndex) 0 0],noiseColor,'FaceAlpha',falseAlarmAlpha,'EdgeAlpha',0);end
end  

% plot the miss and correct-reject filled in regions
if dispSignal,fill(signalStrengths([1 1:criterionIndex criterionIndex 1]),[0 signalPDF([1:criterionIndex]) 0 0],signalColor,'FaceAlpha',missAlpha,'EdgeAlpha',0);end
% again, lowerCriterion is special case in which we have to draw
% either side, usually we just draw the region below criterion
if isempty(lowerCriterion)
  % draw from criterion to the left
  if dispNoise,fill(signalStrengths([1 1:criterionIndex criterionIndex 1]),[0 noisePDF([1:criterionIndex]) 0 0],noiseColor,'FaceAlpha',correctRejectAlpha,'EdgeAlpha',0);end
else
  % draw the two regions for the lowerCriterion special case
  if dispNoise
    fill(signalStrengths([1 1:lowerCriterionIndex lowerCriterionIndex 1]),[0 noisePDF([1:lowerCriterionIndex]) 0 0],noiseColor,'FaceAlpha',correctRejectAlpha,'EdgeAlpha',0);
    fill(signalStrengths([criterionIndex criterionIndex:end end criterionIndex]),[0 noisePDF([criterionIndex:end]) 0 0],noiseColor,'FaceAlpha',correctRejectAlpha,'EdgeAlpha',0);
  end
end  
  

% draw nice axis
drawPublishAxis('whichAxis=horizontal','xAxisMin',min(signalStrengths),'xAxisMax',max(signalStrengths),'xLabel','Signal Strength','yLabel','Probability','xTickLabelHide',true,'labelFontSize',20);

if dispROC
  % plot ROC curve
  subplot(1,2,2);hold on
  % plotROC == 0 means plot everything
  if plotROC==0
    plot(falseAlarmRate,hitRate,'k-','LineWidth',lineWidth*1.5);
    % negative one means plot from the lowest criterion to current criterion
  elseif plotROC == -1
    plot(falseAlarmRate(1:criterionIndex),hitRate(1:criterionIndex),'k-','LineWidth',lineWidth*1.5);
    % positive one means plot from the highest criterion to current criterion
  elseif plotROC == 1
    plot(falseAlarmRate(criterionIndex:end),hitRate(criterionIndex:end),'k-','LineWidth',lineWidth*1.5);
    % if empty, then just plot the one point
  else
    plot(falseAlarmRate,hitRate,'k-');
    plot(falseAlarmRate(criterionIndex),hitRate(criterionIndex),'ko','MarkerFaceColor','k','MarkerSize',12);
  end

  xaxis(0,1);yaxis(0,1);
  axis square;
  % draw diagonal
  dline('k--');
  %draw grid
  for gridPos = 0:0.1:1
    vline(gridPos);hline(gridPos);
  end
  
  % string to display area under ROC and d'
  rocStr = sprintf('Area under ROC = %0.3f (d''=%.3f)',areaUnderROC,(signalMean-noiseMean)/sqrt(0.5*signalSTD^2 + 0.5*noiseSTD^2));
  % put into title
  if isempty(titleStr)
    titleStr = rocStr;
  else
    titleStr = sprintf('%s\n%s',rocStr,titleStr);
  end
  if dispTitle,title(titleStr);end
  drawPublishAxis('xLabel','False Alarms (proportion)','yLabel','Hits (proportion)','labelFontSize',20,'xTick',[0:0.5:1],'yTick',[0:0.5:1])
end

%%%%%%%%%%%%%%%%%%%
%    gaussDist    %
%%%%%%%%%%%%%%%%%%%
function y = gaussDist(x,m,s)

y = (1/(s*sqrt(2*pi)))*exp(-(((x-m).^2)/(2*s^2)));

%%%%%%%%%%%%%%%%%%%%%%%
%    proportionToZ    %
%%%%%%%%%%%%%%%%%%%%%%%
function z = proportionToZ(proportion)

if proportion >= 1
  z = -inf
elseif proportion <= 0
  z = inf;
elseif proportion > 0.5
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    useHitsFalseAlarms    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signalMean noiseMean titleStr] = useHitsFalseAlarms(pHits,pFalseAlarms,signalSTD,noiseSTD,criterion)

signalMean = [];
noiseMean = [];
titleStr = [];

if isempty(pHits) || isempty(pFalseAlarms)
  disp(sprintf('(sdtPlot) Must specify both pHits and pFalseAlarms'));
  return
else
  % make sure we have positive numbers for pHits and pFalseAlarms
  if (pHits < 0) || (pFalseAlarms < 0)
    disp(sprintf('(sdtPlot) pHits and pFalseAlarms cannot be less than 0'));
    return
  end
  % make sure we are given a proportion not a percentile (and convert if we think differently)
  if (pHits>1) || (pFalseAlarms>1)
    if (pHits>100) || (pFalseAlarms>100)
      disp(sprintf('(sdtPlot) pHits and pFalseAlarms cannot be greater than 100'));
      return
    end
    pHits = pHits/100;
    pFalseAlarms = pFalseAlarms/100;
  end
  % ok, now that we have a proportion, convert into z-scores using the inverse complementary error function
  % which is a a way to find the point in a gaussian distribution (in std units) where the area under
  % the curve from there to infinity will be equal to the proportion of hits or false alarms
  zHits = -sqrt(2)*erfcinv(2*pHits);
  zFalseAlarms = -sqrt(2)*erfcinv(2*pFalseAlarms);
  % ok, now that we have zHits and zFalseAlarms, we should be able to figure out where
  % to place the signal and noiseMeans to achieve that
  signalMean = criterion+signalSTD*zHits;
  noiseMean = criterion+noiseSTD*zFalseAlarms;
  % set title string
  titleStr = sprintf('d''=%0.3f (zHits = %0.3f zFalseAlarms = %0.3f)',zHits-zFalseAlarms,zHits,zFalseAlarms);
end

%%%%%%%%%%%%%%%%
%    use2AFC   %
%%%%%%%%%%%%%%%%
function [signalMean noiseMean titleStr] = use2AFC(p2AFC,signalSTD,noiseSTD,criterion)

% default values
signalMean = [];
noiseMean = [];
titleStr = '';
numSignalSteps = 100;

% start parameters
stepSize = 1;
distDiff = 1;

% number of reversals in our search
reversals = 0;
searchDirection = 1;

% arbitrary number controls how good our estimate will be
numReversals = 19;

% we are just going to do a search for the signal mean and noise mean
% that gives us an area under the ROC equal to the desired 2AFC. We
% do this by a staircasing procedure - starting at arbitrary difference
% in the signal and noise means, going in one direction until we over
% step the area under the roc and then switching back (making the
% stepping size smaller) until we converge
disppercent(-inf,sprintf('(sdtPlot) Searching for signal difference that gives area under the ROC of %f',p2AFC));
while (reversals < numReversals)
  % try a signal and noise mean
  signalMean = criterion+distDiff/2;
  noiseMean = criterion-distDiff/2;
  % get signal strengths to compute over
  signalStrengths = getSignalStrengths(signalMean,signalSTD,noiseMean,noiseSTD,numSignalSteps);
  % compute signal and noise distributions
  signalPDF = gaussDist(signalStrengths,signalMean,signalSTD);
  noisePDF = gaussDist(signalStrengths,noiseMean,noiseSTD);
  % compute hit/false alarm rates
  for iSignalStrength = 1:length(signalStrengths)
    hitRate(iSignalStrength) = sum(signalPDF(iSignalStrength:end)/sum(signalPDF));
    falseAlarmRate(iSignalStrength) = sum(noisePDF(iSignalStrength:end)/sum(signalPDF));
  end
  % compute aread under ROC
  areaUnderROC = sum(hitRate(2:end).*-diff(falseAlarmRate));
  % check to see how good we are doing
  if p2AFC > areaUnderROC
    % see if we need to switch directions
    if searchDirection<0
      % switch directions
      searchDirection = 1;
      % increment number of reversals
      reversals = reversals + 1;
      % step size gets smaller
      stepSize = stepSize/2;
      % increase the number of signal steps to compute (so that we 
      % get a more accurate area - the 1.2 is arbitrary here - just
      % a number that allows signals steps to grow but not out of
      % control)
      numSignalSteps = round(numSignalSteps*1.2);
      % update disppercent
      disppercent(reversals/numReversals);
    end
  else
    % see if we need to switch directions
    if searchDirection>0
      % switch directions
      searchDirection = -1;
      % increment number of reversals
      reversals = reversals + 1;
      % step size gets smaller
      stepSize = stepSize/2;
      % increase the number of signal steps to compute (so that we 
      % get a more accurate area - the 1.2 is arbitrary here - just
      % a number that allows signals steps to grow but not out of
      % control)
      numSignalSteps = round(numSignalSteps*1.2);
      % update disppercent
      disppercent(reversals/numReversals);
    end
  end
  % update dist diff
  distDiff = distDiff + stepSize*searchDirection;
end
disppercent(inf);

titleStr = sprintf('d''=%0.3f',signalMean/signalSTD-noiseMean/noiseSTD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getSignalStrengths    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function signalStrengths = getSignalStrengths(signalMean,signalSTD,noiseMean,noiseSTD,numSteps)

% how fine to sample
if nargin < 5,numSteps = 100;end

% get minimum and maximum based on arbitarily going 3.5 std from the min or max
minSignalStrength = min(signalMean,noiseMean)-3.5*max(signalSTD,noiseSTD);
maxSignalStrength = max(signalMean,noiseMean)+3.5*max(signalSTD,noiseSTD);

% get the signal strengths
signalStrengths = minSignalStrength:(maxSignalStrength-minSignalStrength)/10000:maxSignalStrength;



