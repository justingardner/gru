% emriPlot(v,overlayNum,scan,x,y,z)
%
%      usage: emriPlot(v,overlayNum,scan,x,y,z,roi
%         by: justin gardner (based on corAnalPlot)
%       date: 12/16/2022
%    purpose: emriPlot function
%
function emriPlot(v,overlayNum,scan,x,y,z,roi)

% display what we are doing (as this function can take a little bit to put
% up the figure)
if ~isempty(roi) && isfield(roi{1},'name')
  disp(sprintf('(emriPlot) Creating plot for roi: %s',roi{1}.name));
else
  disp(sprintf('(emriPlot) Creating plot for voxel: [%i %i %i]',x,y,z));
end    

% see if the shift key is down on MLR fig
shiftDown = any(strcmp(get(viewGet(v,'figureNumber'),'SelectionType'),'extend'));

% get time series
[t, tSeries,tSeriesSte,headerStr] = getTSeries(v,scan,x,y,z,roi,shiftDown);
if isempty(tSeries),return,end

% calculate fft
ft = fft(tSeries);
absft = abs(ft(1:round((length(ft)/2))));

% plot the basic time series and FFT
nRows = 3; nCols = 3;
emriPlotTSeries(t, tSeries, tSeriesSte, absft, nRows, nCols, headerStr);

% plot auto correlation
emriPlotAutoCorrelation(v, tSeries, nRows, nCols)

% get the coherence/amp/ph overlays with the largest coherence value
[co,amp,ph] = getOverlaysWithMaxCoVal(v,scan,x,y,z,roi);

% if we found one (i.e. corAnal was run)
if ~isempty(co)
  % get parameters for corAnal
  corAnalParams = co.params;curScan = viewGet(v,'curScan');
  nCycles = corAnalParams.ncycles(curScan);
  detrend = corAnalParams.detrend{curScan};
  spatialnorm = corAnalParams.spatialnorm{curScan};
  trigonometricFunction = corAnalParams.trigonometricFunction{curScan};
  
  % run the coranal (note that we do this because the tSeries here may be
  % the average of an ROI for which coranal has not yet been run.
  [coVal, ampVal, phVal, corAnalTSeries] = computeCoranal(tSeries,nCycles,detrend,spatialnorm,trigonometricFunction);

  % plot the cor anal
  emriPlotCorAnal(t, tSeries, tSeriesSte, absft, nCycles, coVal, ampVal, phVal, trigonometricFunction, nRows, nCols, headerStr)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot auto-correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function emriPlotAutoCorrelation(v, tSeries, nRows, nCols)

% plot settings
fontSize = 14;
markerSize = 12;

% normalize time series, so that auto-correlation goes from -1 to 1
tSeries = tSeries-mean(tSeries);
tSeries = tSeries/sqrt(sum(tSeries.^2));
corrSeries = xcorr(tSeries,tSeries);

% compute the permutation analysis for null
nPermutations = 100;
for iRand = 1:nPermutations
  randTSeries = tSeries(randperm(length(tSeries)));
  xCorrRand(iRand,:) = xcorr(randTSeries,randTSeries);
end
xCorrRand = sort(xCorrRand);
minXCorr = xCorrRand(round(0.025*nPermutations),:);
maxXCorr = xCorrRand(round(0.975*nPermutations),:);

% get the shift values for the x axes
shiftVal = -length(tSeries)+1:length(tSeries)-1;
shiftVal = shiftVal * viewGet(v,'framePeriod');

% see how many values are below or above 95 % confidence intreval
outsideConfidenceInterval = 100* sum((corrSeries' > maxXCorr) | (corrSeries' < minXCorr)) / length(corrSeries);

% set the 0 lag to nan (as it is by definition equal to 1
minXCorr(shiftVal == 0) = nan;
maxXCorr(shiftVal == 0) = nan;
corrSeries(shiftVal == 0) = nan;

% plot the auto-correlation function
subplot(nRows,nCols,nCols*2+1:nCols*3);
plot(shiftVal,minXCorr,'r-');hold on
set(gca,'FontSize',fontSize);
plot(shiftVal,maxXCorr,'r-');
plot(shiftVal,corrSeries,'k.-','LineWidth',1,'MarkerSize',markerSize);
title(sprintf('Time series auto-correlation (0 lag correlation set to nan)\n%0.1f%% outside 95%% confidence interval for null distribution',outsideConfidenceInterval));
set(gca,'XLim',[min(shiftVal) max(shiftVal)]);
xlabel('Lag (s)');
ylabel('Correlation (r)');
mylegend({'Auto-correlation','95% confidence interval'},{'k-','r-'});
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getPeakCorrelationValue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [co,amp,ph] = getOverlaysWithMaxCoVal(v,scan,x,y,z,roi)

% default return values
co = [];amp = [];ph = [];

% get overlay names
overlayNames = viewGet(v,'overlayNames');

% look for the co overlay that has the highest coherence value for this
% voxel
maxCo = -inf;maxCoOverlayNum = nan;
for iOverlay = 1:length(overlayNames)
  if ~isempty(strfind(overlayNames{iOverlay},'co'))
    co = viewGet(v,'overlay',overlayNames{iOverlay});
  end
  % if the coherence value is greater than current maximum
  if co.data{scan}(x,y,z) > maxCo
    %then remember this and move to the next one
    maxCoOverlayNum = iOverlay;
    maxCo = co.data{scan}(x,y,z);
  end
end

% if no maximum, there must not have been a coherence overlay
if maxCo == -inf
  disp(sprintf('(emriPlot) Could not find any cohernece overlays'));
  return
end

% load the coherence / amplitude and phase overlay that had the largest
% coherence value
co = viewGet(v,'overlay',overlayNames{maxCoOverlayNum});
amp = viewGet(v,'overlay',overlayNames{maxCoOverlayNum+1});
ph = viewGet(v,'overlay',overlayNames{maxCoOverlayNum+2});
if isempty(co) || isempty(amp) || isempty(ph)
  disp(sprintf('(emriPlot) Could not find correct overlays'));
end

%%%%%%%%%%%%%%%%
%% getTSeries
%%%%%%%%%%%%%%%%
function [t, tSeries,tSeriesSte,headerStr] = getTSeries(v,scan,x,y,z,roi,shiftDown)

% default return parameters
t = [];
tSeries = [];
tSeriesSte = [];
headerStr = '';

% frame clipping parameters
junkframes = viewGet(v,'junkframes',scan);
nframes = viewGet(v,'nframes',scan);
framePeriod = viewGet(v,'framePeriod',scan);

% don't do roi if shift key is down
if shiftDown
  roi = [];
elseif length(roi) > 0
  oneTimeWarning('eventRelatedPlotShiftKey',sprintf('(emriPlot) To avoid showing ROI plots, hold shift down when clicking'),1);
end

% now if we have an roi then load its time series
% and get the mean and plot that instead of the
% time series for the voxel
roiPlot = 0;
if ~isempty(roi)
  roi = loadROITSeries(v,roi{1},viewGet(v,'curScan'),viewGet(v,'curGroup'));
  tSeries = mean(roi.tSeries,1)';
  tSeriesSte = std(100*roi.tSeries/mean(tSeries),1,1)'/sqrt(roi.n);
  headerStr = sprintf('Time series from roi %s (n=%i)\n%s',roi.name,roi.n,viewGet(v,'description',viewGet(v,'curScan')));
  roiPlot = 1;
else
  % Load tseries from file. Error if file doesn't exist.
  pathStr = viewGet(v,'tseriesPathStr',scan);
  if ~exist(pathStr,'file')
    disp(sprintf('(emriPlot) Could not find tSeries data file: %s',pathStr));
    return
  end
  [tSeries,hdr] = mlrImageReadNifti(pathStr,{x,y,z,[]});
  headerStr = sprintf('Time series from voxel [%i %i %i]\n%s',x,y,z,viewGet(v,'description',viewGet(v,'curScan')));
  tSeries = squeeze(tSeries);
  tSeriesSte = nan(1,length(tSeries));
end

% get the time series
t = linspace(framePeriod/2,(length(tSeries)-.5)*framePeriod,length(tSeries))';

% reomove junk frames and adjust length (if needed)
tSeries = tSeries(junkframes+1:junkframes+nframes);
tSeriesSte = tSeriesSte(junkframes+1:junkframes+nframes);
t = t(junkframes+1:junkframes+nframes);

% get percentTSeries
meanTSeries = mean(tSeries);
tSeries = 100*(tSeries-meanTSeries)/meanTSeries;

%%%%%%%%%%%%%%%%%%%%%%%
%% computeCorAnalFit
%%%%%%%%%%%%%%%%%%%%%%%
function [corAnalFit, singleCycle, singleCycleSte] = computeCoranalFit(tSeries, nCycles, coVal, ampVal, phVal,trigonometricFunction)

% default return values
corAnalFit = [];
corAnalSingleCycleFit = [];
singleCycle = [];
singleCycleSte = [];

% get length of tSeries
nFrames = length(tSeries);

% calculate mean cycle
if rem(nFrames,nCycles)
  % chop off the tSeries at the end. Could so something
  % fancier and add the partial cycle back into the man and ajdust
  % the ste for the number of points, but seems like overkill.
  nFrames = floor(nFrames/nCycles)*nCycles;
end

% compute the single cycle and ste
cycleLen = floor(nFrames/nCycles);
singleCycleTSeries = tSeries(1:nFrames);
singleCycle = mean(reshape(singleCycleTSeries,cycleLen,nCycles)');
singleCycleSte = std(reshape(singleCycleTSeries,cycleLen,nCycles)')/sqrt(nCycles);

% Create model fit
switch (trigonometricFunction)
  case 'Cosine'
    corAnalFit = ampVal * cos(2*pi*nCycles/nFrames * [0:nFrames-1]' - phVal);
  case 'Sine'
    corAnalFit = ampVal * sin(2*pi*nCycles/nFrames * [0:nFrames-1]' - phVal);
end

%%%%%%%%%%%%%%%%%%%%%%
% plot the tSeries
%%%%%%%%%%%%%%%%%%%%%%
function emriPlotTSeries(t, tSeries, tSeriesSte, absft, nRows, nCols, headerStr)

% plot settings
markerSize = 12;
fontSize = 14;

% Select window
selectGraphWin;

% Plot timeSeries
subplot(nRows,nCols,1:nCols)
if all(isnan(tSeriesSte))
  plot(t,tSeries,'k.-','LineWidth',1,'MarkerSize',markerSize);hold on
else
  myerrorbar(t,tSeries,'yError',tSeriesSte,'LineWidth',1,'MarkerSize',markerSize);hold on
end
title(headerStr);
xlabel('Time (s)');
ylabel('Response amplitude (% signal change)');
set(gca,'xgrid','on')
set(gca,'XLim',[min(t) max(t)]);
set(gca,'FontSize',fontSize);
grid on

% Plot Fourier Transform
subplot(nRows,nCols,nCols+1:2*nCols)
plot(0:(length(absft)-1),absft,'k.-','MarkerSize',markerSize);hold on
set(gca,'FontSize',fontSize);
ylabel('Magnitude');
xlabel('Fourier component number');
set(gca,'XLim',[0 (length(absft)-1)]);


%%%%%%%%%%%%%%%%%%%%%%
% plot corAnal
%%%%%%%%%%%%%%%%%%%%%%
function emriPlotCorAnal(t, tSeries, tSeriesSte, absft, nCycles, coVal, ampVal, phVal, trigonometricFunction, nRows, nCols, headerStr)

% plot se
markerSize = 12;
fontSize = 14;

% compute cor anal model fits
[corAnalFit, singleCycle, singleCycleSte] = computeCoranalFit(tSeries, nCycles, coVal, ampVal, phVal, trigonometricFunction);

% calculate SNR as the signal at the analysis frequency divided by the
% mean at the top 1/3 of high frequencies (as an estimate of noise)
signalAmp = absft(nCycles+1);
noiseFreq = round(2*length(absft)/3):length(absft);
noiseAmp = mean(absft(noiseFreq));
snr = signalAmp/noiseAmp;

% Select window
selectGraphWin(true);

% Set title for time series
subplot(nRows,nCols,1:nCols);
headerStr = sprintf('%s\nfreq=%i cycles/scan (%0.1f cycles/sec) co=%f amp=%f ph(%s)=%f (%i deg)',...
  headerStr,nCycles,nCycles*(1/t(end)),coVal,ampVal,trigonometricFunction,phVal,round(180*phVal/pi));
title(headerStr,'Interpreter','none');

% put the ticks upon each cycle
cycleLenSec = t(length(singleCycle))-t(1);

% Plot model fit
hold on; plot(t(1:length(corAnalFit)),corAnalFit,'r-','LineWidth',1);

% Plot fft
subplot(nRows,nCols,nCols+2:nCols*2)
plot(0:(length(absft)-1),absft,'k.-','MarkerSize',markerSize);hold on
set(gca,'FontSize',fontSize);
ylabel('Magnitude');
xlabel('Fourier component number');
plot(nCycles,signalAmp,'ro');
plot(noiseFreq-1,absft(noiseFreq),'go');
title(sprintf('Stimulus (red): %f Noise (Mean of green): %f CNR: %f',signalAmp,noiseAmp,snr));
set(gca,'XLim',[0 (length(absft)-1)]);

% if auto-correlation was plotted
if nRows >= 3
  % then draw the cycle length
  subplot(nRows,nCols,nCols*2+1:nCols*3);
  vline([-nCycles*cycleLenSec:cycleLenSec:nCycles*cycleLenSec],'r:');
  mylegend({'Auto-correlation','95% confidence interval'},{'k-','r-'});
end

% Plot single cylce
subplot(nRows,nCols,nCols+1)
if nCycles > 1
  % plot mean and standard error over cycles
  myerrorbar(t(1:length(singleCycle)),singleCycle,'yError',singleCycleSte,'LineWidth',1,'MarkerSize',markerSize);hold on
  plot(t(1:length(singleCycle)),corAnalFit(1:length(singleCycle)),'r-','LineWidth',1.5);
  set(gca,'XLim',[t(1) t(length(singleCycle))]);
else
  % just replot, only one cycle
  plot(t,tSeries,'k.-','LineWidth',1,'MarkerSize',markerSize);hold on
  plot(t,corAnalFit,'r-','LineWidth',1.5);
  cycleLenSec = t(end)-t(1);
end
    
% labels
set(gca,'FontSize',fontSize);
title(sprintf('Single cycle with STE across cycles\n(cycleLen: %0.1fms)',cycleLenSec*1000));
xlabel('Time (s)');
ylabel('Response amplitude (% signal change)');

