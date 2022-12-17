% emriPlot(v,overlayNum,scan,x,y,z)
%
%      usage: emriPlot(v,overlayNum,scan,x,y,z,roi
%         by: justin gardner (based on corAnalPlot)
%       date: 12/16/2022
%    purpose: emriPlot function
%
function emriPlot(v,overlayNum,scan,x,y,z,roi)

overlayNames = viewGet(v,'overlayNames');

% look for the co overlay that has the highest coherence value for this
% voxel
maxCo = -inf;
for iOverlay = 1:length(overlayNames)
  if ~isempty(strfind(overlayNames{iOverlay},'co'))
    co = viewGet(v,'overlay',overlayNames{iOverlay});
  end
  % if the coherence value is greater than current maximum
  if co.data{scan}(x,y,z) > maxCo
    %then remember this and move to the next one
    coOverlayNum = iOverlay;
    maxCo = co.data{scan}(x,y,z);
  end
end

% if no maximum, there must not have been a coherence overlay
if maxCo == -inf
  disp(sprintf('(emriPlot) Could not find any cohernece overlays'));
  return
end

% If corAnal is loaded, then use it. Otherwise, error.
co = viewGet(v,'overlay',overlayNames{coOverlayNum});
amp = viewGet(v,'overlay',overlayNames{coOverlayNum+1});
ph = viewGet(v,'overlay',overlayNames{coOverlayNum+2});
if isempty(co) || isempty(amp) || isempty(ph)
  disp(sprintf('(emriPlot) Could not find correct overlays'));
  return
end

% see if the shift key is down on MLR fig
shiftDown = any(strcmp(get(viewGet(v,'figureNumber'),'SelectionType'),'extend'));

% get the peak frequency (i.e. the one with the highest frequency
ncycles = co.params.ncycles;

% Analysis parameters
detrend = viewGet(v,'detrend',scan);
spatialnorm = viewGet(v,'spatialnorm',scan);
junkframes = viewGet(v,'junkframes',scan);
nframes = viewGet(v,'nframes',scan);
framePeriod = viewGet(v,'framePeriod',scan);
trigonometricFunction = viewGet(v,'trigonometricFunction',scan);

if ncycles == 0
  mrWarnDlg(sprintf('(corAnalPlot) Number of cycles per scan for the corAnal should not be: %i',ncycles));
  return
end

% don't do roi if shift key is down
if shiftDown
  roi = [];
elseif length(roi) > 0
  oneTimeWarning('eventRelatedPlotShiftKey',sprintf('(corAnalPlot) To avoid showing ROI plots, hold shift down when clicking'),1);
end

% now if we have an roi then load its time series
% and get the mean and plot that instead of the
% time series for the voxel
roiPlot = 0;
if ~isempty(roi)
    roi = loadROITSeries(v,roi{1},viewGet(v,'curScan'),viewGet(v,'curGroup'));
    tseries = mean(roi.tSeries,1)';
    ptseriesSte = std(100*roi.tSeries/mean(tseries),1,1)'/sqrt(roi.n);
    ptseriesSte = ptseriesSte(junkframes+1:junkframes+nframes);
    headerStr = sprintf('Times series from roi %s (n=%i)',roi.name,roi.n);
    roiPlot = 1;
else
    % Load tseries from file. Error if file doesn't exist.
    pathStr = viewGet(v,'tseriesPathStr',scan);
    if ~exist(pathStr,'file')
        mrErrorDlg(['File ',pathStr,' not found']);
    end
    [tseries,hdr] = mlrImageReadNifti(pathStr,{x,y,z,[]});
    headerStr = sprintf('Times series from voxel [%i %i %i] ',x,y,z);
    tseries = squeeze(tseries);
end

tseries = tseries(junkframes+1:junkframes+nframes);

%recalculate CorAnal (not necessary for single voxels, but greatly simplifies the code)
[coval, ampval, phval, ptseries] = computeCoranal(tseries,ncycles,detrend,spatialnorm,trigonometricFunction);
 
% calculate mean cycle
if rem(nframes,ncycles)
  disp(sprintf('(corAnalPlot) Scan has %i frames, not evenly divisible by %i cycles',nframes,ncycles));
end
%return
cycleLen = floor(nframes/ncycles);
singleCycleTSeries = ptseries(1:cycleLen*ncycles);
singleCycle = mean(reshape(singleCycleTSeries,cycleLen,ncycles)');
singleCycleSte = std(reshape(singleCycleTSeries,cycleLen,ncycles)')/sqrt(ncycles);

% calculate fft
ft = fft(ptseries);
absft = abs(ft(1:round((length(ft)/2))));
signalAmp = absft(ncycles+1);
noiseFreq = round(2*length(absft)/3):length(absft);
noiseAmp = mean(absft(noiseFreq));
snr = signalAmp/noiseAmp;

% Create model fit
switch (trigonometricFunction)
  case 'Cosine'
    model = ampval * cos(2*pi*ncycles/nframes * [0:nframes-1]' - phval);
  case 'Sine'
    model = ampval * sin(2*pi*ncycles/nframes * [0:nframes-1]' - phval);
end

% Change xlabel to be in seconds
% Change ylabel contingent upon detrend to include units

% Select window
selectGraphWin;

% Plot it
subplot(2,3,1:3)
time = linspace(framePeriod/2,(nframes-.5)*framePeriod,nframes)';
if roiPlot
    errorbar(time,ptseries,ptseriesSte,'k.-','LineWidth',1);hold on
else
    plot(time,ptseries,'k.-','LineWidth',1);hold on
end
plot(time,model,'r-','LineWidth',1.5);

% Ticks and labels
fontSize = 14;
set(gca,'FontSize',fontSize);
headerStr = sprintf('%s freq=%i/scan co=%f amp=%f ph(%s)=%f (%i deg)\n%s',...
  headerStr,ncycles,coval,ampval,trigonometricFunction,phval,round(180*phval/pi),viewGet(v,'description',viewGet(v,'curScan')));
title(headerStr,'Interpreter','none');
xtickStep = nframes*framePeriod/ncycles;
xtick = ceil([0:xtickStep:nframes*framePeriod]);
xtick = sort(unique(xtick));
set(gca,'xtick',xtick);
set(gca,'XLim',ceil([0,nframes*framePeriod]));
set(gca,'xgrid','on')
xlabel('Time (sec)');
switch spatialnorm
    case {'Divide by mean'}
        ylabel('fMRI response (% change in image intensity)');
    otherwise
        ylabel('fMRI response(arbitrary units)');
end


% Plot single cylce
subplot(2,3,4)
time = linspace(framePeriod/2,(nframes/ncycles-.5)*framePeriod,nframes/ncycles)';
errorbar(time,singleCycle,singleCycleSte,'k.-','LineWidth',1);hold on
plot(time,model(1:length(time)),'r-','LineWidth',1.5);

set(gca,'FontSize',fontSize);
title('Single cycle with STE across cycles');
xlabel('Time (sec)');
switch spatialnorm
    case {'Divide by mean'}
        ylabel('fMRI response (% change in image intensity)');
    otherwise
        ylabel('fMRI response(arbitrary units)');
end

% plot fourier amplitude
subplot(2,3,5:6)
plot(0:(length(absft)-1),absft,'k.-');hold on
plot(ncycles,signalAmp,'ro');
plot(noiseFreq-1,absft(noiseFreq),'go');
set(gca,'FontSize',fontSize);
title(sprintf('Stimulus (red): %f Noise (Mean of green): %f CNR: %f',signalAmp,noiseAmp,snr));
ylabel('Magnitude');
xlabel('Fourier component number');

return


