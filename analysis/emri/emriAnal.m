function [v params] = emriAnal(v,params,varargin)
%
%      usage: emriAnal(v,params,varargin)
%         by: justin gardner (Adapted from corAnal and pRF code)
%       date: 12/16/2022
%    purpose: compute various analyses for DIANA sequence data
%
%             if you just want a default parameter structure you
%             can do:
% 
%             v = newView;
%             [v params] = emriAnal(v,[],'justGetParams=1','defaultParams=1','scanList=1')
%

% check arguments
if nargin < 1
  help emri
  return
end

d = [];

% params defaults to empty
if nargin < 2,params =[];end

% other arguments
justGetParams=[];defaultParams=[];scanList=[];groupNum=[];
getArgs(varargin,{'justGetParams=0','defaultParams=0','scanList=[]','groupNum=[]'});

% first get parameters
if isempty(params)
  % get group
  if isempty(groupNum),groupNum = viewGet(v,'curGroup');end
  % put up the gui
  params = emriGUI('v',v,'groupNum',groupNum,'defaultParams',defaultParams,'scanList',scanList);
end

% just return parameters
if justGetParams,d = params;return,end

% Reconcile params with current status of group and ensure that it has
% the required fields. 
params = defaultReconcileParams([],params);

% Abort if params empty
if isempty(params),return,end

% check the params
params = checkEmriParams(params);

% set the group
v = viewSet(v,'curGroup',params.groupName);

% filter the time series, save to use later instead of loadTSeries
    % if params.temporalFiltering = 0, filteredTSeries is same as unfiltered.
    filteredTSeries = temporalFilterTSeries(v,params);
    
    % same thing but with spatial
    filteredTSeries = spatialFilterTSeries(v,params,filteredTSeries)

% run the frequency analysis
if params.frequencyAnalysis
  runCorAnal(v,params,filteredTSeries);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initOverlay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = initOverlay(v,nScans,varargin)

getArgs(varargin,{'name=default','functionName=emriAnal','reconcileFunction=defaultReconcileParams','mergeFunction=defaultMergeParams','colormap',[],'interrogatorFunction=emriPlot','params',[],'overlayRange',[0 1]});

o.name = name;
o.function = functionName;
o.reconcileFunction = reconcileFunction;
o.mergeFunction = mergeFunction;
if isempty(colormap)
  o.colormap = jet(256);
else
  o.colormap = colormap;
end    
o.date = datestr(now);
o.interrogator = interrogatorFunction;
o.groupName = viewGet(v,'groupName');
o.params = params;
o.range = overlayRange;
o.data = cell(1,nScans);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temporalFilterTSeries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filteredTSeries = temporalFilterTSeries(v,params)

% get params and scan list
corAnalParams.recompute = zeros(1,viewGet(v,'nScans'));
corAnalParams.recompute(params.scanNum) = 1;
scanList = find(corAnalParams.recompute(:));

% load the unfiltered time series for all slices and scans
for scanIndex=1:length(scanList)
    scanNum = scanList(scanIndex);
    unfilteredTSeries{scanNum} = loadTSeries(v,scanNum);
end

% if the temporal filtering box was checked, filter the data
if params.temporalFiltering

    % make the filter. add more options here.  
        %box smoothing
        if strcmp(params.temporalFilter,'Box')
            filter = ones(1,params.temporalFilterWidth)/params.temporalFilterWidth; %placeholder box filter.
        end
    
        %gausian smoothing
        if strcmp(params.temporalFilter,'Gaussian')
            filterWidth = params.temporalFilterWidth*2+1;
            filter = gausswin(filterWidth)/sum(gausswin(filterWidth));
        end

        %TODO - add other filters.
    
    % do the filtering
    for scanIndex=1:length(scanList)
    
        scanNum = scanList(scanIndex);
    
        % get scan dimensions
        [scanDim1 scanDim2 nslices nTimepoints] = size(unfilteredTSeries{scanNum});
    
        for dim1 = 1:scanDim1
            for dim2 = 1:scanDim2
                for slice = 1:nslices
    
                    %get the unfiltered tSeries for the voxel at dim1, dim2 in slice+scan
                    tSeriesToFilter = unfilteredTSeries{scanNum}(dim1,dim2,slice,:);
                    %convolve with the filter you made earlier
                    filteredTSeries{scanNum}(dim1,dim2,slice,:) = conv(tSeriesToFilter(:),filter,'same');
    
                end
            end
        end
    end

% if we didn't check the filter box, return the unfiltered TSeries
else
    filteredTSeries = unfilteredTSeries;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spatialFilterTSeries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filteredTSeries = spatialFilterTSeries(v,params,unfilteredTSeries)

% get params and scan list
corAnalParams.recompute = zeros(1,viewGet(v,'nScans'));
corAnalParams.recompute(params.scanNum) = 1;
scanList = find(corAnalParams.recompute(:));

% if the spatial filtering box was checked, filter the data
if params.spatialFiltering

    % make the filter. add more options here.  
        %box smoothing
        if strcmp(params.spatialFilter,'Box')
            filter = ones(1,params.spatialFilterWidth); %placeholder box filter.
            filter = filter'*filter; filter = filter/(params.spatialFilterWidth^2)
        end
    
        %gausian smoothing
        if strcmp(params.spatialFilter,'Gaussian')
            filterWidth = params.spatialFilterWidth*2+1;
            filter = gausswin(filterWidth);
            filter = filter*filter';
        end

        %TODO - add other filters.
    
    % do the filtering
    for scanIndex=1:length(scanList)
    
        scanNum = scanList(scanIndex);
    
        % get scan dimensions
        [scanDim1 scanDim2 nslices nTimepoints] = size(unfilteredTSeries{scanNum});
    
        for slice = 1:nslices
            for timepoint = 1:nTimepoints
    
                    %get the unfiltered tSeries for all voxels in slice, at timepoint
                    tSeriesToFilter = unfilteredTSeries{scanNum}(:,:,slice,timepoint);
                    %convolve with the 2d filter you made earlier
                    filteredTSeries{scanNum}(:,:,slice,timepoint) = conv2(tSeriesToFilter,filter,'same');
                    %%% NOTE -- USING 'same' GIVES YOU EDGE EFFECTS! how do we fix this? -jw 1/09/23
    
            end
        end
    end

% if we didn't check the filter box, return the unfiltered TSeries
else
    filteredTSeries = unfilteredTSeries;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% runCorAnal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runCorAnal(v,params,filteredTSeries)

% get number of scans
nScans = length(params.scanNum);

% set up the params structure as if this were a corAnal
corAnalParams.recompute = zeros(1,viewGet(v,'nScans'));
corAnalParams.recompute(params.scanNum) = 1;
for iScan = 1:viewGet(v,'nScans')
  corAnalParams.ncycles(iScan) = nan;
  corAnalParams.detrend{iScan} = params.detrend;
  if params.divideByMean
    corAnalParams.spatialnorm{iScan} = 'Divide by mean';
  else
    corAnalParams.spatialnorm{iScan} = 'None';
  end
  corAnalParams.trigonometricFunction{iScan} = params.trigonometricFunction;
end

co = {};amp={};ph={};
% run correlation analysis at each frequency 
for iFreq = params.cyclesScanMin:params.cyclesScanMax

  % set the frequency for the analysis
  corAnalParams.ncycles(:) = iFreq;

  % init overlays
  co{end+1} = initOverlay(v,nScans,'name',sprintf('%03i co',iFreq),'params',corAnalParams,'overlayRange',[0 1]);
  amp{end+1} = initOverlay(v,nScans,'name',sprintf('%03i amp',iFreq),'params',corAnalParams,'overlayRange',[0 1]);
  ph{end+1} = initOverlay(v,nScans,'name',sprintf('%03i ph',iFreq),'params',corAnalParams,'overlayRange',[0 2*pi]);

  % Compute it (function is a copy of what is in corAnal
  [co{end},amp{end},ph{end}] = computeCorrelationAnalysis(v,corAnalParams,filteredTSeries,co{end},amp{end},ph{end});

  % Fill range field for amp
  ampMin = realmax; ampMax = 0;
  for iScan=1:nScans
    if ~isempty(amp{end}.data{iScan})
      ampMin = min([ampMin min(amp{end}.data{iScan}(:))]);
      ampMax = max([ampMax max(amp{end}.data{iScan}(:))]);
    end
  end
  if (ampMin <= ampMax)
    amp{end}.range = [ampMin ampMax];
  end
end

% Install emriAnal in the view
emriAnal.name = 'emriAnal';  % This can be reset by editAnalysisGUI
emriAnal.type = 'emriAnal';
emriAnal.groupName = params.groupName;
emriAnal.function = 'emriAnal';
emriAnal.reconcileFunction = 'defaultReconcileParams';
emriAnal.mergeFunction = 'defaultMergeParams';
emriAnal.guiFunction = 'emriGUI';
emriAnal.overlayInterpFunction = 'corAnalInterp';
emriAnal.params = params;
emriAnal.date = datestr(now);
v = viewSet(v,'newanalysis',emriAnal);
for iOverlay = 1:length(co)
  v = viewSet(v,'newoverlay',co{iOverlay});
  v = viewSet(v,'newoverlay',amp{iOverlay});
  v = viewSet(v,'newoverlay',ph{iOverlay});
end
if ~isempty(viewGet(v,'fignum'))
  refreshMLRDisplay(viewGet(v,'viewNum'));
end

% Save it
saveAnalysis(v,emriAnal.name);

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [co,amp,ph] = computeCorrelationAnalysis(view,params,filteredTSeries,co,amp,ph)
% Required fields in params: 'recompute','ncycles','detrend','spatialnorm'

% Get scanList from params.recompute field
scanList = find(params.recompute(:));

disp('Computing corAnal...');
warning('off','MATLAB:divideByZero');
for scanIndex=1:length(scanList)
  scanNum = scanList(scanIndex);
  waitHandle = mrWaitBar(0,['Computing Correlation Analysis for scan ' int2str(scanNum) ':']);

  % sliceDims: [ydim xdim] for single slice
  % volDims; [ydim xdim nslices] for single scan
  sliceDims = viewGet(view,'sliceDims',scanNum);
  volDims = viewGet(view,'dims',scanNum);

  % Initialize data with NaNs
  co.data{scanNum} = NaN*ones(volDims);
  amp.data{scanNum} = NaN*ones(volDims);
  ph.data{scanNum} = NaN*ones(volDims);

  % check for corAnal cycles set to 0
  if params.ncycles(scanList(scanIndex)) == 0
    mrWarnDlg(sprintf('(corAnal:computeCorrelationAnalysis) !!! Scan %i has ncycles set to 0 - this needs to be set to how many cycles of the stimulus you had per scan. Skipping this scan !!!',scanList(scanIndex)));
    continue;
  end
  
  nslices = viewGet(view,'nslices',scanNum);
  for sliceNum = 1:nslices    
    
    % Analysis parameters for this scan
    junkframes = viewGet(view,'junkframes',scanNum);
    nframes = viewGet(view,'nframes',scanNum);
    
    % Load tSeries from filtered TSeries, the way loadTSeries would do it
    tSeries = filteredTSeries{scanNum}(:,:,sliceNum,:);
    
    % Reshape the tSeries
    % ATTN: added reshapeTSeries function, since loadTSeries not longer reshapes when it loads -eli
    tSeries = reshapeTSeries(tSeries);

    % check that junkframes and nframes settings are ok
    if size(tSeries,1) < (junkframes+nframes)
      mrErrorDlg(sprintf('(corAnal) Number of junkframes (%i) plus nframes (%i) should not be larger than number of volumes in scan %i',junkframes,nframes,size(tSeries,1)));
    end
    % Remove junkFrames
    tSeries = tSeries(junkframes+1:junkframes+nframes,:);
    %compute corAnal
    [coSeries,ampSeries,phSeries] = computeCoranal(tSeries,params.ncycles(scanNum),params.detrend{scanNum},params.spatialnorm{scanNum},params.trigonometricFunction{scanNum});

    switch view.viewType
      case {'Volume'}
          co.data{scanNum}(:,:,sliceNum) = reshape(coSeries,sliceDims);
          amp.data{scanNum}(:,:,sliceNum) = reshape(ampSeries,sliceDims);
          ph.data{scanNum}(:,:,sliceNum) = reshape(phSeries,sliceDims);
      case {'Surface'}
          co.data{scanNum} = coSeries;
          amp.data{scanNum} = ampSeries;
          ph.data{scanNum} = phSeries;
      case {'Flat'}
          co.data{scanNum}(:,:,sliceNum) = reshape(coSeries,sliceDims);
          amp.data{scanNum}(:,:,sliceNum) = reshape(ampSeries,sliceDims);
          ph.data{scanNum}(:,:,sliceNum) = reshape(phSeries,sliceDims);
    end
    % Update waitbar
    mrWaitBar(sliceNum/nslices,waitHandle);
  end
  mrCloseDlg(waitHandle);
end
warning('on','MATLAB:divideByZero');

%%%%%%%%%%%%%%%%%%%%%%%%
%    checkEmriParams    %
%%%%%%%%%%%%%%%%%%%%%%%%
function params = checkEmriParams(params)

% just a stub for now in case we have things we need to check
