% pRFRunSplits.m
%
%      Usage: pRFRunSplits(varargin)
%         by: akshay jagadeesh
%       date: 06/26/2017
%
function splits = pRFRunSplits(saveName, whichSplit)

splitFileName = sprintf('%s_split%d.mat', saveName, whichSplit);

% Load the split file
l = load(sprintf('Splits/%s', splitFileName));
s = l.split;
% split contains: nVoxels, scanCoords, tSeries, stim, concatInfo, prefit, pRFFitParams

% load master split file
m = load(sprintf('Splits/%s_master.mat', saveName));

%get the view
v = newView;
v = viewSet(v, 'curGroup', 'Concatenation');
v = viewSet(v, 'curScan', m.scanNum);

% Call pRFFit on the voxels
if exist('Splits/Analysis')~=7
  mkdir('Splits/Analysis');
end

x = s.scanCoords(1,:); y = s.scanCoords(2,:); z = s.scanCoords(3,:);

% Create splits struct
splits.scanCoords = s.scanCoords;

%% Crossval
if m.params.pRFFit.crossval && s.concatInfo.n <= 1
  disp('Cannot run cross-validation with less than 2 runs');
elseif m.params.pRFFit.crossval
  disp('Running leave-one-out crossvalidation...');
  cvt = tic;
  splits.crossval_r2 = nan(s.concatInfo.n, s.nVoxels);
  for ci = 1:s.concatInfo.n
    % Get the held out run's timeseries
    heldoutT = s.concatInfo.runTransition(ci,1):s.concatInfo.runTransition(ci,2);
    heldoutTSeries = s.tSeries(:,heldoutT);

    % Get all the other runs' timeseries
    trainRuns = setdiff(1:s.concatInfo.n, ci);
    trainT = setdiff(1:size(s.tSeries,2), heldoutT);
    train = s.tSeries(:, trainT);

    % Get the stimuli from the training runs
    stimTrain = {s.stim{trainRuns}};
    xPrefit = m.prefit;
    xPrefit.modelResponse = xPrefit.modelResponse(:,trainT);

    % Convert concat info struct for training runs
    trainConcatInfo = setConcatInfo(s.concatInfo,trainRuns, trainT);
    hoConcatInfo = setConcatInfo(s.concatInfo, ci, heldoutT);

    cv_r2 = nan(1,s.nVoxels);
    for vi = 1:s.nVoxels
      tv = tic;
      % Train pRF params on training timeseries
      fit = pRFFit(v, [], x(vi),y(vi),z(vi), 'tSeries', train(vi,:)', 'stim', stimTrain, 'fitTypeParams', s.pRFFitParams, ...
                   'prefit', xPrefit, 'concatInfo', trainConcatInfo);

      if ~isempty(fit)
        thisR2(vi) = fit.r2;
        thisPolarAngle(vi) = fit.polarAngle;
        thisEccentricity(vi) = fit.eccentricity;
        thisRfHalfWidth(vi) = fit.std;
        rawParams(:,vi) = fit.params(:);
        r(vi,:) = fit.r;
      end

      % Test model on heldout time series
      test = pRFFit(v, [], x(vi),y(vi),z(vi), 'tSeries', heldoutTSeries(vi,:)', 'stim', {s.stim{ci}}, 'params', fit.params,...
                    'getModelResponse=1', 'fitTypeParams', s.pRFFitParams, 'concatInfo', hoConcatInfo);

      % Calculate crossval r2.
      cv_r2(vi) = corr(heldoutTSeries(vi,:)', test.modelResponse)^2;
      disp(sprintf('Voxel [%d,%d,%d] crossvalidated r2 = %.02f. Took %.02f seconds.', x(vi),y(vi),z(vi),cv_r2(vi), toc(tv)));

    end
    % Save out crossval params
    splits.crossval_r2(ci,:) = cv_r2;
    splits.trainR2(ci,:) = thisR2;
    splits.rawParams(:,:,ci) = rawParams;
    splits.polarAngle(ci,:) = thisPolarAngle;
    splits.eccentricity(ci,:) = thisEccentricity;
    splits.rfHalfWidth(ci,:) = thisRfHalfWidth;
    splits.r(:,:,ci) = r;
  end
  toc(cvt);
  save(sprintf('Splits/Analysis/%s_Anal.mat', splitFileName(1:end-4)), 'splits');
  return
end

ncv = tic;
for i = 1:s.nVoxels

  fit = pRFFit(v, s.scanNum, x(i),y(i),z(i), 'stim', s.stim, 'concatInfo', s.concatInfo, 'prefit', m.prefit, 'fitTypeParams', s.pRFFitParams, 'tSeries', s.tSeries(i,:)', 'paramsInfo', s.paramsInfo);

  if ~isempty(fit)
    thisR2(i) = fit.r2;
    thisPolarAngle(i) = fit.polarAngle;
    thisEccentricity(i) = fit.eccentricity;
    thisRfHalfWidth(i) = fit.std;
    rawParams(:,i) = fit.params(:);
    r(i,:) = fit.r;
  end

end
toc(ncv)
splits.scanCoords = s.scanCoords;
splits.r2 = thisR2;
splits.polarAngle = thisPolarAngle;
splits.eccentricity = thisEccentricity;
splits.rfHalfWidth = thisRfHalfWidth;
splits.params = rawParams;
splits.r = r;
save(sprintf('Splits/Analysis/%s_Anal.mat', splitFileName(1:end-4)), 'splits');


function thisConcatInfo = setConcatInfo(concatInfo, runs, times)

thisConcatInfo = concatInfo;

thisConcatInfo.n = length(runs);
thisConcatInfo.whichScan = concatInfo.whichScan(times);
thisConcatInfo.whichVolume = concatInfo.whichVolume(times);
thisConcatInfo.nifti = [concatInfo.nifti(runs)];
thisConcatInfo.filename = {concatInfo.filename{runs}};
thisConcatInfo.path = {concatInfo.path{runs}};
thisConcatInfo.junkFrames = concatInfo.junkFrames(runs);
thisConcatInfo.hipassfilter = {concatInfo.hipassfilter{runs}};
thisConcatInfo.runTransition = concatInfo.runTransition(1:length(runs),:);
thisConcatInfo.totalJunkedFrames = concatInfo.totalJunkedFrames(runs);
