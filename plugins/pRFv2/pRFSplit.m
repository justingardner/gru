% pRFSplit.m
%
%       usage: pRF(v, scanNum, params, x,y,z,n, fit)
%          by: akshay jagadeesh
%        date: 06/20/2017
%     purpose: Split pRF data into chunks of voxels, so that they can be more easily computed.
%
%
function splits = pRFSplit(v, scanNum, params, x,y,z, n, fit, overlays)

%% Parameters

vnum = 5; % how many voxels to use to estimate run time
splitTime = params.pRFFit.splitTime; % how many minutes to use per split

%% Remove old splits
% Clean up by deleting the splits folder
disp('Clean up! Deleting splits folder before start..');
system('rm -r Splits');

%% Get current user and current session dir
curPath = pwd;
sherlockSessionPath = ['/share/PI/jlg/' curPath(findstr(curPath, 'data'):end)];
suid = getsuid;

% Set split directory and scripts directory
splitDir = 'Splits';
scriptsDir = 'Splits/Scripts';
% blockSize = ceil(n/numSplits);
whichSplit = 1;
scanDims = viewGet(v, 'scanDims');
pRFAnal = overlays.pRFAnal;

% Make Splits directory if it doesn't exist 
if exist(splitDir) ~= 7
  mkdir(splitDir);
end
if exist(scriptsDir) ~=7
  mkdir(scriptsDir);
end

system('echo "#\!/bin/bash" >! "Splits/Scripts/runAll.sh"');

% Load up the entire timeseries by making a temporary ROI
maxBlockSize = mrGetPref('maxBlockSize');
loadROI = makeEmptyROI(v, 'scanNum', scanNum, 'groupNum', params.groupName);
loadROI.coords(1,:) = x;
loadROI.coords(2,:) = y;
loadROI.coords(3,:) = z;
loadROI = loadROITSeries(v, loadROI, scanNum, params.groupName);

disp(sprintf('Running on %i voxels to calculate runtime',vnum));
% Run on a few voxels to calculate runtime -- use this later to calculate number of splits
tic
for vi = 1:vnum
  xi = loadROI.scanCoords(1,vi);
  yi = loadROI.scanCoords(2,vi);
  zi = loadROI.scanCoords(3,vi);
  pRFFit(v, scanNum, xi, yi, zi, 'stim', fit.stim, 'concatInfo', fit.concatInfo, 'prefit', fit.prefit, 'fitTypeParams', params.pRFFit, 'tSeries', loadROI.tSeries(vi,:)', 'paramsInfo', fit.paramsInfo); 
end
runtime = toc;
blockSize = round((splitTime*60/runtime)*vnum);
disp(sprintf('5 voxel runtime was estimated to be %03.1f seconds: splitting into %i minute chunks using a block size of %1.0f voxels.',runtime,splitTime,blockSize));

for blockStart = 1:blockSize:n
  blockEnd = min(blockStart+blockSize-1,n);
  blockSize = blockEnd-blockStart+1;

  % Make a temporary ROI using the coordinates specified.
%   loadROI = makeEmptyROI(v,'scanNum',scanNum,'groupNum',params.groupName);
%   loadROI.coords(1,1:blockSize) = x(blockStart:blockEnd);
%   loadROI.coords(2,1:blockSize) = y(blockStart:blockEnd);
%   loadROI.coords(3,1:blockSize) = z(blockStart:blockEnd);

  % Load time series
  %loadROI = loadROITSeries(v, loadROI, scanNum, params.groupName);

  % Grab the timeseries
  tSeries = loadROI.tSeries(blockStart:blockEnd,:);

  % reorder x,y,z coordinates since they can get scrambled in loadROITSeries
  x(blockStart:blockEnd) = loadROI.scanCoords(1,blockStart:blockEnd);
  y(blockStart:blockEnd) = loadROI.scanCoords(2,blockStart:blockEnd);
  z(blockStart:blockEnd) = loadROI.scanCoords(3,blockStart:blockEnd);
  % Keep the linear coords
  pRFAnal.d{scanNum}.linearCoords = [pRFAnal.d{scanNum}.linearCoords sub2ind(scanDims,x(blockStart:blockEnd),y(blockStart:blockEnd),z(blockStart:blockEnd))];

  % Add all the needed fields to a split struct
  %split.v = v;
  split.nVoxels = blockSize;
  split.scanCoords = [x; y; z];  
  split.tSeries = loadROI.tSeries;
  split.stim = fit.stim;
  split.concatInfo = fit.concatInfo;
  %split.prefit = fit.prefit;
  split.pRFFitParams = params.pRFFit;
  split.paramsInfo = fit.paramsInfo;
  split.scanNum = scanNum;

  % Save split struct to the specified directory
  split.splitName = sprintf('split%d',whichSplit);
  filename = sprintf('%s_%s', params.saveName, split.splitName);
  saveFile = sprintf('%s/%s.mat', splitDir, filename);
  split.file = saveFile;
  split.num = whichSplit;
  save(saveFile, 'split');
  disp(sprintf('Data split %d saved to %s', whichSplit, saveFile));

  % Call bash script to output a .sbatch file
  %disp('Generating bash scripts');
  %system(sprintf('sh ~/proj/mrTools/mrLoadRet/Plugin/pRF/Sherlock/generateBatchScripts.sh "%s" "%s" "%s" "%d"',params.saveName,sherlockSessionPath, suid, whichSplit));

  splits{whichSplit} = split;
  whichSplit = whichSplit+1;

end

% Save master split struct locally
prefit = fit.prefit;
scanCoords = loadROI.scanCoords;
disp('Saving master struct');
save(sprintf('Splits/%s_master.mat', params.saveName), 'fit', 'x', 'y', 'z', 'scanNum', 'overlays', 'pRFAnal', 'v', 'params', 'sherlockSessionPath', 'suid', 'prefit','scanCoords');
