% t1t2.m
%
%        $Id:$ 
%      usage: t1t2(t1_fidname,t2_fidname,<verbose>,<T1T2alignment>);
%         by: justin gardner
%       date: 02/18/10
%    purpose: does t1/t2 correction saves out a file called t1t2.hdr when everything is done.
%       e.g.: t1t2('t1.fid','t2.fid');
%
%             t1_fidname can be a cell array of fid names (or nifti files) that you want to have averaged
%             t2_fidname can be a cell array of fid names (or nifti files) that you want to have averaged
%
%             This program will align and average the t1 images and the t2 images
%             Note that you will have to set a crop region to do the alignments (crop around the brain).
%
%             Then it will align using AFNI the t1 and t2 images
%
%             Then it will put up a GUI to set the threshold. Click Ok and it will save the t1t2.hdr
% 
% 
%
function retval = t1t2(t1_fidname,t2_fidname,varargin)

% check arguments
if ~any(nargin == [1 2 3 4 5 6 7 8 9 10])
  help t1t2
  return
end

% default is no t2
if nargin == 1
  t2_fidname = {};
end

% get arguments
verbose=[];
T1T2alignment=[];
getArgs(varargin,{'verbose=1','T1T2alignment=1'});

% check for correct commands to run this program
if ~checkCommands(verbose),return,end

% check filenames
[tf t1_fidname] = checkFilenames(t1_fidname,verbose);
if ~tf,return,end
[tf t2_fidname] = checkFilenames(t2_fidname,verbose);
if ~tf,return,end

% load t1 files, align images and avergae
t1 = loadImageFiles(t1_fidname,verbose);
if isempty(t1),return,end
t1 = alignImages(t1,verbose);
t1 = averageImages(t1,verbose);
%t1ic = t1t2IntensityContrastCorrection(t1,verbose);

% load t2 files, align images and average, then align to T1
t2 = loadImageFiles(t2_fidname,verbose);
t2 = alignImages(t2,verbose);
t2 = averageImages(t2,verbose);
if T1T2alignment
  t2 = alignT1toT2(t1,t2,verbose);
end
t2 = blurImage(t2,verbose);

% now set threshold and divide T1 by T2
if ~isempty(t2)
  vol = divideT1byT2(t1,t2);
  if ~isempty(vol),saveNifti(vol,'t1t2');end
end

%%%%%%%%%%%%%%%%%%%%%%%
%    checkFilenames    %
%%%%%%%%%%%%%%%%%%%%%%%
function [tf filenames] = checkFilenames(filenames,verbose)

% default return
tf = 0;

% make sure we have a cell array of names
filenames = cellArray(filenames);

% check for proper filenames
for i = 1:length(filenames)
  %  check if it is an hdr
  filename = setext(filenames{i},'hdr',0);
  if isfile(filename)
    filenames{i} = filename;
    continue
  end
  % check if it is a fid
  filenames{i} = setext(filenames{i},'fid',0);
  if ~isdir(filenames{i})
    disp(sprintf('(t1t2) Could not find t1 fid: %s',filenames{i}));
    return
  end
end

tf = 1;

%%%%%%%%%%%%%%%%%%%%%%%%
%    loadImageFiles    %
%%%%%%%%%%%%%%%%%%%%%%%%
function im = loadImageFiles(filenames,verbose)

% default return
im = [];

% make sure we have a cell array of names
filenames = cellArray(filenames);

% load the fids with fid2nifti
for i = 1:length(filenames)
  % check if it is a nifti file
  if strcmp(getext(filenames{i}),'hdr')
    im{end+1}.filename = filenames{i};
    [im{end}.d im{end}.hdr] = cbiReadNifti(filenames{i});
    continue;
  end
  % otherwise it is a fid
  im{end+1}.filename = filenames{i};
  if verbose,disppercent(-inf,sprintf('(t1t2) Loading %s',filenames{i}));,end
  [im{end}.d im{end}.hdr] = fid2nifti(filenames{i});
  if verbose,disppercent(inf);end
  if isempty(im{end}.d) 
    disp(sprintf('(t1t2) Empty data when loading %s',filenames{i}));
    return
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   intensityCorrection   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function im = t1t2IntensityContrastCorrection(im,verbose)

% get crop region
if ~isfield(im,'crop')
  im.crop = selectCropRegion(im.d);
  drawnow;
end
if verbose,disppercent(-inf,'(t1t2) intensityContrastCorrection');end
im.d = intensityContrastCorrection(im.d,im.crop);
if verbose,disppercent(inf);end


%%%%%%%%%%%%%%%%%%%%%
%    alignImages    %
%%%%%%%%%%%%%%%%%%%%%
function im = alignImages(im,verbose)

% nothing to do if there is not more than 1 image
if length(im) <= 1
  return
end

% which image to align everybody to
baseImage = 1; 
restOfImages = setdiff(1:length(im),baseImage);

% set up parameters for estMotionIter3
niters = 10;
Minitial = eye(4);
rotFlag = 1;
robust = 0;
phaseFlag = 0;
im{baseImage}.crop = selectCropRegion(im{1}.d);
drawnow;

% parameters for warping
badVal = nan;
border = 0;
interpMethod = 'linear';

% now align everything to baseImage
disppercent(-inf,'Aligning volumes');
for i = restOfImages
  % compute motion comp parameters
  M = estMotionIter3(im{baseImage}.d,im{i}.d,niters,Minitial,rotFlag,robust,phaseFlag,im{baseImage}.crop);
  % display parameters
  if verbose
    disp(sprintf('(t1t2) Translating volume %s by [%s] to align to %s',im{i}.filename,mynum2str(M(1:3,4)),im{baseImage}.filename));
    disp(sprintf('(t1t2) Rotating volume %s to align to %s',im{i}.filename,im{baseImage}.filename));
    disp(mynum2str(M(1:3,1:3)));
  end
  % warp the volume
  im{i}.d = warpAffine3(im{i}.d,M,badVal,border,interpMethod);
  disppercent(i/restOfImages);
end
disppercent(inf);

%%%%%%%%%%%%%%%%%%%%%%%
%    averageImages    %
%%%%%%%%%%%%%%%%%%%%%%%
function im = averageImages(im,verbose)

if length(im)<=0,return,end

baseImage = 1;
restOfImages = setdiff(1:length(im),baseImage);

% sum
for i = restOfImages
  im{baseImage}.d = im{baseImage}.d+im{i}.d;
end

% and divide
im{baseImage}.d = im{baseImage}.d/length(im);

% remove other images
im = im{baseImage};

%%%%%%%%%%%%%%%%%%%%
%    alignT1toT2  %%
%%%%%%%%%%%%%%%%%%%%
function t2 = alignT1toT2(t1,t2,verbose)

% AFNI command for aligning
AFNIcommand = '3dAllineate';

% filenames
t1filename = 't1.hdr';
t2filename = 't2_preAlignment.hdr';
t2alignedFilename = 't2.hdr';

% write the files out to disk
saveNifti(t1,t1filename);
saveNifti(t2,t2filename);

% run the afni command
disppercent(-inf,sprintf('(t1t2) Aligning T1 to T2 using AFNI:%s',AFNIcommand));
system(sprintf('%s -verb -warp shift_rotate -base %s -source %s -prefix %s',AFNIcommand,t1filename,t2filename,t2alignedFilename));
disppercent(inf);

% read back in aligned volume
[t2.d t2.hdr] = cbiReadNifti(t2alignedFilename);

% remove temporary files
%  delete(setext(t1filename,'*'));
%  delete(setext(t2filename,'*'));
%  delete(setext(t2alignedFilename,'*'));

%%%%%%%%%%%%%%%%%%%
%    saveNifti    %
%%%%%%%%%%%%%%%%%%%
function saveNifti(im,name,verbose)

name = setext(name,'hdr');
cbiWriteNifti(name,im.d,im.hdr);

%%%%%%%%%%%%%%%%%%%%%%%
%    checkCommands    %
%%%%%%%%%%%%%%%%%%%%%%%
function retval = checkCommands(verbose)

retval = 1;
hailstr = 't1t2';

% commands to check
commandNames = {'3dAllineate'};
helpFlag = {'-help'};
for i = 1:length(commandNames)
  if verbose,disp(sprintf('(%s) Checking for existence of command: %s',hailstr,commandNames{i}));end
  % suse which to tell if we have the command
  [commandStatus commandRetval] = system(sprintf('which %s',commandNames{i}));
  % check for commandStatus error
  if commandStatus~=0
    disp(sprintf('(%s) Could not find command: %s',hailstr,commandNames{i}));
    disp(sprintf('            See http://gru.brain.riken.jp/doku.php?id=gru:segmentationeasy for help setting up your computer'));
    retval = 0;
    return
  end
  % run the command to see what happens
  [commandStatus commandRetval] = system(sprintf('%s %s',commandNames{i},helpFlag{i}));
  % check for commandStatus error
  if commandStatus>1
    disp(commandRetval);
    disp(sprintf('(%s) Found command: %s, but could not run (possibly missing fink library?)',hailstr,commandNames{i}));
    disp(sprintf('            See http://gru.brain.riken.jp/doku.php?id=gru:segmentationeasy for help setting up your computer'));

    retval = 0;
    return
  end
end

%%%%%%%%%%%%%%%%%%%%%%
%    divideT1byT2    %
%%%%%%%%%%%%%%%%%%%%%%
function vol = divideT1byT2(t1,t2,verbose)

vol = [];

% some global variables for displaying
global t1t2fig;

% display figure
t1t2fig.f = smartfig('t1t2');
t1t2fig.f2 = smartfig('t1t2final');


t1t2fig.t1 = t1;
t1t2fig.t2 = t2;

maxSlice = size(t1.d,3);
midSlice = round(maxSlice/2);

% set up params
paramsInfo{1} = {'threshold',1,'numeric=1','incdec=[-0.1 0.1]','minmax=[0 inf]','callback',@displayT1T2,'passParams=1','Controls the threshold of the T2 image that creates the mask. Set this to make as clean a mask as possible'};
paramsInfo{end+1} = {'sliceNum',midSlice,'numeric=1','incdec=[-1 1]',sprintf('minmax=[1 %i]',maxSlice),'callback',@displayT1T2,'passParams=1','Which slice to display'};
paramsInfo{end+1} = {'blurLevel',t2.blurLevel,'numeric=1','incdec=[-1 1]','minmax=[0 inf]','round=1','callback',@displayT1T2,'passParams=1','Which slice to display'};
paramsInfo{end+1} = {'gamma',0.4,'numeric=1','incdec=[-0.1 0.1]','minmax=[0 inf]','callback',@displayT1T2,'passParams=1','Changes the display gamma of the T1T2 image'};

		
% get default params to display initial image
params = mrParamsDefault(paramsInfo);
displayT1T2(params);

% now put up dialog box 
params = mrParamsDialog(paramsInfo,'Set threshold');

% clean up
close(t1t2fig.f);
close(t1t2fig.f2);
drawnow;
clear global t1t2fig;

% user cancel
if isempty(params),return,end

vol = computeT1T2(t1,t2,params);


%%%%%%%%%%%%%%%%%%%%%
%%   computeT1T2   %%
%%%%%%%%%%%%%%%%%%%%%
function [t1t2 mask t2] = computeT1T2(t1,t2,params,sliceNum)

% see if we have to blur T2
if params.blurLevel ~= t2.blurLevel
  t2 = blurImage(t2,1,params.blurLevel);
end

% default to all slices
if ieNotDefined('sliceNum'),sliceNum = 1:size(t1,3);end

% get mask
mask = t2.d(:,:,sliceNum)>params.threshold;
divmask = mask.*t2.blurd(:,:,sliceNum);
divmask(divmask==0) = inf;

% do the division
t1t2.d = t1.d(:,:,sliceNum)./divmask;
t1t2.hdr = t1.hdr;

%%%%%%%%%%%%%%%%%%%%
%    displayT1T2   %
%%%%%%%%%%%%%%%%%%%%
function displayT1T2(params)

global t1t2fig;

% display T1
figure(t1t2fig.f);
subplot(1,3,1);
imagesc(t1t2fig.t1.d(:,:,params.sliceNum));
axis square;
title('T1');
drawnow

% compute division
if ~isfield(t1t2fig,'params') || (t1t2fig.params.threshold ~= params.threshold) || (t1t2fig.params.blurLevel ~= params.blurLevel) || (t1t2fig.params.sliceNum ~= params.sliceNum)
  [t1t2fig.t1t2 t1t2fig.mask t1t2fig.t2] = computeT1T2(t1t2fig.t1,t1t2fig.t2,params,params.sliceNum);
end

% display T2
figure(t1t2fig.f);
subplot(1,3,2);
imagesc(t1t2fig.t2.blurd(:,:,params.sliceNum));
axis square;
title('T2');

% display mask
subplot(1,3,3);
figure(t1t2fig.f);
imagesc(t1t2fig.mask);
axis square;
title('Mask');

% set colormap
colormap(gray);

% draw t1/t2
figure(t1t2fig.f2);
imageg(t1t2fig.t1t2.d,params.gamma);
title('T1/T2');

colormap(gray);

%%%%%%%%%%%%%%%%%%%
%%   blurImage   %%
%%%%%%%%%%%%%%%%%%%
function img = blurImage(img,verbose,nIter)

% nIter and blurLevel are same thing
if ieNotDefined('nIter'),nIter = 1;end
img.blurLevel = nIter;

% no bluring
if nIter == 0
  img.blurd = img.d;
  return
end

% pre allocate memory
if ~isfield(img,'blurd')
  img.blurd = zeros(size(img.d));
end

% do bluring
if verbose,disppercent(-inf,'(t1t2) Blurring T2 image');end
for i = 1:size(img.d,3)
  img.blurd(:,:,i) = blur(img.d(:,:,i),nIter);
  if verbose,disppercent(i/size(img.d,3));end
end
if verbose,disppercent(inf);end

