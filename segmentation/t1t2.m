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
%             e.g.: t1t2({'hires3d01','hires3D03'},'hires3D02')
%             t2_fidname can be a cell array of fid names (or nifti files) that you want to have averaged
%
%             This program will align and average the t1 images and the t2 images
%             Note that you will have to set a crop region to do the alignments (crop around the brain).
%
%             Then it will align using AFNI the t1 and t2 images - if you set the argument
%             t1t2('t1','t2','T1T2alignment=1');
%
%             Then it will put up a GUI to set the threshold and do various bluring or filtering
%             of the T2 image.
%
%             Note that you can also do normalization by dividing by a blurred version 
%             of the same T1. This will save an image called t1t1.hdr
%
%             t1t2('t1');
%
function retval = t1t2(t1_fidname,t2_fidname,varargin)

% check arguments
if nargin == 0
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
getArgs(varargin,{'verbose=1','T1T2alignment=0','roi=[]','ref=[]'});

% load ROI
if ~isempty(roi),loadROI(roi);end

% load ref image
if ~isempty(ref),loadRef(ref);end

% check for correct commands to run this program
if ~isempty(t2_fidname) && T1T2alignment
  if ~checkCommands(verbose),return,end
end

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
% save to disk if this is an averaged one
if isfield(t1,'averaged') && t1.averaged
  saveNifti(t1,'t1.hdr',verbose);
end

% load t2 files, align images and average
t2 = loadImageFiles(t2_fidname,verbose);
t2 = alignImages(t2,verbose);
t2 = averageImages(t2,verbose);
if isfield(t2,'averaged') && t2.averaged
  saveNifti(t2,'t2.hdr',verbose);
end

% now do T1/T2 alignment using AFNI
if T1T2alignment && ~isempty(t2)
  t2 = alignT1toT2(t1,t2,verbose);
end

% now set threshold and divide T1 by T2
if ~isempty(t2)
  [vol processingList] = divideT1byT2(t1,t2,1);
  if ~isempty(vol),saveNifti(vol,'t1t2',verbose,processingList);end
else
  [vol processingList] = divideT1byT2(t1,t1,3);
  if ~isempty(vol),saveNifti(vol,'t1t1',verbose,processingList);end
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
  disp(sprintf('(t1t2:loadImageFiles) Loading %s',filenames{i}));
  im{end+1}.filename = filenames{i};
  [im{end}.d im{end}.hdr] = mlrImageLoad(filenames{i},'orient=LPI','nifti=1');
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
disp(sprintf('(t1t2:alignImages) Select crop region around the brain for alignment. This may take a few moments to put up the figure...'));
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
  [trash indexi] = ismember(i,restOfImages);
  disppercent(indexi/length(restOfImages));
end
disppercent(inf);

%%%%%%%%%%%%%%%%%%%%%%%
%    averageImages    %
%%%%%%%%%%%%%%%%%%%%%%%
function im = averageImages(im,verbose)


nImages = length(im);
if nImages<=0,return,end

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

if nImages > 1
  im.averaged = true;
else
  im.averaged = false;
end

%%%%%%%%%%%%%%%%%%%%
%    alignT1toT2  %%
%%%%%%%%%%%%%%%%%%%%
function t2 = alignT1toT2(t1,t2,verbose)

% AFNI command for aligning
AFNIcommand = '3dAllineate';

% filenames
t1filename = 't1_forAlignment.hdr';
t2filename = 't2_preAlignment.hdr';
t2alignedFilename = 't2.hdr';

% write the files out to disk
disp(sprintf('(t1t2:alignT1toT2) Saving temporary files that are used by AFNI to do alignment'));
saveNifti(t1,t1filename,verbose);
saveNifti(t2,t2filename,verbose);

% run the afni command
disppercent(-inf,sprintf('(t1t2) Aligning T1 to T2 using AFNI:%s',AFNIcommand));
system(sprintf('%s -verb -warp shift_rotate -base %s -source %s -prefix %s',AFNIcommand,t1filename,t2filename,t2alignedFilename));
disppercent(inf);

% read back in aligned volume
[t2.d t2.hdr] = cbiReadNifti(t2alignedFilename);

% remove temporary files
delete(setext(t1filename,'*'));
delete(setext(t2filename,'*'));
%  delete(setext(t2alignedFilename,'*'));

%%%%%%%%%%%%%%%%%%%
%    saveNifti    %
%%%%%%%%%%%%%%%%%%%
function saveNifti(im,name,verbose,txt)

if nargin < 4,txt = -1;end

if verbose,disp(sprintf('(t1t2) Saving %s',name));end
name = setext(name,'hdr');

if isfile(name)
  % ask the user what to do
  saveMethodTypes = {'Rename','Overwrite',};
  paramsInfo = {{'saveMethod',saveMethodTypes,'Choose how you want to save the volume'}};
  params = mrParamsDialog(paramsInfo,sprintf('%s already exists',name));
  if isempty(params),return,end
  if strcmp(params.saveMethod,'Rename')
    % put up a dialog to get new filename
    [name pathStr] = uiputfile({'*.hdr'},'Enter new name to save volume',name);
    % user hit cancel
    if isequal(name,0),return, end
    % otherwise accept the filename, make sure it has an .hdr extension
    name = setext(name,'hdr');
  end
end

% write out image
cbiWriteNifti(name,im.d,im.hdr);

% write out processing list
if ~isequal(txt,-1)
  f = fopen(setext(name,'txt'),'w');
  fprintf(f,txt);
  fclose(f);
end
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
function [vol processingList] = divideT1byT2(t1,t2,blurLevel,verbose)

vol = [];
processingList = '';

% some global variables for displaying
global gt1t2;

% set range - this is what the output will be scaled to
gt1t2.range = 256;

% display figure
gt1t2.f = smartfig('t1t2');
gt1t2.f2 = smartfig('t1t2final');
gt1t2.histfig = smartfig('t1t2histfig');
if isfield(gt1t2,'ref') && ~isempty(gt1t2.ref)
  gt1t2.reffig = smartfig('t1t2reffig');
end
% set fields
gt1t2.t1 = t1;
gt1t2.t2 = t2;
gt1t2.params = [];

% set mouse handler
set(gt1t2.f,'WindowButtonDownFcn',@setMaxFillCenter);

gt1t2.roiPath = '.';

maxSlice = size(t1.d,3);
midSlice = round(maxSlice/2);

initThreshold = 10;%round(median(gt1t2.t2.d(:)));
% set up params
paramsInfo{1} = {'sliceOrientation',{'saggital','coronal','axial'},'callback',@displayT1T2,'passParams=1','Which slice orientation to display'};
paramsInfo{end+1} = {'sliceNum',midSlice,'numeric=1','incdec=[-1 1]',sprintf('minmax=[1 %i]',maxSlice),'callback',@displayT1T2,'passParams=1','Which slice to display'};
paramsInfo{end+1} = {'dispInMlrVol',0,'type=pushbutton','callback',@dispInMlrVol,'passParams=1','buttonString=Display in mlrVol','Display output volume in mlrVol'};
paramsInfo{end+1} = {'dispProcessing',0,'type=pushbutton','callback',@dispProcessing,'passParams=1','buttonString=Display processing list','Display list of all processing done on T2 volume'};
paramsInfo{end+1} = {'setRef',0,'type=pushbutton','callback',@setRefVolume,'passParams=1','buttonString=Set reference','Set the current T1/T2 as a reference which can be viewed when you display in mlrVol. Hold shift down and hit button to clear reference.'};
paramsInfo{end+1} = {'loadROI',0,'type=pushbutton','callback',@loadROI,'passParams=1','buttonString=Load ROI','Load ROI for displaying localized histogram'};
paramsInfo{end+1} = {'revert',0,'type=pushbutton','passParams=1','callback',@revert,'buttonString=Revert T2','Revert to original T1/T2 without any blur or other filtering applied'};
paramsInfo{end+1} = {'threshold',initThreshold,'numeric=1','incdec=[-0.1 0.1]','minmax=[0 inf]','callback',@displayT1T2,'passParams=1','Controls the threshold of the T2 image that creates the mask. Set this to make as clean a mask as possible'};
paramsInfo{end+1} = {'blurLevel',1,'numeric=1','incdec=[-1 1]','minmax=[0 inf]','round=1','How much to blur the image'};
paramsInfo{end+1} = {'blur',0,'type=pushbutton','callback',@doblur,'passParams=1','buttonString=2D blur T2','Blur T2 image, using blurLevel above'};
paramsInfo{end+1} = {'blurWidth',1,'numeric=1','incdec=[-1 1]','minmax=[0 inf]','Standard deviation in pixels of the 3D gaussian to blur with when the button below is pressed'};
paramsInfo{end+1} = {'blurGauss',0,'type=pushbutton','callback',@doGaussianBlur,'passParams=1','buttonString=3D gaussian blur T2','Blur T2 image, using a 3D gaussian with standard deviation as set above'};
paramsInfo{end+1} = {'maxFill',0,'type=pushbutton','callback',@t1t2MaxFill,'passParams=1','buttonString=Max fill T2','Max fills the volume by replacing each voxel by the maximum of its neighbors within a distance as specified by maxfillRadius. If maxfillMask is set then those voxel outside the current mask will get the maximum of its neighbors from the maxfilled volume - which propogats the max value into the unmasked places. The algorithm starts at the maxFillCenter and moves outward. Click on a point in the T2 volume to set this center point'};
paramsInfo{end+1} = {'maxFillCenter',round(size(t1.d)/3),'See above','round=1'};
paramsInfo{end+1} = {'maxFillRadius',1.5,'See above. 1.5 will get all immediate neighbors'};
paramsInfo{end+1} = {'maxFillMask',1,'type=checkbox','See above'};
% get default params to display initial image
params = mrParamsDefault(paramsInfo);
displayT1T2(params);

% now put up dialog box 
params = mrParamsDialog(paramsInfo,'T1T2 controls','modal=1');

% clean up
close(gt1t2.f);
close(gt1t2.f2);
close(gt1t2.histfig);
if isfield(gt1t2,'reffig')
  close(gt1t2.reffig);
end
drawnow;
t2 = gt1t2.t2;

% user cancel
if isempty(params),return,end

% compute T1/T2
vol = computeT1T2(t1,t2,params);

% return processing list
if isfield(gt1t2.t2,'processingList')
  processingList = gt1t2.t2.processingList;
end

clear global gt1t2;

%%%%%%%%%%%%%%%%%%%%%%
%    setRefVolume    %
%%%%%%%%%%%%%%%%%%%%%%
function retval = setRefVolume(params)

retval = 0;
global gt1t2;

if any(strcmp(get(gcf,'CurrentModifier'),'shift'))
  disp(sprintf('(t1t2) Clearing reference volume'));
  if isfield(gt1t2,'refvol')
    gt1t2= rmfield(gt1t2,'refvol');
  end
else
  disp(sprintf('(t1t2) Setting reference volume for viewing in mlrVol'));
  gt1t2.refvol = gt1t2.t1t2;
  gt1t2.refvol.processingList = sprintf('%sThreshold %0.1f\n',gt1t2.t2.processingList,params.threshold);
end
%%%%%%%%%%%%%%%%%%%%%%%%
%    dispProcessing    %
%%%%%%%%%%%%%%%%%%%%%%%%
function retval = dispProcessing(params)

retval = 0;
global gt1t2;
disp(repmat('=',1,60));
if ~isempty(gt1t2.t2.processingList)
  disp(sprintf('%s',gt1t2.t2.processingList(1:end-1)));
else
  disp(sprintf('No processing has been done'));
end
disp(repmat('=',1,60));

if isfield(gt1t2,'refvol')
  disp(sprintf('Reference volume processing'));
  if ~isempty(gt1t2.refvol.processingList)
    disp(sprintf('%s',gt1t2.refvol.processingList(1:end-1)));
  else
    disp(sprintf('No processing has been done'));
  end
  disp(repmat('=',1,60));
end

%%%%%%%%%%%%%%%%%%%%%%%%
%    doGaussianBlur    %
%%%%%%%%%%%%%%%%%%%%%%%%
function retval = doGaussianBlur(params)


% get some parameters
retval = 0;
global gt1t2;
dim = size(gt1t2.t2.blurd);
sigma = params.blurWidth;

disppercent(-inf,sprintf('(t1t2:doGaussianBlur) Computing 3D gaussian blur with standard deviation %0.2f',params.blurWidth));

% compute the gaussian filter
[x y z] = ndgrid(0:dim(1)-1,0:dim(2)-1,0:dim(3)-1);
x = x-dim(1)/2;
y = y-dim(2)/2;
z = z-dim(3)/2;
g = exp(-(x.^2 + y.^2 + z.^2)/(sigma^2));
g = ifftshift(g);

% test
%uhm = zeros(dim);
%uhm(2,1,5) = 1; (gaussian should be centered on this point)
%duh = ifftn(fftn(uhm).*fftn(g));

% compute convolution
gt1t2.t2.blurd = ifftn(fftn(gt1t2.t2.blurd).*fftn(g));

disppercent(inf);

% redisplay
[gt1t2.t1t2 gt1t2.mask gt1t2.t2] = computeT1T2(gt1t2.t1,gt1t2.t2,params);
displayT1T2(params);

% update string of what we have done
gt1t2.t2.processingList = sprintf('%sGaussian blur with standard deviation %0.2f\n',gt1t2.t2.processingList,params.blurWidth);
%%%%%%%%%%%%%%
%    blur    %
%%%%%%%%%%%%%%
function retval = doblur(params)

retval = 0;
global gt1t2;

if params.blurLevel == 0
  disp(sprintf('(t1t2:doblur) Blur level set to 0, not do anything'));
  return
end

% pre allocate memory
if ~isfield(gt1t2.t2,'blurd')
  gt1t2.t2.blurd = zeros(size(gt1t2.t2.d));
end

% do bluring
verbose = 1;
if verbose,disppercent(-inf,'(t1t2) Blurring T2 image');end
for i = 1:size(gt1t2.t2.blurd,3)
  gt1t2.t2.blurd(:,:,i) = blur(gt1t2.t2.blurd(:,:,i),params.blurLevel);
  if verbose,disppercent(i/size(gt1t2.t2.blurd,3));end
end
if verbose,disppercent(inf);end

% redisplay
[gt1t2.t1t2 gt1t2.mask gt1t2.t2] = computeT1T2(gt1t2.t1,gt1t2.t2,params);
displayT1T2(params);

% update string of what we have done
gt1t2.t2.processingList = sprintf('%sBlur level %i\n',gt1t2.t2.processingList,params.blurLevel);

%%%%%%%%%%%%%%%%
%    revert    %
%%%%%%%%%%%%%%%%
function retval = revert(params)

retval = 0;
global gt1t2;
gt1t2.t2.blurd = gt1t2.t2.d;
[gt1t2.t1t2 gt1t2.mask gt1t2.t2] = computeT1T2(gt1t2.t1,gt1t2.t2,params);
displayT1T2(params);

% update string to default
gt1t2.t2.processingList = '';

%%%%%%%%%%%%%%%%%%%%%%
%    dispInMlrVol    %
%%%%%%%%%%%%%%%%%%%%%%
function retval = dispInMlrVol(params)

retval = 0;

global gt1t2;
vol = computeT1T2(gt1t2.t1,gt1t2.t2,params);
if isfield(gt1t2,'refvol')
  mlrVol(vol.d,vol.hdr,gt1t2.refvol.d,gt1t2.refvol.hdr,'toggleOverlay=0','showControls=0');
else
  mlrVol(vol.d,vol.hdr,'showControls=0');
end

%%%%%%%%%%%%%%%%%%%%%%%
%%   paramsChanged   %%
%%%%%%%%%%%%%%%%%%%%%%%
function tf = paramsChanged(params1,params2)

tf = 1;
% if either is empty need to recompute
if isempty(params1),return,end
if isempty(params2),return,end

% if these parameters have changed, need to recompute
if ~isequal(params1.threshold,params2.threshold),return,end
if ~isequal(params1.blurLevel,params2.blurLevel),return,end

% otherwise they are effectively the same
tf = 0;

%%%%%%%%%%%%%%%%%%%%%
%%   computeT1T2   %%
%%%%%%%%%%%%%%%%%%%%%
function [t1t2 mask t2] = computeT1T2(t1,t2,params,sliceNum)

global gt1t2;

% default to all slices
if ieNotDefined('sliceNum'),sliceNum = 1:size(t1.d,3);end

disppercent(-inf,'(t1t2) Computing t1/t2');
% get mask
mask = t2.d(:,:,sliceNum)>params.threshold;
divmask = mask.*t2.blurd(:,:,sliceNum);
divmask(divmask==0) = inf;

% do the division
t1t2.d = t1.d(:,:,sliceNum)./divmask;
t1t2.hdr = t1.hdr;

% now compute normalization factors
t1t2.min = min(t1t2.d(mask));
t1t2.max = max(t1t2.d(mask));

% normalize
t1t2.d = gt1t2.range*(t1t2.d-t1t2.min)/(t1t2.max-t1t2.min);
disppercent(inf);

%%%%%%%%%%%%%%%%%%%%
%    displayT1T2   %
%%%%%%%%%%%%%%%%%%%%
function displayT1T2(params)

global gt1t2;

sliceOrientation = find(strcmp(params.sliceOrientation,{'axial','coronal','saggital'}));

% display T1
figure(gt1t2.f);
subplot(1,3,1);
switch sliceOrientation
  case {1}
   imagesc(flipud(gt1t2.t1.d(:,:,params.sliceNum)'));
  case {2}
   imagesc(flipud(squeeze(gt1t2.t1.d(:,params.sliceNum,:))'));
  case {3}
   imagesc(flipud(squeeze(gt1t2.t1.d(params.sliceNum,:,:))'));
end

axis square;
axis off;
title('T1');
drawnow

% set current t2 image
if ~isfield(gt1t2.t2,'blurd')
  gt1t2.t2.blurd = gt1t2.t2.d;
  gt1t2.t2.processingList = '';
end

% compute division
if paramsChanged(gt1t2.params,params)
  [gt1t2.t1t2 gt1t2.mask gt1t2.t2] = computeT1T2(gt1t2.t1,gt1t2.t2,params);
end
gt1t2.params = params;

% display T2
figure(gt1t2.f);
subplot(1,3,2);
switch sliceOrientation
 case {1}
  imagesc(flipud(gt1t2.t2.blurd(:,:,params.sliceNum)'));
 case {2}
  imagesc(flipud(squeeze(gt1t2.t2.blurd(:,params.sliceNum,:))'));
 case {3}
  imagesc(flipud(squeeze(gt1t2.t2.blurd(params.sliceNum,:,:))'));
end

axis square;
axis off;
title('T2');

% display mask
subplot(1,3,3);
figure(gt1t2.f);
switch sliceOrientation
 case {1}
  maskSlice = flipud(squeeze(gt1t2.mask(:,:,params.sliceNum))');
 case {2}
  maskSlice = flipud(squeeze(gt1t2.mask(:,params.sliceNum,:))');
 case {3}
  maskSlice = flipud(squeeze(gt1t2.mask(params.sliceNum,:,:))');
end
imagesc(maskSlice);

axis square;
axis off;
title('Mask');

% set colormap
colormap(gray);

% get slice
switch sliceOrientation
 case {1}
  imSlice = gt1t2.t1t2.d(:,:,params.sliceNum);
 case {2}
  imSlice = gt1t2.t1t2.d(:,params.sliceNum,:);
 case {3}
  imSlice = squeeze(gt1t2.t1t2.d(params.sliceNum,:,:));
end
imSlice = squeeze(imSlice);

% compute and display histogram
figure(gt1t2.histfig);cla;

if ~isfield(gt1t2,'roi') || ~isroi(gt1t2.roi)
  vals = imSlice(find(maskSlice));
  vals = vals(~isnan(vals));
  vals = vals(~isinf(vals));
  if ~isempty(vals)
    myhist(vals,100);
  end
  imDisplaySlice = imSlice;
else
  % load the voxels 
  % get the linear coords
  linearCoords = sub2ind(size(gt1t2.t1t2.d),gt1t2.roi.coords(1,:),gt1t2.roi.coords(2,:),gt1t2.roi.coords(3,:));
  vals = gt1t2.t1t2.d(linearCoords);
  subplot(2,1,1);
  H = myhist(vals,100);
  [m s] = mixgauss(vals,2);
  title(sprintf('Histogram for ROI: %s',gt1t2.roi.name));
  ylabel('N voxels');
  xlabel('Pixel intensity');
  a = axis;
  subplot(2,1,2);
  myhist(imSlice(:),H.bins);
  % color slice according to mean and standard deviation of gray lump
  imDisplaySlice = colorMatchingVoxels(imSlice,m(1),s(1));
end
% set title for distribution across slice
title(sprintf('Histogram of values across slice in mask'));
ylabel('N voxels');
xlabel('Pixel intensity');

figure(gt1t2.f2);
imagesc(flipud(imDisplaySlice'));
title('T1/T2');
axis square;
axis off;

colormap(gray(gt1t2.range));

% now display ref image if it exists
if isfield(gt1t2,'ref') && ~isempty(gt1t2.ref)
  % get the slice
  switch sliceOrientation
   case {1}
    refSlice = gt1t2.ref(:,:,params.sliceNum);
   case {2}
    refSlice = gt1t2.ref(:,params.sliceNum,:);
   case {3}
    refSlice = squeeze(gt1t2.ref(params.sliceNum,:,:));
  end
  refSlice = squeeze(refSlice);

  
  figure(gt1t2.reffig);
  clf;
  subplot(1,2,1);
  imagesc(flipud(refSlice'));
  title(gt1t2.refname);
  axis square;
  axis off;
  colormap(gray(gt1t2.range));
  subplot(1,2,2);
  plot(imSlice(:),refSlice(:),'k.');
  xlabel('t1t2 pixel intensity');
  ylabel('ref pixel intensity');
  axis square;
end

%%%%%%%%%%%%%%%
%%   check   %%
%%%%%%%%%%%%%%%
function check

% test outtput to compare different methods
d1 = cbiReadNifti('t1t2.hdr');d1 = (d1-min(d1(:)))/(max(d1(:))-min(d1(:)));
d2 = cbiReadNifti('t1t1.hdr');d2 = (d2-min(d2(:)))/(max(d2(:))-min(d2(:)));
d3 = cbiReadNifti('t1.hdr');d3 = (d3-min(d3(:)))/(max(d3(:))-min(d3(:)));

d = d1;
d(:,:,:,2) = d2;
d(:,:,:,3) = d3;
mlrDisplayEPI(d);

%%%%%%%%%%%%%%%%%
%%   loadROI   %%
%%%%%%%%%%%%%%%%%
function val = loadROI(params)

val = [];
global gt1t2;

% see if we are passed in a roi name
justLoad = 0;
if isstr(params)
  justLoad = 1;
  fullPath = params;
else
  % get the path to the roi the user wants to load
  fullPath = getPathStrDialog(gt1t2.roiPath,'Load ROI','*.mat');
end

% if no path, return
if isempty(fullPath),return,end

% check for file
if ~isfile(fullPath)
  disp(sprintf('(t1t2) Could not load roi %s',fullPath));
  return;
end

% load file, and get roi variable
roi = load(fullPath);
roiname = fieldnames(roi);
if length(roiname) > 1
  disp(sprintf('(t1t2) File %s contains more than one variable',fullPath));
  return
end
roi = roi.(roiname{1});

% check if it is a roi
if ~isroi(roi)
  disp(sprintf('(t1t2) File %s is not an ROI',fullPath));
  return
end

% save in global variable and redraw
gt1t2.roi = roi;

disp(sprintf('(t1t2) Loaded ROI %s',fullPath));
if ~justLoad
  % redraw display
  displayT1T2(params);
end

%%%%%%%%%%%%%%%%%
%%   loadRef   %%
%%%%%%%%%%%%%%%%%
function loadRef(filename)

global gt1t2;

% load nifti pair
filename = setext(filename,'hdr');

% check for file
if ~isfile(filename)
  disp(sprintf('(t1t2) Could not load ref image %s',filename));
  return;
end

% load file
disp(sprintf('(t1t2:loadRef) Loading ref image %s',filename));
gt1t2.ref = cbiReadNifti(filename);
gt1t2.refname = filename;

%%%%%%%%%%%%%%%%%%
%%   mixgauss   %%
%%%%%%%%%%%%%%%%%%
function [m s] = mixgauss(x,k)

% make x a row vector and replicate k times for use later
x = x(:)';
xk = ones(k,1)*x;

% get number of elements
n = length(x);

% initialize which distribution each point lies in
whichDist = ceil(rand(1,n)*k);

% initialize mean and stdard deviations
for i = 1:k
  m(i) = mean(x(whichDist==i));
  s(i) = std(x(whichDist==i));
end

% choose how little mean has to change to stop algorithm
epsilon = 0.0000001;

oldm = inf;
disppercent(-inf,'(mixgauss) Doing k-means');
while max(abs(oldm(:)-m(:)))>epsilon
  pWhichDist = [];
  % now compute the probability each point lies in each distribution
  for i = 1:k
    pWhichDist(i,:) = normpdf(x,m(i),s(i));
  end
  % and normalize
  pWhichDist = pWhichDist./(ones(k,1)*sum(pWhichDist,1));
  % uncomment the next line if you want to do
  % a soft k-means (EM)
  pWhichDist = pWhichDist>0.5;
  % make into a weight vector that sums to 1
  w = pWhichDist./(sum(pWhichDist,2)*ones(1,size(pWhichDist,2)));

  % no compute means and standard deviations weighted by the
  % weight vector we just computed
  oldm = m;
  m = mean(w.*xk,2)*n;
  s = sqrt((1./(1-sum(w.^2,2))).*sum(w.*(xk-m*ones(1,n)).^2,2));

  % display
  cla;
  H = myhist(x,100,'k',1);
  hold on
  a = axis;
  x1 = a(1):(a(2)-a(1))/100:a(2);
  for i = 1:k
    plot(H.bins,normpdfbins(H.bins,m(i),s(i)),getcolor(i));
  end
  vline(m);
  drawnow
end

% sort the means
[m i] = sort(m);
s = s(i);
disppercent(inf);

%%%%%%%%%%%%%%%%%%%%%
%%   normpdfbins   %%
%%%%%%%%%%%%%%%%%%%%%
function y = normpdfbins(x,m,s)

% get binsize
binsize = median(diff(x));

x = [x(:)'-binsize/2 x(end)+binsize/2];
y = diff(normcdf(x,m,s));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   colorMatchingVoxels   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imOut = colorMatchingVoxels(im,m,s)

% find all matching voxels
dist = abs((im-m)/s);
dist(dist>2) = nan;
dist = (2-dist)/2;
dist(isnan(dist)) = 0;

match = find((im(:) < (m+s)) & (im(:) > (m-s)));
disp(sprintf('(t1t2:colorMtachingVoxels) Found %i voxels that match within 1 std',length(match)));

% make in image where they are red
im1 = im;
im0 = im;
im0 = im.*(1-dist);
imOut(:,:,1) = im1;
imOut(:,:,2) = im0;
imOut(:,:,3) = im0;


%%%%%%%%%%%%%%%%%%%%%%%%%%
%    setMaxFillCenter    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function setMaxFillCenter(src,eventdata)

global gt1t2;

% figure out which axis we are on
pointerLoc = get(gt1t2.f,'CurrentPoint');
pos = get(gt1t2.f,'Position');
pos = pointerLoc./pos(3:4);
subplotNum = ceil(pos(1)*3);
if subplotNum == 2
  a = subplot(1,3,subplotNum);
else
  return
end

% get the point on that axis
pointerLoc = get(a,'CurrentPoint');
pointerX = round(pointerLoc(1,1));
pointerY = round(pointerLoc(1,2));

% get image dims
dims = size(gt1t2.t2.d);

% get image dimension
switch gt1t2.params.sliceOrientation
  case {'axial'}
   x = pointerX; y = dims(2)-pointerY+1; z = gt1t2.params.sliceNum;
  case {'coronal'}
   x = pointerX; y = gt1t2.params.sliceNum; z= dims(3)-pointerY+1;
  case {'saggital'}
   x = gt1t2.params.sliceNum; y = pointerX; z = dims(3)-pointerY+1;
end

% check dimensions
if (x < 1) || (y < 0) || (z< 0) || (x > dims(1)) || (y > dims(2)) || (z > dims(3))
  return
end

% display coordinates
disp(sprintf('(t1t2:t1t2floodfill) Pointer at [%i,%i,%i]',x,y,z));
params.maxFillCenter = [x y z];
mrParamsSet(params);
return;


%%%%%%%%%%%%%%%%%%%%%
%    t1t2Maxfill    %
%%%%%%%%%%%%%%%%%%%%%
function retval = t1t2MaxFill(params)

retval = 0;

global gt1t2;
dims = size(gt1t2.t2.d);

maxfill = nan(dims);

visitedPoints = [];
currentPoints = params.maxFillCenter';
currentPoints = sub2ind(dims,currentPoints(1,:),currentPoints(2,:),currentPoints(3,:));

% get how the neighbors are spaced in a linear array, using getNeighbors to
% generate the list of neighbors. This is for calculating the next set
% of points to fill
nextNeighborDist = currentPoints-getNeighbors(currentPoints,dims,1.5)';

% now get another list, this is for checking neighbors for max filling
maxFillNeighborDist = currentPoints-getNeighbors(currentPoints,dims,params.maxFillRadius)';

disppercent(-inf,'(t1t2) Filling volume');
while ~isempty(currentPoints)
  % get the neighbors using the fast algorithm (which simply adds the distances
  % between the current points and the neighbors in linear coordinates)
  currentNeighbors = getNeighborsFast(currentPoints(:)',maxFillNeighborDist,dims);
  % get the max across neighbors for each current point in the original
  % and in the current state of the filled image
  maxFromOriginal = max(gt1t2.t2.blurd(currentNeighbors));
  maxFromFilled = max(maxfill(currentNeighbors));
  % find which of the current points are in the mask
  inmask = gt1t2.mask(currentPoints)';
  if ~params.maxFillMask
    inmask(:) = 1;
  end
  % if they are in the mask then take the max from the original
  maxfill(currentPoints(inmask)) = maxFromOriginal(inmask);
  % if not in the mask then use the max from the original and from the filled
  % this should propogate bright values outside of the brain
  maxfill(currentPoints(~inmask)) = max([maxFromOriginal(~inmask);maxFromFilled(~inmask)]);
  % update the visitedPoints list
  visitedPoints = union(visitedPoints,currentPoints);
  % get the neighbors of the current point
  currentPoints = getNeighborsFast(currentPoints(:)',nextNeighborDist,dims);
  % remove any points already visited.
  currentPoints = setdiff(currentPoints(:),visitedPoints);
  % update percent done display
  disppercent(length(visitedPoints)/prod(dims));
end
disppercent(inf);
gt1t2.t2.blurd = maxfill;
displayT1T2(gt1t2.params);

% this will be the starting point of the floodfill
%gt1t2.t2.blurd(x,y,z) = 0;
%displayT1T2(gt1t2.params)
%gt1t2.t2.blurd(neighbors) = 0;

%keyboard

% update string of what we have done
gt1t2.t2.processingList = sprintf('%sMaxfill Center: [%s] Radius: %0.2f Mask: %i\n',gt1t2.t2.processingList,num2str(params.maxFillCenter,'%i '),params.maxFillRadius,params.maxFillMask);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getNeighborsFast    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function neighbors = getNeighborsFast(inlist,neighborDist,dims)

% convert the inlist of voxels into a matrix where each
% column contains the voxel number repeated by the 
% number of voxels that are its neighbors
inlist = repmat(inlist,length(neighborDist),1);
% then add the distance to the neighbors, thus we 
% will have a matrix in which each column contains
% all the neighbors of each voxel
neighbors = inlist + repmat(neighborDist,1,size(inlist,2));
% do a boundary check to make sure that we don't go past the
% first or last voxel. Note that this does NOT check for
% whether the neighbor has wrapped around an edge of the image
boundaryValues = (neighbors < 1) | (neighbors > prod(dims));
neighbors(boundaryValues) = inlist(boundaryValues);

%%%%%%%%%%%%%%%%%%%%%%
%    getNeighbors    %
%%%%%%%%%%%%%%%%%%%%%%
function outlist = getNeighbors(inlist,dims,dist)

% convert from linear coordinates
[x y z] = ind2sub(dims,inlist);

% how many voxels we have
nVoxels = size(inlist,2);

% add or subtract 1 in every combination of dimensions to get 
% all neighbors (will also include voxel itself)
outlist = [];
for xOffset = -floor(dist):ceil(dist)
  for yOffset = -floor(dist):ceil(dist)
    for zOffset = -floor(dist):ceil(dist)
      if sqrt(xOffset^2+yOffset^2+zOffset^2) < dist
	outlist = [outlist [x;y;z]+repmat([xOffset yOffset zOffset]',1,nVoxels)];
      end
    end
  end
end

% convert to linear coordinates
outlist = mrSub2ind(dims,outlist(1,:),outlist(2,:),outlist(3,:));

% remove nans (which will be ones outside the volume as computed
% by mrSub2ind
outlist = outlist(~isnan(outlist));
outlist = unique(outlist);
