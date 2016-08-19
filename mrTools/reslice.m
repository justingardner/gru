% reslice.m
%
%        $Id: reslice.m,v 1.19 2008/08/19 12:00:21 justin Exp $
%      usage: reslice(volume)
%         by: justin gardner
%       date: 05/04/07
%    purpose: 
%             program to reslice a volume, call with filename
%             reslice('jg041001');
%
function retval = reslice(event,viewNum)

% check arguments
if ~any(nargin == [1 2 3])
  help reslice
  return
end

% init arguments
if nargin == 1
  % if we are passed in a structure then this is a button callback
  if isstruct(event)
    viewNum = event.viewNum;
    event = event.event;
    retval = [];
    % otherwise it is init event
  elseif isstr(event)
    filename = sprintf('%s.img',stripext(event));
    if ~isfile(filename)
      filename = setext(filename,'nii');
    end
    event = 'init';
    % check for file
    if isfile(filename)
      % read header
      header = mlrImageHeaderLoad(filename);
      % if this is a 4D volume then we have to take a particular volume or
      % or the mean. Ask the user what to do.
      doMean = 0;volNum = [];
      if (header.nDim >= 4) && (header.dim(4) > 1)
	paramsInfo = {{'volNum',0,'incdec=[-1 1]','round=1',sprintf('minmax=[0 %i]',header.dim(4)),'Choose volume number to display (0 for mean)'}};
	params = mrParamsDialog(paramsInfo,'Choose volume to use (0 for mean)');
	if isempty(params),return,end
	if params.volNum ~= 0
	  volNum = params.volNum;
	else
	  doMean = 1;
	end
      end
      % read it
      disppercent(-inf,sprintf('(reslice) Loading %s',filename));
      [vol header] = mlrImageLoad(filename,'volNum',volNum,'orient=LPI');
      if doMean,vol = mean(vol,4);end
      disppercent(inf);
    else
      disp(sprintf('(reslice) Could not open file %s',filename));
      return
    end
  % passed in volume
  elseif isnumeric(event)
    filename = '';
    vol = event;
    event = 'init';
    hdr = [];
  else
    help reslice;
    return
  end
end

switch (event)
 case 'init'
  initHandler(filename,vol,header);
 case 'end'
  endHandler(viewNum);
 case 'mouseMove'
  %  mouseMoveHandler(viewNum);
 case 'mouseUp'
  mouseUpHandler(viewNum);
 case 'save'
  saveHandler(viewNum);
 case 'clearROI'
  clearROIHandler(viewNum);
 case 'loadROI'
  loadROIHandler(viewNum);
 case 'print'
  printHandler(viewNum);
 case 'exportXform'
  exportXformHandler(viewNum);
 case 'import'
  importHandler(viewNum);
 case 'quit'
  endHandler(viewNum);
 case 'mouseDown'
  mouseDownHandler(viewNum);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mousemove
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseMoveHandler(viewNum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mouseup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseUpHandler(viewNum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mouseup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseDownHandler(viewNum)

global gReslice;

% get the pointer location
pointerLoc = get(gca,'CurrentPoint');
% get mouse, remembering that we have swapped y/x
% (see getImageSlice) 
mouseX = round(pointerLoc(1,1));
mouseY = round(pointerLoc(1,2));

% which figure we are on
if gcf == gReslice{viewNum}.fig(1)
  a = [1 2 3];
  %  mouseX = size(gReslice{viewNum}.vol,2)-mouseX+1;
  mouseY = size(gReslice{viewNum}.vol,3)-mouseY+1;
  disp(sprintf('(reslice) z:%i y:%i',mouseY,mouseX));
elseif gcf == gReslice{viewNum}.fig(2)
  a = [2 1 3];
  %  mouseX = size(gReslice{viewNum}.vol,1)-mouseX+1;
  mouseY = size(gReslice{viewNum}.vol,3)-mouseY+1;
  disp(sprintf('(reslice) z:%i x:%i',mouseY,mouseX));
elseif gcf == gReslice{viewNum}.fig(3)
  a = [3 1 2];
  %  mouseX = size(gReslice{viewNum}.vol,1)-mouseX+1;
  mouseY = size(gReslice{viewNum}.vol,2)-mouseY+1;
  disp(sprintf('(reslice) y:%i x:%i',mouseY,mouseX));
else
  return
end

if ((mouseX > 0) && (mouseX < gReslice{viewNum}.dim(a(2))) && ...
    (mouseY > 0) && (mouseY < gReslice{viewNum}.dim(a(3))))
  figure(gReslice{viewNum}.fig(a(2)));
  dispVolumeSlice(viewNum,a(2),mouseX);
  gReslice{viewNum}.pos(a(2)) = mouseX;

  figure(gReslice{viewNum}.fig(a(3)));
  dispVolumeSlice(viewNum,a(3),mouseY);
  gReslice{viewNum}.pos(a(3)) = mouseY;
  set(gReslice{viewNum}.fig(a(1)),'pointer','fullcrosshair');
else
  set(gReslice{viewNum}.fig(a(1)),'pointer','arrow');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get current reslice rotation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rotmatrix = getRotMatrix(viewNum,params,offset)


% get the transformation marix
% rotxy
c = cos(d2r(params.xyRot));
s = sin(d2r(params.xyRot));
rotxy = [c -s 0 0;s  c 0 0;0  0 1 0;0  0 0 1];

% rotyz
c = cos(d2r(params.yzRot));
s = sin(d2r(params.yzRot));
rotyz = [1  0  0 0;0  c -s 0;0  s  c 0;0  0  0 1];

% rot xz
c = cos(d2r(params.xzRot));
s = sin(d2r(params.xzRot));
rotxz = [c  0 -s 0;0  1  0 0;s  0  c 0;0  0  0  1];

% if viewNum is empty, then just make rot matrix
if isempty(viewNum)
  offset = eye(4);
  sliceOffset = eye(4);
  flipMatrix = eye(4);
% otherwise composit with offset matrices
else
  global gReslice;
  dim = gReslice{viewNum}.dim;

  % offset
  if ieNotDefined('offset')
    offset = [1  0  0 params.xCenter+round(dim(1)/2);
	      0  1  0 params.yCenter+round(dim(2)/2);
	      0  0  1 params.zCenter+round(dim(3)/2);
	      0  0  0    1];
  end

  sliceOffset = [1 0 0 -params.width/2;
		 0 1 0 -params.height/2;
		 0 0 1 0;
		 0 0 0 1];
  
  flipMatrix = gReslice{viewNum}.flipMatrix;
end
  
% full rotation matrix (with possible flip matrix - which is used
% only to make importHeader work)
rotmatrix = offset*rotxy*rotyz*rotxz*flipMatrix*sliceOffset;

% testing
r = rotxy*rotyz*rotxz;

a = atan2(-r(2,3),r(2,2));
b = atan2(-r(3,1),r(1,1));
g = asin(r(2,1));
disp(sprintf('[%0.1f %0.1f %0.1f]',r2d(a),r2d(b),r2d(g)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clearROIHandler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clearROIHandler(viewNum)

global gReslice;
gReslice{viewNum}.roiCoords = [];
% refresh the display
refreshResliceDisplay(viewNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exportXform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function exportXformHandler(viewNum)

global gReslice;
params = gReslice{viewNum}.params;

% get vol2mag
if ~isempty(gReslice{viewNum}.header.sform)
  vol2mag = gReslice{viewNum}.header.sform;
else
  disp(sprintf('(reslice) Volume %s does not have sform set. Using qform instead',getLastDir(gReslice{viewNum}.filename)));
  vol2mag = gReslice{viewNum}.header.qform;
end  

% get the rotation matrix that transforms to the
% resliced coordinates
reslice2vol = getReslice2vol(viewNum);

% comput inp2mag -- this can be used as a qform44 or an sform44
inp2mag = vol2mag*reslice2vol;

% put up dialog to select how to export
nparamsInfo = {};
vol2mag
inp2mag
saveType = questdlg('Export transformation to Nifti file as','Export xform','qform','sform','cancel','cancel');
if strcmp(saveType,'cancel'),return,end

% pick file
[filename, pathname, filterindex] = uigetfile('*.hdr', 'Pick a Nifti file');

if isequal(filename,0),return,end

% read nifti file
filename = fullfile(pathname,filename);
if ~isfile(filename),disp(sprintf('(reslice) Could not open file %s',filename));return,end
[d header] = mlrImageLoad(filename);
if isempty(d),return,end

% get voxel sizes
pixdim = header.pixdim(1:3);
disp(sprintf('(reslice) Pixdims for %s are: [%0.2f %0.2f %0.2f]',getLastDir(filename),pixdim(1),pixdim(2),pixdim(3)));

% scale the qform appropriately, according to the pixdims -- so, make pixdims into
% a 4x4 matrix and multiply with the inp2mag (which is for a 1x1x1 voxel size)
pixdim = [[diag(pixdim); 0 0 0] [0 0 0 1]'];
inp2mag = inp2mag*pixdim;

% set the approprate field
if strcmp(saveType,'qform')
  header.qform = inp2mag;
  % display qform
  qform = inp2mag
else
  h.sform = inp2mag;
  % display sform
  sform = inp2mag
end

% write out the file
mlrImageSave(filename,d,header);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loadROIHandler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadROIHandler(viewNum)

global gReslice;

% get the roi to load
% now get the filename to save to
thispath = pwd;
if isfield(gReslice{viewNum},'loadROIPath')
  cd(gReslice{viewNum}.loadROIPath);
end
[loadnames loadpath] = uigetfile({'*.mat','mrLoadRet 4 ROI file (*.mat)'},'Load ROI','MultiSelect','on');
cd(thispath);

% check for user cancel
if isequal(loadnames,0) || isequal(loadpath,0)
  return
end

% make sure that loadname is a cell array
loadnames = cellArray(loadnames);

% otherwise set the default loadROIPath
gReslice{viewNum}.loadROIPath = loadpath;

% set proper extension
for j = 1:length(loadnames)
  loadname = loadnames{j};
  loadname = fullfile(loadpath,sprintf('%s.mat',stripext(loadname)));
  disp(sprintf('Loading ROI %s',loadname));

  ROI = load(loadname);
  roiName = fieldnames(ROI);
  if ~isfield(gReslice{viewNum},'header')
    disp(sprintf('(reslice) Need to run reslice with a file name to have nifti hdr'));
    return
  else
    % loop over each roi in file
    for i = 1:length(roiName)
      if ~isfield(ROI.(roiName{i}),'xform')
	disp(sprintf('(reslice) ROI %s is not a mrLoadRet-4 ROI. It has no xform.',roiName{i}));
      else
	% get xforms and voxel sizes
	volXform = gReslice{viewNum}.header.sform;
	volVoxelSize = gReslice{viewNum}.header.pixdim(1:3)';
	roiXform = ROI.(roiName{i}).xform;
	roiVoxelSize = ROI.(roiName{i}).voxelSize;
	roiCoords = ROI.(roiName{i}).coords;
	% and concatenate on to out roiCoords
	gReslice{viewNum}.roiCoords = [gReslice{viewNum}.roiCoords round(xformROIcoords(roiCoords,inv(volXform)*roiXform,roiVoxelSize,volVoxelSize))];
      end
    end
  end
end

% refresh the display
refreshResliceDisplay(viewNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% printHandler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printHandler(viewNum)

global gReslice;
params = gReslice{viewNum}.params;

% choose which view to print out
sliceDim = find(strcmp(params.printView,{'saggital','coronal','axial'}));

% get current slice num
sliceNum = gReslice{viewNum}.sliceNum(sliceDim);

% get the other dimensions
otherDim = setdiff([1 2 3],sliceDim);

% slices to print
printSlices = sliceNum+params.printSlices;

% slice title
sliceTitle = sprintf('Rot (yz=%0.1f xz=%0.1f xy=%0.1f) Pos (y=%0.1f x=%0.1f z=%0.1f) %0.1fx%0.1f',params.yzRot,params.xzRot,params.xyRot,params.yCenter,params.xCenter,params.zCenter,params.width,params.height);

% cycle through each figure and print
for i = 1:length(printSlices)
  figure(gReslice{viewNum}.fig(sliceDim));
  dispVolumeSlice(viewNum,sliceDim,printSlices(i));
  % set the title
  title(sprintf('%s (slice=%i)\n%s',gReslice{viewNum}.filename,printSlices(i),sliceTitle));
  drawnow
  printdlg(gReslice{viewNum}.fig(sliceDim));
end

% redisplay original 
figure(gReslice{viewNum}.fig(sliceDim));
dispVolumeSlice(viewNum,sliceDim,sliceNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getReslice2vol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function reslice2vol = getReslice2vol(viewNum)

global gReslice;
reslice = [];
dim = size(gReslice{viewNum}.vol);
params = gReslice{viewNum}.params;

% get the rotation matrix (do not offset to center of volume)
rotmatrix = getRotMatrix(viewNum, gReslice{viewNum}.params, eye(4));
% offset to zCenter
zOffset = eye(4);
zOffset(3,4) = -params.slices(1);

offsetSliceCenter = [1  0  0 params.xCenter;
  		     0  1  0 params.yCenter;
		     0  0  1 params.zCenter;
		     0  0  0    1];

offsetVolumeCenter = [1  0  0 round(dim(1)/2);
		      0  1  0 round(dim(2)/2);
		      0  0  1 round(dim(3)/2);
		      0  0  0    1];

% get slice thickness
zThick = eye(4);
if ~isempty(first(diff(params.slices)))
  zThick(3,3) = first(diff(params.slices));
end

% and compute the transform that goes from the volume coordinates
% to the resliced coordiantes
reslice2vol = shiftOriginXform*offsetVolumeCenter*offsetSliceCenter*rotmatrix*inv(zOffset)*zThick;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% importHandler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function importHandler(viewNum)

% pick file
[filename, pathname, filterindex] = uigetfile('*.hdr', 'Pick a Nifti file to import reslice settings from');
if isequal(filename,0),return,end

% read nifti file
filename = fullfile(pathname,filename);
if ~isfile(filename),disp(sprintf('(reslice) Could not open file %s',filename));return,end
[header] = mlrImageHeaderLoad(filename);
if isempty(header),return,end

global gReslice;
dim = gReslice{viewNum}.dim;

% check for sform being set on loaded volume
volheader = gReslice{viewNum}.header;
if isempty(volheader.sform)
  disp(sprintf('(reslice) *** Loaded volume does not have sform set! Using qform ***'));
  vol2mag = volheader.qform;
else
  vol2mag = volheader.sform;
end

disp(sprintf('(reslice) volume pixdims: [%0.1f %0.1f %01.f]',volheader.pixdim(1),volheader.pixdim(2),volheader.pixdim(3)));
disp(sprintf('(reslice) import pixdims: [%0.1f %0.1f %01.f]',header.pixdim(1),header.pixdim(2),header.pixdim(3)));
disp(sprintf('(reslice) import dims: [%i %i %i]',header.dim(1),header.dim(2),header.dim(3)));

disp(sprintf('(reslice) volume qform\n%s',mlrnum2str(volheader.qform,'compact=0')));
disp(sprintf('(reslice) volume sform\n%s',mlrnum2str(volheader.sform,'compact=0')));
disp(sprintf('(reslice) import qform\n%s',mlrnum2str(header.qform,'compact=0')));
disp(sprintf('(reslice) import sform\n%s',mlrnum2str(header.sform,'compcat=0')));

% check for sform being set on import file
if isempty(header.qform)
  disp(sprintf('(reslice) *** Import file %s does not have sform set - has it been aligned to volume with mrAlign? Using qform ***',filename));
  inplane2mag = header.qform;
else
  inplane2mag = header.sform;
end

% get reslice2vol - transform from reslice to the base being displayed
% so we invert the vol2mag of the base being displayed and multiply
% it by the sform44
reslice2vol = inv(vol2mag)*inplane2mag;

% undo pixel scaling
reslice2vol = reslice2vol*inv(diag([header.pixdim(1:2) 1 1]));

% get the width and height
params.width = round(header.pixdim(1)*header.dim(1));
params.height = round(header.pixdim(2)*header.dim(2));

% get slices
params.slices = sprintf('0:%0.1f:%0.1f',header.pixdim(3),(header.dim(3)-1)*header.pixdim(3));

offset = eye(4);

% get center of slice
sliceOffset = [1 0 0 -params.width/2;
	       0 1 0 -params.height/2;
	       0 0 1 0;
	       0 0 0 1];

% get the slice thickness
zThick = eye(4);
zThick(3,3) = header.pixdim(3);

% get the offset to the first slice, this can always be identity
% since we will set the "slices" params to start at 0.
zOffset = eye(4);

% this is the offset to the middle of the loaded volume
offsetVolumeCenter = [1  0  0 round(dim(1)/2);
		      0  1  0 round(dim(2)/2);
		      0  0  1 round(dim(3)/2);
		      0  0  0    1];

% get the rotation+scaling matrix
% These are the two transformations that need to be inverted. 
%
% reslice2vol = shiftOriginXform*offsetVolumeCenter*offsetSliceCenter*rotmatrix*inv(zOffset)*zThick;
% rotmatrix = offset*rotxy*rotyz*rotxz*sliceOffset;
%
% Note that offsetSliceCenter and rotmatrix are the transformations we want.
% From these we can extract the parameters needed for the gui.
r = inv(shiftOriginXform)*inv(offsetVolumeCenter)*inv(offset)*reslice2vol*inv(sliceOffset)*zOffset*inv(zThick);

% grab the offset of the slice center
params.xCenter = r(1,4);
params.yCenter = r(2,4);
params.zCenter = r(3,4);

% extract the Euler angles. This comes from proper inversion of terms in
% our rotation matrix (i.e. the one in getRotMatrix). Note, that there
% are probably two cases, in which the following equations have a singularity
% need to special case those....
%[c1c3+s1s2s3 -s1c2 -c1s3+s1s2c3 0]
%[s1c3-c1s2s3  c1c2 -s1s3-c1s2c3 0]
%[c2s3           s2         c2c3 0]
%[   0            0            0 1]
%
% 1 = xy
% 2 = yz
% 3 = xz
%
% Note, that there are multiple ambiguities with this scheme. For example,
% tan(theta) = tan(180+theta) and sin(theta) = sin(-theta)
% so, we add those possibilities in. Also, this scheme does not allow for
% axis flips. So we also add the possibilities for axis flips in which gives
% us 4 possible angles for each rotation.
xyRot = [r2d(atan2(-r(1,2),r(2,2))) 180+r2d(atan2(-r(1,2),r(2,2))) r2d(atan2(r(1,2),r(2,2))) 180+r2d(atan2(r(1,2),r(2,2)))];
yzRot = [r2d(asin(r(3,2))) 180-r2d(asin(r(3,2))) r2d(-asin(r(3,2))) 180-r2d(-asin(r(3,2)))];
xzRot = [r2d(atan2(r(3,1),r(3,3))) 180+r2d(atan2(r(3,1),r(3,3))) r2d(-atan2(r(3,1),r(3,3))) 180+r2d(atan2(-r(3,1),r(3,3)))];

% now we check each possibility to see which gives us a match to the desired rotation matrix
matched = 0;
% try all combination of angles from above
for xyi = 1:length(xyRot)
  if matched, break, end
  for yzi = 1:length(yzRot)
    if matched, break, end
    for xzi = 1:length(xzRot)
      if matched, break, end
      % get the angles for this try
      params.xyRot = xyRot(xyi);
      params.yzRot = yzRot(yzi);
      params.xzRot = xzRot(xzi);
      % create the rotation matrix
      r2 = getRotMatrix([],params);
      % see if the rotation matrix is the same in absolute value
      % this means that we may have the right rotations, but be
      % off by a few axis flips
      if nearlyequal(abs(r2(1:3,1:3)),abs(r(1:3,1:3)),2)
	for i1 = [1 -1]
	  if matched, break, end
	  for j1 = [1 -1]
	    if matched, break, end
	    for k1 = [1 -1]
	      % check if we now match with these axis flips
	      if nearlyequal(r2(1:3,1:3)*diag([i1 j1 k1]),r(1:3,1:3),2)
		disp(sprintf('(reslice) Matched rotation with: [%i %i %i],[%i %i %i]\n%s',i1,j1,k1,xyi,yzi,xzi,mynum2str(r2(1:3,1:3)*diag([i1 j1 k1]),'sigfigs=3')));
		matched = 1;
		break;
	      end
	    end
	  end
	end
      end
    end
  end
end

% make sure we have matched, otherwise something is broken in the code above
% and we need to fix it since the settings won't be correct.
if ~matched
  disp(sprintf('(reslice) %s',repmat('*',1,40)));
  disp(sprintf('(reslice) *** Could not find angles to exactly match rotation matrix ****'));
  disp(sprintf('(reslice) *** reslice:importHandler slice placement may not be correct ****'));
  params.xyRot = xyRot(1);
  params.yzRot = yzRot(1);
  params.xzRot = xzRot(1);
  i1 = 1;j1 = 1;k1 = 1;
  disp(sprintf('(reslice) Rotation matrix should be\n%s',mynum2str(r(1:3,1:3),'sigfigs=3')));
  disp(sprintf('(reslice) Euler angle reconstructed rotation matrix\n%s',mynum2str(getRotMatrix([],params),'sigfigs=3')));
  disp(sprintf('(reslice) %s',repmat('*',1,40)));
end

% set the flip matrix
gReslice{viewNum}.flipMatrix = diag([i1 j1 k1 1]);

% display settings
disp(sprintf('(reslice) rotation: [%0.1f %0.1f %0.1f]\n(reslice) center: [%0.1f %0.1f %0.1f]\n(reslice) width x height: [%0.1fx%0.1f]\n(reslice) slices: %s',params.yzRot,params.xzRot,params.xyRot,params.yCenter,params.xCenter,params.zCenter,params.width,params.height,params.slices));

% set the params and refresh display
params = mrParamsSet(params);
resliceControls(params,viewNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveHandler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveHandler(viewNum)

global gReslice;
reslice = [];
dim = size(gReslice{viewNum}.vol);
params = gReslice{viewNum}.params;

for sliceNum = 1:length(gReslice{viewNum}.params.slices)
  reslice(:,:,sliceNum) = getReslice(viewNum,sliceNum)';
end

if ~isempty(gReslice{viewNum}.header)
  reslice2vol = getReslice2vol(viewNum)
  % get the original header
  header = gReslice{viewNum}.header;
  % and set the qform/sform appropriately - that is, they should be the reslice2mag xforms, so we
  % need to composite reslice2vol with the qform44/sform44 (vol2mag). 
  header.qform = header.qform*reslice2vol;
  if ~isempty(header.sform)
    header.sform = header.sform*reslice2vol;
  end
else
  disp(sprintf('(reslice) No hdr, if you need a valid qform, run reslice(''filename.img'')'));
  header = [];
end

% now get the filename to save to
thispath = pwd;
if isfield(gReslice{viewNum},'savepath')
  cd(gReslice{viewNum}.savepath);
end
%if ~isfield(gReslice{viewNum},'savepath')
[savename savepath] = uiputfile({'*.hdr','Nifti files (*.hdr)'},'Save as');
%else
%  savepath = gReslice{viewNum}.savepath;savename = 'reslice.hdr';
%end
cd(thispath);

% check for user cancel
if isequal(savename,0) || isequal(savepath,0)
  return
end

% otherwise set the default savepath
gReslice{viewNum}.savepath = savepath;

% set proper extension
savename = fullfile(savepath,sprintf('%s.hdr',stripext(savename)));
disp(sprintf('Saving %s',savename));

%save it
mlrImageSave(savename,reslice,header);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end the mrInterrogator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function endHandler(viewNum)

global gReslice;
% only run this if it is not being called by someonw else
if ~gReslice{viewNum}.shutdown
  gReslice{viewNum}.shutdown = 1;
  if ishandle(gReslice{viewNum}.controlFig)
    close(gReslice{viewNum}.controlFig);
  end
  if ishandle(gReslice{viewNum}.resliceFignum)
    close(gReslice{viewNum}.resliceFignum);
  end
  for i = 1:3
    if ishandle(gReslice{viewNum}.fig(i))
      close(gReslice{viewNum}.fig(i));
    end
  end
  gReslice{viewNum}.shutdown = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init the interrogator handler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initHandler(filename,vol,header)

global gReslice;

% get figure handles
viewNum = mlrSmartfig('reslice1');
gReslice{viewNum} = [];
gReslice{viewNum}.fig(1) = viewNum;
gReslice{viewNum}.fig(2) = mlrSmartfig('reslice2');
gReslice{viewNum}.fig(3) = mlrSmartfig('reslice3');
gReslice{viewNum}.resliceFignum = mlrSmartfig('reslice4');
gReslice{viewNum}.roiCoords = [];
gReslice{viewNum}.shutdown = 0;
gReslice{viewNum}.filename = filename;
gReslice{viewNum}.flipMatrix = eye(4);
% set the callbacks appropriately
for i= 1:3
  %  set(gReslice{viewNum}.fig(i),'MenuBar','none');
  set(gReslice{viewNum}.fig(i),'WindowButtonMotionFcn',sprintf('reslice(''mouseMove'',%i)',viewNum));
  set(gReslice{viewNum}.fig(i),'WindowButtonDownFcn',sprintf('reslice(''mouseDown'',%i)',viewNum));
  set(gReslice{viewNum}.fig(i),'WindowButtonUpFcn',sprintf('reslice(''mouseUp'',%i)',viewNum));
  set(gReslice{viewNum}.fig(i),'DeleteFcn',sprintf('reslice(''end'',%i)',viewNum));
end
set(gReslice{viewNum}.resliceFignum,'DeleteFcn',sprintf('reslice(''end'',%i)',viewNum));

% set pointer to crosshairs
set(gReslice{viewNum}.fig(1),'pointer','fullcrosshair');
set(gReslice{viewNum}.fig(2),'pointer','fullcrosshair');
set(gReslice{viewNum}.fig(3),'pointer','fullcrosshair');

% set info for callback
gReslice{viewNum}.viewNum = viewNum;

% set volume
gReslice{viewNum}.vol = vol;
gReslice{viewNum}.header = header;
gReslice{viewNum}.dim = size(vol);
gReslice{viewNum}.gamma = 0.4;

gReslice{viewNum}.pos(1) = round(size(vol,1)/2);
gReslice{viewNum}.pos(2) = round(size(vol,2)/2);
gReslice{viewNum}.pos(3) = round(size(vol,3)/2);

% compute color map
g = gray(256);
y = g;y(:,3) = 0;
c = g;c(:,1) = 0;
m = g;m(:,2) = 0;
r = g;r(:,2:3) = 0;
myColormap = [g;y;c;r;m;m];

% set up controls
paramsInfo = {};
paramsInfo{end+1} = {'gamma',gReslice{viewNum}.gamma,'incdec=[-0.1 0.1]','Gamma of display'};
paramsInfo{end+1} = {'yzRot',0,'incdec=[-1 1]','Rotation around yz axis in degrees'};
paramsInfo{end+1} = {'xzRot',0,'incdec=[-1 1]','Rotation around xz axis in degrees'};
paramsInfo{end+1} = {'xyRot',0,'incdec=[-1 1]','Rotation around xy axis in degrees'};
paramsInfo{end+1} = {'yCenter',0,'incdec=[-1 1]','y Center position in voxels (for standard 1x1x1 volumes this will be equivalent to mm)'};
paramsInfo{end+1} = {'xCenter',0,'incdec=[-1 1]','x Center position in voxels (for standard 1x1x1 volumes this will be equivalent to mm)'};
paramsInfo{end+1} = {'zCenter',0,'incdec=[-1 1]','z Center position in voxels (for standard 1x1x1 volumes this will be equivalent to mm)'};
paramsInfo{end+1} = {'width',192,'incdec=[-1 1]','width of slice in voxels (assuming volume is 1x1x1)'};
paramsInfo{end+1} = {'height',192,'incdec=[-1 1]','height of slice in voxels (assuming volume is 1x1x1)'};
paramsInfo{end+1} = {'slices',0,'slices to create, put in a matlab style array like -3:1.5:3. This will make slices starting -3 voxels (mm for a standard 1x1x1 volume) away from the center specified with the above controls, with each slice placed in 1.5 voxel increments up to +3 voxels away from the center. Note that each slice will be displayed with a 1 voxel thickness, and that round off error may make it show less than the number of slices you have set, e.g. slice at 0.5 and 1 will show up in the same place because it rounds both of these to 1','type=string'};
paramsInfo{end+1} = {'viewSlice',1,'Slice to view','round=1','minmax=[1 inf]','incdec=[-1 1]'};
% set callback argument for save button
callbackArg.viewNum = viewNum;
callbackArg.event = 'save';
paramsInfo{end+1} = {'save',0,'type=pushbutton','callback=reslice','buttonString=Save resliced volume','callbackArg',callbackArg,'Save out the resliced volume to a nifti file. This will give a correct qform and sform for the resliced volume.'};
callbackArg.viewNum = viewNum;
callbackArg.event = 'loadROI';
paramsInfo{end+1} = {'loadROI',0,'type=pushbutton','callback=reslice','buttonString=Load ROI','callbackArg',callbackArg,'Load a mrLoadRet 4 ROI, all ROIs will be displayed in red, you can load multiple ROIs but they all show up in the same color'};
callbackArg.viewNum = viewNum;
callbackArg.event = 'clearROI';
paramsInfo{end+1} = {'clearROI',0,'type=pushbutton','callback=reslice','buttonString=Clear ROIs','callbackArg',callbackArg,'Clear all the ROIs that you have loaded'};
paramsInfo{end+1} = {'dispROI',1,'type=checkbox','Display the ROI in the resliced image or not'};

callbackArg.viewNum = viewNum;
callbackArg.event = 'print';
paramsInfo{end+1} = {'printView',{'saggital','corona','axial'},'Which view to print out when you press the print button'};
paramsInfo{end+1} = {'printSlices',0,'Set a matlab array to specify what slices to print out (e.g.) -5:5 will print out slices from -5 to 5 mm'};
paramsInfo{end+1} = {'print',0,'type=pushbutton','callback=reslice','buttonString=Print','callbackArg',callbackArg,'Print'};
callbackArg.viewNum = viewNum;
callbackArg.event = 'exportXform';
paramsInfo{end+1} = {'exportXform',0,'type=pushbutton','callback=reslice','buttonString=Export xform','callbackArg',callbackArg,'Click this button to export a transformation matrix that can be used as the qform or sform in a nifti header.'};
callbackArg.viewNum = viewNum;
callbackArg.event = 'import';
paramsInfo{end+1} = {'import',0,'type=pushbutton','callback=reslice','buttonString=Import nifti','callbackArg',callbackArg,'Import slice settings from a nifti file. For example, show the slice locations of an inplane that has been aligned to the volume with mrAlign'};
callbackArg.viewNum = viewNum;
callbackArg.event = 'quit';
paramsInfo{end+1} = {'quit',0,'type=pushbutton','callback=reslice','buttonString=Quit','callbackArg',callbackArg,'quit'};

[gReslice{viewNum}.controlFig gReslice{viewNum}.params] = mrParamsDialog(paramsInfo,'Reslice controls',[],@resliceControls,viewNum);

% compute slice coordinates
gReslice{viewNum}.params.slices = eval(gReslice{viewNum}.params.slices);
gReslice{viewNum}.resliceCoords = calcResliceCoords(viewNum,gReslice{viewNum}.params);

% compute each dimension coordinates
[x y] = meshgrid(1:gReslice{viewNum}.dim(2),1:gReslice{viewNum}.dim(3));
gReslice{viewNum}.sliceCoords{1}(1,1:length(x(:))) = gReslice{viewNum}.pos(1);
gReslice{viewNum}.sliceCoords{1}(2,:) = x(:);
gReslice{viewNum}.sliceCoords{1}(3,:) = y(:);
gReslice{viewNum}.sliceMin{1}(1:3) = [gReslice{viewNum}.pos(1) min(x(:)) min(y(:))];
gReslice{viewNum}.sliceMax{1}(1:3) = [gReslice{viewNum}.pos(1) max(x(:)) max(y(:))];
[x y] = meshgrid(1:gReslice{viewNum}.dim(1),1:gReslice{viewNum}.dim(3));
gReslice{viewNum}.sliceCoords{2}(1,:) = x(:);
gReslice{viewNum}.sliceCoords{2}(2,:) = gReslice{viewNum}.pos(2);
gReslice{viewNum}.sliceCoords{2}(3,:) = y(:);
gReslice{viewNum}.sliceMin{2}(1:3) = [min(x(:)) gReslice{viewNum}.pos(1) min(y(:))];
gReslice{viewNum}.sliceMax{2}(1:3) = [max(x(:)) gReslice{viewNum}.pos(1) max(y(:))];
[x y] = meshgrid(1:gReslice{viewNum}.dim(1),1:gReslice{viewNum}.dim(2));
gReslice{viewNum}.sliceCoords{3}(1,:) = x(:);
gReslice{viewNum}.sliceCoords{3}(2,:) = y(:);
gReslice{viewNum}.sliceCoords{3}(3,:) = gReslice{viewNum}.pos(3);
gReslice{viewNum}.sliceMin{3}(1:3) = [min(x(:)) min(y(:)) gReslice{viewNum}.pos(1)];
gReslice{viewNum}.sliceMax{3}(1:3) = [max(x(:)) max(y(:)) gReslice{viewNum}.pos(1)];

% set up cache
gReslice{viewNum}.c = mrCache('init',50);

% display three windows
figure(gReslice{viewNum}.fig(1));
dispVolumeSlice(viewNum,1,gReslice{viewNum}.pos(1));
colormap(myColormap);axis equal;axis tight;axis on

figure(gReslice{viewNum}.fig(2));
dispVolumeSlice(viewNum, 2,gReslice{viewNum}.pos(2));
colormap(myColormap);axis equal;axis tight;axis on

figure(gReslice{viewNum}.fig(3));
dispVolumeSlice(viewNum,3,gReslice{viewNum}.pos(3));
colormap(myColormap);axis equal;axis tight;axis on

figure(gReslice{viewNum}.resliceFignum);
dispReslice(viewNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% callback for handling controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resliceControls(params,viewNum)

global gReslice

params.slices = eval(params.slices);
params.viewSlice = min(params.viewSlice,length(params.slices));
% compute reslice coordinates after transformation
gReslice{viewNum}.resliceCoords = calcResliceCoords(viewNum,params);

% set the global
gReslice{viewNum}.gamma = params.gamma;
gReslice{viewNum}.params = params;

refreshResliceDisplay(viewNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% redisplay everything
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refreshResliceDisplay(viewNum)

global gReslice;

% and redisplay
for i = 1:3
  figure(gReslice{viewNum}.fig(i));
  dispVolumeSlice(viewNum,i,gReslice{viewNum}.pos(i));
end

figure(gReslice{viewNum}.resliceFignum);
dispReslice(viewNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate reslice coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resliceCoords = calcResliceCoords(viewNum,params)

global gReslice;

% get the rotation matrix that transforms to the
% desired resliced coordinates
rotmatrix = getRotMatrix(viewNum,params);

% get the initial coords
% [x y] = meshgrid(-params.width/2:params.width/2,-params.height/2:params.height/2);
[x y] = meshgrid(0:(round(params.width)-1),0:(round(params.height)-1));
for i = 1:length(params.slices)
  initCoords(1,:) = x(:);
  initCoords(2,:) = y(:);
  initCoords(3,:) = params.slices(i);
  initCoords(4,:) = 1;

  % compute coordinates after transformation
  resliceCoords{i} = round(rotmatrix*initCoords);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% display slice of volume
%%%%%%%%%%%%%%%%%%%%%%%%%
function slice = dispVolumeSlice(viewNum,sliceDim,sliceNum)

global gReslice

% save slice num
gReslice{viewNum}.sliceNum(sliceDim) = sliceNum;

% get the other dimensions
otherDim = setdiff([1 2 3],sliceDim);

% clear axis so we don't keep drawing over old
cla

% get slice out of volume
switch (sliceDim)
 case {1}
  slice = squeeze(gReslice{viewNum}.vol(sliceNum,:,:));
 case {2}
  slice = squeeze(gReslice{viewNum}.vol(:,sliceNum,:));
 case {3}
  slice = squeeze(gReslice{viewNum}.vol(:,:,sliceNum));
end

% update slice coords
gReslice{viewNum}.sliceCoords{sliceDim}(sliceDim,:) = sliceNum;
gReslice{viewNum}.sliceMin{sliceDim}(sliceDim) = sliceNum;
gReslice{viewNum}.sliceMax{sliceDim}(sliceDim) = sliceNum;

% gamma correct
smax = max(slice(:));smin = min(slice(:));
slice = (slice-smin)./(smax-smin);
slice = 256*(slice.^gReslice{viewNum}.gamma);

% now check to see if there is any overlap between this slice
% and the reslice
intersectCoords = [];displayIntersectCoords = [];
sMin = gReslice{viewNum}.sliceMin{sliceDim};
sMax = gReslice{viewNum}.sliceMax{sliceDim};
for i = 1:length(gReslice{viewNum}.params.slices)
  % this code used to use intersect, but that is hoeplessly
  % slow. Because we are always displaying cardinal views
  % we can just do min/max checking
  resliceCoords = gReslice{viewNum}.resliceCoords{i}(1:3,:)';
  thisSliceCoords = ...
      ((resliceCoords(:,1) >= sMin(1)) & ...
       (resliceCoords(:,1) <= sMax(1)) & ...
       (resliceCoords(:,2) >= sMin(2)) & ...
       (resliceCoords(:,2) <= sMax(2)) & ...
       (resliceCoords(:,3) >= sMin(3)) & ...
       (resliceCoords(:,3) <= sMax(3)));
  thisIntersectCoords = resliceCoords(thisSliceCoords,:);
  %  thisIntersectCoords = intersect(resliceCoords,gReslice{viewNum}.sliceCoords{sliceDim}(1:3,:)','rows');
  if i ~= gReslice{viewNum}.params.viewSlice
    intersectCoords = [intersectCoords;thisIntersectCoords];
  else
    displayIntersectCoords = thisIntersectCoords;
  end
end

% get the indexes in the current slice for those intersection coordinates
if ~isempty(intersectCoords)
  intersectIndexes = sub2ind(gReslice{viewNum}.dim(otherDim),intersectCoords(:,otherDim(1)),intersectCoords(:,otherDim(2)));
else
  intersectIndexes = [];
end

if ~isempty(displayIntersectCoords)
  displayIntersectIndexes = sub2ind(gReslice{viewNum}.dim(otherDim),displayIntersectCoords(:,otherDim(1)),displayIntersectCoords(:,otherDim(2)));
else
  displayIntersectIndexes = [];
end

% and set all the intersection coordinates to yellow
slice(intersectIndexes) = slice(intersectIndexes)/2+384;
% and those for the current slice to be cyan
slice(displayIntersectIndexes) = slice(displayIntersectIndexes)/2+640;

if gReslice{viewNum}.params.dispROI
  % now get all the roi coords that intersect this slice
  roiCoords = gReslice{viewNum}.roiCoords';
  if ~isempty(roiCoords)
    thisROICoords = ...
	((roiCoords(:,1) >= sMin(1)) & ...
	 (roiCoords(:,1) <= sMax(1)) & ...
	 (roiCoords(:,2) >= sMin(2)) & ...
	 (roiCoords(:,2) <= sMax(2)) & ...
	 (roiCoords(:,3) >= sMin(3)) & ...
	 (roiCoords(:,3) <= sMax(3)));
    thisROICoords = roiCoords(thisROICoords,:);
    if ~isempty(thisROICoords)
      % get the indexes
      roiIndexes = sub2ind(gReslice{viewNum}.dim(otherDim),thisROICoords(:,otherDim(1)),thisROICoords(:,otherDim(2)));
      % and set them to red/magenta
      slice(roiIndexes) = slice(roiIndexes)+768;
    end
  end
end

% flipud and transpose
slice = flipud(slice');

% display
image(slice);
axis equal
axis tight
axis on

titles = {sprintf('%s\nSaggital',gReslice{viewNum}.filename),sprintf('%s\nCoronal',gReslice{viewNum}.filename),sprintf('%s\nAxial',gReslice{viewNum}.filename)};
title(titles{sliceDim},'Interpreter','none');
dimLabels = {'left <- x -> right','posterior <- y -> anterior','inferior <- z -> superior'};
xlabel(dimLabels{otherDim(1)});
ylabel(dimLabels{otherDim(2)});
axis off
% flip the y-axis up-down
set(gca,'YTickLabel',flipud(get(gca,'YTickLabel')));
YLim = get(gca,'YLim');
set(gca,'YTick',YLim(2)-fliplr(get(gca,'YTick')));
axis on
%%%%%%%%%%%%%%%%%%%%%%%%%
% display slice of volume
%%%%%%%%%%%%%%%%%%%%%%%%%
function dispReslice(viewNum);

global gReslice
params = gReslice{viewNum}.params;
% get the slice
reslice = getReslice(viewNum,gReslice{viewNum}.params.viewSlice);
reslice = 256*(((reslice-min(reslice(:)))/(max(reslice(:))-min(reslice(:)))).^gReslice{viewNum}.gamma);

if gReslice{viewNum}.params.dispROI
  % get the coords for this slice
  resliceCoords = gReslice{viewNum}.resliceCoords{gReslice{viewNum}.params.viewSlice};

  % get ROI coords
  roiCoords = gReslice{viewNum}.roiCoords;

  if ~isempty(roiCoords)
    % nan out entries that are not in the volume
    for i = 1:3
      resliceCoords(i,resliceCoords(i,:) < 1) = nan;
      resliceCoords(i,resliceCoords(i,:) > gReslice{viewNum}.dim(i)) = nan;
      roiCoords(i,roiCoords(i,:) < 1) = nan;
      roiCoords(i,roiCoords(i,:) > gReslice{viewNum}.dim(i)) = nan;
    end
    
    % and convert to linear coordinates
    resliceCoords = sub2ind(gReslice{viewNum}.dim,resliceCoords(1,:),resliceCoords(2,:),resliceCoords(3,:));
    roiCoords = sub2ind(gReslice{viewNum}.dim,roiCoords(1,:),roiCoords(2,:),roiCoords(3,:));

    % find the roi coords that exist on this slice
    [matchCoords matchROICoords matchResliceCoords] = intersect(roiCoords,resliceCoords);

    % and set all those coords to red
    reslice(matchResliceCoords) = reslice(matchResliceCoords)+256;
  end
end

% clear axis so we don't keep drawing over old
cla

% draw it
image(flipud(reslice));
g = gray(256);
r = g;r(:,2:3) = 0;
colormap([g;r]);
axis equal;axis tight;axis on;
paramsTitle = sprintf('Rot (yz=%0.1f xz=%0.1f xy=%0.1f) Pos (y=%0.1f x=%0.1f z=%0.1f) %0.1fx%0.1f',params.yzRot,params.xzRot,params.xyRot,params.yCenter,params.xCenter,params.zCenter,params.width,params.height);
title(sprintf('%s\nResliced image slice %i\n%s',gReslice{viewNum}.filename,gReslice{viewNum}.params.viewSlice,paramsTitle),'Interpreter','none');

%%%%%%%%%%%%%%%%%%%%%%%%%
% get slice of volume
%%%%%%%%%%%%%%%%%%%%%%%%%
function reslice = getReslice(viewNum,sliceNum);

global gReslice
dim = gReslice{viewNum}.dim;
params = gReslice{viewNum}.params;

% get reslice
[reslice gReslice{viewNum}.c] = mrCache('find',gReslice{viewNum}.c,getCacheID(viewNum,sliceNum));
if isempty(reslice)
  disppercent(-inf,'(reslice) Computing reslice coordinates');
  xi = gReslice{viewNum}.resliceCoords{sliceNum}(1,:);
  yi = gReslice{viewNum}.resliceCoords{sliceNum}(2,:);
  zi = gReslice{viewNum}.resliceCoords{sliceNum}(3,:);
  reslice = interp3(gReslice{viewNum}.vol,yi,xi,zi,'cubic');
  reslice = reshape(reslice,params.height,params.width);
  % add to cache
  gReslice{viewNum}.c = mrCache('add',gReslice{viewNum}.c,getCacheID(viewNum,sliceNum),reslice);
  disppercent(inf);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
% get reslice cache val
%%%%%%%%%%%%%%%%%%%%%%%%%%
function cacheID = getCacheID(viewNum,sliceNum)

global gReslice;
p = gReslice{viewNum}.params;

cacheID = sprintf('%0.1f_%0.1f_%0.1f_%0.1f_%0.1f_%0.1f_%i_%i_%0.1f',p.yzRot,p.xzRot,p.xyRot,p.xCenter,p.yCenter,p.zCenter,p.width,p.height,sliceNum);


% nearlyequal.m
%
%      usage: nearlyequal(a,b,sigfigs)
%         by: justin gardner
%       date: 12/12/03
%       e.g.: nearlyequal(0.123,0.124,2);
%    purpose: compare numbers up to sigfigs to see if they
%             are the same
%
function retval = nearlyequal(a,b,sigfigs)

if (nargin == 2)
  sigfigs = 0;
elseif (nargin ~= 3)
  help nearlyequal
  return
end

retval = round(a * 10^sigfigs) == round(b * 10^sigfigs);
retval = all(retval(:));

% mynum2str.m
%
%        $Id:$ 
%      usage: mynum2str(num,<sigfigs=2>,<doFixBadChars=false>,<tabs=false>,'compact=true')
%         by: justin gardner
%       date: 09/07/09
%    purpose: num2str that allows setting # of significant figures and also doesn't make so many spaces
%             but still will align numbers from line to line.
%
%             For example, if you set sigfigs=-1 then the program will figure out how many sigfigs are needed
%             to display all numbers and make them align across lines. For example, the following
%             line will produce:
%       e.g.: mynum2str([12.1 -0.001 10.1;-1.3 -12.4 30.01],'sigfigs=-1','compact',false)
%              12.100  -0.001  10.100
%              -1.300 -12.400  30.010
%
%             set sigfigs to the number of significant digits you want to show (numbers are rounded to
%               show the appropriate number of sigfigs).
%             set tabs to true if you want to have each number followed by a tab
%             set compact to false if you want numbers to align from row to row (like in the above example)
%                 the default gives the most compact string possible
%
function s = mynum2str(num,varargin)

% check arguments
if nargin == 0
  help mynum2str
  return
end

% check for empty num, in which case return empty
s = '';
if isempty(num),return,end

% evaluate arguments
sigfigs = [];
doFixBadChars = [];
tabs = [];
compact = [];
getArgs(varargin,{'sigfigs=2','doFixBadChars',false,'tabs',false,'compact',true});

% automatic sigfig
sigFigsEachNum = [];

if sigfigs == -1
  % get how many sigfigs each number needs
  sigFigsEachNum = getSigFigs(num);
  % and get the maximum needed sigfigs
  sigfigs = max(sigFigsEachNum(:));
end

% need sigFigsEachNum if we are doing a compact display
if compact && isempty(sigFigsEachNum)
  sigFigsEachNum = getSigFigs(num,sigfigs);
end

% intialize return
s = '';

% check for non 1d or 2d array
if length(size(num))>2
  disp(sprintf('(mynum2str) Can not display %i dimensional matrix',length(size(s))));
  return
elseif length(size(num)) == 2
  num = num';
  % compute maxnumdigits on the left side of decimal excluding inf and nan
  maxnumdigits = length(sprintf('%i',floor(nanmax(abs(num(~isinf(num(:))))))));
  % if there is an inf or nan, account for that
  if any(isnan(num(:))) | any(isinf(num(:)))
    % if sigfigs is zero then Inf and Nan cannot line up with decimal, so shift over by a space
    if sigfigs == 0
      maxnumdigits = max(maxnumdigits,3);
      if any(num(isinf(num(:)))<0)
	maxnumdigits = max(maxnumdigits,4);
      end
    else
      maxnumdigits = max(maxnumdigits,2);
      if any(num(isinf(num(:)))<0)
	maxnumdigits = max(maxnumdigits,4);
      end
    end
  end
end

% make the string
for j = 1:size(num,2)
  for i = 1:size(num,1)
    % create the formatting string
    if compact
      % create string
      formatString = sprintf('%%s%%0.0%if ',sigFigsEachNum(j,i));
    else
      % figure out how many spaces to add
      if num(i,j) < 0,addspace = '';else addspace = ' ';end
      % get the number of digits to the left of decimal point
      if isnan(num(i,j)) | isinf(num(i,j))
	% for nan or inf, then set spaces so that the Inf or NaN string lines up with decimal
	if sigfigs > 0
	  numDigitsToLeftOfDecimal = 2;
	else
	  numDigitsToLeftOfDecimal = 3;
	end
      else
	% for all other numbers count how many places it takes to represent 
	numDigitsToLeftOfDecimal = length(sprintf('%i',floor(abs(num(i,j)))));
      end
      addspace = [addspace repmat(' ',1,maxnumdigits-numDigitsToLeftOfDecimal)];
      % create string
      formatString = sprintf('%%s%s%%0.0%if ',addspace,sigfigs);
      % add spaces after a nan or inf
      if isnan(num(i,j)) | isinf(num(i,j))
	formatString = sprintf('%s%s',formatString,repmat(' ',1,sigfigs));
      end
    end
    % add tabs if called for
    if tabs
      formatString = sprintf('%s\t',formatString(1:end-1));
    end
    % and update the full string
     s = sprintf(formatString,s,round(num(i,j)*(10^sigfigs))/(10^sigfigs));
  end
  % strip off last space
  s = s(1:end-1);
  % add new line character
  if j ~= size(num,2),s = sprintf('%s\n',s);end
end

% fix bad chars
if doFixBadChars
  s = fixBadChars(s,[],{'.','p'});
  s = s(2:end);
end

%%%%%%%%%%%%%%%%%%%%
%    getSigFigs    %
%%%%%%%%%%%%%%%%%%%%
function retval = getSigFigs(num,maxSigFigs)

% if we want at most 0 sigfigs then just return zero for all elements
if (nargin>=2) && (maxSigFigs == 0)
  retval = zeros(size(num));
  return
end

if nargin == 1
  % maximum number of sigfigs
  maxSigFigs = 6;
  minDiff = 1e-10;
else
  % get the minimum difference that is still considered the same number
  % this is a bit of a hack since at some point due to numerical round
  % off you might be smaller than this limit, but still be different numbers
  % but, this should only mean that you will get less digits than asked
  % for only in the compact case
  minDiff = 10^(-2*maxSigFigs);
end

for i = 1:size(num,1)
  for j = 1:size(num,2)
    % for nan and inf the sigfigs needed are 0
    if isnan(num(i,j)) | isinf(num(i,j))
      retval(i,j) = 0;
    else
      % otherwise find out how many sigfigs are needed
      for iSigfig =  1:maxSigFigs
	% check to see if this number is evenly roundable by this
	% many sig digits
	if (abs(round(num(i,j)*(10^iSigfig))/(10^iSigfig) - num(i,j))) < minDiff
	  break
	end
      end
      retval(i,j) = iSigfig;
    end
  end
end

