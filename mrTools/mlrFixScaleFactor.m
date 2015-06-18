% mlrFixScaleFactor.m
%
%        $Id: mlrFixScaleFactor.m 965 2008-02-21 00:03:20Z justin $
%      usage: mlrFixScaleFactor(scaleFactor)
%         by: justin gardner
%       date: 02/04/2010
%    purpose: Changes the scale factor of all the epis, base
%             anatomies and rois in a session
%             
%
function retval = mlrFixScaleFactor(scaleFactor)

defaultScaleFactor = [1 1 1];

% keep backup headers with the following tagged on 
backupHeaderPostfix = 'preScaleFactorFix';

% check arguments
if nargin == 0
  scaleFactor = defaultScaleFactor;
elseif nargin ~= 1
  help mlrFixScaleFactor
  return
end

if length(scaleFactor) ~= 3
  disp(sprintf('(mlrFixScaleFactor) Scale factor must by [x y z]'));
  return
end

% make into a row vector
scaleFactor = scaleFactor(:)';

% get scale factor matrix
scaleFactorMatrix = diag([scaleFactor 1]);

% first check this directory for mrSession.m
path = '';
if isfile('mrSession.mat')
  % if there is a session then ask if the user wants to export to this directory
  answer = questdlg(sprintf('Fix scale factor in %s?',getLastDir(pwd)),'Export');
  if strcmp(answer,'Cancel'),return,end
  if strcmp(answer,'Yes')
    path = pwd;
  end
end

% if path is not set then ask the user to choose a path
if isempty(path)
  [filename path] = uigetfile('*.mat','Choose a mrSession.mat file to fix scale factor in');
  if filename == 0,return,end
  if ~strcmp(filename,'mrSession.mat')
    mrWarnDlg(sprintf('(mlrFixScaleFactor) %s is not an mrSession file',filename));
    return
  end
end

% start a view in the corresponding location
cd(path);
v = newView;
mrGlobals
% just make sure that the home dir matches
if ~strcmp(MLR.homeDir,path)
  answer = questdlg(sprintf('mrLoadRet is open on session %s? Ok to close and open on %s',getLastDir(MLR.homeDir),getLastDir(path)));
  if ~strcmp(answer,'Yes')
    mrWarnDlg(sprintf('(mlrFixScaleFactor) Could not open a view to %s',getLastDir(path)));
    deleteView(v);
    return
  end
  % clear MLR and start over
  deleteView(v);
  clear global MLR;
  v = newView;
end

% ask the user which groups to export to
groupNames = viewGet(v,'groupNames');
for i = 1:length(groupNames)
  paramsInfo{i}{1} = groupNames{i};
  paramsInfo{i}{2} = 1;
  paramsInfo{i}{3} = 'type=checkbox';
  paramsInfo{i}{4} = sprintf('Fix scale factor in group %s',groupNames{i});
end
params = mrParamsDialog(paramsInfo,'Choose groups to fix scale factors in');
if isempty(params),return,end

% now go through and update the headers
for iGroup = 1:viewGet(v, 'numberofGroups')
  % see if we are supposed to update the group
  if params.(viewGet(v,'groupName',iGroup))
    for iScan = 1:viewGet(v, 'nScans', iGroup);
      % load the nifti header from the mrSession file
      curhdr = viewGet(v, 'niftiHdr', iScan, iGroup);
      % get the nifti filename
      filename = setext(viewGet(v, 'tseriesPath', iScan, iGroup),'hdr');
      % get the nifti backup filenmae
      backupFilename = setext(sprintf('%s_%s',stripext(filename),backupHeaderPostfix),'hdr');
      % check if it is there
      if isfile(filename)
	if isfile(backupFilename)
	  if ~askuser(sprintf('(mlrFixScaleFactor) %s has already been scaled. Rescale from original',backupFilename))
	    disp(sprintf('(mlrFixScaleFactor) Skipping %s',filename));
	    continue
	  end
	  hdr = cbiReadNiftiHeader(backupFilename);
	else
	  % load the header
	  hdr = cbiReadNiftiHeader(filename);
	end
	% check the current sformcode
	if hdr.sform_code ~= 1
	  disp(sprintf('(mlrFixScaleFactor) %s does not have sform set to 1',filename));
	  keyboard
	end
	% set the sform
	hdr = cbiSetNiftiSform(hdr,scaleFactorMatrix*hdr.sform44);
	% move the header to a backup filename
	system(sprintf('mv %s %s',filename,backupFilename));
	% and write it back
	hdr = cbiWriteNiftiHeader(hdr,filename);
	% read it back, (I think there is a slight numerical
	% difference in the sform44 from when it is
	% written to when it is read). This doesn't
	% affect anything, but gives better consistency
	% checking for mrUpdateNiftiHeader
	hdr = cbiReadNiftiHeader(filename);
	% now save it in the session
	v = viewSet(v,'niftiHdr',hdr,iScan,iGroup);
        v = viewSet(v,'scanVoxelSize',scaleFactor.*hdr.pixdim(2:4),iScan,iGroup);
	
	vol2mag = viewGet(v,'scanVol2mag',iScan,iGroup);
	v = viewSet(v,'scanvol2mag',scaleFactorMatrix*vol2mag,iScan,iGroup);
	
      else
	disp(sprintf('(mlrFixScaleFactor) Could not open file %s',filename));
      end
    end
  end
end

% save the session
saveSession
% also remove any base anatomies from mrLastView if it is
% there since those might have a different sform
if isfile('mrLastView.mat')
  load mrLastView
  if ~isempty(view.baseVolumes)
    if ~askuser('(mlrFixScaleFactor) Ok to dump baseVolumes loaded into mrLoadRet viewer? If they are not saved, you may want to save them. Otherwise, they will not get scaled')
      disp(sprintf('(mlrFixScaleFactor) Keeping base anatomies in mrLastView. Note that these loaded bases will not be fixed'));
    else
      view.baseVolumes = [];
    end
  end
  if ~isempty(view.ROIs)
    if ~askuser('(mlrFixScaleFactor) Ok to dump ROIs loaded into mrLoadRet viewer? If they are not saved, you may want to save them. Otherwise, they will not get scaled')
      disp(sprintf('(mlrFixScaleFactor) Keeping ROIs in mrLastView. Note that these loaded ROIs will not be fixed'));
    else
      view.ROIs = [];
    end
  end
  save mrLastView view viewSettings
end

v = loadROI(v);
for i = 1:viewGet(v,'nROIs')
  % get roi name
  roiname = viewGet(v,'roiname',i);
  % get roi xform
  xform = viewGet(v,'roixform',i);
  % get roi voxel size
  voxelSize = viewGet(v,'roiVoxelSize',i);
  % get notes
  notes = viewGet(v,'roinotes',i);
  % check for scaleFactorFix note
  scaleKeyStr = 'scaleFactorFix: ';
  scaleKeyLoc = findstr(scaleKeyStr,notes);
  if ~isempty(scaleKeyLoc)
    [oldScaleFactor count errmsg nextindex] = sscanf(notes(scaleKeyLoc+length(scaleKeyStr):end),'%f %f %f');
    if count ~= 3
      disp(sprintf('(mlrFixScaleFactor) Error reading old scale factor from roi %s. Skipping.',roiname));
      continue
    end
    if ~askuser(sprintf('(mlrFixScaleFactor) Roi %s has already been scaled, rescale based on original xform',roiname))
      disp(sprintf('(mlrFixScaleFactor) Skipping ROI %s',roiname));
      continue
    end
    oldScaleFactor = 1./oldScaleFactor';
    % reset xform
    xform = diag([oldScaleFactor 1])*xform;
    % rest voxelSize
    voxelSize = oldScaleFactor.*voxelSize;
    % reset notes
    if scaleKeyLoc > 1
      notes = sprintf('%s%s',notes(1:scaleKeyLoc),notes(scaleKeyLoc+length(scaleKeyStr)+nextindex:end));
    else
      notes = sprintf('%s',notes(scaleKeyLoc+length(scaleKeyStr)+nextindex:end));
    end
  end
  % fix scale factor
  v = viewSet(v,'roixform',scaleFactorMatrix*xform,i);
  v = viewSet(v,'roiVoxelSize',scaleFactor.*voxelSize,i);
  % set notes
  notes = sprintf('%s%s %f %f %f',notes,scaleKeyStr,scaleFactor(1),scaleFactor(2),scaleFactor(3));
  v = viewSet(v,'roinotes',notes,i);
  % save changed roi
  saveROI(v,i,0);
end


v = loadAnat(v);
anatomyDir = viewGet(v,'anatomyDir');
for i = 1:viewGet(v,'numbase')
  % get base name
  basename = viewGet(v,'basename',i);
  % get old hdr name
  currentBasename = fullfile(anatomyDir,sprintf('%s.hdr',basename));
  backupBasename = fullfile(anatomyDir,sprintf('%s_%s.hdr',basename,backupHeaderPostfix));
  if isfile(backupBasename)
    hdr = cbiReadNiftiHeader(backupBasename);
    if ~askuser(sprintf('(mlrFixScaleFactor) %s has aleady been scaled, fix from original xform',basename))
      disp(sprintf('(mlrFixScaleFactor) Skipping %s',basename));
      continue
    end
    xform = hdr.sform44;
    load(sprintf('%s.mat',stripext(backupBasename)));
    baseVol2mag = base.vol2mag;
  else
    % make a backup header
    system(sprintf('mv %s %s',currentBasename,backupBasename));
    system(sprintf('mv %s.mat %s.mat',stripext(currentBasename),stripext(backupBasename)));
    % get base xform
    xform = viewGet(v,'basexform',i);
    baseVol2mag = viewGet(v,'baseVol2mag',i);
  end
  % fix scale factor
  v = viewSet(v,'basexform',scaleFactorMatrix*xform,i);
  v = viewSet(v,'baseVol2mag',scaleFactorMatrix*baseVol2mag,i);
  % save changed base
  saveAnat(v,i,1);
end

deleteView(v);
