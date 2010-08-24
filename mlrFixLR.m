% mlrFixLR.m
%
%        $Id: mlrFixLR.m 965 2008-02-21 00:03:20Z justin $
%      usage: mlrFixLR(<anatdir>)
%         by: justin gardner
%       date: 02/04/2010
%    purpose: composes all the xforms found in a directory with 
%             anatomies and rois in a session           
%
function retval = mlrFixLR(type)

% keep backup headers with the following tagged on 
backupHeaderPostfix = 'preLRFix';

% check arguments
if ~any(nargin == [0 1])
  help mlrFixLR
  return
end

if nargin == 1
  fixAnatDirLR
end

% no scaling of voxels
scaleFactor = [1 1 1];

% get scale factor matrix
fixMatrix = diag([-1 1 1 1]);

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
    mrWarnDlg(sprintf('(mlrFixLR) %s is not an mrSession file',filename));
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
    mrWarnDlg(sprintf('(mlrFixLR) Could not open a view to %s',getLastDir(path)));
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
	  if ~askuser(sprintf('(mlrFixLR) %s has already been scaled. Rescale from original',backupFilename))
	    disp(sprintf('(mlrFixLR) Skipping %s',filename));
	    continue
	  end
	  hdr = cbiReadNiftiHeader(backupFilename);
	else
	  % load the header
	  hdr = cbiReadNiftiHeader(filename);
	end
	% check the current sformcode
	if hdr.sform_code ~= 1
	  disp(sprintf('(mlrFixLR) %s does not have sform set to 1',filename));
	  keyboard
	end
	% set the sform
	hdr = cbiSetNiftiSform(hdr,fixMatrix*hdr.sform44);
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
	v = viewSet(v,'scanvol2mag',fixMatrix*vol2mag,iScan,iGroup);
	
      else
	disp(sprintf('(mlrFixLR) Could not open file %s',filename));
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
    if ~askuser('(mlrFixLR) Ok to dump baseVolumes loaded into mrLoadRet viewer? If they are not saved, you may want to save them. Otherwise, they will not get scaled')
      disp(sprintf('(mlrFixLR) Keeping base anatomies in mrLastView. Note that these loaded bases will not be fixed'));
    else
      view.baseVolumes = [];
    end
  end
  if ~isempty(view.ROIs)
    if ~askuser('(mlrFixLR) Ok to dump ROIs loaded into mrLoadRet viewer? If they are not saved, you may want to save them. Otherwise, they will not get scaled')
      disp(sprintf('(mlrFixLR) Keeping ROIs in mrLastView. Note that these loaded ROIs will not be fixed'));
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
  % check for fix note
  scaleKeyStr = 'LRFix';
  scaleKeyLoc = findstr(scaleKeyStr,notes);
  if ~isempty(scaleKeyLoc)
    disp(sprintf('(mlrFixLR) Roi %s has already been LR fixed, skipping...',roiname));
    continue
  end
  % fix LR
  v = viewSet(v,'roixform',fixMatrix*xform,i);
  v = viewSet(v,'roiVoxelSize',scaleFactor.*voxelSize,i);
  % set notes
  notes = sprintf('%s%s',notes,scaleKeyStr);
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
    if ~askuser(sprintf('(mlrFixLR) %s has aleady been scaled, fix from original xform',basename))
      disp(sprintf('(mlrFixLR) Skipping %s',basename));
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
  v = viewSet(v,'basexform',fixMatrix*xform,i);
  v = viewSet(v,'baseVol2mag',fixMatrix*baseVol2mag,i);
  % save changed base
  saveAnat(v,i,1);
end

deleteView(v);

%%%%%%%%%%%%%%%%%%%%%%
%    fixAnatDirLR    %
%%%%%%%%%%%%%%%%%%%%%%
function fixAnatDirLR


