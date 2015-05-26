% dofmricni.m
%
%        $Id:$ 
%      usage: dofmricni
%         by: justin gardner
%       date: 03/17/2015
%    purpose: do initial processing stream for Stanford data brought down
%             from NIMS. Based on dofmrigru code.
%
%             To use this, just call dofmricni and select the session you
%             want to bring to your local computer. It will connect to the
%             cniComputer (you can change computer name with cniComputerName)
%             You must put in your password. You will see a list of sessions
%             When you have selected which one you want to load then
%             you will need to put your password in again so that it can 
%             copy the files locally. It will put files into the correct directories
%             for mlr (and decompress the niftis). It will also put the dicomInfo
%             into the scan
%
%                'numMotionComp=1': Set to 0 if you don't want to run MLR motion comp. Set to > 1 if you want
%                    to set multiple motionComp parameters (e.g. for motionComping two sets of scans taken at
%                    different resolutions)
%                'unwarp=1': Set to 0 if you do not want to run FSL distortion correction (which looks at
%                    images created with different phase encode directions - you need a calibration image
%                    with phase encode taken in the opposite direction as in your main data set, see 
%                    calibrationNameStrings below for how it should be named. Then FSL uses the different
%                    images to decide how how to stretch and compress image apropriately to undo distortions
%                'calibrationNameStrings={'CAL','pe0'}: If these strings are in the filename, then it will
%                   assume these are calibrations for use with the FSL unwarp correction
%                'minVolumes=10': Ignores scans that have less than this number of volumes
%                'removeInitialVols=2': Will remove this many initial volumes (steady-states) and correct stimfile
%                   to match. Note that setting to 0 will correct stimfile. Set to [] if you do not want to correct
%                   stimfile either.
%                'stimfileRemoveInitialVols=[]': Typically just set to empty and the program will calculate
%                   the right number of vols to remove from the stimfiles (this is needed because there are
%                   some initial vols that are recorded that are just calibrations). Override the setting
%                   by putting the number of vols you want to actually remove
%                'cleanUp=1': Sets whether to keep the temporary staging directory or not.
%                'useLocalData=0': Set this to 1 if you want to use already downloaded data
%                   in your temp direcotry. Or set it to a string containing the path which
%                   you want to use to get data from.
%
function retval = dofmricni(varargin)

% todo: Also would be nice to default motionComp parameters to
% what we are using these days

% Default arguments
getArgs(varargin,{'stimfileDir=[]','numMotionComp=1','cniComputerName=cnic7.stanford.edu','localDataDir=~/data','stimComputerName=oban','stimComputerUserName=gru','username=[]','unwarp=1','minVolumes=10','removeInitialVols=2','stimfileRemoveInitialVols=[]','calibrationNameStrings',{'CAL','pe0'},'cleanUp=1','useLocalData=0','spoofTriggers=1','fixMuxXform=2'});

% clear screen
clc;

% set up system variable (which gets passed around with important system info) - this
% means we have to copy these variables into s.
sParams = {'cniComputerName','username','stimComputerUserName','stimComputerName','stimfileDir','numMotionComp','minVolumes','removeInitialVols','stimfileRemoveInitialVols','calibrationNameStrings','cleanUp','useLocalData','spoofTriggers','unwarp','fixMuxXform'};
for iParam = 1:length(sParams)
  s.(sParams{iParam}) = eval(sParams{iParam});
end
s.localDataDir = mlrReplaceTilde(localDataDir);

% range for te to be considered a BOLD scan
s.teLower = 20;
s.teHigher = 40;

% check to make sure we have the computer setup correctly to run
% gunzip, FSL and other unix commands
[tf s] = checkCommands(s);
if ~tf,return,end

% check to see if we are to use downloaded data or not
% if this is set (usually not) then we will check
% the local temp directory, otherwise we will connect
% to he CNI computer
if isequal(s.useLocalData,0)
  % choose which directory to download from cni
  s = getCNIDir(s);
  if isempty(s.cniDir),return,end

  % now move data into temporary directory on local machine so that we can analyze it
  [tf s] = getCNIData(s);
  if ~tf,return,end
else
  % get downloaded data
  [tf s] = getLocalData(s);
  if ~tf,return,end
end

% make sure that MLR is not running
mlrPath mrTools;
mrQuit;
mlrPath mrTools;

% check the data
[tf s] = examineData(s);
if ~tf,return,end

% get stimfiles
[tf s] = getStimfiles(s);
if ~tf,return,end

% match stimfiles
[tf s] = matchStimfiles(s);
if ~tf,return,end

% setup FSL distortion correction
if s.unwarp
  [tf s] = doFSLunwarp(s);
  if ~tf,return,end
end

% now propose a move
[tf s] = moveAndPreProcessData(true,s);
if ~tf,return,end

% and do it
[tf s] = moveAndPreProcessData(false,s);
if ~tf,return,end

% now run mrInit
runMrInit(s);

%%%%%%%%%%%%%%%%%%%%%%
%    getLocalData    %
%%%%%%%%%%%%%%%%%%%%%%
function [tf s] = getLocalData(s)

tf = false;

% if useLocalData is string then it means to check in that directory
% for data
if isstr(s.useLocalData)
  s.localDir = mlrReplaceTilde(s.useLocalData);
else
  % list all directroies in temp and let user select
  tempDir = mlrReplaceTilde(fullfile(s.localDataDir,'temp/dofmricni'));
  tempDirListing = dir(tempDir);
  % cycle through the tempDirListing looking for downloadedData
  downloadedData = {};
  for i = 1:length(tempDirListing)
    if tempDirListing(i).isdir && ~isequal(tempDirListing(i).name(1),'.')
      downloadedData{end+1} = tempDirListing(i).name;
    end
  end
  % if now directories, then warn and return
  if isempty(downloadedData)
    disp(sprintf('(dofmricni) No downloaded data found in %s. Consider setting useLocalData equal to the directory that you want to use that contains directory',tempDir));
    return
  elseif length(downloadedData)==1
    % just one directory, so use that
    s.localDir = fullfile(tempDir,downloadedData{1});
  else
    % multiple directories, so ask the user to choose which one
    paramsInfo = {{'localDir',downloadedData,sprintf('Directories in your %s folder which may contain downloaded data.',tempDir)}};
    params = mrParamsDialog(paramsInfo,'Choose which directory you want to use');
    if isempty(params),return,end
    s.localDir = fullfile(tempDir,params.localDir);
  end
end

tf = true;

%%%%%%%%%%%%%%%%%%%%%
%%   examineData   %%
%%%%%%%%%%%%%%%%%%%%%
function [tf s] = examineData(s)

tf = false;
% get the list of filest that we have
fileList = getFileList(s.localDir);

% get dicom info
disppercent(-inf,'(dofmricni) Getting dicom info');
s.subjectID = [];
s.magnet = [];
s.operatorName = [];
s.receiveCoilName = [];
s.studyDate = '';
for i = 1:length(fileList)
  fileList(i).dicomInfo = getDicomInfo(fileList(i).dicom,s);
  fileList(i).tr = nan;
  fileList(i).te = nan;
  % pull out info from dicom header
  if ~isempty(fileList(i).dicomInfo) 
    % pull out tr
    if isfield(fileList(i).dicomInfo,'RepetitionTime')
      fileList(i).tr = fileList(i).dicomInfo.RepetitionTime;
    end
    % pull out TE
    if isfield(fileList(i).dicomInfo,'EchoTime')
      fileList(i).te = fileList(i).dicomInfo.EchoTime;
    end
    % pull out subjectID
    if isfield(fileList(i).dicomInfo,'subjectID') && ~isempty(fileList(i).dicomInfo.subjectID)
      s.subjectID = fileList(i).dicomInfo.subjectID;
    end
    % pull out date
    if isfield(fileList(i).dicomInfo,'StudyDate')
      s.studyDate = fileList(i).dicomInfo.StudyDate;
    end
    % pull out magnet
    if isempty(s.magnet)
      if isfield(fileList(i).dicomInfo,'Manufacturer')
	s.magnet = fileList(i).dicomInfo.Manufacturer;
      end
      if isfield(fileList(i).dicomInfo,'ManufacturerModelName')
	s.magnet = sprintf('%s %s',s.magnet,fileList(i).dicomInfo.ManufacturerModelName);
      end
      if isfield(fileList(i).dicomInfo,'MagneticFieldStrength')
	s.magnet = sprintf('%s %iT',s.magnet,fileList(i).dicomInfo.MagneticFieldStrength);
      end
      if isfield(fileList(i).dicomInfo,'ImagingFrequency')
	s.magnet = sprintf('%s %sMhz',s.magnet,mlrnum2str(fileList(i).dicomInfo.ImagingFrequency,'compact=1'));
      end
      if isfield(fileList(i).dicomInfo,'InstitutionName')
	s.magnet = sprintf('%s %s',s.magnet,fileList(i).dicomInfo.InstitutionName);
      end
    end
    % pull out operator
    if isfield(fileList(i).dicomInfo,'OperatorName')
      s.operatorName = fileList(i).dicomInfo.OperatorName;
      if ~isempty(s.operatorName)
	if isfield(s.operatorName,'FamilyName') && isfield(s.operatorName,'GivenName')
	  s.operatorName = strtrim(sprintf('%s %s',s.operatorName.GivenName,s.operatorName.FamilyName));
	elseif isfield(s.operatorName,'FamilyName')
	  s.operatorName = strtrim(sprintf('%s',s.operatorName.FamilyName));
	elseif isfield(s.operatorName,'GivenName')
	  s.operatorName = strtrim(sprintf('%s',s.operatorName.GivenName));
	else
	  s.operatorName = '';
	end
      end
    end
    % pull out coil
    if isfield(fileList(i).dicomInfo,'ReceiveCoilName')
      if ~isempty(fileList(i).dicomInfo.ReceiveCoilName)
	s.receiveCoilName = fileList(i).dicomInfo.ReceiveCoilName;
      end
    end
  end
  disppercent(i/length(fileList));
end
disppercent(inf);

% get scan start time
for i = 1:length(fileList)
  fileList(i).startDate = [];
  if all(isfield(fileList(i).dicomInfo,{'AcquisitionTime','AcquisitionDate'}))
    % get start time
    fileList(i).startTime = str2num(fileList(i).dicomInfo.AcquisitionTime);
    % try to get datenum
    try
      fileList(i).startDate = datestr(datenum([fileList(i).dicomInfo.AcquisitionDate fileList(i).dicomInfo.AcquisitionTime],'yyyymmddHHMMSS'));
    catch
      disp(sprintf('(dofmricni) Unrecognized date format in dicom: %s %s',fileList(i).dicomInfo.AcquisitionDate,fileList(i).dicomInfo.AcquisitionTime));
    end
      
  else
    fileList(i).startTime = inf;
  end
  fileList(i).startHour = floor(fileList(i).startTime/10000);
  fileList(i).startMin = floor(fileList(i).startTime/100)-fileList(i).startHour*100;
end

% sort by start time
fileList = sortFileList(fileList);

% get list of bold scans
boldNum = 0;
s.seriesDescription = [];
s.boldScans = [];
for i = 1:length(fileList)
  % check by whether name contains BOLD or TE is in range 
  if ~isempty(findstr('bold',lower(fileList(i).filename))) || (~isempty(fileList(i).te) && (fileList(i).te >= s.teLower) && fileList(i).te <= s.teHigher)
    fileList(i).bold = true;
    fileList(i).flipAngle = nan;
    % also get receiverCoilName and sequence type info
    if isfield(fileList(i).dicomInfo,'SeriesDescription')
      % keep the first one in the list as the series description
      if isempty(s.seriesDescription)
	s.seriesDescription = fileList(i).dicomInfo.SeriesDescription;
      end
    end
    % get mux factor from name
    muxloc = findstr('mux',lower(fileList(i).filename));
    fileList(i).mux = [];
    if ~isempty(muxloc)
      fileList(i).mux = str2num(strtok(fileList(i).filename(muxloc(1)+3:end),'_ abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'));
    end
    % update bold count
    s.boldScans(end+1) = i;
    boldNum = boldNum+1;
  else
    fileList(i).bold = false;
  end
end

% get list of anat scans
anatNum = 0;
for i = 1:length(fileList)
  % check by whether name contains t1w
  if ~isempty(findstr('t1w',lower(fileList(i).filename))) 
    fileList(i).anat = true;
    anatNum = anatNum+1;
    s.anatScans(anatNum) = i;
    fileList(i).toName = setext(sprintf('anat%02i_%s',anatNum,fileList(i).filename),fileList(i).niftiExt);
  else
    fileList(i).anat = false;
  end
end

% uncompress nifti files
disppercent(-inf,'(dofmricni) Uncompressing nifti files');
for i = 1:length(fileList)
  % BOLD scan or anat scans may need to be uncompressed
  if (fileList(i).bold || fileList(i).anat) && ~fileList(i).ignore
    % check if compressed
    if (length(fileList(i).niftiExt)>2) && strcmp(fileList(i).niftiExt(end-1:end),'gz')
      % change mode just in case
      system(sprintf('chmod 644 %s',fileList(i).nifti),'-echo');
      % make command to gunzip
      uncompressedFilename = stripext(fileList(i).nifti);
      if ~isfile(uncompressedFilename)
	system(sprintf('%s -c %s > %s',s.commands.gunzip,fileList(i).nifti,uncompressedFilename),'-echo');
      end
      % make sure that the mode is set correctly
      system(sprintf('chmod 644 %s',uncompressedFilename),'-echo');
      % change nifti name
      fileList(i).nifti = uncompressedFilename;
      fileList(i).niftiExt = getext(uncompressedFilename);
    end
  end
  disppercent(i/length(fileList));
end
disppercent(inf);

% read nifti headers
% initialze to know calibration files
s.calibrationFile = [];
disppercent(-inf,'(dofmricni) Reading nifti headers');
for i = 1:length(fileList)
  if ~isempty(fileList(i).nifti)
    system(sprintf('chmod 644 %s',fileList(i).nifti));
    fileList(i).h = mlrImageHeaderLoad(fileList(i).nifti);
  else
    fileList(i).h = [];
  end
  % get other info from description field
  if isfield(fileList(i).h,'hdr') && isfield(fileList(i).h.hdr,'descrip')
    descrip = fileList(i).h.hdr.descrip;
    % parse it
    while ~isempty(descrip)
      [thisVar descrip] = strtok(descrip,';');
      [varName varValue] = strtok(thisVar,'=');
      if ~isempty(deblank(varName))
	fileList(i).descrip.(varName) = str2num(varValue(2:end));
      end
    end
  end
  % get flip angle
  if isfield(fileList(i),'descrip') && isfield(fileList(i).descrip,'fa')
    fileList(i).flipAngle = fileList(i).descrip.fa;
  end
  % see if this is a calibration scan
  calibrationFile = false;
  for iString = 1:length(s.calibrationNameStrings)
    if ~isempty(findstr(fileList(i).filename,s.calibrationNameStrings{iString}))
      calibrationFile = true;
      s.calibrationFile(end+1) = i;
      % remove it from the bold list
      if fileList(i).bold
	fileList(i).bold = false;
	s.boldScans = setdiff(s.boldScans,i);
      end
      % set its too name
      fileList(i).toName = sprintf('%02i_CAL_%s',length(s.calibrationFile),getLastDir(fileList(i).nifti));
    end
  end
  % check if this is a bold
  if fileList(i).bold
    % if it is and we do not have a nifti header then remove from list
    % or if we have too few volumes
    if isempty(fileList(i).h) || size(fileList(i).h.dim,2) < 4 || (fileList(i).h.dim(4) < s.minVolumes)
      if isempty(fileList(i).h)
	disp(sprintf('(dofmricni) !!! Scan %s has no nifti header, removing from list of BOLD scans !!!',fileList(i).filename));
      else
	% tell  user we are dropping
	if size(fileList(i).h.dim,2) < 4
	  disp(sprintf('\n(dofmricni) !!! Scan %s has 0 volumes, removing from list of BOLD scans.',fileList(i).filename));
	else
	  disp(sprintf('\n(dofmricni) !!! Scan %s has only %i volumes, removing from list of BOLD scans. Decrease minVolumes setting in dofmricni to keep. !!!',fileList(i).filename,fileList(i).h.dim(4)));
	end
      end
      % remove from list
      fileList(i).bold = false;
      s.boldScans = setdiff(s.boldScans,i);
    end
  end
  disppercent(i/length(fileList));
end
disppercent(inf);

% set the bold scans toName
for iBOLD = 1:length(s.boldScans)
  fileList(s.boldScans(iBOLD)).toName = sprintf('%03i_BOLD_%s',iBOLD,getLastDir(fileList(s.boldScans(iBOLD)).nifti));
end

% get the subject id
if isempty(s.subjectID)
  mrParams = {{'subjectID',0,'incdec=[-1 1]','minmax=[0 inf]','Subject ID'}};
  params = mrParamsDialog(mrParams,'Set subject ID');
  if isempty(params),return,end
  s.subjectID = gruSubjectNum2ID(params.subjectID);
end

% set experiment name
dofmricniFilename = fullfile(s.localDir,'dofmricni.mat');
if isfield(s,'cniDir') && ~isempty(s.cniDir)
  % get experiment name from cniDir
  s.experimentName = fileparts(s.cniDir);
  % save the s variable, so next time we can read this value
  save(dofmricniFilename,'s');
else
  % we don't know the experiment name
  s.experimentName = 'unknown';
  % try to load it
  if isfile(dofmricniFilename)
    sold = load(dofmricniFilename);
    if isfield(sold,'s') && isfield(sold.s,'experimentName')
      s.experimentName = sold.s.experimentName;
    end
  end
  % if still unknown then ask
  if isempty(s.experimentName) || strcmp(s.experimentName,'unknown')
    paramsInfo = {{'experimentName',s.experimentName,'Name of the experiment you are running - this will be used as a directory name under data to store all of your files'}};
    params = mrParamsDialog(paramsInfo,'Set experiment name');
    if isempty(params),return,end
    s.experimentName = params.experimentName;
  end
  % save the s variable, so next time we can read this value
  save(dofmricniFilename,'s');
end

% get name of directory to copy
s.localSessionDir = fullfile(s.localDataDir,s.experimentName,sprintf('%s%s',s.subjectID,s.studyDate));

% confirm with user
mrParams = {{'localSessionDir',s.localSessionDir,'Where your data will get stored'}};
params = mrParamsDialog(mrParams,'Confirm where you want the data stored on your local computer');
if isempty(params),return,end
s.localSessionDir = params.localSessionDir;

% return true
tf = true;
s.fileList = fileList;

%%%%%%%%%%%%%%%%%%%
%%   runMrInit   %%
%%%%%%%%%%%%%%%%%%%
function tf = runMrInit(s)

tf = false;
curpwd = pwd;
cd(s.localSessionDir);

% create stimfile list to pass into mrInit
for i = 1:length(s.stimfileMatch)
  if s.stimfileMatch(i)~=0
    stimfileMatchList{i} = s.stimfileInfo(s.stimfileMatch(i)).name;
  else
    % no match
    stimfileMatchList{i} = 'None';
  end
end

% initialize parameters
[sessionParams groupParams] = mrInit([],[],'justGetParams=1','magnet',s.magnet,'operator',s.operatorName,'subject',s.subjectID,'coil',s.receiveCoilName,'pulseSequence',s.seriesDescription,'stimfileMatchList',stimfileMatchList);
if isempty(sessionParams),return,end

% now run mrInit
disp(sprintf('(dofmricni1) Setup mrInit for your directory'));
mrInit(sessionParams,groupParams,'makeReadme=0');

% now set the dicom info
v = newView;
nScans = viewGet(v,'nScans');
if ~isempty(v)
  for iScan = 1:nScans
    scanName = viewGet(v,'tseriesFile',iScan);
    % look for matching scan in fileList
    fileListNum = find(strcmp(scanName,{s.fileList(:).toName}));
    if ~isempty(fileListNum)
      % set the dicom information
      v = viewSet(v,'auxParam','dicomInfo',s.fileList(fileListNum).dicomInfo,iScan);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set framePeriod as recorded in stimfile
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stimfile = viewGet(v,'stimfile',iScan);
    if ~isempty(stimfile)
      if strcmp(stimfile{1}.filetype,'mgl')
	% get all the volume events
	volEvents = find(stimfile{1}.myscreen.events.tracenum == 1);
	if length(volEvents) > 1
	  % get the framePeriod
	  framePeriod = median(diff(stimfile{1}.myscreen.events.time(volEvents)));
	  % round to nearest 1/1000 of a second
	  framePeriod = round(framePeriod*1000)/1000;
	  disp(sprintf('(dofmricni) Frame period as recorded in stimfile is: %0.3f',framePeriod));
	end
      end
    end
  end
  saveSession(0);
  deleteView(v);
end

% run motion comp
for i = 1:s.numMotionComp
  v = newView;
  if ~isempty(v)
    [v params] = motionComp(v,[],'justGetParams=1');
    v = motionComp(v,params);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%    removeTempFiles    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function removeTempFiles(fidList)

% find all temp files to delete
deleteList = {};
fileExt = {'hdr','img','sdt','spr','edt','epr'};
for i = 1:length(fidList)
  filestem = stripext(fidList{i}.fullfile);
  for j = 1:length(fileExt)
    % check for files with that extension
    filename = setext(filestem,fileExt{j});
    if isfile(filename)
      deleteList{end+1} = filename;
    end
  end
end

% check for some other files
otherFiles = {'dofmricni2.log','logbook'};
for i = 1:length(otherFiles)
  filename = fullfile('Pre',otherFiles{i});
  if isfile(filename)
    deleteList{end+1} = filename;
  end
end

if ~isempty(deleteList)
  disp(sprintf('=============================================='));
  disp(sprintf('Temporary files'));
  disp(sprintf('=============================================='));
  for i = 1:length(deleteList)
    disp(sprintf('%i: %s',i,deleteList{i}));
  end
  disp(sprintf('=============================================='));

  if askuser('Found temporary files. Ok to remove them and place them in the directory Pre/deleteme')
    % make the directory if necessary
    if ~isdir('Pre/deleteme')
      mkdir(fullfile('Pre/deleteme'));
    end
    % move the files
    for i = 1:length(deleteList)
      moveToFilename = fullfile('Pre','deleteme',getLastDir(deleteList{i}));
      disp(sprintf('Moving %s -> %s',deleteList{i},moveToFilename));
      movefile(deleteList{i},moveToFilename);
    end
  end
end

%%%%%%%%%%%%
%% myeval %%
%%%%%%%%%%%%
function myeval(command,justDisplay)

if justDisplay
  disp(command);
else
  eval(command);
  dispConOrLog(command,justDisplay,true);
end  

%%%%%%%%%%%%%%%%%%%%%
%%   doFSLunwarp   %%
%%%%%%%%%%%%%%%%%%%%%
function [tf s] = doFSLunwarp(s,doit)

tf = true;
if nargin < 2, doit = false;end

% set unwarp to the files the calibration file we found, and 
if isempty(s.calibrationFile)
  dispConOrLog(sprintf('(dofmricni:doFSLunwarp) !!! No calibration file found. Skipping unwarping !!!!',~doit));
  return
else
  if ~isstruct(s.unwarp)
    % put the calibration files in pe1
    for i = 1:length(s.calibrationFile)
      s.unwarp.calfiles{i} = s.fileList(s.calibrationFile(i)).toName;
    end
    % check length
    if length(s.unwarp.calfiles) > 1
      dispHeader('Choose calibration scan for FSL unwarping');
      % display list
      for i = 1:length(s.unwarp.calfiles)
	calibFile = s.fileList(s.calibrationFile(i));
	disp(sprintf('%i: %s (dims=[%s] tr:%0.1f te:%0.1f start: %s)',i,calibFile.filename,mlrnum2str(calibFile.h.dim,'sigfigs=0'),calibFile.tr,calibFile.te,calibFile.startDate));
      end
      calnum = getnum('(dofmricni) Which scan do you want to use for fsl unwarping calibration scan: ',1:length(s.unwarp.calfiles));
      s.unwarp.calfiles = {s.unwarp.calfiles{calnum}};
    end
  end
end

% put bold scans into pe1
if isempty(s.boldScans)
  s = rmfield(s,'unwarp');
  return
else
  for i = 1:length(s.boldScans)
    s.unwarp.EPIfiles{i} = s.fileList(s.boldScans(i)).toName;
  end
end

% now actually do it
if doit && isfield(s,'unwarp')
  retval = fsl_pe0pe1(fullfile(s.localSessionDir,'Pre'),s.unwarp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    doRemoveInitialsVols    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = doRemoveInitialVols(s,justDisplay)

for iBOLD = 1:length(s.boldScans)
  % get scan info
  boldScan = s.fileList(s.boldScans(iBOLD));
  % get vols before and after
  nVols = boldScan.h.dim(4);
  nVolsAfter = nVols-s.removeInitialVols;
  % display which one it is
  dispConOrLog(sprintf('%i: %s (%i->%i vols)',iBOLD,boldScan.filename,nVols,nVolsAfter),justDisplay);
  % we also fix the muxXform here (since we are already saving and loading)
  if s.fixMuxXform && ~isempty(boldScan.mux) && (boldScan.mux > 1)
    % get slices/mux factor
    dicomSlices = boldScan.dicomInfo.numSlices;
    muxSlices = boldScan.h.dim(3);
    mux = boldScan.mux;
    % compute number of slices to shift
    shiftSlice = -(muxSlices-dicomSlices)/2;
    % computed shifted xform
    xform = boldScan.dicomInfo.xform;
    s.fileList(s.boldScans(iBOLD)).muxXform = xform*[1 0 0 0;0 1 0 0;0 0 1 shiftSlice;0 0 0 1];
    % display what we are doing
    dispConOrLog(sprintf('(dofmricni) Fixing xform (mux=%i nDicomSlices=%i muxSlices: %i) Shift by %0.2f slices',mux,dicomSlices,muxSlices,shiftSlice),justDisplay);
    % check xform
    if (s.fixMuxXform > 1) && (iBOLD == 1) && justDisplay
      % load epi and anat
      disppercent(-inf,sprintf('Loading epi: %s and anat: %s to check fixed header',getLastDir(boldScan.nifti),getLastDir(s.fileList(s.anatScans(1)).nifti)));
      [epid epih] = mlrImageLoad(boldScan.nifti);
      [anatd anath] = mlrImageLoad(s.fileList(s.anatScans(1)).nifti);
      disppercent(inf);
      epih.sform = s.fileList(s.boldScans(iBOLD)).muxXform;
      mlrVol(anatd,anath,epid,epih);
    end
  end
  if ~justDisplay
    % go ahead and remove them, first load the file
    filename = fullfile(s.localSessionDir,'Raw','TSeries',boldScan.toName);
    [d h] = mlrImageLoad(filename);
    % remove the appropriate number of voulems
    d = d(:,:,:,1+s.removeInitialVols:end);
    % fix header if necessary
    if isfield(boldScan,'muxXform')
      h.qform = boldScan.muxXform;
      h.sform = h.qform;
    end
    % write it back
    mlrImageSave(filename,d,h);
  end
  % if there is a matching stimfile, then what are we to do ### DAN FIX:
  % added a check for zero to avoid stimfileMatch failing to get
  % stimfileInfo ####
  if length(s.stimfileMatch) >= iBOLD && s.stimfileMatch(iBOLD)~=0
    % get stimfileInfo for this trial
    stimfileInfo = s.stimfileInfo(s.stimfileMatch(iBOLD));
    % check for missing igorded volumes
    if stimfileInfo.missingIgnoredVols
      dispConOrLog(sprintf('  Fixing missing ignored triggers by %i for associated stimfile: %s (ignoredInitialVols: %i nVols: %i)',stimfileInfo.missingIgnoredVols,stimfileInfo.name,stimfileInfo.ignoredInitialVols,stimfileInfo.numVols),justDisplay);
    end
    % it should have already been called when we matched stimfile
    % to check if we need to fixAcq, so display that info here
    if isfield(stimfileInfo,'fixAcq') && ~isempty(stimfileInfo.fixAcq) && stimfileInfo.fixAcq
      dispConOrLog(sprintf('  Fixing acq triggers = %i for associated stimfile: %s (old values->ignoredInitialVols: %i nVols: %i)',nVolsAfter,stimfileInfo.name,stimfileInfo.ignoredInitialVols,stimfileInfo.numVols),justDisplay);
    end
    % actually do the change
    if ~justDisplay 
      % fix stimfile
      fixStimfileTriggers(s,s.stimfileMatch(iBOLD),boldScan,justDisplay);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    fixStimfileTriggers    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = fixStimfileTriggers(s,stimfileNum,boldScan,justDisplay)

% get the stimfile info
stimfileInfo = s.stimfileInfo(stimfileNum);

% number of volumes we should remove
removeVols = s.stimfileRemoveInitialVols;

% get the mux factor
if isfield(boldScan,'mux') && ~isempty(boldScan.mux)
  mux = boldScan.mux;
else
  mux = 1;
end

% info about scan
nVols = boldScan.h.dim(4);
nSlices = boldScan.h.dim(3);

% check to see if we are closer to the number of volumes
% or number of volumes * slices/mux 
if abs(stimfileInfo.numVols-nVols) < abs(stimfileInfo.numVols-nVols*nSlices/mux)
  % we had a trigger for every slice
  triggerEverySlice = false;
else
  triggerEverySlice = true;
end

% check whether the correct number of initial volumes were received
% We expect mux * calibration volumes (which is usually 2 - passed in argument)
if triggerEverySlice
  calibrationTriggers = nSlices*s.removeInitialVols;
else
  calibrationTriggers = boldScan.mux*s.removeInitialVols;
end

% now check to see if it matches
missingIgnoredVols = stimfileInfo.ignoredInitialVols - calibrationTriggers;
s.stimfileInfo(stimfileNum).missingIgnoredVols = missingIgnoredVols;
if missingIgnoredVols ~= 0
  % print warning message
  if triggerEverySlice
    dispConOrLog(sprintf('(dofmricni) !!! ignoredInitialVols should have been set to nSlices*initialVolumes %ix%i=%i but was set to %i',nSlices,s.removeInitialVols,calibrationTriggers,stimfileInfo.ignoredInitialVols),justDisplay);
  else
    dispConOrLog(sprintf('(dofmricni) !!! ignoredInitialVols should have been set to mux*initialVolumes %ix%i=%i, but was set to %i',boldScan.mux,s.removeInitialVols,calibrationTriggers,stimfileInfo.ignoredInitialVols),justDisplay);
  end
end

% calculate how many acq pulses we expect
acqTriggers = nVols-s.removeInitialVols;

% and see if we have a match. If not, we will fix them
spoofTriggers = false;
if s.spoofTriggers && ((acqTriggers ~= (stimfileInfo.numVols+missingIgnoredVols)) || triggerEverySlice)
  if (stimfileInfo.ignoredInitialVols+stimfileInfo.numVols) <= 0
    % no volumes - totally busted
    dispConOrLog(sprintf('(dofmricni) !!! No volumes found in stimfile !!!'),justDisplay);
    askuser('   This stimfile has no usable information about the scan');
  else
    % display warning
    dispConOrLog(sprintf('%s (%i vols)',boldScan.filename,boldScan.h.dim(4)-s.removeInitialVols));
    dispConOrLog(sprintf('(dofmricni) !!! Stimfile expected to have %i acquisiton triggers but had %i !!!',acqTriggers,stimfileInfo.numVols+missingIgnoredVols));
    if justDisplay 
      if askuser('Do you want to fix the acq triggers')
	stimfileInfo.fixAcq = true;
	s.stimfileInfo(stimfileNum).fixAcq = 1;
      else
	stimfileInfo.fixAcq = false;
	s.stimfileInfo(stimfileNum).fixAcq = 0;
      end
    end
    % triggers needed to be added
    spoofTriggers = true;
  end
end

% this is the code that actually adjusts the stimfiles
if missingIgnoredVols || spoofTriggers 
  % get the stimfile directory
  if justDisplay
    stimfileDir = s.localDir;
  else
    stimfileDir = fullfile(s.localSessionDir,'Etc');
  end
  % and load the stimfile
  stimfileName = fullfile(stimfileDir,stimfileInfo.name);
  stimfile = load(stimfileName);
  % save original
  if ~justDisplay
    % save original
    save(sprintf('%s_original.mat',stripext(stimfileName)),'-struct','stimfile');
  end
  % fix any missing ignored vols
  if missingIgnoredVols && ~justDisplay
    if missingIgnoredVols < 0
      stimfile = removeTriggers(stimfile,1:-missingIgnoredVols);
    else
      stimfile = unignoreTriggers(stimfile,calibrationTriggers+1:stimfileInfo.ignoredInitialVols);
    end
    % update stimfileInfo
    stimfileInfo.numVols = stimfile.myscreen.volnum;
    stimfileInfo.ignoredInitialVols = stimfile.myscreen.ignoredInitialVols;
  end
  % get volume events
  e = stimfile.myscreen.events;
  volTrace = find(strcmp('volume',stimfile.myscreen.traceNames));
  volEvents = find(e.tracenum==volTrace);
  volTimes = e.time(volEvents);
  % fix the triggers by spoofing
  if spoofTriggers %&& ~justDisplay
    % fix (note there is code for actually putting back missing
    % triggers, but removed it - should be in the git repo b0ce70f on May 18, 2015)
    if stimfileInfo.fixAcq
      % remove all volumes, but the first
      stimfile = removeTriggers(stimfile,2:stimfileInfo.numVols);
      % now add volumes for where events should have happened
      framePeriod = boldScan.tr/1000;
      dispConOrLog(sprintf('(dofmricni) Spoofing volumes every %ss in %s',num2str(framePeriod),stimfileName));
      % get all the volumes we need to add, and add them
      addTimes = (volTimes(1)+framePeriod):framePeriod:volTimes(1)+(nVols-s.removeInitialVols-1)*framePeriod;
      stimfile = addVolEvents(stimfile,addTimes);
      stimfile.myscreen.modifiedDate = datestr(now);
    end
  end
  % save
  if ~justDisplay
    % save the stimfile back
    save(stimfileName,'-struct','stimfile');
  end
end

%%%%%%%%%%%%%%%%%%%%%%
%    addVolEvents    %
%%%%%%%%%%%%%%%%%%%%%%
function stimfile = addVolEvents(stimfile,addVolTimes)

% shortcut
e = stimfile.myscreen.events;

% cycle over all vol times to be added
for i = 1:length(addVolTimes)
  volTime = addVolTimes(i);
  % get the event number where we want to insert the event
  addNum = last(find(e.time(1:e.n)<volTime))+1;
  lastEntry = addNum-1;
  nextEntry = min(addNum,e.n);
  if isempty(addNum),addNum=1;lastEntry=1;end
  % get time of this event relative to the ones before it
  % and after it (i.e. 0 if same as last event, 1 if same
  % as next event.
  thisFraction = (volTime-e.time(lastEntry))/(e.time(nextEntry)-e.time(lastEntry));
  % add the event
  e.n = e.n+1;
  e.tracenum = [e.tracenum(1:addNum-1) 1 e.tracenum(addNum:end)];
  e.data = [e.data(1:addNum-1) 1 e.data(addNum:end)];
  % get interpolated ticknum
  ticknum = round(e.ticknum(lastEntry)+thisFraction*(e.ticknum(nextEntry)-e.ticknum(lastEntry)));
  e.ticknum = [e.ticknum(1:addNum-1) ticknum e.ticknum(addNum:end)];
  e.volnum = [e.volnum(1:addNum-1) e.volnum(lastEntry) e.volnum(addNum:end)+1];
  e.time = [e.time(1:addNum-1) volTime e.time(addNum:end)];
  e.force = [e.force(1:addNum-1) 1 e.force(addNum:end)];
end

% fix the volume number of all the volume events (other events are fixed by code above)
% but volume events get updated afterwords - i.e. they are logged starting at 0. 
e.volnum(e.tracenum==1) = 0:(sum(e.tracenum==1)-1);

% put update events back into structure
stimfile.myscreen.events = e;

% record times of events that have been added
if ~isfield(stimfile.myscreen,'addVolTimes')
  stimfile.myscreen.addVolTimes = [];
end
stimfile.myscreen.addVolTimes(end+1:end+length(addVolTimes)) = addVolTimes;

% update number of volumes
stimfile.myscreen.volnum = stimfile.myscreen.volnum+length(addVolTimes);

%%%%%%%%%%%%%%%%%%
%%   moveData   %%
%%%%%%%%%%%%%%%%%%
function [tf s] = moveAndPreProcessData(justDisplay,s)

clc;
curpwd = pwd;

% open the logfile
if ~justDisplay
  % make the local session direcotry
  if ~isdir(s.localSessionDir)
    mkdir(s.localSessionDir);
  end
  % change to that directory
  cd(s.localSessionDir);
  % and start a log there
  openLogfile('dofmricni.log');
end

% now display all the files in the order they were acquired
dispConOrLog(sprintf('=============================================='),justDisplay);
dispConOrLog(sprintf('File list'),justDisplay);
dispConOrLog(sprintf('=============================================='),justDisplay);

dispList(s,justDisplay);

dispConOrLog(sprintf('=============================================='),justDisplay,true);
dispConOrLog(sprintf('Make Directories'),justDisplay,true);
dispConOrLog(sprintf('=============================================='),justDisplay,true);

% list of directories to make
dirList = {'Etc','Pre','Raw','Raw/TSeries','Anatomy'};

% make them
for i = 1:length(dirList)
  if ~isdir(dirList{i})
    command = sprintf('mkdir(''%s'');',fullfile(s.localSessionDir,dirList{i}));
    myeval(command,justDisplay);
  end
end

dispConOrLog(sprintf('=============================================='),justDisplay,true);
dispConOrLog(sprintf('Copy files from staging area'),justDisplay,true);
dispConOrLog(sprintf('=============================================='),justDisplay,true);

commandNum = 0;
for i = 1:length(s.fileList)
  % BOLD scan
  if s.fileList(i).bold || any(i==s.calibrationFile)
    % check for valid nifti
    if isempty(s.fileList(i).nifti)
      dispConOrLog(sprintf('********************************************'),justDisplay,true);
      dispConOrLog(sprintf('(dofmricni) !!!! BOLD nifti file for %s is missing !!!!',s.fileList(i).filename),justDisplay,true);
      % ignore it for later
      s.fileList(i).ignore = true;
    else
      % make full path
      s.fileList(i).toFullfile = fullfile(s.localSessionDir,'Pre',s.fileList(i).toName);
      % make command to copy
      command = sprintf('copyfile %s %s f',s.fileList(i).nifti,s.fileList(i).toFullfile);
      if justDisplay,commandNum=commandNum+1;disp(sprintf('%i: %s',commandNum,command)),else,myeval(command,justDisplay);,end
    end
  % anat scan
  elseif s.fileList(i).anat
    if isempty(s.fileList(i).nifti)
      dispConOrLog(sprintf('********************************************'),justDisplay,true);
      dispConOrLog(sprintf('(dofmricni) !!!! Anatomy nifti file for %s is missing !!!!',s.fileList(i).filename),justDisplay,true);
      % ignore it for later
      s.fileList(i).ignore = true;
    else
      s.fileList(i).toFullfile = fullfile(s.localSessionDir,'Anatomy',s.fileList(i).toName);
      % make command to copy
      command = sprintf('copyfile %s %s f',s.fileList(i).nifti,s.fileList(i).toFullfile);
      if justDisplay,commandNum=commandNum+1;disp(sprintf('%i: %s',commandNum,command)),else,myeval(command,justDisplay);,end
    end
  end
end

% distortion correction with fsl
if isfield(s,'unwarp') && isfield(s.unwarp,'calfiles')
  dispConOrLog(sprintf('=============================================='),justDisplay,true);
  dispConOrLog(sprintf('FSL distortion correction'),justDisplay,true);
  dispConOrLog(sprintf('=============================================='),justDisplay,true)
  % first time, call fsl_pe0pe1 to see what it wants to do
  if justDisplay
    if isfield(s,'unwarp')
      disp(sprintf('Calibration file: %s',s.unwarp.calfiles{1}));
      for i = 1:length(s.unwarp.EPIfiles)
	disp(sprintf('%i: %s',i,s.unwarp.EPIfiles{i}));
      end
    end
  else
    % just call it
    doFSLunwarp(s,true);
  end
end
  
dispConOrLog(sprintf('=============================================='),justDisplay,true);
dispConOrLog(sprintf('Move Files from Pre to Raw'),justDisplay,true);
dispConOrLog(sprintf('=============================================='),justDisplay,true);

commandNum = 0;
for i = 1:length(s.fileList)
  % BOLD scans
  if s.fileList(i).bold
    % make full path
    s.fileList(i).preFullfile = fullfile(s.localSessionDir,'Pre',s.fileList(i).toName);
    s.fileList(i).toFullfile = fullfile(s.localSessionDir,'Raw/TSeries',s.fileList(i).toName);
    % make command to copy
    command = sprintf('movefile %s %s f',s.fileList(i).preFullfile,s.fileList(i).toFullfile);
    if justDisplay,commandNum=commandNum+1;disp(sprintf('%i: %s',commandNum,command)),else,myeval(command,justDisplay);,end
  end
end

% stimfile move to Etc directory
if ~isempty(s.stimfileInfo)
  dispConOrLog(sprintf('=============================================='),justDisplay,true);
  dispConOrLog(sprintf('Copy stimfiles'),justDisplay,true);
  dispConOrLog(sprintf('=============================================='),justDisplay,true)
  commandNum = 0;
  for i = 1:length(s.stimfileInfo)
    % make command to move
    command = sprintf('copyfile %s %s f',fullfile(s.localDir,s.stimfileInfo(i).name),fullfile(s.localSessionDir,'Etc'));
    if justDisplay,commandNum=commandNum+1;disp(sprintf('%i: %s',commandNum,command)),else,myeval(command,justDisplay);,end
  end
end

% disp the stimfile match
if ~isempty(s.stimfileMatch)
  dispConOrLog(sprintf('=============================================='),justDisplay,true);
  dispConOrLog(sprintf('Stimfile match'),justDisplay,true);
  dispConOrLog(sprintf('=============================================='),justDisplay,true)
  dispStimfileMatch(s,s.stimfileMatch,false);
end

% remove initial volumes
if ~isempty(s.removeInitialVols)
  dispConOrLog(sprintf('=============================================='),justDisplay,true);
  dispConOrLog(sprintf('Remove %i initial (steady-state) volumes',s.removeInitialVols),justDisplay,true);
  dispConOrLog(sprintf('=============================================='),justDisplay,true)
  s = doRemoveInitialVols(s,justDisplay);
end

% clean up
dispConOrLog(sprintf('=============================================='),justDisplay,true);
dispConOrLog(sprintf('Clean up'),justDisplay,true);
dispConOrLog(sprintf('=============================================='),justDisplay,true);

if s.cleanUp
  command = sprintf('rm -rf %s',s.localDir);
  if justDisplay,commandNum=commandNum+1;disp(sprintf('%i: %s',commandNum,command)),else,mysystem(command);,end
else
  dispConOrLog(sprintf('Keeping temporary files in %s',s.localDir));
end

dispConOrLog(sprintf('=============================================='),justDisplay,true);
dispConOrLog(sprintf('Done'),justDisplay,true);
dispConOrLog(sprintf('=============================================='),justDisplay,true);

if ~justDisplay
  closeLogfile
  if ~isequal(s.localDir,curpwd)
    cd(curpwd);
  end
end

% ask user if this is ok
if justDisplay
  tf = askuser('(dofmricni) Ok to do the above');
else
  tf = true;
end

%%%%%%%%%%%%%%%%%%
%%   dispList   %%
%%%%%%%%%%%%%%%%%%
function dispList(s,justDisplay)

if nargin < 2,justDisplay = true;end

dispStr = {};
for i = 1:length(s.fileList)
  if ~isinf(s.fileList(i).startTime)
    if ~isempty(s.fileList(i).h)
      pixdim = s.fileList(i).h.pixdim;
      dim = s.fileList(i).h.dim;
      flipAngle = s.fileList(i).flipAngle;
      dispConOrLog(sprintf('%i) %02i:%02i %s [%s] [%s] TR: %0.2f TE: %s flipAngle: %s -> %s',i,s.fileList(i).startHour,s.fileList(i).startMin,s.fileList(i).filename,mlrnum2str(pixdim,'compact=1'),mlrnum2str(dim,'compact=1','sigfigs=0'),s.fileList(i).tr,mlrnum2str(s.fileList(i).te,'compact=1'),mlrnum2str(s.fileList(i).flipAngle,'compact=1'),s.fileList(i).toName),justDisplay);
    else
      dispConOrLog(sprintf('%i) %02i:%02i %s TR: %0.2f TE: %s -> %s',i,s.fileList(i).startHour,s.fileList(i).startMin,s.fileList(i).filename,s.fileList(i).tr,mlrnum2str(s.fileList(i).te,'compact=1'),s.fileList(i).toName),justDisplay);
      
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%
%%   sortFileList  %%
%%%%%%%%%%%%%%%%%%%%%
function fileList = sortFileList(fileList)

% sort by time stamp 
for i = 1:length(fileList)
  for j = 1:length(fileList)-1
    if fileList(j).startTime > fileList(j+1).startTime
      temp = fileList(j);
      fileList(j) = fileList(j+1);
      fileList(j+1) = temp;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%
%%   getFileList   %%
%%%%%%%%%%%%%%%%%%%%%
function fileList = getFileList(dirname)

fileList = [];

% open the directory
dirList = dir(dirname);
if isempty(dirList),return,end

% now go through the directory and fill in some information about what we have
disppercent(-inf,'(dofmricni) Getting file list');
for i = 1:length(dirList)
  match = 0;
  % skip all . files
  if dirList(i).name(1) == '.',continue,end
  % keep the name and date of the file
  fileList(end+1).filename = dirList(i).name;
  fileList(end).fullfile = fullfile(dirname,dirList(i).name);
  fileList(end).date = dirList(i).date;
  fileList(end).datenum = dirList(i).datenum;
  % set ignore to false (this gets set if something goes wrong)
  fileList(end).ignore = false;
  % get name of nifti
  % first look for uncompressed nifti
  fileList(end).nifti = dir(sprintf('%s/*.nii',fullfile(dirname,dirList(i).name)));
  if isempty(fileList(end).nifti)
    % now look for compressed
    fileList(end).nifti = dir(sprintf('%s/*.nii.gz',fullfile(dirname,dirList(i).name)));
  end
  if ~isempty(fileList(end).nifti)
    fileList(end).nifti = strtrim(fullfile(dirname,dirList(i).name,fileList(end).nifti(1).name));
  end
  % get the nifti extension
  if ~isempty(fileList(end).nifti)
    ext = getext(fileList(end).nifti);
    if strcmp(lower(ext),'gz')
      ext = sprintf('%s.gz',getext(stripext(fileList(end).nifti)));
    end
    fileList(end).niftiExt = ext;
  end
  % get name of dicom
  try
    % look for an uncompress dicom directory
    dicomDir = dir(sprintf('%s/*_dicoms',fullfile(dirname,dirList(i).name)));
    if ~isempty(dicomDir)
      fileList(end).dicom = fullfile(dirname,dirList(i).name,dicomDir(1).name);
    else
      % if not looked for a compressed zip file
      fileList(end).dicom = dir(sprintf('%s/*_dicoms.tgz',fullfile(dirname,dirList(i).name)));
      if ~isempty(fileList(end).dicom)
	fileList(end).dicom = fullfile(dirname,dirList(i).name,fileList(end).dicom.name);
      end
    end
  catch
    fileList(end).dicom = [];
  end
  disppercent(i/length(dirList));
end
disppercent(inf);

%%%%%%%%%%%%%%%%%%%%%%
%%   getDicomInfo   %%
%%%%%%%%%%%%%%%%%%%%%%
function info = getDicomInfo(filename,s)

info = [];
if isempty(filename),return,end

%unzip if necessary
if getext(filename,'tgz')
  % make sure the directory is writable
  dirName = fileparts(filename);
  system(sprintf('chmod 755 %s',fileparts(filename)));

  % change to direcotry
  curpwd = pwd;
  cd(dirName);

  % unzip the dicom
  system(sprintf('%s xf %s',s.commands.tar,getLastDir(filename)));

  % change path back
  cd(curpwd);
end

% now do a dir to look at all the files and select the first dicom
d = dir(fullfile(stripext(filename),'*.dcm'));

% if we got one, then load it
if length(d) >= 1
  info = dicominfo(fullfile(stripext(filename),d(1).name));
end

% get the subjectID
subjectID = [];
if isfield(info,'PatientName') && isfield(info.PatientName,'FamilyName')
  subjectID = info.PatientName.FamilyName;
  if ~isempty(subjectID) && (~any(length(subjectID) == [4 5]) || ~isequal(lower(subjectID(1)),'s'))
    mrWarnDlg('(dofmricni) PatientName should always be set to a subjectID (not the real name)!!!');
    subjectID = [];
  end
end
if isempty(subjectID)
  if isfield(info,'PatientName') && isfield(info.PatientName,'GivenName')
    subjectID = info.PatientName.GivenName;
    if ~isempty(subjectID) && (~any(length(subjectID) == [4 5]) || ~isequal(lower(subjectID(1)),'s'))
      mrWarnDlg('(dofmricni) PatientName should always be set to a subjectID (not the real name)!!!');
    end
  end
end
info.subjectID = subjectID;

% scrub all fields that begin with Patient so that
% we don't keep around any identifiers
infoFieldNames = fieldnames(info);
for i = 1:length(infoFieldNames)
  if strncmp(lower('Patient'),lower(infoFieldNames{i}),7)
    info = rmfield(info,infoFieldNames{i});
  end
end

% if we are trying to fix the mux x form, then load a sequence, making sure we get the one with the top slice number
if s.fixMuxXform && (length(d)>1)
  % cycle over reading of dicom headers until we get one as being the top slice
  % (this way we do not assume anything about slice order)
  foundFirstSlice = 0;
  % get the first xform
  refDicomNum = 1;
  refXform = dicom2xform(fullfile(stripext(filename),d(refDicomNum).name));
  % start by assuming this is the top slice
  topSliceXform = refXform;
  minSliceNum = 0;
  sliceNum = [];
  % start by loading each xform
  iDicom = refDicomNum+1;
  while foundFirstSlice < 2
    % get the slice x form for this dicom
    sliceXform = dicom2xform(fullfile(stripext(filename),d(iDicom).name));
    % find out the slice number
    slice2slice = inv(refXform)*sliceXform;
    sliceNum(end+1) = round(slice2slice(3,4));
    % if slice number is the lowest, then replace
    if sliceNum(end) < minSliceNum
      topSliceXform = sliceXform;
      minSliceNum = sliceNum(end);
    end
    iDicom = iDicom+1;
    % stop when we cycle back to the same reference slice
    if (sliceNum(end) == 0) 
      % increment how many times we found first slice
      foundFirstSlice = foundFirstSlice+1;
    end
    % check that we have not run out
    if iDicom>length(d),foundFirstSlice = inf;end
  end
  % set it in the return
  info.xform = topSliceXform;
  info.numSlices = length(unique(sliceNum));
end

%%%%%%%%%%%%%%%%%%
%%   mysystem   %%
%%%%%%%%%%%%%%%%%%
function status = mysystem(commandName)

% display to buffer
disp('=================================================================');
disp(sprintf('%s',datestr(now)));
disp(commandName);
disp('=================================================================');
% run the command
[status result] = system(commandName);
disp(result);

% write into the logfile
writeLogFile('\n=================================================================\n');
writeLogFile(sprintf('%s\n',datestr(now)));
writeLogFile(sprintf('%s\n',commandName));
writeLogFile('=================================================================\n');
writeLogFile(result);

%%%%%%%%%%%%%%%%%%%%%
%%   openLogfile   %%
%%%%%%%%%%%%%%%%%%%%%
function openLogfile(filename)

global gLogfile;

% open logfile
gLogfile.fid = fopen(filename,'w');
if gLogfile.fid == -1
  disp(sprintf('(dofmricni) Could not open logfile %s',filename))
  return
end

% remember filename
gLogfile.filename = filename;

%%%%%%%%%%%%%%%%%%%%%%
%%   closeLogfile   %%
%%%%%%%%%%%%%%%%%%%%%%
function closeLogfile

global gLogfile;

% close logfile
if isfield(gLogfile,'fid') && (gLogfile.fid ~= -1)
  fclose(gLogfile.fid);
end

%%%%%%%%%%%%%%%%%%%%%%
%%   writeLogFile   %%
%%%%%%%%%%%%%%%%%%%%%%
function writeLogFile(text)

global gLogfile;

if isfield(gLogfile,'fid') && (gLogfile.fid ~= -1)
  fprintf(gLogfile.fid,text);
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispConeOrLog    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispConOrLog(textStr,justDisplay,alsoDisplay)

if nargin < 2,justDisplay = true;end
if nargin < 3,alsoDisplay = false;end
if justDisplay
  disp(textStr);
else
  writeLogFile(sprintf('%s\n',textStr));
  if alsoDisplay
    disp(textStr);
  end
end
%%%%%%%%%%%%%%%%%%%%%%%
%    checkCommands    %
%%%%%%%%%%%%%%%%%%%%%%%
function [retval s] = checkCommands(s)

retval = true;

% check mgl
%if isempty(which('mglOpen')) 
%  disp(sprintf('(dofmrcni) You need to install mgl'));
%  retval = false;
%  return
%end

% check mrTools
if isempty(which('mlrVol')) 
  disp(sprintf('(dofmrcni) You need to install mrTools'));
  retval = false;
  return
end

% commands to check
preferredCommandNames = {'/usr/bin/tar','/usr/bin/gunzip'};
commandNames = {'tar','gunzip'};
helpFlag = {'-h','-h'};
[retval s] = checkShellCommands(s,commandNames,preferredCommandNames,helpFlag);
if ~retval,return,end

% needs fsl
if s.unwarp
  preferredCommandNames = {'fslroi','fslmerge','topup','applytopup'};
  commandNames = {'fslroi','fslmerge','topup','applytopup'};
  helpFlag = {'-h','-h','-h','-h'};
  [retval s] = checkShellCommands(s,commandNames,preferredCommandNames,helpFlag);
  if ~retval
    disp(sprintf('(dofmricni) !!! FSL does not appear to be installed correctly !!!'));
    if askuser('Do you want to continue to run w/out FSL? This just means it will skip the distortion correction')
      retval = true;
      s.unwarp = false;
    else
      return
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    checkShellCommands    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [retval s] = checkShellCommands(s,commandNames,preferredCommandNames,helpFlag)

retval = true;

for i = 1:length(commandNames)
  % check if the preferred command exists
  [commandStatus commandRetval] = system(sprintf('which %s',preferredCommandNames{i}));
  if commandStatus==0
    s.commands.(commandNames{i}) = preferredCommandNames{i};
  else
    % not found, so just look for any old command
    [commandStatus commandRetval] = system(sprintf('which %s',commandNames{i}));
    if commandStatus==0
      s.commands.(commandNames{i}) = commandNames{i};
    else
      % could not find anything. Error and give up
      disp(sprintf('(dofmricni) Could not find command: %s',commandNames{i}));
      disp(sprintf('            See http://gru.stanford.edu/doku.php/gruprivate/stanford#computer_setup for help setting up your computer'));
      retval = 0;
      return
    end
  end

  % run the command to see what happens
  [commandStatus commandRetval] = system(sprintf('%s %s',s.commands.(commandNames{i}),helpFlag{i}));
  % check for commandStatus error
  if commandStatus>1
    disp(commandRetval);
    disp(sprintf('(dofmricni) Found command: %s, but running gave an error',s.commands.(commandNames{i})));
    disp(sprintf('            See http://gru.brain.riken.jp/doku.php?id=grupub:dofmricni for help setting up your computer'));
    retval = 0;
    return
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    setStimfileListDispStr    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimfileList = setStimfileListDispStr(stimfileList)

for i = 1:length(stimfileList);
  % try to load it
  if isfile(stimfileList{i}.fullfile)
    stimfile = load(stimfileList{i}.fullfile);
    if ~isfield(stimfile,'myscreen')
      stimfileList{i}.dispstr = sprintf('%s (!!!No myscreen variable!!!)',stimfileList{i}.filename);
      continue;
    end
    myscreen = stimfile.myscreen;
    % make a string of some info myscreen
    stimfileStr = stimfileList{i}.filename;
    if isfield(myscreen,'volnum')
      stimfileStr = sprintf('%s [%i vols] ',stimfileStr,myscreen.volnum);
    end
    if isfield(myscreen,'starttime')
      stimfileStr = sprintf('%s%s ',stimfileStr,myscreen.starttime);
    end
    if isfield(myscreen,'endtime')
      stimfileStr = sprintf('%s(End: %s) ',stimfileStr,myscreen.endtime);
    end
    stimfileList{i}.dispstr = stimfileStr;
  end
end

%%%%%%%%%%%%%%%%%%%
%%   getCNIDir   %%
%%%%%%%%%%%%%%%%%%%
function s = getCNIDir(s)

s.cniDir = [];

% get the username
if isempty(s.username)
  s.username = getusername;
end

% default sunetID to be username
s.sunetID = mglGetParam('sunetID');
if isempty(s.sunetID),s.sunetID = s.username;,end

% put up dialog making sure info is correct
mrParams = {{'cniComputerName',s.cniComputerName,'The name of the computer to ssh into'},...
	    {'sunetID',s.sunetID,'Your sunet ID'}};
params = mrParamsDialog(mrParams,'Login information');
if isempty(params),return,end

% save sunetID
if ~isempty(params.sunetID) mglSetParam('sunetID',params.sunetID,1);end

% get some variables into system variable
s.sunetID = params.sunetID;
s.cniComputerName = params.cniComputerName;
  
% get the list of directoris that live on the cni computer
result = doRemoteCommand(s.sunetID,s.cniComputerName,'/home/jlg/bin/gruDispData');
if isempty(result),return,end

% parse the results
cniDir = [];
while ~isempty(result)
  % get one line
  [thisLine result] = strtok(result,10);
  % try to get the dirname. Should be "dirname,scan:scan:" etc.
  [thisDirName scanNames] = strtok(thisLine,',');
  % if we got something followed by scan names, keep going
  if ~isempty(scanNames) && isempty(strfind(thisDirName,' '))
    % strip comma from scanNames
    if length(scanNames) > 1
      scanNames = scanNames(2:end);
    end
    % get scanNames
    scanNames = textscan(scanNames,'%s','Delimiter',':');
    scanNames = scanNames{1};
    % if we have a dir and scan names keep it
    if ~isempty(thisDirName) && ~isempty(scanNames)
      cniDir(end+1).dirName = thisDirName;
      cniDir(end).scanNames = scanNames;
    end
  end
end

% check that we found something
if isempty(cniDir)
  disp(sprintf('(dofmricni) Could not find any studies'));
  return
end

% get list of all studies
allStudies = {};
for iDir = 1:length(cniDir)
  allStudies = {allStudies{:} cniDir(iDir).scanNames{:}};
end
% sort into reverse cronological order
allStudies = fliplr(sort(allStudies));

% limit to last 25 studies (so we do not get too long a list)
maxAllStudies = min(25,length(allStudies));
allStudies = {allStudies{1:min(maxAllStudies,end)}};

% Now set up variables to have the default list be from all studies
dirNames = {sprintf('From any of last %i studies',maxAllStudies),cniDir(:).dirName};
scanNames = {allStudies cniDir(:).scanNames};
mrParams = {{'chooseNum',1,'minmax',[1 length(dirNames)],'incdec=[-1 1]'},...
	    {'studyName',dirNames,'type=string','Name of studies','group=chooseNum','editable=0'},...
	    {'scanName',scanNames,'Name of scan','group=chooseNum'}};
params = mrParamsDialog(mrParams);
if isempty(params),return,end

% get the scan name at top of list
scanName = params.scanName{params.chooseNum};

% get the directory
if params.chooseNum == 1
  % got to find dir name if they selected from the all studies list
  for iDir = 1:length(cniDir)
    if any(strcmp(scanName,cniDir(iDir).scanNames))
      dirName = cniDir(iDir).dirName;
    end
  end
else
  % otherwise it just the parm
  dirName = params.studyName{params.chooseNum};
end

% set cniDir in system variable
s.cniDir = fullfile(dirName,scanName);
disp(sprintf('(dofmricni:getCNIDir) Directory chosen is: %s',s.cniDir))

% set the directory to which we resync data
toDir = mlrReplaceTilde(fullfile(s.localDataDir,'temp/dofmricni'));
if ~isdir(toDir)
  try
    mkdir(toDir);
  catch
    mrWarnDlg(sprintf('(dofmricni) Cannot make directory %s. Either you do not have permissions, or perhaps you have a line to a drive that is not currently mounted?',toDir));
    return
  end
end
s.localDir = fullfile(toDir,getLastDir(s.cniDir));

%%%%%%%%%%%%%%%%%%%%
%%   getCNIData   %%
%%%%%%%%%%%%%%%%%%%%
function [tf s] = getCNIData(s)

tf = false;

% Tell user what is going on
dispHeader;
disp(sprintf('Copying files from %s to %s',s.cniDir,s.localDir));
disp(sprintf('This could take some time. Using rsync, so that if you quit in the middle'))
disp(sprintf('You can continue where you left off by running dofmricni again'))
dispHeader;

% get dicoms
fromDir = fullfile('/nimsfs/raw/jlg',s.cniDir);
disp(sprintf('(dofmricni) Get files'));
command = sprintf('rsync -rtv --progress --size-only --exclude ''*Screen_Save'' --exclude ''*_pfile*'' --exclude ''*.pyrdb'' --exclude ''*.json'' --exclude ''*.png'' %s@%s:/%s/ %s',s.sunetID,s.cniComputerName,fromDir,s.localDir);
disp(command);
system(command,'-echo');

% got here, so everything is good
tf = true;

%%%%%%%%%%%%%%%%%%%%%
%%   getStimfiles  %%
%%%%%%%%%%%%%%%%%%%%%
function [tf s] = getStimfiles(s)

tf = false;

% stimfile info
s.stimfileInfo = {};

% get stimfile stem
s.stimfileStem = s.studyDate(3:end);

% if not using local directory, go find the stimfiles
if ~s.useLocalData
  % first get experiment folders on stimulus computer
  dataListing = doRemoteCommand(s.stimComputerUserName,s.stimComputerName,sprintf('find data ''-type'' d ''-print'''));

  % look for correct directory
  stimfileListing = [];
  while (~isempty(dataListing))
    [thisListing,dataListing] = strtok(dataListing);
    % find data directory
    [dataDir,thisListing] = strtok(thisListing,filesep);
    [dataDir,thisListing] = strtok(thisListing,filesep);
    if ~isempty(dataDir)
      % get subject dir
      [subjectDir,thisListing] = strtok(thisListing,filesep);
      if ~isempty(subjectDir)
	% check for match
	if isequal(gruSubjectID2num(subjectDir),gruSubjectID2num(s.subjectID)) && strcmp(lower(dataDir),lower(s.experimentName))
	  % now try to get the stimfiles
	  s.stimDataDir = fullfile('data',dataDir,subjectDir,s.stimfileStem);
	  stimfileListing = getRemoteListing(s.stimComputerUserName,s.stimComputerName,sprintf('%s*.mat',s.stimDataDir));
	  break;
	end
      end
    end
  end

  if ~isempty(stimfileListing)
    % ask user if these are correct
    paramsInfo = {};
    if strcmp(questdlg(sprintf('Found %i stimfiles in directory: %s. Use these?',length(stimfileListing),s.stimDataDir),'Confirm stimfile directory','Yes','No','Yes'),'Yes')
								 % get the files
								 getRemoteFiles(s.stimComputerUserName,s.stimComputerName,sprintf('%s*.mat',s.stimDataDir),s.localDir);
    else
      stimfileListing = [];
    end
  end

  % if could not find from above, then ask user to input a name
  while isempty(stimfileListing)
    paramsInfo = {};
    paramsInfo{end+1} = {'computerName',s.stimComputerName,'Name of computer where files are located'};
    paramsInfo{end+1} = {'computerUserName',s.stimComputerUserName,'Name of user on computer for ssh login'};
    paramsInfo{end+1} = {'stimfileStem',s.stimfileStem,'The base of the stimfile names'};
    paramsInfo{end+1} = {'dataDir',fullfile('data',s.experimentName,s.subjectID),'Name of directory on %s where stimfiles are located'};
    params = mrParamsDialog(paramsInfo,sprintf('Indicate location of stimfiles'));
    if isempty(params),break,end
    stimfileListing = getRemoteListing(params.computerUserName,params.computerName,fullfile(params.dataDir,sprintf('%s*.mat',params.stimfileStem)));
    if ~isempty(stimfileListing)
      % get the listing
      stimfileListing = getRemoteFiles(params.computerUserName,params.computerName,fullfile(params.dataDir,sprintf('%s*.mat',params.stimfileStem)),s.localDir);
      s.stimComputerName = params.computerName;
      s.stimComputerUserName = params.computerUserName;
      s.stimDataDir = params.dataDir;
      s.stimfileStem = params.stimfileStem;
    end
  end
else
  % local files then just get a listing
  stimfileDir = dir(fullfile(s.localDir,sprintf('%s_stim*.mat',s.stimfileStem)));
  stimfileListing = {stimfileDir(:).name};
end

% examine the stimfiles that we have
if ~isempty(stimfileListing)
  for i = 1:length(stimfileListing)
      % Another bugfix by dan: ###
      if ~isempty(strfind(stimfileListing{i},'.mat'))
        % remember name
        s.stimfileInfo(end+1).name = stimfileListing{i};
        % load the file
        stimfile = load(fullfile(s.localDir,stimfileListing{i}));
        % figure out tr from stimfile
        volTrace = find(strcmp('volume',stimfile.myscreen.traceNames));
        e = stimfile.myscreen.events;
        s.stimfileInfo(end).tr = median(diff(e.time(e.tracenum==volTrace)));
        % get some other info
        s.stimfileInfo(end).startTime = stimfile.myscreen.starttime;
        s.stimfileInfo(end).endTime = stimfile.myscreen.endtime;
        s.stimfileInfo(end).numVols = stimfile.myscreen.volnum;
        if isfield(stimfile.myscreen,'ignoredInitialVols')
          s.stimfileInfo(end).ignoredInitialVols = (stimfile.myscreen.ignoredInitialVols-stimfile.myscreen.ignoreInitialVols);
        else
          s.stimfileInfo(end).ignoredInitialVols = 0;
        end
      end
  end
end

% basically, always return true even if no stimfiles found - since
% we can always get them later
tf = true;

%%%%%%%%%%%%%%%%%%%%%%%
%%   matchStimfiles  %%
%%%%%%%%%%%%%%%%%%%%%%%
function [tf s] = matchStimfiles(s)

tf = false;
clc;
% if we have non-zero stimfiles and non-zero bold
stimfileMatch = [];
if length(s.boldScans) && length(s.stimfileInfo)
  % get stimfileInfo
  stimfileInfo = s.stimfileInfo;
  availableStimfiles = 1:length(stimfileInfo);
  % cycle over bold scans
  for iBOLD = 1:length(s.boldScans)
    % try to match to stimfile which is closest in time
    if ~isempty(availableStimfiles)
      timeDiff = inf(1,length(stimfileInfo));
      for iStimfile = availableStimfiles
	if ~isempty(s.fileList(s.boldScans(iBOLD)).startDate)
	  timeDiff(iStimfile) = abs(datenum(stimfileInfo(iStimfile).startTime)-datenum(s.fileList(s.boldScans(iBOLD)).startDate));
	else
	  timeDiff(iStimfile) = inf;
	end
      end
      % find closest match in time
      [minTime stimfileMatch(iBOLD)] = min(timeDiff);
      % remove that from the available list
      availableStimfiles = setdiff(availableStimfiles,stimfileMatch(iBOLD));
    else
      stimfileMatch(iBOLD) = 0;
    end
  end
end

if ~isempty(s.stimfileInfo)
  isMatched = false;
  while ~isMatched
    % show match
    dispStimfileMatch(s,stimfileMatch,true)
    % ask user if this is ok
    if askuser(sprintf('(dofmricni) Do you accept this matching of bold files to stimfiles?'))==1
      break;
    else
      stimfileMatch = [];
      while length(stimfileMatch) ~= length(s.boldScans)
	% get what the user wants
	stimfileMatch = getnum(sprintf('(dofmricni) Enter a numeric array for the stimfiles you want to match to each of these bold scans. Set entries to 0 for bold scans you do not want to match to any stimfile. This must be an array of length %i (set to -1 if you do not want to match): ',length(s.boldScans)));
	if isequal(stimfileMatch,-1),isMatched=true;stimfileMatch=[];break,end
	% check length
	if any(stimfileMatch>length(s.stimfileInfo)) || any(stimfileMatch<0)
	  disp(sprintf('(dofmricni) !!! Each element must be 0 or a number between 1:%i !!!',length(s.stimfileInfo)))
	  stimfileMatch = [];
	end
      end
    end
  end
end

% set in structure
s.stimfileMatch = stimfileMatch;
tf = true;

% check for correct number of volumes
if ~isempty(stimfileMatch)
  clc
  for iBOLD = 1:length(s.boldScans)
    % get scan info
    boldScan = s.fileList(s.boldScans(iBOLD));
    % if there is a matching stimfile, then check whether it needs to have
    % its volumes fixed
    if length(s.stimfileMatch) >= iBOLD && s.stimfileMatch(iBOLD)~=0
      % get the stimfile info
      s = fixStimfileTriggers(s,s.stimfileMatch(iBOLD),boldScan,1);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dispStimfileMatch    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispStimfileMatch(s,stimfileMatch,dispAvailable)

if isempty(stimfileMatch),return,end

if nargin == 2,dispAvailable=true;end
if dispAvailable
  % show available stimfiles
  dispHeader('Available stimfiles');
  for iStimfile = 1:length(s.stimfileInfo)
    stimfile = s.stimfileInfo(iStimfile);
    disp(sprintf('%i: %s (%i vols, %i ignored vols %s %s)',iStimfile,stimfile.name,stimfile.numVols,stimfile.ignoredInitialVols,stimfile.startTime,stimfile.endTime));
  end
end

% display header only if we are showing dispAvailable
if dispAvailable
  dispHeader('Proposed match');
end

% propose the match to the user
for iBOLD = 1:length(s.boldScans)
  % get the scan info
  bold = s.fileList(s.boldScans(iBOLD));
  % if we have a match, display it
  if stimfileMatch(iBOLD) ~= 0
    % stimfile info for this stimfile
    stimfile = s.stimfileInfo(stimfileMatch(iBOLD));
    disp(sprintf('%i->%i: %s (%i vols %s) -> %s (%i vols %s %s)',iBOLD,stimfileMatch(iBOLD),bold.filename,bold.h.dim(4),bold.startDate,stimfile.name,stimfile.numVols,stimfile.startTime,stimfile.endTime));
  else
    % otherwise display no match
    disp(sprintf('%i->X: %s (%i vols %s) -> None',iBOLD,bold.filename,bold.h.dim(4),bold.startDate));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%
%    getRemoteFiles  %%
%%%%%%%%%%%%%%%%%%%%%%%
function retval = getRemoteFiles(username,computerName,fromFiles,toDir)

% copy them
[status,retval] = system(sprintf('scp ''%s@%s:%s'' %s',username,computerName,fromFiles,toDir),'-echo');

% bad status, return
if status
  retval = [];
else
  % get listings (the above gets the output of scp which is not as useful)
  retval = getRemoteListing(username,computerName,fromFiles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%    getRemoteListing  %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = getRemoteListing(username,computerName,listDir)

retval = [];
lsRetval = doRemoteCommand(username,computerName,sprintf('ls ''%s''',listDir));
% check for empty
badMatchStr = 'ls: No match.';
if strncmp(retval,badMatchStr,length(badMatchStr))
  retval = [];
else
  % make into a cell array
  while ~isempty(lsRetval)
    [thisFilename lsRetval] = strtok(lsRetval);
    if ~isempty(thisFilename)
      retval{end+1} = getLastDir(thisFilename);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%    doRemoteCommand    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = doRemoteCommand(username,computerName,commandName)

retval = [];
command = sprintf('ssh %s@%s %s',username,computerName,commandName);
disp(sprintf('(dofrmicni) Doing remote command: %s',command));
disp(sprintf('If you have not yet set passwordless ssh (see: http://gru.stanford.edu/doku.php/gruprivate/sshpassless) then enter your password here: ',computerName));
[status,retval] = system(command,'-echo');
if status~=0
  disp(sprintf('(dofmricni) Could not ssh in to do remote command on: %s@%s',username,computerName));
  return
end
disp(sprintf('(dofmricni) Remote command on %s successful',computerName));


%%%%%%%%%%%%%%%%%%%%%
%%   getusername   %%
%%%%%%%%%%%%%%%%%%%%%
% getusername.m
%
%      usage: getusername.m()
%         by: justin gardner
%       date: 09/07/05
%
function username = getusername()

[retval username] = system('whoami');
% sometimes there is more than one line (errors from csh startup)
% so need to strip those off
username = strread(username,'%s','delimiter','\n');
username = username{end};
if (retval == 0)
  % get it again
  [retval username2] = system('whoami');
  username2 = strread(username2,'%s','delimiter','\n');
  username2 = username2{end};
  if (retval == 0)
    % find the matching last characers
    % this is necessary, because matlab's system command
    % picks up stray key strokes being written into
    % the terminal but puts those at the beginning of
    % what is returned by stysem. so we run the
    % command twice and find the matching end part of
    % the string to get the username
    minlen = min(length(username),length(username2));
    for k = 0:minlen
      if (k < minlen)
	if username(length(username)-k) ~= username2(length(username2)-k)
	  break
	end
      end
    end
    if (k > 0)
      username = username(length(username)-k+1:length(username));
      username = lower(username);
      username = username(find((username <= 'z') & (username >= 'a')));
    else
      username = 'unknown';
    end
  else
    username = 'unknown';
  end
else
  username = 'unknown';
end

