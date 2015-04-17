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
%                'stimfileDir=/usr1/justin/data/s00620101001/stimfile': Set this if you want to load the first pass
%                    stimfiles form a specific directory.
%                'numMotionComp=1': Set to 0 if you don't want to run MLR motion comp. Set to > 1 if you want
%                    to set multiple motionComp parameters (e.g. for motionComping two sets of scans taken at
%                    different resolutions)
%
%
function retval = dofmricni(varargin)

% todo: stimfile processing. Also would be nice to default motionComp parameters to
% what we are using these days

% Default arguments
getArgs(varargin,{'stimfileDir=[]','numMotionComp=1','cniComputerName=cnic7.stanford.edu','localDataDir=~/data'});

clc;

% set up system variable (which gets passed around with important system info)
s.cniComputerName = cniComputerName;
s.localDataDir = mlrReplaceTilde(localDataDir);
s.stimfileDir = stimfileDir;
s.numMotionComp = numMotionComp;
s.dispNiftiHeaderInfo = false;
% range for te to be considered a BODL scan
s.teLower = 25;
s.teHigher = 35;

% check to make sure we have the computer setup correctly to run epibsi, postproc and sense
[tf s] = checkCommands(s);
if ~tf,return,end

% choose which directory to download from cni
s = getCNIDir(s);
if isempty(s.cniDir),return,end

% now move data into temporary directory on local machine so that we can analyze it
[tf s] = getCNIData(s);
if ~tf,return,end

% make sure that MLR is not running
mlrPath mrTools;
mrQuit;
mlrPath mrTools;

% check the data
[tf s] = examineData(s);
if ~tf,return,end

% now propose a move
[tf s] = moveData(true,s);
if ~tf,return,end

% and do it
[tf s] = moveData(false,s);
if ~tf,return,end

% now run mrInit
runMrInit(s);

%%%%%%%%%%%%%%%%%%%%%
%%   examineData   %%
%%%%%%%%%%%%%%%%%%%%%
function [tf s] = examineData(s);

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
  if all(isfield(fileList(i).dicomInfo,{'AcquisitionTime','AcquisitionDate'}))
    fileList(i).startTime = str2num(fileList(i).dicomInfo.AcquisitionTime);
  else
    fileList(i).startTime = inf;
  end
  fileList(i).startHour = floor(fileList(i).startTime/10000);
  fileList(i).startMin = floor(fileList(i).startTime/100)-fileList(i).startHour*100;
end

% sort by start time
fileList = sortFileList(fileList);

% read nifti headers
if s.dispNiftiHeaderInfo
  disppercent(-inf,'(dofmricni) Reading nifti headers');
  for i = 1:length(fileList)
    if ~isempty(fileList(i).nifti)
      system(sprintf('chmod 644 %s',fileList(i).nifti));
      fileList(i).h = mlrImageHeaderLoad(fileList(i).nifti);
    else
      fileList(i).h = [];
    end
    disppercent(i/length(fileList));
  end
  disppercent(inf);
else
  for i = 1:length(fileList)
    fileList(i).h = [];
  end
end

% get list of bold scans
boldNum = 0;
s.seriesDescription = [];
for i = 1:length(fileList)
  % check by whether name contains BOLD or TR is in range
  if ~isempty(findstr('bold',lower(fileList(i).filename))) || (~isempty(fileList(i).te) && (fileList(i).te >= s.teLower) && fileList(i).te <= s.teHigher)
    fileList(i).bold = true;
    boldNum = boldNum+1;
    % name to be copied to
    fileList(i).toName = setext(sprintf('bold%02i_%s',boldNum,fileList(i).filename),fileList(i).niftiExt);
    % also get receiverCoilName and sequence type info
    if isfield(fileList(i).dicomInfo,'SeriesDescription')
      if isempty(s.seriesDescription)
	s.seriesDescription = fileList(i).dicomInfo.SeriesDescription;
      end
    end
    
  else
    fileList(i).bold = false;
    fileList(i).toName = '';
  end
end

% get list of anat scans
anatNum = 0;
for i = 1:length(fileList)
  % check by whether name contains BOLD or TR is in range
  if ~isempty(findstr('t1w',lower(fileList(i).filename))) 
    fileList(i).anat = true;
    anatNum = anatNum+1;
    fileList(i).toName = setext(sprintf('anat%02i_%s',anatNum,fileList(i).filename),fileList(i).niftiExt);
  else
    fileList(i).anat = false;
  end
end

% get the subject id
if isempty(s.subjectID)
  mrParams = {{'subjectID',0,'incdec=[-1 1]','minmax=[0 inf]','Subject ID'}};
  params = mrParamsDialog(mrParams,'Set subject ID');
  if isempty(params),return,end
  s.subjectID = sprintf('s%04i',params.subjectID);
end

% get name of directory to copy
s.localSessionDir = fullfile(s.localDataDir,fileparts(s.cniDir),sprintf('%s%s',s.subjectID,s.studyDate));

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
[sessionParams groupParams] = mrInit([],[],'justGetParams=1','magnet',s.magnet,'operator',s.operatorName,'subject',s.subjectID,'coil',s.receiveCoilName,'pulseSequence',s.seriesDescription);
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
    fileListNum = find(strcmp(scanName,{s.fileList(:).toUncompressedName}));
    if ~isempty(fileListNum)
      % and store it
      v = viewSet(v,'auxParam','dicomInfo',s.fileList(fileListNum).dicomInfo,iScan);
    end
    % set framePeriod as recorded in stimfile
    stimfile = viewGet(v,'stimfile',iScan);
    if ~isempty(stimfile)
      if strcmp(stimfile{1}.filetype,'mgl')
	% get all the volume events
	volEvents = find(stimfile{1}.myscreen.events.tracenum == 1);
	if length(volEvents) > 1
	  % get the framePeriod
	  framePeriod = median(diff(stimfile{1}.myscreen.events.time(volEvents)));
	  % now see if there are more volumes than acquisition triggers
	  if volTrigRatio(iScan)>1
	    framePeriod = framePeriod/volTrigRatio(iScan);
	  end
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

% set up motion comp parameters
%for i = 1:numMotionComp
%  v = newView;
%  if ~isempty(v)
%    [v params] = motionComp(v,[],'justGetParams=1');
%    deleteView(v);
%    eval(sprintf('save motionCompParams%i params',i));
%  end
%end

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

%%%%%%%%%%%%%%%%%%
%%   moveData   %%
%%%%%%%%%%%%%%%%%%
function [tf s] = moveData(justDisplay,s)

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
dirList = {'Etc','Raw','Raw/TSeries','Anatomy'};

% make them
for i = 1:length(dirList)
  if ~isdir(dirList{i})
    command = sprintf('mkdir(''%s'');',fullfile(s.localSessionDir,dirList{i}));
    myeval(command,justDisplay);
  end
end

dispConOrLog(sprintf('=============================================='),justDisplay,true);
dispConOrLog(sprintf('Move Files'),justDisplay,true);
dispConOrLog(sprintf('=============================================='),justDisplay,true);

commandNum = 0;
for i = 1:length(s.fileList)
  % BOLD scan
  if s.fileList(i).bold
    % check for valid nifti
    if isempty(s.fileList(i).nifti)
      dispConOrLog(sprintf('********************************************'),justDisplay,true);
      dispConOrLog(sprintf('(dofmricni) !!!! BOLD nifti file for %s is missing !!!!',s.fileList(i).filename),justDisplay,true);
      % ignore it for later
      s.fileList(i).ignore = true;
    else
      % make full path
      s.fileList(i).toFullfile = fullfile(s.localSessionDir,'Raw/TSeries',s.fileList(i).toName);
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

dispConOrLog(sprintf('=============================================='),justDisplay,true);
dispConOrLog(sprintf('Uncompress and set permission of files'),justDisplay,true);
dispConOrLog(sprintf('=============================================='),justDisplay,true);

commandNum = 0;
for i = 1:length(s.fileList)
  % BOLD scan or anat scans may need to be uncompressed
  if (s.fileList(i).bold || s.fileList(i).anat) && ~s.fileList(i).ignore
    % check if compressed
    if (length(s.fileList(i).niftiExt)>2) && strcmp(s.fileList(i).niftiExt(end-1:end),'gz')
      % make command to gunzip
      command = sprintf('%s -f %s',s.commands.gunzip,s.fileList(i).toFullfile);
      s.fileList(i).toUncompressedName = stripext(s.fileList(i).toName);
      if justDisplay,commandNum=commandNum+1;disp(sprintf('%i: %s',commandNum,command)),else,mysystem(command);,end
      command = sprintf('chmod 644 %s',fullfile(fileparts(s.fileList(i).toFullfile),s.fileList(i).toUncompressedName));
      if justDisplay,commandNum=commandNum+1;disp(sprintf('%i: %s',commandNum,command)),else,mysystem(command);,end
    else
      s.fileList(i).toUncompressedName = s.fileList(i).toName;
    end
  end
end
dispConOrLog(sprintf('=============================================='),justDisplay,true);
dispConOrLog(sprintf('Clean up'),justDisplay,true);
dispConOrLog(sprintf('=============================================='),justDisplay,true);

command = sprintf('rm -rf %s',s.localDir);
if justDisplay,commandNum=commandNum+1;disp(sprintf('%i: %s',commandNum,command)),else,mysystem(command);,end

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
    if s.dispNiftiHeaderInfo && ~isempty(s.fileList(i).h)
      pixdim = s.fileList(i).h.pixdim;
      dim = s.fileList(i).h.dim;
      dispConOrLog(sprintf('%i) %02i:%02i %s [%s] [%s] TR: %0.2f TE: %s -> %s',i,s.fileList(i).startHour,s.fileList(i).startMin,s.fileList(i).filename,mlrnum2str(pixdim,'compact=1'),mlrnum2str(dim,'compact=1','sigfigs=0'),s.fileList(i).tr,mlrnum2str(s.fileList(i).te,'compact=1'),s.fileList(i).toName),justDisplay);
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
  if ~isempty(subjectID) || ~any(length(subjectID) == [4 5]) || isequal(lower(subjectID(1)),'s')
    mrWarnDlg('(dofmricni) PatientName should always be set to a subjectID (not the real name)!!!');
    subjectID = [];
  end
end
if isempty(subjectID)
  if ~isfield(info,'PatientName') && isfield(info.PatientName,'GivenName')
    subjectID = info.PatientName.GivenName
    if ~isempty(subjectID) || ~any(length(subjectID) == [4 5]) || isequal(lower(subjectID(1)),'s')
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
helpFlag = {'-h','-h','-h','-h'};


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

s.sunetID = mglGetParam('sunetID');
if isempty(s.sunetID),s.sunetID = getusername;,end

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
command = sprintf('ssh %s@%s /home/jlg/bin/gruDispData',s.sunetID,s.cniComputerName);
disp(command);
disp('Enter password: ');
[status,result] = system(command);

% check return
if ~isequal(status,0) || isempty(result)
  disp(sprintf('(dofmricni:cniDir) Error getting data directories from cni'));
  return
end

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

%%%%%%%%%%%%%%%%%%%%
%%   getCNIData   %%
%%%%%%%%%%%%%%%%%%%%
function [tf s] = getCNIData(s)

tf = false;

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

% Tell user what is going on
dispHeader;
disp(sprintf('Copying files from %s to %s',s.cniDir,toDir));
disp(sprintf('This could take some time. Using rsync, so that if you quit in the middle'))
disp(sprintf('You can continue where you left off by running dofmricni again'))
dispHeader;

% get dicoms
fromDir = fullfile('/nimsfs/raw/jlg',s.cniDir);
disp(sprintf('(dofmricni) Get files'));
command = sprintf('rsync -rtv --progress --size-only --exclude ''*Screen_Save'' --exclude ''*_pfile*'' --exclude ''*.pyrdb'' --exclude ''*.json'' --exclude ''*.png'' %s@%s:/%s %s',s.sunetID,s.cniComputerName,fromDir,toDir);
disp(command);
system(command);

% got here, so everything is good
tf = true;
s.localDir = fullfile(toDir,getLastDir(s.cniDir));

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
