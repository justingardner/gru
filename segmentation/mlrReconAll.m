function mlrReconAll()
%   rP = mlrReconAll
%
% GOAL:
%   mlrReconAll is designed to simplify your life by dealing with the
%   entire freesurfer recon process in the background.
%
% USAGE:
%   mlrReconAll See here for usage notes:
%   http://gru.stanford.edu/doku.php/runningfmri/tutorial_mlrreconall
%
% project       The current project, should be: /retinotopy/ examNum
% Exam number, you can find this on NIMS (cni.stanford.edu/nims/)
% subjectName   who is this person, e.g. s300 or dan username      Your
% Stanford suid, e.g. dbirman password      Your Stanford suid password
%
% MANUAL
%   The command we use is: recon-all -subject SID -i FILE1 -i FILE2 ... -all -sd /data/freesurfer/subjects
%
%   If you copy a file in yourname@cnic7.stanford.edu:/data/... and then
%   specify that you can run by hand

%%
disp('If this code hangs, ctrl+c out and try again.');

% Default arguments
cniComputerName = 'cnic7.stanford.edu';
localDataDir = '~/data';

% set up system variable (which gets passed around with important system
% info)
s.sunetID = mglGetParam('sunetID');
s.cniComputerName = cniComputerName;
s.localDataDir = mlrReplaceTilde(localDataDir);
s.dispNiftiHeaderInfo = false;
s.PI = 'jlg';
% range for te to be considered a BODL scan
s.teLower = 25;
s.teHigher = 35;

%%
% choose which directory to download from cni
%s = getCNIDir(s);
s = getFlywheelDir(s);
if isempty(s.cniDir),return,end

%% Check if we are at Stanford in GRU lab (171.64.40.***)
% disp('Checking ip address...'); ipPlus =
% urlread('http://checkip.dyndns.org/'); if
% isempty(strfind(ipPlus,'171.64.40'))
%     error('mlrReconAll and mlrGetSurf are only intended for use at
%     Stanford, with the NIMS database!');
% end

%% Fill out s
s.cniDirFull = fullfile('jlg',s.cniDir);

% get the subject id
if ~isfield(s,'subjectID') || isempty(s.subjectID)
  mrParams = {{'subjectID',0,'incdec=[-1 1]','minmax=[0 inf]','Subject ID'}};
  params = mrParamsDialog(mrParams,'Set subject ID (Numbers Only, s0021 = 21)');
  if isempty(params),return,end
  s.subjectID = sprintf('s%04i',params.subjectID);
end

s.tempPath = fullfile('/data/temp/',s.subjectID);
s.fstempPath = fullfile('/data/freesurfer/subjects/',s.subjectID);



%% test zone

try
    command = sprintf('ssh %s@%s ls %s',s.sunetID,s.cniComputerName,s.cniDirFull);
    command = sprintf('%s ls %s', s.commands.fw, s.cniDirFull);
    disp(command);
    disp('Enter password: ');
    [status,result] = system(command);
catch e
    warning('Something went wrong with ssh.');
    error(e);
end

%% Check that our folder is correct
if isempty(result)
    error('Failed to find your folder, are you sure everything is set up correctly?');
end

%% Replace result with a cell array
newLines = [1 strfind(result,char(10))];

resultc = {};
for i = 1:length(newLines)-1
    resultc{end+1} = strrep(result(newLines(i):newLines(i+1)),char(10),'');
end

%% Make temp directories

command = sprintf('mkdir %s','/data/temp');
result = doRemoteCommand(s.sunetID,s.cniComputerName,command);
command = sprintf('chmod 777 %s','/data/temp');
result = doRemoteCommand(s.sunetID,s.cniComputerName,command);
command = sprintf('mkdir %s',s.tempPath);
result = doRemoteCommand(s.sunetID,s.cniComputerName,command);

%% Get the local directory
%s.localSessionDir = fullfile(s.localDataDir,s.experimentName,sprintf('%s%s',s.subjectID,s.studyDate));
mrParams = {{'localSessionDir',s.localDataDir,'Where is your data locally stored?'}};
params = mrParamsDialog(mrParams,'Enter the full path to your local session directory');
if isempty(params),return,end
s.localSessionDir = params.localSessionDir;

a = dir([s.localSessionDir '/Anatomy/*T1w*9mm*.nii']);
anat_filename = a(1).name;

%% Figure out which FILE is ours, first try finding a single T1
corrFile = cellfun(@strfind,resultc,repmat({'T1w .9mm'},size(resultc)),'UniformOutput',false);
pos = logical(1-cellfun(@isempty,corrFile));
filePos = find(pos==1);

reconStr = '';
for i = 1:length(filePos)
    pos = filePos(i);
    cFile = resultc{pos};
    if isempty(cFile)
        error('Failed to find your sequence...');
    else
      cFile = cFile(strfind(cFile, 'T1w') : end);
    end
    tempWithFile = fullfile(s.cniDirFull,cFile);
    
    % copy files
    curFilePath = fullfile(s.tempPath,strcat(s.subjectID,'_',num2str(i),'_','c.nii'));
    fileStem = getLastDir(tempWithFile);
    endLoc = findstr(fileStem,'_T1w');
    fileStem = fileStem(1:endLoc-1);
    curFileName = getLastDir(s.tempPath);
    %command = sprintf('cp %s %s',fullfile(tempWithFile,sprintf('%s.nii.gz',fileStem)),curFilePath);
    command = sprintf('rsync -av %s/Anatomy/%s "%s@%s:%s"', s.localSessionDir, anat_filename, s.sunetID, s.cniComputerName, curFilePath);
    result = system(command);
    % in case it doesn't work with no quotation
    %we add double quotation in the command: seems to work
    %commandCheckExist = ['ls ' strcat([s.tempPath '/' s.subjectID,'_',num2str(i),'_','c.nii'])];
    %[~,status] = doRemoteCommand(s.sunetID,s.cniComputerName,commandCheckExist);        
    %if status~=0
    %    command = ['"' command '"'];
    %end
   
    %result = doRemoteCommand(s.sunetID,s.cniComputerName,command);
    % gunzip
    pause(.1)
    %command = sprintf('gunzip -d %s',fullfile(curFilePath));
    %result = doRemoteCommand(s.sunetID,s.cniComputerName,command);
    curFilePath = fullfile(s.tempPath,strcat(s.subjectID,'_',num2str(i),'_','c.nii'));
    
    pause(.1)
    % Build reconStr
    reconStr = sprintf('%s -i %s',reconStr,curFilePath);
end

%% Make freesurfer folder
command = sprintf('mkdir %s','/data/freesurfer');
result = doRemoteCommand(s.sunetID,s.cniComputerName,command);
command = sprintf('chmod 777 %s','/data/freesurfer');
result = doRemoteCommand(s.sunetID,s.cniComputerName,command);


command = sprintf('mkdir %s','/data/freesurfer/subjects');
result = doRemoteCommand(s.sunetID,s.cniComputerName,command);
command = sprintf('chmod 777 %s','/data/freesurfer/subjects');
result = doRemoteCommand(s.sunetID,s.cniComputerName,command);
%% Recon-All
reconCommand = sprintf('recon-all -subject %s %s -all -sd /data/freesurfer/subjects',s.subjectID,reconStr);

% You have to do this into a new terminal, it won't run through MATLAB
% directly.
disp(sprintf('\n\nAll of your files are now organized. Open a new terminal window.\nLeave the terminal window open until these commands complete (6-10 hours):\n\nssh -XY %s@cnic7.stanford.edu\n%s\n\nThat''s it!',s.sunetID,reconCommand));


%%%%%%%%%%%%%%%%%%%%%%%%
%    getFlywheelDir    %
%%%%%%%%%%%%%%%%%%%%%%%%
function s = getFlywheelDir(s)

% tell user what is going on
dispHeader;
dispHeader('Querying Flywheel using fw CLI');
dispHeader;
s.commands.fw = 'fw';
% user CLI to get list of experiments
[status fwListing] = system(sprintf('%s ls %s',s.commands.fw,s.PI));
if ~isequal(status,0)
  disp(sprintf('(mlrReconAll:getFlywheelDir) fw ls command returned status: %i',status));
  disp(sprintf('(mlrReconAll:getFlywheelDir) You may need to login to fw'));
  return
end

% parse the results
studyNames = textscan(fwListing,'%s %s');
if length(studyNames) == 2,studyNames = studyNames{2};else,studyNames = {};end

% check that we found something
if isempty(studyNames)
  disp(sprintf('(mlrReconAll) Could not find any studies'));
  return
end

% now go through and get each of the scans for each studies
for iStudy = 1:length(studyNames)
  % remove annoying escape codes
  studyNames{iStudy} = removeEscapeCodes(studyNames{iStudy});
  % get list
  [status fwListing] = system(sprintf('%s ls %s/%s',s.commands.fw,s.PI,studyNames{iStudy}));
  % parse the results
  expNamesParse = textscan(fwListing,'%s %s %s %s %s %s');
  if length(expNamesParse) == 6
    for iExp = 1:length(expNamesParse{6})
      expNames{iStudy}{iExp} = removeEscapeCodes(expNamesParse{6}{iExp});
      expFullNames{iStudy}{iExp} = sprintf('%s_%s_%s_%s',removeEscapeCodes(expNamesParse{2}{iExp}),removeEscapeCodes(expNamesParse{3}{iExp}),removeEscapeCodes(expNamesParse{5}{iExp}),removeEscapeCodes(expNamesParse{6}{iExp}));
    end
  end
end

% Now set up variables to have the default list be from all studies
mrParams = {{'chooseNum',1,'minmax',[1 length(studyNames)],'incdec=[-1 1]'},...
	    {'studyName',studyNames,'type=string','Name of studies','group=chooseNum','editable=0'},...
	    {'expName',expFullNames,'Name of scan','group=chooseNum'}};
params = mrParamsDialog(mrParams);
if isempty(params),return,end
studyName = params.studyName{params.chooseNum};
expName = params.expName{params.chooseNum};

% get subjectID
underscores = strfind(expName,'_');
s.subjectID = gruSubjectNum2ID(expName(underscores(2)+1:underscores(3)-1));

% convert back expname to not have data and subjectID
expName = expName(last(strfind(expName,'_'))+1:end);

% set cniDir
s.cniDir = fullfile(studyName,expName);
disp(sprintf('(mlrReconAll:getFlywheelDir) Directory chosen is: %s',s.cniDir))

% set the directory to which we resync data
toDir = mlrReplaceTilde(fullfile(s.localDataDir,'temp/dofmricni'));
if ~isdir(toDir)
  try
    mkdir(toDir);
  catch
    mrWarnDlg(sprintf('(mlrReconAll) Cannot make directory %s. Either you do not have permissions, or perhaps you have a line to a drive that is not currently mounted?',toDir));
    return
  end
end
s.localDir = fullfile(toDir,getLastDir(s.cniDir));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    removeEscapeCodes    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = removeEscapeCodes(s)

% strip off weird pre-code escape characters
if s(1) == 27
  % remove pre escape code
  s = s(8:end);
end

if ~isempty(find(s==27))
  s = s(1:first(find(s==27)-1));
end

%%%%%%%%%%%%%%%%%%%
%%   getCNIDir   %%
%%%%%%%%%%%%%%%%%%%
function s = getCNIDir(s)

s.cniDir = [];

% get the username
if isempty(s.sunetID)
  s.sunetID = getusername;
end

% default sunetID to be username
s.sunetID = mglGetParam('sunetID');
if isempty(s.sunetID),s.sunetID = s.sunetID;end

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
result = doRemoteCommand(s.sunetID,s.cniComputerName,sprintf('/home/%s/bin/gruDispData -pi %s',s.sunetID,s.PI));
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
  disp(sprintf('(mlrReconAll) Could not find any studies'));
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
disp(sprintf('(mlrReconAll:getCNIDir) Directory chosen is: %s',s.cniDir))

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


%%%%%%%%%%%%%%%%%%%%%%%%%
%    doRemoteCommand    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function [retval,status] = doRemoteCommand(username,computerName,commandName)

retval = [];
command = sprintf('ssh %s@%s %s',username,computerName,commandName);
disp(sprintf('(mlrReconAll) Doing remote command: %s',command));
disp(sprintf('If you have not yet set passwordless ssh (see: http://gru.stanford.edu/doku.php/gruprivate/sshpassless) then enter your password here: ',computerName));
[status,retval] = system(command,'-echo');
if status~=0
  disp(sprintf('(mlrReconAll) Could not ssh in to do remote command on: %s@%s',username,computerName));
  return
end
disp(sprintf('(mlrReconAll) Remote command on %s successful',computerName));


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
% sometimes there is more than one line (errors from csh startup) so need
% to strip those off
username = strread(username,'%s','delimiter','\n');
username = username{end};
if (retval == 0)
  % get it again
  [retval username2] = system('whoami');
  username2 = strread(username2,'%s','delimiter','\n');
  username2 = username2{end};
  if (retval == 0)
    % find the matching last characers this is necessary, because matlab's
    % system command picks up stray key strokes being written into the
    % terminal but puts those at the beginning of what is returned by
    % stysem. so we run the command twice and find the matching end part of
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
