function mlrReconAll()
%   rP = mlrReconAll(proj,exam,subj,user,pass)
%
% GOAL:
%   mlrReconAll is designed to simplify your life by dealing with the
%   entire freesurfer recon process in the background.
%
% USAGE:
%   mlrReconAll('retinotopy','9077','s300','dbirman','PASSW');
%
% project       The current project, should be: /retinotopy/
% examNum       Exam number, you can find this on NIMS (cni.stanford.edu/nims/)
% subjectName   who is this person, e.g. s300 or dan
% username      Your Stanford suid, e.g. dbirman
% password      Your Stanford suid password

%%
disp('If this code hangs, ctrl+c out and try again.');

% Default arguments
cniComputerName = 'cnic7.stanford.edu';
localDataDir = '~/data';

% set up system variable (which gets passed around with important system info)
s.cniComputerName = cniComputerName;
s.localDataDir = mlrReplaceTilde(localDataDir);
s.dispNiftiHeaderInfo = false;
% range for te to be considered a BODL scan
s.teLower = 25;
s.teHigher = 35;

%%
% choose which directory to download from cni
s = getCNIDir(s);
if isempty(s.cniDir),return,end


%% Check if we are at Stanford in GRU lab (171.64.40.***)
% disp('Checking ip address...');
% ipPlus =  urlread('http://checkip.dyndns.org/');
% if isempty(strfind(ipPlus,'171.64.40'))
%     error('mlrReconAll and mlrGetSurf are only intended for use at Stanford, with the NIMS database!');
% end

%% Fill out s

s.cniDirFull = fullfile('/nimsfs','jlg',s.cniDir);

% get the subject id
if ~isfield(s,'subjectID') || isempty(s.subjectID)
  mrParams = {{'subjectID',0,'incdec=[-1 1]','minmax=[0 inf]','Subject ID'}};
  params = mrParamsDialog(mrParams,'Set subject ID');
  if isempty(params),return,end
  s.subjectID = sprintf('s%04i',params.subjectID);
end

s.tempPath = fullfile('/data/temp/',s.subjectID);
s.fstempPath = fullfile('/data/freesurfer/subjects/',s.subjectID);

%% Connect to LXC Server
try
    curDir = pwd;
    cd('~/proj/gru');
    addpath(genpath('~/proj/gru/ssh2_v2_m1_r6'));
    ssh2_conn = ssh2_config(s.cniComputerName,s.sunetID,input('Password: ','s'));
    cd(curDir);
catch e
    warning('Something went wrong with ssh2, are you sure you have the ssh2 folder on your path?');
    error(e);
end
s.conn = ssh2_conn;

%% Check that our folder is correct
ssh2_conn = ssh2_command(ssh2_conn,sprintf('ls %s',s.cniDirFull));
if isempty(ssh2_conn.command_result)
    error('Failed to find your folder, are you sure everything is set up correctly?');
end

rs = ssh2_conn.command_result;
%% Figure out which FILE is ours, first try finding a single T1
corrFile = cellfun(@strfind,ssh2_conn.command_result,repmat({'T1w_9mm'},size(ssh2_conn.command_result)),'UniformOutput',false);
pos = logical(1-cellfun(@isempty,corrFile));

filePos = find(pos==1);


reconStr = '';
for i = 1:length(filePos)
    pos = filePos(i);
    cFile = rs{pos};
    tempWithFile = fullfile(s.cniDirFull,cFile);
    
    if isempty(rs{i})
        error('Failed to find your sequence...');
    end
    % copy files
    ssh2_conn = ssh2_command(ssh2_conn,sprintf('mkdir %s','/data/temp'));
    ssh2_conn = ssh2_command(ssh2_conn,sprintf('mkdir %s',s.tempPath));
    curFilePath = fullfile(s.tempPath,strcat(s.subjectID,'_',num2str(i),'_','c.nii.gz'));
    ssh2_conn = ssh2_command(ssh2_conn,sprintf('cp %s %s',fullfile(tempWithFile,'*.nii.gz'),curFilePath));
    % gunzip
    pause(.1)
    ssh2_conn = ssh2_command(ssh2_conn,sprintf('gunzip -d %s',fullfile(curFilePath)));
    curFilePath = fullfile(s.tempPath,strcat(s.subjectID,'_',num2str(i),'_','c.nii'));
    
    pause(.1)
    % Build reconStr
    reconStr = sprintf('%s -i %s',reconStr,curFilePath);
end

%% Recon-All
reconCommand = sprintf('recon-all -subject %s %s -all',s.subjectID,reconStr);

% I don't have this working, so for now you have to do the next step by
% hand.
disp(sprintf('\n\nAll of your files are now organized. Open a new terminal window.\nLeave the terminal window open until these commands complete:\n\nssh -XY %s@cnic7.stanford.edu\n%s\n\nThat''s it!',s.sunetID,reconCommand));

%% Close the LXC Server
ssh2_conn = ssh2_close(ssh2_conn);

clear ssh2_conn

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
scanName = params.scanName{1};

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
