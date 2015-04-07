function mlrGetSurf()
%   mlrGetSurf(project)
%
% GOAL:
%   You should have already run mlrReconAll for your volumes. Now we're
%   going to recover the actual data with the same rP settings. Note that
%   the data will go into the folder on your computer: ~/data/retinotopy/
%   which you should have already initialized with mrInit.
%
% USAGE:
%   success = mlrGetSurf(project,password)
%
% project       Current project, e.g. 'retinotopy'
% password      Your suid password
%
% RETURNS:
%   success     1 when your files were successfully moved onto this local
%               computer. 0 otherwise.

%% Search the LXC server for subjects
s.localDataDir= '~/data/mlrAnatDB';
s.cniComputerName = 'cnic7.stanford.edu';

s = getInfo(s);

command = sprintf('ssh %s@%s ls /data/freesurfer/subjects/',s.sunetID,s.cniComputerName);
disp(command);
disp('Enter password: ');
[status,result] = system(command);

%% Parse result for subject folders
allSubjPos = regexp(result,'s\d\d\d\d');
subjFolders = {};
for i = 1:length(allSubjPos)
    subjFolders{end+1} = result(allSubjPos(i):allSubjPos(i)+4);
end

%% Check what folder we want to move
disp('Which folder(s) do you want to copy locally: e.g. 1, or [1 2]');
for i = 1:length(subjFolders);
    disp(sprintf('%i: %s',i,subjFolders{i}));
end

values = [];
while isempty(values)
    values = input('Folders: ');
    if ~isnumeric(values), values = []; end
end

%% Run mlrGetSurfHelper

for i = 1:length(values)
    mlrGetSurfHelper(subjFolders{values(i)},s);
end

%% Helpers

function success = mlrGetSurfHelper(folder,s)
%% Setup

s.aDBServ = fullfile('https://gru.stanford.edu/hg/mlrAnatDB',folder);
s.aDBLocal = fullfile('~/data/mlrAnatDB',folder);
s.subjectID = folder;

%% 
s.fstempPath = fullfile('/data/freesurfer/subjects/',folder);

%% Check if we already ran this...
if checkForSurfFiles(s.localDataDir)
    warning('It looks like you already ran mlrReconAll for this subject, this function will overwrite the files.');
    keyboard
end

%% Check if we are at Stanford in GRU lab (171.64.40.***)
% ipPlus =  urlread('http://checkip.dyndns.org/');
% if isempty(strfind(ipPlus,'171.64.40'))
%     error('mlrReconAll and mlrGetSurf are only intended for use at Stanford!');
% end

%% Folder Structure
if ~isdir(s.aDBLocal)
    disp(sprintf('Please copy into a terminal, then dbcont:'));
    disp(sprintf('hg clone %s %s',s.aDBServ,s.aDBLocal));
    keyboard
else
    disp('Anatomy Database exists for subject. ');
    cur = pwd;
    cd(s.aDBLocal);
    system('hg pull');
    cd(cur);
end
base = {'3D','mlrBaseAnatomies','niftiROIS','localizers','mlrROIs'};
    
succ = 1;

if ~isdir(s.aDBLocal)            
    succ = 0;
else
    succ = 1;
    for i = 1:length(base)
        if ~isdir(fullfile(s.aDBLocal,base{i}))
            succ = 0;
        end
    end
end
if ~succ
    error(sprintf('Anatomy Database folder structure for %s is not correct.',s.subjectID));
end

%% Connect to LXC Server
try
    curDir = pwd;
    cd('~/proj/gru');
    addpath(genpath('~/proj/gru/ssh2_v2_m1_r6'));
    s.conn = ssh2_config(s.cniComputerName,s.sunetID,input('Password: ','s'));
    cd(curDir);
catch e
    warning('Something went wrong with ssh2, are you sure you have the ssh2 folder on your path?');
    error(e);
end
s.conn = s.conn;

%% Check if FreeSurfer is finished

s.fstempPath = fullfile('/data/freesurfer/subjects/',folder);
s.fstempPathMRI = fullfile(s.fstempPath,'mri');
s.fstempPathSurf = fullfile(s.fstempPath,'surf');

% Check if the T1 is finished
command = sprintf('ls %s',s.fstempPathMRI);
s.conn = ssh2_command(s.conn,command);

pos = findFile(s.conn,'T1.mgz');

exist_t1 = 1;
if ~pos
    exist_t1 = 0;
else
    disp('Found T1.mgz file');
    while ~findFile(s.conn,'T1.nii')
        disp('Converting T1.mgz to T1.nii');
        command = sprintf('ssh %s@%s\nmri_convert %s %s',s.sunetID,s.cniComputerName,fullfile(s.fstempPathMRI,'T1.mgz'),fullfile(s.fstempPathMRI,'T1.nii'));
        disp('Enter the following commands into any terminal [dbcont when complete]:');
        disp(command);
        keyboard
        command = sprintf('ls %s',s.fstempPathMRI);
        s.conn = ssh2_command(s.conn,command);
    end
end
if exist_t1
    disp('T1.nii file exists, will copy as canonical.');
end

%% Check for Surfaces
command = sprintf('ls %s',s.fstempPathSurf);
s.conn = ssh2_command(s.conn,command);


files = {{'lh.pial' 'rh.pial' 'lh.smoothwm' 'rh.smoothwm' 'lh.inflated' 'rh.inflated'}, {'T1.mgz', 'T1.nii'}};

exist_surf = 1;
if ~findFile(s.conn,'lh.pial') || ~findFile(s.conn,'rh.pial') || ~findFile(s.conn,'lh.smoothwm') || ~findFile(s.conn,'rh.smoothwm')
    warning('FreeSurfer didn''t finish running. Re-run this later to copy the surfaces.');
    exist_surf = 0;
end

if exist_surf
    disp('FreeSurfer appears to have completed, will copy surfaces.');
end

%% SCP Files to /data/
i = 1;
if ~isdir(fullfile(s.aDBLocal,'mlrBaseAnatomies','FS'))
    mkdir(fullfile(s.aDBLocal,'mlrBaseAnatomies','FS'));
end
s.fullLocal = fullfile(s.aDBLocal,'mlrBaseAnatomies','FS',num2str(i));
while isdir(s.fullLocal)
    warning('Found existing directory %s.',s.fullLocal);
    in = input('Do you want to overwrite? [y/n]: ','s');
    if strcat(in,'y')
        break
    end
    i = i+1;
    s.fullLocal = fullfile(s.aDBLocal,'mlrBaseAnatomies','FS',num2str(i));
end

scpCommand = sprintf('scp -r %s@%s:%s %s',s.sunetID,s.cniComputerName,s.fstempPath,s.fullLocal);

disp('Copying files locally.');
disp(scpCommand);
system(scpCommand);

%% Check that SCP succeeded
[s_t1 s_surf] = checkForSurfFiles(s.fullLocal);
if ~s_surf
    warning('No surface files.');
end
if ~s_t1
    warning('No canonical files.');
end

%% Run mlrImportFreeSurfer
if s_surf
    cDir = pwd;
    cd (s.fullLocal)
    mlrImportFreeSurfer
    cd(cDir);
    warning('mlrImportFreeSurfer was not run.');
end

%% Copy canonical T1
if s_t1
    file = fullfile(s.fullLocal,'mri','T1.nii');
    nfile = fullfile(s.aDBLocal,'3D','T1w_Canonical.nii');
    copyfile(file,nfile);
end

%% Push to hg server
cur = pwd;
cd(s.aDBLocal);
system(sprintf('hg commit -m "Added segmentation files for %s"',s.subjectID));
system(sprintf('hg push'));
cd(cur);

%% Cleanup the server, get rid of the /data/temp/s#### and /data/freesurfer/subjects/s####
if s_t1 && s_surf
    if strcmp(input('Do you want to remove the folders on the LXC server? y/n: ','s'),'y')
        s.conn = ssh2_command(s.conn,sprintf('rm -rf %s',fullfile('/data/temp/',s.subjectID)));
        s.conn = ssh2_command(s.conn,sprintf('rm -rf %s',fullfile(s.fstempPath)));
        disp('Files removed...');
    end
end

%% Close lxc
ssh2_close(s.conn);
s.conn = struct;

%% 
disp('mlrGetSurf for subject %s is complete.',s.subjectID);

function [s_t1, s_surf] = checkForSurfFiles(dir)
dirs = {'surf','mri'};
files = {{'lh.pial' 'rh.pial' 'lh.smoothwm' 'rh.smoothwm' 'lh.inflated' 'rh.inflated'}, {'T1.mgz', 'T1.nii'}};
success = [1 1];
for di = 1:length(dirs)
    cdir = dirs{di};
    cfiles = files{di};
    for fi = 1:length(cfiles)
        cfile = cfiles{fi};
        % check for file
        if ~isfile(fullfile(dir,cdir,cfile))
            success(di) = 0;
        end
    end
end
s_t1 = success(2);
s_surf = success(1);
return



function s = getInfo(s)

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

function pos = findFile(conn,file)
findFile = cellfun(@strfind,conn.command_result,repmat({file},size(conn.command_result)),'UniformOutput',false);
if isempty(findFile)
    pos = 0;
else
    pos = find(logical(1-cellfun(@isempty,findFile))==1);
    if isempty(pos)
        pos = 0;
    end
end