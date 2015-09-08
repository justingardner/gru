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
%   mlrGetSurf()

%% Search the LXC server for subjects
s.localDataDir= '~/data/mlrAnatDB';
s.cniComputerName = 'cnic7.stanford.edu';
s.sunetID = mglGetParam('sunetID');

disp(sprintf('\n\n\n\n'));
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

%% Check if FreeSurfer is finished

s.fstempPath = fullfile('/data/freesurfer/subjects/',folder);
s.fstempPathMRI = fullfile(s.fstempPath,'mri');
s.fstempPathSurf = fullfile(s.fstempPath,'surf');

%% Check for Surfaces 
disp(sprintf('\n\n\n\n'));
command = sprintf('ssh %s@%s ls %s',s.sunetID,s.cniComputerName,s.fstempPathSurf);
disp(command);
disp('Enter password: ');
[status,result] = system(command);

if isempty(strfind(result,'lh.pial')) || isempty(strfind(result,'rh.pial')) || isempty(strfind(result,'lh.smoothwm')) || isempty(strfind(result,'rh.smoothwm'))
    warning('FreeSurfer didn''t finish running. Re-run this later to copy the surfaces.');
    return
end

disp('FreeSurfer appears to have completed, copying surfaces.');

%% SCP Files to /temp/
s.tempPath = fullfile('~/data/temp/mlrGetSurf/');
s.fullLocal = fullfile(s.tempPath,s.subjectID);
if ~isdir('~/data/temp')
    mkdir('~/data/temp');
end
if ~isdir(s.tempPath)
    mkdir(s.tempPath);
end

disp(sprintf('\n\n\n\n'));
command = sprintf('scp -r %s@%s:%s %s',s.sunetID,s.cniComputerName,s.fstempPath,s.tempPath);
disp(command);
disp('This could take up to a minute...');
disp('Enter password: ');
[status,result] = system(command);

%% Check that SCP succeeded
[s_t1, s_surf] = checkForSurfFiles(s.fullLocal);
if ~s_surf
    warning('No surface files.');
    keyboard
end
if ~s_t1
    warning('No canonical files.');
    keyboard
end

%% Run mlrImportFreeSurfer
if s_surf
    cDir = pwd;
    cd(s.fullLocal);    
    mrSetPref('niftiFileExtension','.nii');
    mlrImportFreeSurfer('baseName',s.subjectID);
    cd(cDir);
end

%% Run Justin's import code

s.subjNum = mrStr2num(strrep(s.subjectID,'s',''));

mlrAnatDBPut(s.subjNum,s.fullLocal,'freesurfer');

%% Cleanup the server, get rid of the /data/temp/s#### and /data/freesurfer/subjects/s####
disp(sprintf('\n\n\n\n'));
if s_t1 && s_surf
    if strcmp(input('Do you want to remove the folders on the LXC server? y/n: ','s'),'y')

        command1 = sprintf('ssh %s@%s rm -rf %s',s.sunetID,s.cniComputerName,s.fstempPath);
        command2 = sprintf('ssh %s@%s rm -rf %s',s.sunetID,s.cniComputerName,fullfile('/data/temp/',s.subjectID));
        disp(command1);
        disp(command2);
        disp('Enter password: ');
        [status,result] = system(command1);
    end
end

%% 
disp(sprintf('mlrGetSurf for subject %s is complete.',s.subjectID));

function [s_t1, s_surf] = checkForSurfFiles(dir)
dirs = {'surf','mri'};
files = {{'lh.pial' 'rh.pial' 'lh.smoothwm' 'rh.smoothwm' 'lh.inflated' 'rh.inflated'}, {'T1.mgz'}};
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
if pos > 0
    pos = 1;
end