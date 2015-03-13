function reconParams = mlrReconAll(project,examNum,subjectName,username,password)
%   rP = mlrReconAll(proj,exam,subj,user,pass)
%
% GOAL:
%   mlrReconAll is designed to simplify your life by dealing with the
%   entire freesurfer recon process in the background.
%
% USAGE:
%   rP = mlrReconAll('retinotopy','9077','s300','dbirman','PASSW');
%   succ = mlrGetSurf(rP);
%
% project       The current project, should be: /retinotopy/
% examNum       Exam number, you can find this on NIMS (cni.stanford.edu/nims/)
% subjectName   who is this person, e.g. s300 or dan
% username      Your Stanford suid, dbirman
% password      Your Stanford suid password

%% Check if we are at Stanford in GRU lab (171.64.40.***)
ipPlus =  urlread('http://checkip.dyndns.org/');
if isempty(strfind(ipPlus,'171.64.40'))
    error('mlrReconAll and mlrGetSurf are only intended for use at Stanford, with the NIMS database!');
end

%% Fill out reconParams

reconParams = struct;

reconParams.LXCServer = 'cnic7.stanford.edu';
reconParams.nimsDataPath = fullfile('/nimsfs','jlg',project);
reconParams.tempPath = fullfile('/data/freesurfer/subjects/',subjectName);
reconParams.localDataPath = fullfile('/data',project,'Anatomy',subjectName);
reconParams.subj = subjectName;
reconParams.fs_subj = strcat('fs_',subjectName);
reconParams.fstempPath = fullfile('/data/freesurfer/subjects/',reconParams.fs_subj);
reconParams.user = username;
reconParams.pass = password;

%% Connect to LXC Server
try
    ssh2_conn = ssh2_config(reconParams.LXCServer,reconParams.user,reconParams.pass);
catch e
    warning('Something went wrong with ssh2, are you sure you have the ssh2 folder on your path?');
    error(e);
end
reconParams.conn = ssh2_conn;

%% Figure out the full folder name of our exam
ssh2_conn = ssh2_command(ssh2_conn,sprintf('ls %s',reconParams.nimsDataPath));
corrFolder = cellfun(@strfind,ssh2_conn.command_result,repmat({num2str(examNum)},size(ssh2_conn.command_result)),'UniformOutput',false);
pos = logical(1-cellfun(@isempty,corrFolder));
tempWithFolder = fullfile(reconParams.nimsDataPath,ssh2_conn.command_result{pos});

if isempty(ssh2_conn.command_result{pos})
    error('Failed to find your folder, are you sure everything is set up correctly?');
end

%% Figure out which FILE is ours, first try finding a single T1
ssh2_conn = ssh2_command(ssh2_conn,sprintf('ls %s',tempWithFolder));
corrFile = cellfun(@strfind,ssh2_conn.command_result,repmat({'T1w_9mm'},size(ssh2_conn.command_result)),'UniformOutput',false);
pos = logical(1-cellfun(@isempty,corrFile));
if sum(pos) > 1 || sum(pos) == 0
    warning('Can''t handle multiple T1 files in the same folder... We will just use the first');
    allOnes = find(pos==1);
    pos(allOnes(1)+1:end) = 0;
end
tempWithFile = fullfile(tempWithFolder,ssh2_conn.command_result{pos});

if isempty(ssh2_conn.command_result{pos})
    error('Failed to find your sequence...');
end

%% There should be just one file here, the .nii.gz, copy this to /data/
ssh2_conn = ssh2_command(ssh2_conn,sprintf('ls %s',reconParams.tempPath));
if isempty(ssh2_conn.command_result{1})
    ssh2_conn = ssh2_command(ssh2_conn,sprintf('mkdir %s',reconParams.tempPath));
    ssh2_conn = ssh2_command(ssh2_conn,sprintf('cp %s %s',fullfile(tempWithFile,'*.nii.gz'),reconParams.tempPath));
    %% Gunzip
    ssh2_conn = ssh2_command(ssh2_conn,sprintf('gunzip -d %s',fullfile(reconParams.tempPath,'*')));
end

%% And get the actual file name
ssh2_conn = ssh2_command(ssh2_conn,sprintf('ls %s',reconParams.tempPath));
reconParams.tempFilePath = fullfile(reconParams.tempPath,ssh2_conn.command_result{1});

%% Recon-All
reconCommand = sprintf('recon-all -subject %s -i %s -all',reconParams.fs_subj,reconParams.tempFilePath);

% I don't have this working, so for now you have to do the next step by
% hand.
disp(sprintf('\n\nAll of your files are now organized.\nPlease execute the following commands in the terminal:\n\nssh -XY %s@cnic7.stanford.edu\n%s\n\nThat''s it! When FreeSurfer finishes call\n\nmlrGetSurf(reconParams);',username,reconCommand));

%% Close the LXC Server
ssh2_conn = ssh2_close(ssh2_conn);