function success = mlrGetSurf(reconParams)

% You should have already run rP = mlrReconAll for your volumes. Now we're
% going to recover the actual data with the same rP settings. The data will
% get saved into the same directory structure that it was pulled out of
% NIMS on, in this subject's folder.

% TODO:
%   Add code that connects to server
%   Check whether freesurfer is done:
%       You need *h.pial, *h.smoothwm, *h.inflated, and mri/T1.mgz for both
%       left and right hemispheres
%   If freesurfer is done, scp the files locally
%   Check that all relevant files exist
%   Check if mrInit has been run
%   If yes, run mlrImportFreeSurfer
%   Exit success

%% Check if we are at Stanford in GRU lab (171.64.40.***)
ipPlus =  urlread('http://checkip.dyndns.org/');
if isempty(strfind(ipPlus,'171.64.40'))
    error('mlrReconAll and mlrGetSurf are only intended for use at Stanford, with the NIMS database!');
end

%% Check that reconParams has everything necessary
if ~isfield(reconParams,'conn') || ~isfield(reconParams,'LXCServer') || ~isfield(reconParams,'fstempPath') || ~isfield(reconParams,'localDataPath') || ~isfield(reconParams,'user') || ~isfield(reconParams,'pass')
    error('reconParams was not built correctly... cannot continue');
end

%% Connect to the LXC Server
ssh2_conn = reconParams.conn;

%% Check if FreeSurfer is finished
ssh2_conn = ssh2_command(ssh2_conn,sprintf('ls %s',fullfile(reconParams.fstempPath,'surf')));
if isempty(ssh2_conn.command_result{1})
    warning('FreeSurfer didn''t finish running. Be patient!');
    success = 0;
    return
end
disp('FreeSurfer appears to have completed, checking for specific files.');

files = ssh2_con.command_result;
leftFiles = [0 0 0];
rightFiles = [0 0 0];
T1 = [0 0 0];
for i = 1:length(files)
    cfile = files{i};
    
end

%% SCP Files to /data/

%% Check that SCP succeeded

%% Check if mrInit was run

%% Run mlrImportFreeSurfer

%% Success

success = 1;