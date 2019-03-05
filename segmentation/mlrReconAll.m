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

%% Get the subject ID
s.sid = mglGetSID;
if isempty(s.sid)
    disp('*****');
    disp('(mlrReconAll) Please set the subject ID by calling mglSetSID(#)');
    disp('*****');
    return
end

disp('*******************************************************');
disp(sprintf('*** mlrReconAll RUNNING FOR SUBJECT %s ***',s.sid));
disp('*******************************************************');

% get the numerical version and the s0### version
s.nsid = gruSubjectID2num(s.sid);
s.sid0 = gruSubjectNum2ID(s.nsid);

%% Search for anatomy files using Flywheel CLI

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

datas = {};

% find studies that match the right subject
for si = 1:length(studyNames)
  studyName = studyNames{si};
  disp(studyName);
  % remove annoying escape codes
  studyName = removeEscapeCodes(studyName);
  % get list
  [status fwListing] = system(sprintf('%s ls %s/%s',s.commands.fw,s.PI,studyName));
  % parse the results
  expNamesParse = textscan(fwListing,'%s %s %s %s %s %s');
  if length(expNamesParse) == 6
    for iExp = 1:length(expNamesParse{6})
      expNames{si}{iExp} = removeEscapeCodes(expNamesParse{6}{iExp});
      expFullNames{si}{iExp} = sprintf('%s_%s_%s_%s',removeEscapeCodes(expNamesParse{2}{iExp}),removeEscapeCodes(expNamesParse{3}{iExp}),removeEscapeCodes(expNamesParse{5}{iExp}),removeEscapeCodes(expNamesParse{6}{iExp}));
    end
  end
  
  % keep studies + sessions which match the subject id (with or without the
  % zero)
  for ei = 1:length(expNamesParse{1})
    subj = expNamesParse{5}{ei};
    subj = removeEscapeCodes(subj);
    
    if strcmp(s.sid,subj) || strcmp(s.sid0,subj) || strcmp(s.nsid,subj)
        sess = expNamesParse{6}{ei};
        sess = removeEscapeCodes(sess);
        % find the anatomy files
        [status fwListing] = system(sprintf('%s ls %s/%s/%s',s.commands.fw,s.PI,studyName,sess));
        % because the file names have spaces you can't use textscan, so
        % parse by the position of the word 'admin' instead
        aIdx = strfind(fwListing,'admin');
        aIdx(end+1) = length(fwListing);
        for ai = 1:(length(aIdx)-1)
            str = fwListing(aIdx(ai):(aIdx(ai+1)));
            % chop everything except the file name
            scanType = removeEscapeCodes(str(27:end));
            scanType_ = lower(scanType);
            % check if the file type is something about anatomy
            anatomy = false;
            if ~isempty(strfind(scanType_,'anat')), anatomy = true; end
            if ~isempty(strfind(scanType_,'t1w .9')), anatomy = true; end
            
            if anatomy
                disp('**This file is an anatomy')
                datas{end+1,1} = s.PI;
                datas{end,2} = studyName;
                datas{end,3} = sess;
                datas{end,4} = scanType;
                disp(sprintf('Downloading: %s', datas{end}));
            end
        end
    end
  end
end

%% Create the temp download directory
s.tempDir = fullfile('~/data/temp');
if ~isdir(s.tempDir), mkdir(s.tempDir); end
s.tempDir = fullfile(s.tempDir,'mra');
if ~isdir(s.tempDir), mkdir(s.tempDir); end

%% Confirm the files that we will use, drop files if needed

if isempty(datas)
    disp('(mlrReconAll) No anatomy files were found');
    return
end

picked = false;
while ~picked
    disp('*******************************************************');
    disp('******* mlrReconAll Using the following files: ********');
    disp('*******************************************************');
    for di = 1:size(datas,1)
        disp(sprintf('File %i: %s/%s/%s/%s',di,datas{di,1},datas{di,2},datas{di,3},datas{di,4}));
    end
    res = input('(mlrReconAll) Enter to continue, or choose files to keep (e.g. [1 3]): ');
    if isempty(res)
        picked = true;
    else
        datas = datas(res);
    end
end

if isempty(datas)
    disp('(mlrReconAll) No anatomy files were found');
    return
end


disp('*******************************************************');

%% Download the anatomy files
cd(s.tempDir);
files = {};
for di = 1:size(datas,1)
    command = sprintf('%s download -f "%s"',s.commands.fw,fullfile(datas{di,1},datas{di,2},datas{di,3},datas{di,4}));
    disp(command);
    system(command);
    
    files{end+1} = fullfile(s.tempDir,sprintf('%s.%s',strrep(lower(datas{di,4}),' ','-'),'tar'));
    
    % check to see if the tar file downloaded properly
    if ~isfile(files{end})
      disp(sprintf('(dofmricni:getFlywheelData) Tarfile %s did not get downloaded properly from fw',files{end}));
      return
    end
end

%% Check if we are at Stanford in GRU lab (171.64.40.***)
% disp('Checking ip address...'); 
%ipPlus = urlread('http://checkip.dyndns.org/'); if
% isempty(strfind(ipPlus,'171.64.40'))
%     error('mlrReconAll and mlrGetSurf are only intended for use at
%     Stanford, with the NIMS database!');
% end

%% Fill out s
s.tempPath = fullfile('/data/temp/',s.sid0);
s.fstempPath = fullfile('/data/freesurfer/subjects/',s.sid0);

%% Untar and copy files to CNI server
% gzFiles = {};
for fi = 1:length(files)
    [status, result] = system(sprintf('%s xfv %s','tar',files{fi}));
    if ~isequal(status,0)
      disp(sprintf('(dofmricni:getFlywheelData) Could not untar: %s',tarfileName));
      return
    end

    % this would work, but if you're stupid and put spaces in your scan
    % sequences then it can't handle it:
%     niiFiles = textscan(result,'%s %s');
%     niiFiles = niiFiles{2};
%     for ni = 1:length(niiFiles)
%         file = niiFiles{ni};
%         if strfind(file,'nii')
%             gzFiles{end+1} = file;
%         end
%     end
end

% the files are now in the folder
% ~/data/temp/mra/jlg/studyName/sid/session#/fileType/*.gz

% because wtf flywheel seriously we asked for a single file not a folder
% structure... also why the f does it add the subject into the folder
% structure? That wasn't even there originally... 

% warning: the sid could be anything, since we're not always consistent
% about using the 0 or not, so instead of pulling the file directly let's
% do a depth search of the folder structure and just grab anything that is
% .nii.gz, gunzip it, and pop it into the top level. Then we can check if 
% we are missing any files.

gzFiles = {};
search = {s.tempDir};

T1 = zeros(1,size(datas,1));
T2 = zeros(1,size(datas,1));

while ~isempty(search)
    curPath = search{1};
    search = search(2:end);
    
    files = dir(curPath);
    for fi = 1:length(files)
        if ~(strcmp(files(fi).name,'.') || strcmp(files(fi).name,'..'))
            filePath = fullfile(curPath,files(fi).name);
            disp(filePath);
            if isfile(filePath)
                if strfind(filePath,'nii')
                    % move to top level folder
                    newPath = fullfile(s.tempDir,files(fi).name);
                    movefile(filePath,newPath);
                    gzFiles{end+1} = newPath;
                    if ~isempty(strfind(lower(curPath),'t1'))
                        T1(length(gzFiles))= 1;
                    elseif ~isempty(strfind(lower(curPath),'t2'))
                        T2(length(gzFiles)) = 1;
                    end
                    disp('Copying file');
                end
            elseif isdir(filePath)
                search{end+1} = filePath;
                disp('Directory');
            end
        end
    end
end

% gunzip files
for gi = 1:length(gzFiles)
    system(sprintf('gunzip %s',gzFiles{gi}));
    gzFiles{gi} = gzFiles{gi}(1:end-3);
end

% check that this all worked
if size(datas,1)~=length(gzFiles)
    disp('(mlrReconAll) A file got lost during data copying.');
    return
end

for gi = 1:length(gzFiles)
    if ~isfile(gzFiles{gi})
        disp(sprintf('(mlrReconAll) During gunzip a file disappeared: %s',gzFiles{gi}));
    end
end

if any(T1==T2)
    disp('(mlrReconAll) A file was marked as both a T1 and a T2, that''s not good');
end

%% report results back to user

disp('*******************************************************');
disp('********** mlrReconAll: DOWNLOAD SUCCESSFUL ***********');
disp('*******************************************************');

for gi = 1:length(gzFiles)
    if T1(gi)
        disp(sprintf('File %i is a T1 anatomy: %s',gzFiles{gi}));
    else
        disp(sprintf('File %i is a T2 anatomy: %s',gzFiles{gi}));
    end
end


%% Make temp directories

command = sprintf('mkdir %s','/data/temp');
result = doRemoteCommand(s.sunetID,s.cniComputerName,command);
command = sprintf('chmod 777 %s','/data/temp');
result = doRemoteCommand(s.sunetID,s.cniComputerName,command);
command = sprintf('mkdir %s',s.tempPath);
result = doRemoteCommand(s.sunetID,s.cniComputerName,command);

%% copy files to CNIC server
cniFiles = cell(size(gzFiles));
fullcniFiles = cell(size(gzFiles));
disp('(mlrReconAll) Copying files to CNI server');
for gi = 1:length(gzFiles)
    if T1(gi)
        cniFiles{gi} = sprintf('anat%i_T1.nii',gi);
    else
        cniFiles{gi} = sprintf('anat%i_T2.nii',gi);
    end
    fullcniFiles{gi} = fullfile(fullfile(s.tempPath,cniFiles{gi}));
    command = sprintf('scp %s %s@%s:%s',gzFiles{gi},s.sunetID,s.cniComputerName,fullcniFiles{gi});
    disp(command);
    system(command);
end

%% delete local files
disp('(mlrReconAll) removing local files');
system(sprintf('rm -rf %s',s.tempDir));

%% test zone

try
    command = sprintf('ssh %s@%s ls %s',s.sunetID,s.cniComputerName,s.tempPath);
    disp(command);
    disp('Testing SSH. Enter password (if needed): ');
    [status,result] = system(command);
catch e
    warning('Something went wrong with ssh.');
    error(e);
end

%% Check that our folder is correct
if isempty(result)
    error('Failed to find your folder, are you sure everything is set up correctly?');
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

%% Build reconStr
reconStr = '';
for ci = 1:length(fullcniFiles)
    if T1(ci)
        reconStr = [reconStr,sprintf('-i %s ',fullcniFiles{ci})];
    end
end
for ci = 1:length(fullcniFiles)
    if T2(ci)
        reconStr = [reconStr,sprintf('-T2 %s',fullcniFiles{ci})];
    end
end
if any(T2)
    reconStr = strcat(reconStr,' -T2pial');
else
    reconStr = reconStr(1:end-1);
end

%% Recon-All
reconCommand = sprintf('recon-all -subject %s %s -all -sd /data/freesurfer/subjects',s.sid0,reconStr);

disp('*******************************************************');
disp('***** mlrReconAll: COMPLETE -- SEE INSTRUCTIONS: ******');
disp('*******************************************************');

% You have to do this into a new terminal, it won't run through MATLAB
% directly.
disp(sprintf('\n\nAll of your files are now organized. Open a new terminal window.\nLeave the terminal window open until these commands complete (6-10 hours):\n\nssh -XY %s@cnic7.stanford.edu\n%s\n\nThat''s it!',s.sunetID,reconCommand));

%%%%%%%%%%%%%%%%%%%%%%%%
%    getFlywheelDir    %
%%%%%%%%%%%%%%%%%%%%%%%%
function s = getFlywheelDir(s)



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

if isempty(s)
    return
end

% strip off weird pre-code escape characters
if s(1) == 27
  % remove pre escape code
  s = s(8:end);
end

if ~isempty(find(s==27))
  s = s(1:first(find(s==27)-1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%    doRemoteCommand    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function [retval,status] = doRemoteCommand(username,computerName,commandName)

retval = [];
command = sprintf('ssh %s@%s %s',username,computerName,commandName);
disp(sprintf('(mlrReconAll) Doing remote command: %s',command));
[status,retval] = system(command,'-echo');
if status~=0
  disp(sprintf('(mlrReconAll) Could not ssh in to do remote command on: %s@%s',username,computerName));
  disp(sprintf('If you have not yet set passwordless ssh (see: http://gru.stanford.edu/doku.php/gruprivate/sshpassless) then enter your password here: ',computerName));
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
