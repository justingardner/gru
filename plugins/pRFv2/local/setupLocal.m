function setupSherlock()

curPath = pwd;
sherlockSessionPath = ['/share/PI/jlg/' curPath(findstr(curPath, 'data'):end)];
suid = getsuid;

% Check if session directory exists on Sherlock - and make it otherwise.
[~,out] = system(sprintf('ssh %s@sherlock.stanford.edu "[ -d %s ] && echo exists || echo does not exist"', suid, sherlockSessionPath));
if ~strcmp(deblank(out), 'exists')
  disp('Session directory does not exist on Sherlock. Transferring session dir to Sherlock');
  system(sprintf('ssh %s@sherlock.stanford.edu "mkdir -p %s"', suid, sherlockSessionPath));
  system(sprintf('rsync -q %s/Anatomy/* %s@sherlock.stanford.edu:%s/Anatomy/', curPath, suid, sherlockSessionPath));
  system(sprintf('rsync -q %s/Etc/* %s@sherlock.stanford.edu:%s/Etc/', curPath, suid, sherlockSessionPath));
  system(sprintf('rsync -q %s/%s/* %s@sherlock.stanford.edu:%s/%s/', curPath, params.groupName, suid, sherlockSessionPath, params.groupName)); 
  system(sprintf('rsync -q %s/mrSession.mat %s@sherlock.stanford.edu:%s/.', curPath, suid, sherlockSessionPath));
end

% Use rsync to transfer split structs to Sherlock
disp('Copying split structs (.mat) and batch scripts to sherlock server');
system(sprintf('rsync -q Splits/%s*.mat %s@sherlock.stanford.edu:%s/Splits/', params.saveName, suid, sherlockSessionPath));
system(sprintf('rsync -q %s/* %s@sherlock.stanford.edu:%s/%s', scriptsDir, suid, sherlockSessionPath, scriptsDir));



%% OLD CODE
% Call batch submission scripts on Sherlock
%disp('Submitting batch scripts...');
%system(sprintf('ssh %s@sherlock.stanford.edu "cd %s/%s/; sh runAll.sh"', suid, sherlockSessionPath, scriptsDir));
