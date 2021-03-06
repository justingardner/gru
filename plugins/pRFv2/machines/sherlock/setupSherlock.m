function sherlockParams = setupSherlock(params)
global controller

% clear the splits folder on sherlock
system(sprintf('ssh %s@sherlock.stanford.edu "rm -r %s/Splits"',controller.suid,controller.sherlockSessionPath));

% Check if session directory exists on Sherlock - and make it otherwise.
[~,out] = system(sprintf('ssh %s@sherlock.stanford.edu "[ -d %s ] && echo exists || echo does not exist"', controller.suid, controller.sherlockSessionPath));
if ~strcmp(deblank(out), 'exists')
  disp('Session directory does not exist on Sherlock. Transferring session dir to Sherlock');
  system(sprintf('ssh %s@sherlock.stanford.edu "mkdir -m 777 -p %s"', controller.suid, controller.sherlockSessionPath));
  system(sprintf('rsync -q %s/Anatomy/* %s@sherlock.stanford.edu:%s/Anatomy/', controller.curPath, controller.suid, controller.sherlockSessionPath));
  system(sprintf('rsync -q %s/Etc/* %s@sherlock.stanford.edu:%s/Etc/', controller.curPath, controller.suid, controller.sherlockSessionPath));
  system(sprintf('rsync -q %s/%s/* %s@sherlock.stanford.edu:%s/%s/', controller.curPath, params.groupName, controller.suid, controller.sherlockSessionPath, params.groupName)); 
  system(sprintf('rsync -q %s/mrSession.mat %s@sherlock.stanford.edu:%s/.', controller.curPath, controller.suid, controller.sherlockSessionPath));
  system(sprintf('ssh %s@sherlock.stanford.edu "chmod -R 777 %s"', controller.suid, controller.sherlockSessionPath));
else
  % If it does exist, delete the Splits folder
  disp('Session directory found on Sherlock. Cleaning up previous Splits...');
  system(sprintf('ssh %s@sherlock.stanford.edu "rm -r %s/Splits"', controller.suid, controller.sherlockSessionPath));
end

% Use rsync to transfer split structs to Sherlock
% disp('Copying split structs (.mat) to sherlock server');
% system(sprintf('rsync -q Splits/%s*.mat %s@sherlock.stanford.edu:%s/Splits/', params.saveName, controller.suid, controller.sherlockSessionPath));

% make a Scripts folder
system(sprintf('ssh %s@sherlock.stanford.edu "mkdir -m 777 %s/%s"',controller.suid,controller.sherlockSessionPath,controller.scriptsDir));

%% Return 10 bins
sherlockParams = struct;
sherlockParams.bins = 10;

%% OLD CODE
% Call batch submission scripts on Sherlock
%disp('Submitting batch scripts...');
%system(sprintf('ssh %s@sherlock.stanford.edu "cd %s/%s/; sh runAll.sh"', controller.suid, controller.sherlockSessionPath, controller.scriptsDir));
