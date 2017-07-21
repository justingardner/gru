function gruParams = setupGRU(params)

warning('This is an example function: replace GRU with the machine name and save into the appropriate folders');

machine = 'gru';

global controller

% clear the splits folder on sherlock
system(sprintf('ssh gru@%s.stanford.edu "rm -r %s/Splits"',machine,controller.sherlockSessionPath));

% Check if session directory exists on Sherlock - and make it otherwise.
disp(sprintf('Transferring session dir to %s',machine));
system(sprintf('ssh gru@%s.stanford.edu "mkdir -p %s"', machine, controller.curPath));
system(sprintf('rsync -q %s/Anatomy/* gru@%s.stanford.edu:%s/Anatomy/', controller.curPath, machine, controller.curPath));
system(sprintf('rsync -q %s/Etc/* gru@%s.stanford.edu:%s/Etc/', controller.curPath, machine, controller.curPath));
system(sprintf('rsync -q %s/%s/* gru@%s.stanford.edu:%s/%s/', controller.curPath, params.groupName, machine, controller.curPath, params.groupName)); 
system(sprintf('rsync -q %s/mrSession.mat gru@%s.stanford.edu:%s/.', controller.curPath, machine, controller.curPath));

% Use rsync to transfer split structs to Sherlock
disp('Copying split structs (.mat) to sherlock server');
system(sprintf('rsync -q Splits/%s*.mat gru@%s.stanford.edu:%s/Splits/', params.saveName, machine, controller.curPath));

%% Count bins, or default to 10

try
    [~,bins] = system(sprintf('ssh gru@%s.stanford.edu "matlab -nodesktop;cd ~/proj/gru;startup;n = mlrNumWorkers(16);disp(sprintf(''bins=%i'',n));exit',machine));
    gruParams = struct;
    nbins = bins(strfind(bins,'bins=')+4:end);
    if nbins == -1
        nbins = 1;
    end
    gruParams.bins = nbins;
catch
    gruParams.bins = 10;
end