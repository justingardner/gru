
pathNames = {'~/proj/att_awe','~/proj/unlearning','~/Box Sync/COHCON_DATA'};

addpath('~/proj/gru/mac_helpers'); % explicitly do this so we can use genpath_exclude to remove .git/.svn
for i = 1:length(pathNames)
  if isdir(pathNames{i})
    addpath(genpath_exclude(pathNames{i},{'.git','.svn'}));
  end
end

mglSetParam('sunetID','dbirman');
mrSetPref('maxBlockSize',4000000000);