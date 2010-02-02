function y = applyScalingFactorToEveryHdr

disp('Scaling factor not specified. Will use default values finalized 20jan10.');
disp('sFacs = [1.08468528829887 1.0531954310169 1.05883615483564]');
sFacs = [1.08468528829887 1.0531954310169 1.05883615483564];

fileList = getFileList(pwd,'hdr');

for i = 1:length(fileList)
    baseImg = fileList{i}.filename;
    applyScalingFactor(baseImg,sFacs);
end

y = 0;
disp('Finished scaling every header in this directory. Exiting.');



%%%%%%%%%%%%%%%%%%%%%
%%   getFileList   %%
%%%%%%%%%%%%%%%%%%%%%
function [fileList fileListArray] = getFileList(dirname,extList,matchExtList,filenameMatch)

fileList = {};fileListArray = {};

% make into a cell array
extList = cellArray(extList);

if ~ieNotDefined('matchExtList')
  matchExtList = cellArray(matchExtList);
else
  matchExtList = [];
end

% default to no filename matching 
if ieNotDefined('filenameMatch')
  filenameMatch = '';
else
  filenameMatch = cellArray(filenameMatch);
end

% open the directory
dirList = dir(dirname);
if isempty(dirList),return,end

% now go through the directory looking for matches
for i = 1:length(dirList)
  match = 0;
  % skip all . files
  if dirList(i).name(1) == '.',continue,end
  % skip all files that don't match the filename match if specified
  if ~isempty(filenameMatch)
    noFilenameMatch = 0;
    for j = 1:length(filenameMatch)
      if isempty(strfind(lower(dirList(i).name),lower(filenameMatch{j}))),noFilenameMatch=1;,end
    end
    if noFilenameMatch,continue,end
  end
  % get the file extension
  thisExt = getext(dirList(i).name);
  % check for match
   for j = 1:length(extList)
    if strcmp(extList{j},thisExt)
      % found the match
      if ~match
	% keep the name of the file
	fileList{end+1}.filename = dirList(i).name;
	fileListArray{end+1} = dirList(i).name;
	fileList{end}.fullfile = fullfile(dirname,dirList(i).name);
	fileList{end}.date = dirList(i).date;
	fileList{end}.datenum = dirList(i).datenum;
	% if we need to find matches
	if ~isempty(matchExtList)
	  fileList{end}.filename = dirList(i).name;
	  fileListArray{end+1} = dirList(i).name;
	  fileList{end}.fullfile = fullfile(dirname,dirList(i).name);
	  fileList{end}.path = dirname;
	  fileList{end}.(sprintf('%sfilename',extList{j})) = dirList(i).name;
	  % check for matching file
	  stemName = stripext(fullfile(dirname,dirList(i).name));
	  for k = 1:length(matchExtList)
	    if ~isfile(setext(stemName,matchExtList{j}))
	      disp(sprintf('(dofmrigru1:getFileList) No matching %s file for %s',matchExtList{j},dirList(i).name));
	    else
	      fileList{end}.(sprintf('%sfilename',matchExtList{j})) = setext(dirList(i).name,matchExtList{j});
	    end
	  end
	  
	end
      end
      match = 1;
    end
  end
end
end
end
