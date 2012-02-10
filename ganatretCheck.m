% checkFlats.m
%
%        $Id:$ 
%      usage: checkFlats()
%         by: justin gardner
%       date: 02/10/12
%    purpose: Checks whether the flats exist in the right place 
%             for the surfaces in ganatret. S
%
function retval = checkFlats(subjectDir)

% check arguments
if ~any(nargin == [0 1])
  help checkFlats
  return
end

% volume directory
volumeDir = mrGetPref('volumeDirectory');
% check for volume directory
if ~isdir(volumeDir)
  disp(sprintf('(checkFlats) Could not find volumeDir: %s',volumeDir));
  return
end

% with no input arguments, check each one of the subject directories
if nargin == 0
  % get all sub directories
  volumeDirList = dir(fullfile(volumeDir,'s*'));
  % check each directory
  for i = 1:length(volumeDirList)
    if volumeDirList(i).isdir
      checkFlats(volumeDirList(i).name);
    end
  end
  return
end

% check the called for subjetDir
subjectDir = fullfile(volumeDir,subjectDir);
if ~isdir(subjectDir)
  disp(sprintf('(checkFlats) Could not find subjectDir: %s',subjectDir));
  return
end

% get baseAnatomy
baseAnatDir = getSubdir(subjectDir,{'mlrBaseAnatomies','baseAnatomies'});
if isempty(baseAnatDir)
  disp(sprintf('(checkFlats) No base anatomies found in %s',getLastDir(subjectDir)));
  return
end

% get surface directory
surfaceDir = getSubdir(subjectDir,{'surface','surfaces','surfRelax'});
if isempty(surfaceDir)
  % make a surface direcotry
  surfaceDir = fullfile(subjectDir,'surfaces');
  mkdir(surfaceDir);
end

% make empty view for loading the anatomies
v = makeEmptyView;

% load all the anatomies and get information about them
extList = {'hdr','nii'};
base = [];
copyFrom = {};copyTo = {};surfaceFiles = {};canonical = {};
for iExt = 1:length(extList)
  baseAnatDirList = dir(fullfile(baseAnatDir,sprintf('*.%s',extList{iExt})));
  for iFile = 1:length(baseAnatDirList)
    % load the anatomy
    v = loadAnat(v,baseAnatDirList(iFile).name,baseAnatDir);
    iBase = viewGet(v,'curBase');
    % get info
    base(iBase).name = viewGet(v,'baseName');
    base(iBase).type = viewGet(v,'baseType');
    base(iBase).dir = viewGet(v,'baseCoordMapPath');
    coordMap = viewGet(v,'baseCoordMap');
    base(iBase).filenameFrom = {};
    base(iBase).filenameTo = {};
    base(iBase).fileExists = [];
    base(iBase).fileExistsSurfaceDir = [];
    base(iBase).resave = false;
    % get the names of the files that are needed
    if base(iBase).type == 1
      base(iBase).filenameFrom = {'flatFileName','innerCoordsFileName','outerCoordsFileName','curvFileName','anatFileName'};
      base(iBase).filenameTo = {setext(base(iBase).name,'off')};
      % check whether flat name and off match
      if ~strcmp(stripext(base(iBase).name),stripext(coordMap.flatFileName))
	base(iBase).resave = true;
	% change the coord map and put it back in the view
	coordMap.flatFileName = setext(base(iBase).name,'off');
	v = viewSet(v,'baseCoordMap',coordMap);
      end
    elseif base(iBase).type == 2
      base(iBase).filenameFrom = {'innerSurfaceFileName','innerCoordsFileName','outerSurfaceFileName','outerCoordsFileName','curvFileName','anatFileName'};
      base(iBase).filenameTo = {};
    end
    % check for the files
    for iFilename = 1:length(base(iBase).filenameFrom)
      % if its canonical keep in special list
      if strcmp(base(iBase).filenameFrom{iFilename},'anatFileName')
	canonical{end+1} = stripext(coordMap.(base(iBase).filenameFrom{iFilename}));
      end
      % get filename
      base(iBase).filenameFrom{iFilename} = coordMap.(base(iBase).filenameFrom{iFilename});
      % check for file existence
      if isfile(fullfile(base(iBase).dir,base(iBase).filenameFrom{iFilename}))
	base(iBase).fileExists(iFilename) = true;
      else
	base(iBase).fileExists(iFilename) = false;
      end
      % check for file existence
      if isfile(fullfile(surfaceDir,base(iBase).filenameFrom{iFilename}))
	base(iBase).fileExistsSurfaceDir(iFilename) = true;
      else
	base(iBase).fileExistsSurfaceDir(iFilename) = false;
      end
      % set the to filename
      if length(base(iBase).filenameTo) < iFilename
	base(iBase).filenameTo{iFilename} = base(iBase).filenameFrom{iFilename};
      end
      % if the filename to does not exist, then figure out if we can copy it
      % from someplace
      if ~isfile(fullfile(surfaceDir,base(iBase).filenameTo{iFilename}))
	if base(iBase).fileExists(iFilename)
	  copyFrom{end+1} = fullfile(base(iBase).dir,base(iBase).filenameFrom{iFilename});
	  copyTo{end+1} = fullfile(surfaceDir,base(iBase).filenameTo{iFilename});
	elseif base(iBase).fileExistsSurfaceDir(iFilename)
	  copyFrom{end+1} = fullfile(surfaceDir,base(iBase).filenameFrom{iFilename});
	  copyTo{end+1} = fullfile(surfaceDir,base(iBase).filenameTo{iFilename});
	end
      end
      surfaceFiles = union(surfaceFiles,{base(iBase).filenameTo{iFilename}});
    end
  end
end
canonical = unique(canonical);

% check to see which files aren't needed in surfacedir
surfaceDirList = dir(surfaceDir);
% get only filenames
surfaceDirList = {surfaceDirList(find(~[surfaceDirList(:).isdir])).name};
for i = 1:length(surfaceDirList)
  if any(strcmp(stripext(surfaceDirList{i}),canonical))
    isCanonical(i) = true;
  else
    isCanonical(i) = false;
  end
end
% find which ones are not in the surfaceFiles list from above
isSurfaceFile = ismember(surfaceDirList,surfaceFiles);
% delete files are neither surface or canonical
deleteSurfaceFiles = {surfaceDirList{~(isCanonical|isSurfaceFile)}};

% display info about what was found
dispHeader(getLastDir(subjectDir));
for iBase = 1:viewGet(v,'numBase')
  if any(base(iBase).type == [1 2])
    if base(iBase).type == 1
      baseType = 'flat';
    else
      baseType = 'surface';
    end
    fileExists = or(base(iBase).fileExists,base(iBase).fileExistsSurfaceDir);
    disp(sprintf('(checkFlats) Found %s: %s and %i/%i supporting files',baseType,base(iBase).name,sum(fileExists),length(base(iBase).fileExists)));
    % display what's missing
    if any(~fileExists)
      missingFiles = '';
      for i = 1:length(fileExists)
	if fileExists(i) == 0
	  missingFiles = sprintf('%s%s, ',missingFiles,base(iBase).filenameFrom{i});
	end
      end
      disp(sprintf('(checkFlats) !!! Missing files: %s !!!',missingFiles(1:end-2)));
    end
  end
end

dispHeader('todo');
for iBase = 1:viewGet(v,'numBase')
  if base(iBase).resave
    disp(sprintf('(checkFlats) Resave base anatomy: %s',base(iBase).name));
  end
end
for iCopy = 1:length(copyFrom)
  disp(sprintf('(checkFlats) Copy from %s to %s',copyFrom{iCopy},copyTo{iCopy}));
end
for iDelete = 1:length(deleteSurfaceFiles)
  disp(sprintf('(checkFlats) Delete %s',fullfile(surfaceDir,deleteSurfaceFiles{iDelete})));
end
keyboard

% To do
%  make list of things to do and display
% change flat name in surface file if necessary
% copy any files necessary to surface file
% change filename of flat if necessary
% give list of unaccounted for files
deleteView(v);

%%%%%%%%%%%%%%%%%%%
%    getSubdir    %
%%%%%%%%%%%%%%%%%%%
function subdirName = getSubdir(dirname,possibleSubdirs)

for i = 1:length(possibleSubdirs)
  subdirName = fullfile(dirname,possibleSubdirs{i});
  if isdir(subdirName),return,end
end
subdirName = [];

