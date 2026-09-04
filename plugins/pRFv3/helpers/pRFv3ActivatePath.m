function pluginRoot = pRFv3ActivatePath(pluginLocation)
% pRFv3ActivatePath
%
% Put the complete pRFv3 runtime tree ahead of older pRF versions and
% verify that every shared runtime name resolves back inside this copy of
% pRFv3. pluginLocation may point to the plugin root, any subdirectory, or
% a file within the plugin. With no input, the root is found by walking up
% from this function until pRFv3Plugin.m and pRF.m are found.

if nargin < 1 || isempty(pluginLocation)
  pluginLocation = mfilename('fullpath');
end
pluginRoot = findPluginRoot(pluginLocation);

% Use the tree that actually exists rather than prescribing which helper
% must live in which subdirectory. Add in reverse so the plugin root remains
% first while every current and future runtime subdirectory is available.
runtimePaths = strsplit(genpath(pluginRoot),pathsep);
runtimePaths = runtimePaths(~cellfun(@isempty,runtimePaths));
runtimePaths = unique(runtimePaths,'stable');
for iPath = numel(runtimePaths):-1:1
  addpath(runtimePaths{iPath},'-begin');
end
rehash path

% These are public or cross-file names used by pRFv3. Their precise
% subdirectories are intentionally not encoded here, so reorganizing files
% within pRFv3 cannot by itself break activation. The ownership check still
% catches missing functions and accidental resolution to pRFv2.
requiredFunctions = { ...
  'pRF','pRFFit','pRFGUI','pRFPlot','pRFv3Plugin','mlrRunPRF', ...
  'pRFFitPrepared','pRFFormatVoxelFit','pRFGetStimImageFromStimfile', ...
  'pRFInstallGetptsCompatibility','pRFIsShiftClick', ...
  'pRFProgressReporter','pRFResolveNumWorkers','pRFSaveAnalysis', ...
  'pRFv3ActivatePath','pRFv3Startup','pRFMergeParams', ...
  'pRFMakeGaussian','pRFPrepareStimulusProjection','pRFProjectStimulus', ...
  'pRFStimulusProjectionBank','pRF_DoG_CSS','pRF_css', ...
  'pRF_diffGaussian','pRF_divGaussian','pRF_divNorm', ...
  'pRF_gaussian','pRF_gaussianhdr','fminsearchbnd'};

for iFunction = 1:numel(requiredFunctions)
  functionName = requiredFunctions{iFunction};
  expectedPaths = {};
  for iPath = 1:numel(runtimePaths)
    candidatePath = fullfile(runtimePaths{iPath},[functionName '.m']);
    if isfile(candidatePath)
      expectedPaths{end+1} = canonicalPath(candidatePath);
    end
  end
  expectedPaths = unique(expectedPaths,'stable');
  if isempty(expectedPaths)
    error('pRFv3ActivatePath:MissingFunction', ...
      'Could not find %s.m anywhere inside %s.',functionName,pluginRoot);
  elseif numel(expectedPaths) > 1
    error('pRFv3ActivatePath:DuplicateFunction', ...
      'Found multiple pRFv3 copies of %s.m: %s', ...
      functionName,strjoin(expectedPaths,', '));
  end

  resolvedPath = which(functionName);
  if isempty(resolvedPath) || ~isfile(resolvedPath) || ...
      ~strcmpi(canonicalPath(resolvedPath),expectedPaths{1})
    if isempty(resolvedPath)
      resolvedPath = '<not found>';
    end
    error('pRFv3ActivatePath:ResolutionFailed', ...
      '%s resolved to %s instead of %s.', ...
      functionName,resolvedPath,expectedPaths{1});
  end
end

%%%%%%%%%%%%%%%%%%%%%%
%   findPluginRoot   %
%%%%%%%%%%%%%%%%%%%%%%
function pluginRoot = findPluginRoot(pluginLocation)

if ~(ischar(pluginLocation) || ...
    (isstring(pluginLocation) && isscalar(pluginLocation)))
  error('pRFv3ActivatePath:InvalidPluginLocation', ...
    'The pRFv3 location must be a character vector or scalar string.');
end
pluginLocation = char(pluginLocation);

if isfile(pluginLocation)
  searchDir = fileparts(pluginLocation);
elseif isfile([pluginLocation '.m'])
  % MFILE('fullpath') reports a function file without its .m extension.
  searchDir = fileparts([pluginLocation '.m']);
elseif isfolder(pluginLocation)
  searchDir = pluginLocation;
else
  error('pRFv3ActivatePath:MissingPluginLocation', ...
    'Could not find the supplied pRFv3 location: %s',pluginLocation);
end

% Convert relative input to the absolute spelling MATLAB reports from
% WHICH. This also makes the containment check independent of Current Folder.
[locationExists,locationInfo] = fileattrib(searchDir);
if locationExists
  searchDir = locationInfo.Name;
end

while true
  if isfile(fullfile(searchDir,'pRFv3Plugin.m')) && ...
      isfile(fullfile(searchDir,'pRF.m'))
    pluginRoot = searchDir;
    return
  end
  parentDir = fileparts(searchDir);
  if isempty(parentDir) || strcmp(parentDir,searchDir)
    break
  end
  searchDir = parentDir;
end

error('pRFv3ActivatePath:PluginRootNotFound', ...
  ['Could not locate the pRFv3 root above %s. A valid root contains both ' ...
   'pRFv3Plugin.m and pRF.m.'],pluginLocation);

%%%%%%%%%%%%%%%%%%%%%%
%   canonicalPath    %
%%%%%%%%%%%%%%%%%%%%%%
function pathname = canonicalPath(pathname)

[pathExists,pathInfo] = fileattrib(pathname);
if pathExists
  pathname = pathInfo.Name;
end
