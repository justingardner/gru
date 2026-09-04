function [v,savedAnalysisName] = pRFSaveAnalysis(v,analysisName,forceMerge,varargin)
% pRFSaveAnalysis
%
% Save a pRF analysis without opening an overwrite or rename dialog.
%
% [v savedAnalysisName] = pRFSaveAnalysis(v,analysisName,forceMerge)
%
% If the target already exists, the current mrTools overwrite preference is
% honored. Ask is replaced with a Command Window choice, while Rename asks
% for a new analysis name in the Command Window. forceMerge preserves the
% automatic merge behavior used when continuing an existing pRF analysis.

if nargin < 3 || isempty(forceMerge)
  forceMerge = false;
end
if ~((islogical(forceMerge) || isnumeric(forceMerge)) && ...
    isreal(forceMerge) && isscalar(forceMerge) && isfinite(forceMerge))
  error('pRFSaveAnalysis:InvalidForceMerge', ...
    'forceMerge must be a logical or numeric scalar.');
end
forceMerge = logical(forceMerge);

parser = inputParser;
parser.addParameter('Dependencies',struct,@(value) isstruct(value) && isscalar(value));
parser.parse(varargin{:});
dependencies = makeDependencies(parser.Results.Dependencies);

if isstring(analysisName) && isscalar(analysisName)
  analysisName = char(analysisName);
end
if ~ischar(analysisName) || ~isrow(analysisName) || isempty(strtrim(analysisName))
  error('pRFSaveAnalysis:InvalidAnalysisName', ...
    'analysisName must be a nonempty character vector or string scalar.');
end

savedAnalysisName = analysisName;
analysisNum = dependencies.viewGetFcn(v,'analysisNum',analysisName);
if isempty(analysisNum)
  error('pRFSaveAnalysis:AnalysisNotFound', ...
    'Could not find the in-memory analysis named %s.',analysisName);
end
analysisDir = dependencies.viewGetFcn(v,'analysisdir',[],analysisNum);
targetFilename = fullfile(analysisDir,[analysisName '.mat']);

% The ordinary no-conflict path remains exactly the mrTools save path.
if ~dependencies.fileExistsFcn(targetFilename)
  dependencies.saveFcn(v,analysisName);
  return
end

originalPolicy = dependencies.getPrefFcn('overwritePolicy');
if forceMerge
  saveAction = 'Merge';
else
  saveAction = actionFromPolicy(originalPolicy);
  if isempty(saveAction)
    saveAction = promptForAction(analysisName,targetFilename,dependencies.inputFcn);
  end
end

if strcmp(saveAction,'Cancel')
  fprintf('(pRF) Analysis remains loaded in MATLAB but was not saved.\n');
  return
elseif strcmp(saveAction,'Copy')
  existingAnalysisNames = dependencies.viewGetFcn(v,'analysisNames');
  [copyName,cancelCopy] = promptForCopyName(analysisName,analysisDir, ...
    existingAnalysisNames,dependencies);
  if cancelCopy
    fprintf('(pRF) Analysis remains loaded in MATLAB but was not saved.\n');
    return
  end

  % Match saveAnalysis's Rename behavior by renaming the newly installed
  % in-memory analysis before saving. Roll the name back if saving fails.
  v = dependencies.viewSetFcn(v,'analysisName',copyName,analysisNum);
  try
    dependencies.saveFcn(v,copyName);
  catch saveException
    try
      dependencies.viewSetFcn(v,'analysisName',analysisName,analysisNum);
    catch rollbackException
      warning('pRFSaveAnalysis:RenameRollbackFailed', ...
        'Could not restore the analysis name after save failure: %s', ...
        rollbackException.message);
    end
    rethrow(saveException);
  end
  savedAnalysisName = copyName;
  return
end

% saveAnalysis already implements the established merge and overwrite
% semantics. Select that policy only for this call, restoring the user's
% preference even if saving throws or the run is interrupted.
saveWithTemporaryPolicy(v,analysisName,saveAction,originalPolicy,dependencies);


function dependencies = makeDependencies(overrides)

dependencies.inputFcn = @input;
dependencies.saveFcn = @saveAnalysis;
dependencies.getPrefFcn = @mrGetPref;
% saveAnalysis reads this preference through the global mrDEFAULTS value.
% Change only that in-memory value for the duration of this save instead of
% writing a temporary policy to ~/.mrDefaults twice. This also prevents a
% hard MATLAB crash from persisting the temporary choice to the next run.
dependencies.setPrefFcn = @setPreferenceInMemory;
dependencies.fileExistsFcn = @(filename) exist(filename,'file') == 2;
dependencies.viewGetFcn = @viewGet;
dependencies.viewSetFcn = @viewSet;

overrideNames = fieldnames(overrides);
dependencyNames = fieldnames(dependencies);
for iOverride = 1:numel(overrideNames)
  overrideName = overrideNames{iOverride};
  if ~ismember(overrideName,dependencyNames)
    error('pRFSaveAnalysis:UnknownDependency', ...
      'Unknown dependency override: %s',overrideName);
  end
  if ~isa(overrides.(overrideName),'function_handle')
    error('pRFSaveAnalysis:InvalidDependency', ...
      'Dependency %s must be a function handle.',overrideName);
  end
  dependencies.(overrideName) = overrides.(overrideName);
end


function action = actionFromPolicy(policy)

action = '';
if isstring(policy) && isscalar(policy)
  policy = char(policy);
end
if ~ischar(policy),return,end

switch lower(strtrim(policy))
  case 'merge'
    action = 'Merge';
  case 'overwrite'
    action = 'Overwrite';
  case 'rename'
    action = 'Copy';
end


function action = promptForAction(analysisName,targetFilename,inputFcn)

fprintf('\n(pRF) An analysis named "%s" already exists:\n  %s\n', ...
  analysisName,targetFilename);
fprintf('  [M] Merge with the existing analysis\n');
fprintf('  [O] Overwrite the existing analysis\n');
fprintf('  [C] Save this result as a separate copy\n');
fprintf('  [Q] Cancel without saving\n');

action = '';
while isempty(action)
  [response,inputAvailable] = readCommandLine( ...
    inputFcn,'Choose M, O, C, or Q [default: M]: ');
  if ~inputAvailable
    action = 'Cancel';
    return
  end
  switch lower(strtrim(response))
    case {'','m','merge','1'}
      action = 'Merge';
    case {'o','overwrite','2'}
      action = 'Overwrite';
    case {'c','copy','r','rename','s','separate','3'}
      action = 'Copy';
    case {'q','quit','cancel','4'}
      action = 'Cancel';
    otherwise
      fprintf('(pRF) Invalid choice "%s". Enter M, O, C, or Q.\n',response);
  end
end


function [copyName,cancelCopy] = promptForCopyName(analysisName,analysisDir, ...
    existingAnalysisNames,dependencies)

if isempty(existingAnalysisNames)
  existingAnalysisNames = {};
elseif ischar(existingAnalysisNames)
  existingAnalysisNames = cellstr(existingAnalysisNames);
elseif isstring(existingAnalysisNames)
  existingAnalysisNames = cellstr(existingAnalysisNames);
end

defaultCopyName = findDefaultCopyName( ...
  analysisName,analysisDir,existingAnalysisNames,dependencies.fileExistsFcn);
copyName = '';
cancelCopy = false;
while isempty(copyName)
  prompt = sprintf( ...
    'New analysis name [default: %s; type Q to cancel]: ',defaultCopyName);
  [response,inputAvailable] = readCommandLine(dependencies.inputFcn,prompt);
  if ~inputAvailable
    cancelCopy = true;
    return
  end
  response = strtrim(response);
  if isempty(response)
    candidate = defaultCopyName;
  elseif any(strcmpi(response,{'q','quit','cancel'}))
    cancelCopy = true;
    return
  else
    candidate = removeMatExtension(response);
  end

  [isValid,rejectionReason] = isAvailableCopyName( ...
    candidate,analysisDir,existingAnalysisNames,dependencies.fileExistsFcn);
  if isValid
    copyName = candidate;
  else
    fprintf('(pRF) %s\n',rejectionReason);
  end
end


function defaultName = findDefaultCopyName(analysisName,analysisDir, ...
    existingAnalysisNames,fileExistsFcn)

copyNumber = 1;
defaultName = [analysisName '_copy'];
while ~isAvailableCopyName(defaultName,analysisDir,existingAnalysisNames,fileExistsFcn)
  copyNumber = copyNumber+1;
  defaultName = sprintf('%s_copy%i',analysisName,copyNumber);
end


function [isAvailable,rejectionReason] = isAvailableCopyName(candidate, ...
    analysisDir,existingAnalysisNames,fileExistsFcn)

candidate = strtrim(candidate);
isAvailable = false;
rejectionReason = '';
if isempty(candidate)
  rejectionReason = 'The new analysis name cannot be empty.';
  return
end
if any(candidate == '/') || any(candidate == '\') || any(candidate == ':') || ...
    strcmp(candidate,'.') || strcmp(candidate,'..')
  rejectionReason = 'Enter an analysis name, not a path.';
  return
end
if any(strcmpi(candidate,existingAnalysisNames))
  rejectionReason = sprintf('An in-memory analysis named "%s" already exists.',candidate);
  return
end
candidateFilename = fullfile(analysisDir,[candidate '.mat']);
if fileExistsFcn(candidateFilename)
  rejectionReason = sprintf('The analysis file %s already exists.',candidateFilename);
  return
end
isAvailable = true;


function name = removeMatExtension(name)

name = strtrim(name);
if length(name) >= 4 && strcmpi(name(end-3:end),'.mat')
  name = strtrim(name(1:end-4));
end


function [response,inputAvailable] = readCommandLine(inputFcn,prompt)

response = '';
inputAvailable = true;
try
  response = inputFcn(prompt,'s');
  if isstring(response) && isscalar(response)
    response = char(response);
  end
  if ~ischar(response)
    response = '';
  end
catch inputException
  warning('pRFSaveAnalysis:CommandLineInputUnavailable', ...
    ['Could not read a Command Window save choice (%s). The analysis ' ...
     'will remain loaded but will not be saved.'],inputException.message);
  inputAvailable = false;
end


function saveWithTemporaryPolicy(v,analysisName,saveAction,originalPolicy,dependencies)

if ischar(originalPolicy) && strcmpi(originalPolicy,saveAction)
  dependencies.saveFcn(v,analysisName);
  return
end

restoreCleanup = onCleanup(@() restoreOverwritePolicy(originalPolicy,dependencies));
dependencies.setPrefFcn('overwritePolicy',saveAction);
dependencies.saveFcn(v,analysisName);
clear restoreCleanup


function restoreOverwritePolicy(originalPolicy,dependencies)

try
  dependencies.setPrefFcn('overwritePolicy',originalPolicy);
catch restoreException
  warning('pRFSaveAnalysis:PreferenceRestoreFailed', ...
    'Could not restore the overwritePolicy preference: %s', ...
    restoreException.message);
end


function setPreferenceInMemory(preferenceName,value)

global mrDEFAULTS %#ok<GVMIS> mrTools stores live preferences in this global.
if isempty(mrDEFAULTS) || ~isfield(mrDEFAULTS,'prefs')
  % Initialize the normal mrTools defaults without changing them on disk.
  mrGetPref(preferenceName);
end
mrDEFAULTS.prefs.(preferenceName) = value;
