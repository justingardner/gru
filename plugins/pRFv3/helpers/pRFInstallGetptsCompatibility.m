% pRFInstallGetptsCompatibility
%
%      usage: pRFInstallGetptsCompatibility(<figureHandle>)
%             pRFInstallGetptsCompatibility('install',<figureHandle>)
%             pRFInstallGetptsCompatibility('teardown')
%         by: Austin Kuo with Codex
%       date: 09/02/2026
%    purpose: Guard legacy getpts modifier-key handling while pRF is active.
%
function pRFInstallGetptsCompatibility(action,figureHandle)

if nargin < 1
  action = 'install';
  figureHandle = [];
elseif ~(ischar(action) || (isstring(action) && isscalar(action)))
  % Preserve the original one-argument calling convention, in which the
  % first argument is the figure to guard.
  figureHandle = action;
  action = 'install';
elseif nargin < 2
  figureHandle = [];
end
action = lower(char(action));

rootListenerAppData = 'pRFGetptsCompatibilityRootListener';
figureListenerAppData = 'pRFGetptsCompatibilityListener';
compatibilityVersion = 3;
compatibilityOwner = mfilename('fullpath');

switch action
  case {'teardown','remove','uninstall'}
    teardownCompatibility(rootListenerAppData,figureListenerAppData);
    return
  case {'install','i'}
    % Continue below.
  otherwise
    error('pRFInstallGetptsCompatibility:UnknownAction', ...
      'Unknown compatibility action ''%s''.',action);
end

% mrTools sometimes calls GETPTS without an explicit figure, so it binds to
% whichever figure is current. Watch figure switches and install the guard
% before GETPTS assigns its legacy KeyPressFcn on the selected figure.
rootRecord = getappdata(groot,rootListenerAppData);
if listenerRecordMatches(rootRecord,compatibilityOwner,compatibilityVersion)
  installationToken = rootRecord.token;
else
  removeListenerAppData(groot,rootListenerAppData);
  installationToken = tempname;
  propertyListener = addlistener(groot,'CurrentFigure','PostSet', ...
    @(~,~) installFigureGuard(get(groot,'CurrentFigure'),installationToken, ...
      compatibilityOwner,compatibilityVersion,figureListenerAppData));
  rootRecord = makeListenerRecord(propertyListener,installationToken, ...
    compatibilityOwner,compatibilityVersion);
  setappdata(groot,rootListenerAppData,rootRecord);
end

% Cover the requested view, any currently active point-selection figure,
% and figures that existed before the pRF plugin was installed or reloaded.
installFigureGuard(figureHandle,installationToken,compatibilityOwner, ...
  compatibilityVersion,figureListenerAppData);
figureHandles = findall(groot,'Type','figure');
for iFigure = 1:numel(figureHandles)
  installFigureGuard(figureHandles(iFigure),installationToken, ...
    compatibilityOwner,compatibilityVersion,figureListenerAppData);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeListenerRecord     %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function listenerRecord = makeListenerRecord(propertyListener,token,owner,version)

listenerRecord = struct('listener',propertyListener,'token',token, ...
  'owner',owner,'version',version);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% listenerRecordMatches  %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = listenerRecordMatches(listenerRecord,owner,version,token)

tf = isstruct(listenerRecord) && isscalar(listenerRecord) && ...
  all(isfield(listenerRecord,{'listener','token','owner','version'})) && ...
  isequal(listenerRecord.owner,owner) && ...
  isequal(listenerRecord.version,version) && ...
  isListenerValid(listenerRecord.listener);
if tf && nargin >= 4
  tf = isequal(listenerRecord.token,token);
end

%%%%%%%%%%%%%%%%%%%%%%%%
% isListenerValid      %
%%%%%%%%%%%%%%%%%%%%%%%%
function tf = isListenerValid(propertyListener)

tf = false;
if isempty(propertyListener)
  return
end
try
  tf = isvalid(propertyListener);
catch
  tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%
% installFigureGuard   %
%%%%%%%%%%%%%%%%%%%%%%%%
function installFigureGuard(figureHandle,token,owner,version,listenerAppData)

if isempty(figureHandle) || ~isgraphics(figureHandle,'figure')
  return
end

listenerRecord = getappdata(figureHandle,listenerAppData);
if ~listenerRecordMatches(listenerRecord,owner,version,token)
  removeListenerAppData(figureHandle,listenerAppData);
  propertyListener = addlistener(figureHandle,'KeyPressFcn','PostSet', ...
    @(~,~) wrapLegacyGetptsKeyPress(figureHandle));
  listenerRecord = makeListenerRecord(propertyListener,token,owner,version);
  setappdata(figureHandle,listenerAppData,listenerRecord);
end

% Also repair a callback that was already installed before this helper ran.
wrapLegacyGetptsKeyPress(figureHandle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wrapLegacyGetptsKeyPress   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wrapLegacyGetptsKeyPress(figureHandle)

if ~isgraphics(figureHandle,'figure')
  return
end

keyPressFcn = get(figureHandle,'KeyPressFcn');
currentCompatibilityCallback = @getptsKeyPressCompat;
if isequal(keyPressFcn,currentCompatibilityCallback)
  return
end

replaceCallback = false;
if ischar(keyPressFcn) || (isstring(keyPressFcn) && isscalar(keyPressFcn))
  callbackText = regexprep(char(keyPressFcn),'\s','');
  replaceCallback = any(strcmp(callbackText, ...
    {'getpts(''KeyPress'')','getpts(''KeyPress'');'}));
elseif isa(keyPressFcn,'function_handle')
  % Function handles retain the source file in which they were created.
  % Replace a compatibility callback from an older version/location before
  % that removed file can be invoked by a figure event.
  replaceCallback = isGetptsCompatibilityCallback(keyPressFcn);
end

if replaceCallback
  set(figureHandle,'KeyPressFcn',currentCompatibilityCallback);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% isGetptsCompatibilityCallback          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = isGetptsCompatibilityCallback(callback)

tf = false;
if ~isa(callback,'function_handle')
  return
end
try
  callbackName = func2str(callback);
  tf = ~isempty(regexp(callbackName, ...
    '(^|[>/])getptsKeyPressCompat$','once'));
catch
  tf = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% getptsKeyPressCompat    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function getptsKeyPressCompat(source,~)

% GETPTS uses these globals for the active point-selection state. If the
% callback was left behind after an interruption, quietly ignore all keys.
global GETPTS_FIG GETPTS_H1 GETPTS_H2
if isempty(GETPTS_FIG) || ~isequal(source,GETPTS_FIG) || ...
    isempty(GETPTS_H1) || ~ishghandle(GETPTS_H1) || ...
    isempty(GETPTS_H2) || ~ishghandle(GETPTS_H2)
  return
end

key = get(source,'CurrentCharacter');
if ~ischar(key) || numel(key) ~= 1
  return
end

% Forward supported character keys to MATLAB's original implementation.
% Modifier-only keys, including Shift, return above. The subsequent mouse
% callback still sees SelectionType='extend', preserving Shift-click.
getpts('KeyPress');

%%%%%%%%%%%%%%%%%%%%%%%%%%
% teardownCompatibility  %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function teardownCompatibility(rootListenerAppData,figureListenerAppData)

% Remove the root listener first so changing or deleting figure state cannot
% cause it to reinstall a figure listener during teardown.
removeListenerAppData(groot,rootListenerAppData);
figureHandles = findall(groot,'Type','figure');
for iFigure = 1:numel(figureHandles)
  removeListenerAppData(figureHandles(iFigure),figureListenerAppData);
  restoreLegacyGetptsKeyPress(figureHandles(iFigure));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restoreLegacyGetptsKeyPress        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function restoreLegacyGetptsKeyPress(figureHandle)

if ~isgraphics(figureHandle,'figure')
  return
end
keyPressFcn = get(figureHandle,'KeyPressFcn');
if isGetptsCompatibilityCallback(keyPressFcn)
  set(figureHandle,'KeyPressFcn','getpts(''KeyPress'');');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% removeListenerAppData  %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function removeListenerAppData(graphicsHandle,appDataName)

if ~isappdata(graphicsHandle,appDataName)
  return
end

listenerRecord = getappdata(graphicsHandle,appDataName);
rmappdata(graphicsHandle,appDataName);
if isstruct(listenerRecord) && isscalar(listenerRecord) && ...
    isfield(listenerRecord,'listener')
  propertyListener = listenerRecord.listener;
else
  % Compatibility with the unversioned implementation, which stored the
  % listener handle directly in appdata.
  propertyListener = listenerRecord;
end
if ~isempty(propertyListener)
  try
    delete(propertyListener);
  catch
    % A stale listener whose defining file was removed may already be
    % unusable. Removing its appdata is still sufficient to replace it.
  end
end
