function [stim,summary] = pRFPrepareStimulusProjection(stim)
% pRFPrepareStimulusProjection  Prepare exact duplicate-frame projection.
%
%   preparedStim = pRFPrepareStimulusProjection(stim)
%   [preparedStim,summary] = pRFPrepareStimulusProjection(stim)
%
% STIM can be the usual cell array of stimulus structs (each containing an
% .im movie), one stimulus struct, or one numeric movie. Compatible runs are
% combined before exact duplicate detection, so identical frames and whole
% runs share one immutable bank. The original .im fields are retained.
%
% Prepared stimuli should be treated as immutable. Re-run this function
% after changing any .im movie.

[stimCells,containerType,containerSize] = normalizeStimulusContainer(stim);
nRuns = numel(stimCells);

summary.nRuns = nRuns;
summary.nPreparedRuns = 0;
summary.nGroups = 0;
summary.nFrames = 0;
summary.nUniqueFrames = 0;
summary.nDuplicateFrames = 0;
summary.nDuplicateRuns = 0;

% Remove stale metadata before rebuilding it. This makes repeated calls on
% an already-prepared context safe and gives every preparation a new bank
% identity without depending on counters, hashes, or random identifiers.
for iRun = 1:nRuns
  if isstruct(stimCells{iRun}) && isscalar(stimCells{iRun}) && ...
      isfield(stimCells{iRun},'projection')
    stimCells{iRun} = rmfield(stimCells{iRun},'projection');
  end
end

groupKeys = {};
groupRuns = {};
for iRun = 1:nRuns
  images = getStimulusImages(stimCells{iRun});
  if ~isSupportedMovie(images)
    continue
  end
  key = sprintf('%s|%d|%d|%d',class(images),size(images,1),size(images,2),isreal(images));
  iGroup = find(strcmp(groupKeys,key),1);
  if isempty(iGroup)
    groupKeys{end+1} = key;
    groupRuns{end+1} = iRun;
  else
    groupRuns{iGroup}(end+1) = iRun;
  end
end

for iGroup = 1:numel(groupRuns)
  runIndices = groupRuns{iGroup};
  runFrameCounts = zeros(1,numel(runIndices));
  movies = cell(1,numel(runIndices));
  for iMember = 1:numel(runIndices)
    movies{iMember} = getStimulusImages(stimCells{runIndices(iMember)});
    runFrameCounts(iMember) = size(movies{iMember},3);
  end
  nFrames = sum(runFrameCounts);
  if nFrames == 0
    continue
  end

  allImages = cat(3,movies{:});
  flattenedFrames = reshape(allImages,[],nFrames).';
  [~,firstOccurrences,frameToUnique] = unique(flattenedFrames,'rows','stable');
  uniqueImages = allImages(:,:,firstOccurrences);
  nUnique = numel(firstOccurrences);

  summary.nGroups = summary.nGroups+1;
  summary.nFrames = summary.nFrames+nFrames;
  summary.nUniqueFrames = summary.nUniqueFrames+nUnique;
  summary.nDuplicateFrames = summary.nDuplicateFrames+(nFrames-nUnique);

  % With no duplicate in the group, retaining a second copy of the entire
  % movie would add memory without avoiding a projection.
  if nUnique == nFrames
    continue
  end

  bank = pRFStimulusProjectionBank(uniqueImages);
  firstFrame = 1;
  preparedMaps = cell(1,numel(runIndices));
  for iMember = 1:numel(runIndices)
    nRunFrames = runFrameCounts(iMember);
    lastFrame = firstFrame+nRunFrames-1;
    thisMap = reshape(frameToUnique(firstFrame:lastFrame),1,[]);
    if nUnique <= double(intmax('uint32'))
      thisMap = uint32(thisMap);
    end
    preparedMaps{iMember} = thisMap;

    projection.version = 1;
    projection.bank = bank;
    projection.frameToUnique = thisMap;
    stimCells{runIndices(iMember)} = setProjection( ...
      stimCells{runIndices(iMember)},movies{iMember},projection);
    firstFrame = lastFrame+1;
  end
  summary.nPreparedRuns = summary.nPreparedRuns+numel(runIndices);

  for iMember = 2:numel(preparedMaps)
    for iEarlier = 1:iMember-1
      if isequal(preparedMaps{iMember},preparedMaps{iEarlier})
        summary.nDuplicateRuns = summary.nDuplicateRuns+1;
        break
      end
    end
  end
end

stim = restoreStimulusContainer(stimCells,containerType,containerSize);


function [stimCells,containerType,containerSize] = normalizeStimulusContainer(stim)

containerSize = size(stim);
if iscell(stim)
  containerType = 'cell';
  stimCells = reshape(stim,1,[]);
elseif isstruct(stim)
  containerType = 'struct';
  stimCells = num2cell(reshape(stim,1,[]));
else
  containerType = 'numeric';
  stimCells = {struct('im',stim)};
end


function stim = restoreStimulusContainer(stimCells,containerType,containerSize)

switch containerType
  case 'cell'
    stim = reshape(stimCells,containerSize);
  case 'struct'
    % A MATLAB struct array must have the same fields in every element. A
    % compatibility group with duplicates receives projection metadata,
    % while an all-unique group intentionally does not. Add an empty field
    % to the latter before concatenating the scalar structs again.
    for iStim = 1:numel(stimCells)
      if ~isfield(stimCells{iStim},'projection')
        stimCells{iStim}.projection = [];
      end
    end
    stim = reshape([stimCells{:}],containerSize);
  otherwise
    stim = stimCells{1};
end


function images = getStimulusImages(stimulus)

if isstruct(stimulus) && isscalar(stimulus) && isfield(stimulus,'im')
  images = stimulus.im;
elseif isnumeric(stimulus) || islogical(stimulus)
  images = stimulus;
else
  images = [];
end


function tf = isSupportedMovie(images)

tf = (isa(images,'double') || isa(images,'single') || islogical(images)) && ...
     isreal(images) && ~issparse(images) && ndims(images) <= 3;


function stimulus = setProjection(stimulus,images,projection)

if isstruct(stimulus) && isscalar(stimulus)
  stimulus.projection = projection;
else
  stimulus = struct('im',images,'projection',projection);
end
