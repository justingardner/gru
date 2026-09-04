% pRFProjectStimulus
%
% Project a spatial receptive-field model onto every stimulus frame.
%
% Dense double-precision inputs are flattened and projected with one matrix
% multiplication. Mixed-class, complex, sparse, and otherwise unsupported
% inputs retain the general frame-by-frame calculation as a safe fallback.
%
% The legacy calculation is
%
%   sum(sum(rfModel .* stimulusImages(:,:,frameNum)))
%
function modelResponse = pRFProjectStimulus(rfModel,stimulusImages,nOutputFrames)

% A prepared stimulus carries an immutable shared bank of unique frames and
% a run-specific expansion map. The projection cache is local to each
% MATLAB process (and therefore to each parallel worker). It is deliberately
% bounded; the conservative byte accounting includes the stimulus bank for
% every entry so retained memory cannot grow with the number of fits.
persistent projectionCache cacheClock cacheStats cacheInitialized
if isempty(cacheInitialized)
  projectionCache = emptyProjectionCache;
  cacheClock = uint64(0);
  cacheStats = emptyCacheStats;
  cacheInitialized = true;
end

% Explicit cache control is useful when a long-lived worker releases a scan
% context. It is not required for correctness: cache entries compare the
% identity of immutable bank handle objects rather than a reusable numeric
% identifier or a probabilistic hash.
if nargin == 1 && (ischar(rfModel) || (isstring(rfModel) && isscalar(rfModel)))
  switch lower(char(rfModel))
    case {'reset','resetcache'}
      projectionCache = emptyProjectionCache;
      cacheClock = uint64(0);
      cacheStats = emptyCacheStats;
      cacheInitialized = true;
      modelResponse = [];
      return
    case {'cacheinfo','cachestats'}
      modelResponse = cacheStats;
      modelResponse.entries = numel(projectionCache);
      modelResponse.bytes = sum([projectionCache.bytes]);
      return
    otherwise
      error('pRFProjectStimulus:UnknownCommand','Unknown command: %s',char(rfModel));
  end
end

if nargin < 2
  error('pRFProjectStimulus:NotEnoughInputs','RF model and stimulus are required.');
end

% Model functions pass the complete stimulus struct so prepared metadata is
% available. Direct callers can continue to pass a numeric movie exactly as
% before.
preparedProjection = [];
if isstruct(stimulusImages) && isscalar(stimulusImages) && isfield(stimulusImages,'im')
  stimulusStruct = stimulusImages;
  stimulusImages = stimulusStruct.im;
  if isfield(stimulusStruct,'projection')
    preparedProjection = stimulusStruct.projection;
  end
end

nStimFrames = size(stimulusImages,3);
if nargin < 3 || isempty(nOutputFrames)
  nOutputFrames = nStimFrames;
end

if isValidPreparedProjection(preparedProjection,stimulusImages,nStimFrames)
  bank = preparedProjection.bank;
  [uniqueResponse,projectionCache,cacheClock,cacheStats] = ...
    projectPreparedBank(rfModel,bank,projectionCache,cacheClock,cacheStats);

  % Match legacy preallocation and expansion behavior, including the case in
  % which nOutputFrames is shorter or longer than the stimulus movie.
  modelResponse = zeros(1,nOutputFrames);
  if nStimFrames ~= 0
    frameToUnique = double(preparedProjection.frameToUnique);
    modelResponse(1:nStimFrames) = uniqueResponse(frameToUnique);
  end
  return
end

modelResponse = projectNumericFrames(rfModel,stimulusImages,nOutputFrames);


function isValid = isValidPreparedProjection(projection,stimulusImages,nStimFrames)

isValid = isstruct(projection) && isscalar(projection) && ...
          isfield(projection,'bank') && ...
          isa(projection.bank,'pRFStimulusProjectionBank') && ...
          isscalar(projection.bank) && isvalid(projection.bank) && ...
          isfield(projection,'frameToUnique') && ...
          numel(projection.frameToUnique) == nStimFrames;
if ~isValid
  return
end

uniqueImages = projection.bank.Images;
frameToUnique = double(projection.frameToUnique(:));
isValid = (size(uniqueImages,1) == size(stimulusImages,1)) && ...
          (size(uniqueImages,2) == size(stimulusImages,2)) && ...
          strcmp(class(uniqueImages),class(stimulusImages)) && ...
          isreal(uniqueImages) == isreal(stimulusImages) && ...
          all(isfinite(frameToUnique)) && ...
          all(frameToUnique == fix(frameToUnique)) && ...
          all(frameToUnique >= 1) && ...
          all(frameToUnique <= size(uniqueImages,3));


function [response,cache,clock,stats] = projectPreparedBank( ...
  rfModel,bank,cache,clock,stats)

clock = clock+1;
for iEntry = 1:numel(cache)
  if cache(iEntry).bank == bank && exactArrayMatch(cache(iEntry).rfModel,rfModel)
    response = cache(iEntry).response;
    cache(iEntry).lastUsed = clock;
    stats.hits = stats.hits+1;
    return
  end
end

stats.misses = stats.misses+1;
uniqueImages = bank.Images;
response = projectNumericFrames(rfModel,uniqueImages,size(uniqueImages,3));

% Avoid retaining pathological inputs. Byte accounting is intentionally
% conservative: charging the complete bank to each RF entry overestimates
% memory when entries share the same immutable handle.
maxCacheBytes = 128*2^20;
maxCacheEntries = 64;
entryBytes = double(bank.NumBytes)+arrayBytes(rfModel)+arrayBytes(response);
if entryBytes > maxCacheBytes
  stats.bypasses = stats.bypasses+1;
  return
end

while ~isempty(cache) && ((sum([cache.bytes])+entryBytes > maxCacheBytes) || ...
    (numel(cache) >= maxCacheEntries))
  [~,leastRecentlyUsed] = min([cache.lastUsed]);
  cache(leastRecentlyUsed) = [];
  stats.evictions = stats.evictions+1;
end

newEntry.bank = bank;
newEntry.rfModel = rfModel;
newEntry.response = response;
newEntry.bytes = entryBytes;
newEntry.lastUsed = clock;
cache(end+1) = newEntry;
stats.stores = stats.stores+1;


function tf = exactArrayMatch(firstValue,secondValue)

% ISEQUALN intentionally regards equal-valued numeric classes as equal in
% MATLAB. Cache identity must be stricter: arithmetic class, complexness,
% and sparse storage can select different projection code paths or results.
tf = strcmp(class(firstValue),class(secondValue)) && ...
     (isreal(firstValue) == isreal(secondValue)) && ...
     (issparse(firstValue) == issparse(secondValue)) && ...
     isequaln(firstValue,secondValue);


function nBytes = arrayBytes(value)

valueInfo = whos('value');
nBytes = double(valueInfo.bytes);


function cache = emptyProjectionCache

cache = struct('bank',{},'rfModel',{},'response',{},'bytes',{},'lastUsed',{});


function stats = emptyCacheStats

stats = struct('hits',0,'misses',0,'stores',0,'evictions',0,'bypasses',0);


function modelResponse = projectNumericFrames(rfModel,stimulusImages,nOutputFrames)

nStimFrames = size(stimulusImages,3);

% Match the legacy functions, which always preallocated a double row.
modelResponse = zeros(1,nOutputFrames);
if nStimFrames == 0
  return
end

% Limit the optimized path to the built-in, real array types exercised by
% the pRF models. Keep the original frame loop as a compatibility fallback.
isBuiltinReal = (isa(rfModel,'double') || isa(rfModel,'single') || islogical(rfModel)) && ...
                (isa(stimulusImages,'double') || isa(stimulusImages,'single') || islogical(stimulusImages)) && ...
                isreal(rfModel) && isreal(stimulusImages) && ...
                ~issparse(rfModel) && ~issparse(stimulusImages);
hasExpectedShape = ismatrix(rfModel) && ...
                   (ndims(stimulusImages) <= 3) && ...
                   (size(rfModel,1) == size(stimulusImages,1)) && ...
                   (size(rfModel,2) == size(stimulusImages,2));

% MATLAB's matrix multiplication uses optimized BLAS for dense doubles.
% Reshape is copy-on-write, so this does not require a second stimulus
% movie. Its reduction order can differ from the general fallback by a few
% floating-point bits.
canUseMatrixMultiply = isa(rfModel,'double') && ...
                       isa(stimulusImages,'double') && ...
                       isBuiltinReal && hasExpectedShape;
if canUseMatrixMultiply
  flattenedStimulus = reshape(stimulusImages,numel(rfModel),nStimFrames);
  modelResponse(1:nStimFrames) = rfModel(:).' * flattenedStimulus;
  return
end

if ~(isBuiltinReal && hasExpectedShape)
  modelResponse = projectFramesLegacy(rfModel,stimulusImages,modelResponse,nStimFrames);
  return
end

maxTemporaryElements = 2^22;
nSpatialElements = max(1,numel(rfModel));
nFramesPerChunk = max(1,floor(maxTemporaryElements/nSpatialElements));

for firstFrame = 1:nFramesPerChunk:nStimFrames
  lastFrame = min(firstFrame+nFramesPerChunk-1,nStimFrames);
  frameRange = firstFrame:lastFrame;

  weightedStimulus = rfModel .* stimulusImages(:,:,frameRange);
  chunkResponse = sum(sum(weightedStimulus,1),2);
  modelResponse(frameRange) = reshape(chunkResponse,1,[]);
end


function modelResponse = projectFramesLegacy(rfModel,stimulusImages,modelResponse,nStimFrames)

for frameNum = 1:nStimFrames
  modelResponse(frameNum) = sum(sum(rfModel .* stimulusImages(:,:,frameNum)));
end
