function nWorkers = pRFResolveNumWorkers(requestedWorkers)
% pRFResolveNumWorkers
%
% nWorkers = pRFResolveNumWorkers
%   Return the number of workers available in the local profile without
%   starting or changing a parallel pool.
%
% nWorkers = pRFResolveNumWorkers(requestedWorkers)
%   Resolve a positive integer worker request noninteractively. A request of
%   one is true serial execution and deliberately does not inspect, start, or
%   close an existing pool. Requests greater than one reuse an exact-size
%   pool or replace a mismatched pool. Pool failures fall back to serial.

if nargin < 1 || isempty(requestedWorkers)
  nWorkers = getAvailableLocalWorkers;
  return
end

if ~(isnumeric(requestedWorkers) && isreal(requestedWorkers) && ...
    isscalar(requestedWorkers) && isfinite(requestedWorkers) && ...
    requestedWorkers >= 1 && requestedWorkers == fix(requestedWorkers))
  error('pRFResolveNumWorkers:InvalidNumWorkers', ...
    'numWorkers must be a positive integer.');
end
requestedWorkers = double(requestedWorkers);

% Serial execution must not call gcp, parpool, or mlrNumWorkers. In
% particular, leave a caller-owned pool open but unused.
if requestedWorkers == 1
  nWorkers = 1;
  return
end

nWorkers = 1;
try
  if exist('gcp','file') ~= 2 || exist('parpool','file') ~= 2
    warning('pRFResolveNumWorkers:ParallelUnavailable', ...
      ['Parallel Computing Toolbox is unavailable. Continuing with ' ...
       'serial pRF execution.']);
    return
  end

  existingPool = gcp('nocreate');
  if ~isempty(existingPool)
    if existingPool.NumWorkers == requestedWorkers
      [poolIsReady,activationError] = activatePRFv3Workers(existingPool);
      if poolIsReady
        nWorkers = requestedWorkers;
      else
        warning('pRFResolveNumWorkers:WorkerPathActivationFailed', ...
          ['Could not activate pRFv3 on the existing parallel workers: %s. ' ...
           'Continuing with serial pRF execution.'],activationError);
      end
      return
    end
    delete(existingPool);
  end

  requestedPool = parpool('local',requestedWorkers);
  if isempty(requestedPool) || requestedPool.NumWorkers ~= requestedWorkers
    if ~isempty(requestedPool)
      delete(requestedPool);
    end
    warning('pRFResolveNumWorkers:WrongPoolSize', ...
      ['Could not create the requested %i-worker pool. Continuing with ' ...
       'serial pRF execution.'],requestedWorkers);
    return
  end
  [poolIsReady,activationError] = activatePRFv3Workers(requestedPool);
  if ~poolIsReady
    delete(requestedPool);
    warning('pRFResolveNumWorkers:WorkerPathActivationFailed', ...
      ['Could not activate pRFv3 on the requested parallel workers: %s. ' ...
       'Continuing with serial pRF execution.'],activationError);
    return
  end
  nWorkers = requestedPool.NumWorkers;
catch poolException
  warning('pRFResolveNumWorkers:PoolStartupFailed', ...
    ['Could not start the requested %i-worker pool: %s. Continuing with ' ...
     'serial pRF execution.'],requestedWorkers,poolException.message);
  nWorkers = 1;
end


function nWorkers = getAvailableLocalWorkers

nWorkers = 1;
try
  if exist('parcluster','file') == 2 && exist('parpool','file') == 2
    localCluster = parcluster('local');
    nWorkers = localCluster.NumWorkers;
  end
catch
  % Parallel Computing Toolbox is optional; retain the serial default.
end

if ~(isnumeric(nWorkers) && isscalar(nWorkers) && isfinite(nWorkers))
  nWorkers = 1;
end
nWorkers = max(1,floor(double(nWorkers)));


function [poolIsReady,errorMessage] = activatePRFv3Workers(pool)

poolIsReady = false;
errorMessage = '';
try
  % Reasserting the path after a pool exists also updates process workers.
  % The explicit worker call makes that synchronization deterministic and
  % validates all shared fitter/model/helper names before any prefit PARFOR.
  pluginRoot = pRFv3ActivatePath;
  if exist('parfevalOnAll','file') ~= 2
    error('pRFResolveNumWorkers:MissingParfevalOnAll', ...
      'PARFEVALONALL is unavailable.');
  end
  activationFuture = parfevalOnAll(pool,@pRFv3ActivatePath,0,pluginRoot);
  fetchOutputs(activationFuture);
  poolIsReady = true;
catch activationException
  errorMessage = activationException.message;
end
