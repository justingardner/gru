function output = pRFFitPrepared(action,varargin)
% pRFFitPrepared
%
% Prepare scan-invariant pRF fitting state once and reuse it for voxels.
%
%   context = pRFFitPrepared('prepare',v,scanNum,<pRFFit name/value args>)
%   fit = pRFFitPrepared('fit',context,tSeries,x,y,z,<dispIndex>,<dispN>)
%
% The fit operation bypasses pRFFit's generic argument parser, viewGet
% queries, stimulus setup, parameter setup, and prefit-bank construction.
% It deliberately calls the same objective and solvers as pRFFit.

% Process workers start from MATLAB's ordinary startup path and therefore
% may initially prefer pRFv2 even when the client has installed pRFv3. A
% unique bootstrap name lets each process correct that ordering once before
% dispatching to the shared public pRFFit name.
persistent pRFv3PathIsActive
if isempty(pRFv3PathIsActive)
  pRFv3ActivatePath;
  pRFv3PathIsActive = true;
end

if nargin < 1 || ~(ischar(action) || (isstring(action) && isscalar(action)))
  error('pRFFitPrepared:InvalidAction','The first input must be ''prepare'' or ''fit''.');
end

switch lower(char(action))
  case 'prepare'
    if length(varargin) < 2
      error('pRFFitPrepared:MissingInputs','Prepare requires a view and scan number.');
    end
    output = pRFFit('__prepareContext__',varargin{:});

  case 'fit'
    if length(varargin) < 5
      error('pRFFitPrepared:MissingInputs', ...
        'Fit requires a prepared context, time series, and x/y/z coordinates.');
    end
    output = pRFFit('__fitPrepared__',varargin{:});

  otherwise
    error('pRFFitPrepared:InvalidAction','Unknown action: %s',char(action));
end
