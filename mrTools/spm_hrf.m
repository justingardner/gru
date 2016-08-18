function [hrf,p] = spm_hrf(RT,P)
% returns a hemodynamic response function
% FORMAT [hrf,p] = spm_hrf(RT,[p]);
% RT   - scan repeat time
% p    - parameters of the response function (two gamma functions)
%
%							defaults
%							(seconds)
%	p(1) - delay of response (relative to onset)	   6
%	p(2) - delay of undershoot (relative to onset)    16
%	p(3) - dispersion of response			   1
%	p(4) - dispersion of undershoot			   1
%	p(5) - ratio of response to undershoot		   6
%	p(6) - onset (seconds)				   0
%	p(7) - length of kernel (seconds)		  32
%
% hrf  - hemodynamic response function
% p    - parameters of the response function
%_______________________________________________________________________
% @(#)spm_hrf.m	2.8 Karl Friston 02/07/31

% global parameter
%-----------------------------------------------------------------------
global defaults
if ~isempty(defaults),
	fMRI_T = defaults.stats.fmri.t;
else,
	fMRI_T = 16;
end;

% default parameters
%-----------------------------------------------------------------------
p     = [6 16 1 1 6 0 32];
if nargin > 1
      p(1:length(P)) = P;
end

% modelled hemodynamic response function - {mixture of Gammas}
%-----------------------------------------------------------------------
dt    = RT/fMRI_T;
u     = [0:(p(7)/dt)] - p(6)/dt;
hrf   = spm_Gpdf(u,p(1)/p(3),dt/p(3)) - spm_Gpdf(u,p(2)/p(4),dt/p(4))/p(5);
hrf   = hrf([0:(p(7)/RT)]*fMRI_T + 1);
hrf   = hrf'/sum(hrf);

function f = spm_Gpdf(x,h,l)
% Probability Density Function (PDF) of Gamma distribution
% FORMAT f = spm_Gpdf(x,h,l)
%
% x - Gamma-variate   (Gamma has range [0,Inf) )
% h - Shape parameter (h>0)
% l - Scale parameter (l>0)
% f - PDF of Gamma-distribution with shape & scale parameters h & l
%__________________________________________________________________________
%
% spm_Gpdf implements the Probability Density Function of the Gamma
% distribution.
%
% Definition:
%--------------------------------------------------------------------------
% The PDF of the Gamma distribution with shape parameter h and scale l
% is defined for h>0 & l>0 and for x in [0,Inf) by: (See Evans et al.,
% Ch18, but note that this reference uses the alternative
% parameterisation of the Gamma with scale parameter c=1/l)
%
%           l^h * x^(h-1) exp(-lx)
%    f(x) = ----------------------
%                   gamma(h)
%
% Variate relationships: (Evans et al., Ch18 & Ch8)
%--------------------------------------------------------------------------
% For natural (strictly +ve integer) shape h this is an Erlang distribution.
%
% The Standard Gamma distribution has a single parameter, the shape h.
% The scale taken as l=1.
%
% The Chi-squared distribution with v degrees of freedom is equivalent
% to the Gamma distribution with scale parameter 1/2 and shape parameter v/2.
%
% Algorithm:
%--------------------------------------------------------------------------
% Direct computation using logs to avoid roundoff errors.
%
% References:
%--------------------------------------------------------------------------
% Evans M, Hastings N, Peacock B (1993)
%       "Statistical Distributions"
%        2nd Ed. Wiley, New York
%
% Abramowitz M, Stegun IA, (1964)
%       "Handbook of Mathematical Functions"
%        US Government Printing Office
%
% Press WH, Teukolsky SA, Vetterling AT, Flannery BP (1992)
%       "Numerical Recipes in C"
%        Cambridge
%__________________________________________________________________________
% Copyright (C) 1993-2011 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_Gpdf.m 4182 2011-02-01 12:29:09Z guillaume $


%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<3, error('Insufficient arguments'), end

ad = [ndims(x);ndims(h);ndims(l)];
rd = max(ad);
as = [[size(x),ones(1,rd-ad(1))];...
      [size(h),ones(1,rd-ad(2))];...
      [size(l),ones(1,rd-ad(3))]];
rs = max(as);
xa = prod(as,2)>1;
if sum(xa)>1 && any(any(diff(as(xa,:)),1))
    error('non-scalar args must match in size');
end

%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
f = zeros(rs);

%-Only defined for strictly positive h & l. Return NaN if undefined.
md = ( ones(size(x))  &  h>0  &  l>0 );
if any(~md(:))
    f(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Degenerate cases at x==0: h<1 => f=Inf; h==1 => f=l; h>1 => f=0
ml = ( md  &  x==0  &  h<1 );
f(ml) = Inf;
ml = ( md  &  x==0  &  h==1 ); if xa(3), mll=ml; else mll=1; end
f(ml) = l(mll);

%-Compute where defined and x>0
Q  = find( md  &  x>0 );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qh=Q; else Qh=1; end
if xa(3), Ql=Q; else Ql=1; end

%-Compute
f(Q) = exp( (h(Qh)-1).*log(x(Qx)) +h(Qh).*log(l(Ql)) - l(Ql).*x(Qx)...
-gammaln(h(Qh)) );