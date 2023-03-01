function [xl, yl, p] = gtm_ppd(t, Y, beta, X, xDim, yDim)
% Latent space posterior probability distribution for a given data point. 
%
%		This function calculates the posterior probability 
%		distribution induced in the latent space of a 
%		trained GTM model for a given data point, and
%		returns it in a format suitable for MATLAB's 2D or
%		3D graphic plotting routines, depending on the 
%		latent space dimensionality.
%
% Synopsis:	[xl, yl, p]  = gtm_ppd(t, Y, beta, X, xDim, yDim)
%		[xl, p]  = gtm_ppd(t, Y, beta, X)
%
% Arguments:	t -		a point in the data space; 1-by-D
%
%		Y -		centres of the Gaussian mixture generated by 
%				the GTM  in the data space, Y = FI*W; K-by-D
%
%		beta -		variance of Gaussian mixture; scalar
%
%		X -		latent sample
%
%		xDim, yDim -	number of points along the 2 dimensions
%				of the latent space meshgrid sample
%
% Return:	xl, yl -	latent sample; if the latent space is 
%				2D, xl and yl are mesh matrices; if it is
%				1D, xl is identical to X
%
%		p -		posterior distribution over latent space
%				given data point t; if the latent space is
%				2D, p is a mesh matrix, otherwise it is
%				a vector of same length as xl
%				
% Notes:	If the latent sample X is 2 dimensional, it is assumed
%		to have been constructed from a mesh-grid, e.g. as if
%		generated by gtm_stp2
%
% See also:	gtm_pmn, gtm_pmd, gtm_stp2
%


[N, D] = size(t);
[K, L] = size(X);

[DIST, minDist, maxDist] = gtm_dist(t, Y, 1);
[err, R] = gtm_resp(DIST, minDist, maxDist, beta, D, 1);

if (L==1 & nargout==2)
  xl = X;
  yl = R;
elseif (L==2 & nargout==3)
  [xl, yl, p] = gtm_r2m(X(:,1), X(:,2), R, xDim, yDim);
else
  error('Mismatch between input and output arguments');
end
