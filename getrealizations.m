function fr = getrealizations(xp, p)
%
% Get realizations for the PC experiments
% Usage is
%   fr = getrealizations(xp, p)
% where p is vector of uncertain coefficients
%       xp is the (vector) of points in physical space
%       fr is the realization at all points in physical space
%
fs = inline('(exp(P1*x/P2) - 1)/(exp(P1/P2)-1)',2); % 1D advection-diffusion
fr = fs(xp, p(1), p(2));
