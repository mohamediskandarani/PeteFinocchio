params=load('pcdims.dat');
kind = params(1); % type of quadrature
nord = params(2); % order of PC expansion in 1D, 0:nord
ndim = params(3); % number of stochastic dimensions, 1:ndim
nup  = params(4); % number of multiD polynomials, 0:nup
nr   = params(5); % number of realizations

xw   = load('pcquad.dat'); % pc quadrature information
xq   = xw(:,1);            % quadrature points
wq   = xw(:,2);            % quadrature weights
nc   = length(xq);         % number of quadrature points

ind  = load('pcind.dat');    % multi-index of realization
ind  = reshape(ind,nr,ndim); % ind(r,d) = index of quad point along dimension d

pcnpt = load('pcnpt.dat');         % multi-index of polynomial basis
pcnpt = reshape(pcnpt,ndim,nup+1)'; % pcnpt(p,d) = degree of p-th
                                   % multi-dimensional basis function
                                   % along dimension d

Anisp  = load('pcanisp.dat');      % projection matrix of realization
Anisp  = reshape(Anisp,nup+1,nr);

psinor = load('pcpsinor.dat');     % weighed 2-norm of multi-dimensional basis
