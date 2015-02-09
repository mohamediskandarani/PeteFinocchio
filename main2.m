% First initialize the PC-information
pcstuff;

pmin=[1.0;0.5];   % minima of uncertain parameters
pmax=[3.0;1.0];   % maxima of uncertain parameters
% Read-in the realizations
% Make it up for now
fs = inline('(exp(P1*x/P2) - 1)/(exp(P1/P2)-1)',2); % 1D advection-diffusion
nx = 41;               % number of spatial points
xp = [0:1/(nx-1):1]';  % location of spatial points
fr = zeros([length(xp) nr]);
for ir = 1:nr
  iq = ind(ir,:);                 % realization ir corresponds to quad
  p = (pmax-pmin) .* (xq(iq)+1)/2 + pmin;
  fr(:,ir) = fs(xp, p(1), p(2));  % calculate realization ir
end

% Compute PC-coefficients
fh = (fr*Anisp');  % produce the PC coefficients \hat{f}_n
                   % fh(:,1) is the 0-th mode or the mean
                   % fh(:,m) is the (m-1)-th multi-dimensional mode,
                   % m=1,...nup+1

% Compute the standard deviation
fstdev = sqrt( (fh(:,2:nup+1).^2) * psinor(2:nup+1) );

% Plot mean and plus or minus a standard deviation
plot(xp, fh(:,1),'k',...
     xp, fh(:,1)+fstdev,'r--',...
     xp, fh(:,1)-fstdev,'r--');
 xlabel('x');
 legend('Mean','Mean+\sigma','Mean-\sigma');

% Validation points
nv1d = 11;
xv1d = [-1:2/(nv1d-1):1]';
nvtotal = nv1d^ndim;
pvind = setndrlz(ndim,nv1d);        % multi-dimensional index of realizations
fv = zeros([nx nvtotal]);           % validation realizations
for iv = 1:nvtotal
  iq = pvind(iv,:);                 % coordinate of validation point in 
                                    % parameter space
  p = (pmax-pmin) .* (xv1d(iq)+1)/2 + pmin;
  fv(:,iv) = fs(xp, p(1), p(2));    % validation realization
end

% Building the surrogate estimate on validation points
Pv1d = legendrepols(nord,xv1d);     % 1D polynomials
Pvnd = zeros([nup+1 1]);            % nD polynomials
for iv = 1:nvtotal
  iq = pvind(iv,:);                 % iv corresponds to points iq
  xiv = xv1d(iq);                   % multiD coordinates of validation point
  for n = 0:nup
    kord = pcnpt(n+1,:);
    Pvnd(n+1) = prod(prod(Pv1d(iq,kord+1)));
  end
  fvpc(:,iv) = fh * Pvnd;
end
