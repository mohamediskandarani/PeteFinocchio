function psinormnd = normnd(multiindex,psinorm1d)
%
% get the 2-norm of the multi-dimensional polynomials
% by tensor product construction of 1D basis norms.
% Usage is
%  psind = normnd(multiindex,psin1d)
%     multiindex is the (nup+1) by ndim multi-index set
%     psinorm1d = the 1D norm of the (isotropic) basis
%     multiindex(k,idim) = degree of the k-th basis function along
%                          dimension i
%
psinormnd = ones([nup+1 1]);
for m = 1:nup+1
  psinormnd(:) = psinormnd(:).*psinorm1d(multiindex(:,m));
end
