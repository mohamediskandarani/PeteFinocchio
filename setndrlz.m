function ind = setndrlz(ndim,nq1d);
%
% Returns the matrix mapping the nr realizations into
% a multi-index coordinates
% usage is ind = setndrlz(ndim,nq1d);
%   nq1d is the number of 1D quadrature points
%   ndim is the number of uncertain parameters
%   ind(r, d) is the coordinate index of realization r along dimension d
%   ind is of size nr by ndim
% We need to map a collocation point in ndim-space
% to a single index. For tensor-based realization this is fairly simple.
% Let the multi-dimensional index of the cubature be
% (m_1, m_2, m_3,..., m_D), 
% m_d is the m-th quadrature point along the d-th direction.
%  A simple mapping to a single index notation is
% r = (m_D-1) * N^{d-1} + (m_{D-1}-1) * N^{d-2} +...+ (m_3-1)*N^2 + (m_2-1)*N +m_1
%
%  The individual coordinate index range is 1<= m_d <= N
%
nr = nq1d^ndim;    % total number of realizations
ind = zeros([nr ndim]);
for r = 1:nr
  count = r;
  m = zeros([ndim 1]);
  for d = ndim:-1:1
    m(d) = floor((count-1)/nq1d^(d-1)) + 1;
    count = count - (m(d)-1) * nq1d^(d-1);
    ind(r,d) = m(d);
  end
end
