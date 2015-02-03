fs = inline('sin(pi*x)+0.3*cos(2*pi*x)'); % make up a function
%fs = inline('tanh(pi/2*cos(pi*x)).*(x-2)');
%fs = inline('max(0,1-4*(x-1/4).^2)'); % make up a function
xp =[-1:2/400:1]';        % some plotting and validation points
fp = fs(xp);              % values of the reference function at validation points
plot(xp,fp);

nordTests = 2:15;
for itest = 1:length(nordTests)
  nord = nordTests(itest);  % degree of 1D expansion
  nup  = nord + 1;          % number of basis functions
  nr   = nord + 1;          % number of quadrature points
                            % or (polynomial degree+1)

% 1D setup
  ab=r_jacobi(nord+1);      % recurrence coefficients for orthogonal polynomials (legendre)
  xw=gauss(nord+1,ab);      % corresponding Gauss roots and Gauss weights
  xw(:,2)  = 0.5*xw(:,2);   % scale quadrature weights to have sum=1

   % Build the 1D vectors and matrices
  P = legendrepols(nord,xw(:,1));   % sample basis function on quadrature points
                                    % P(i,n) is the basis function of degree (n-1) at quadrature point i
  psinor = (P.^2)'*xw(:,2);         % norm of 1D basis functions

  % Build the projection matrix
  Prj = P;
  for n = 1:nup
     Prj(:,n) = P(:,n).*xw(:,2);
  end
  for i = 1:nr
     Prj(i,:) = Prj(i,:)./psinor(:)';
  end
  
  fq = fs(xw(:,1));                 %  sample function at the quadrature points

  fhproj = Prj'*fq;                 %  coefficients by projection

% validation points
  Pp = legendrepols(nord,xp);       % sample basis  at validation points
  faproj = Pp*fhproj;               % approximation at validation points

  erms(itest) = sqrt(sum((faproj-fp).^2)/length(fp)); % approximation error on validation points
  [nord erms(itest)]


  plot(xp,fp,'k',xp,faproj,'b', xw(:,1),fq,'rx');
  title(strcat('nord=',int2str(nord)));
  legend('Exact','Projection','Samples');
  pause
end
semilogy(nordTests-1, erms,'r');
xlabel('Polynomial degree $P$','Interpret','latex');
ylabel('$\|\epsilon\|_2$','Interpret','latex');
legend('Projection');
title('Convergence Rates');
print('-dpdf', 'ConvergenceCurve.pdf');
