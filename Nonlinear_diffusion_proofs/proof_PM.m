function [rmin,rmax,Abar,w]=proof_PM(u,alpha,beta,g,nu,Abar,w)

%The last 2 arguments are optionals, and are computed inside this function
%if not given as input.

if exist('intval','file') && isintval(u(1))
    fprintf('\nRigorous validation with interval arithmetic\n')
    ipi = intval('pi');
else
    fprintf('\nPrevalidation without interval arithmetic\n')
    ipi = pi;
end

%% Init
N = length(u);
e = zeros(N,1);
e(1) = 1;

xi_nu = [1, 2*nu.^(1:4*N-4)];

u_ext = [u;zeros(N-1,1)];
[F_ext,DF_ext] = F_DF_PM(u_ext,alpha,beta,g);

if nargin<6
    Mu = convomat(u);
    w = 1/2*(Mu\e);
    Abar = inv(DF_ext);
    Abar = Abar(1:N,1:N);
end

A = convomat(w,4*N-3);
A = A * diag([0,-1./(ipi^2*(1:(4*N-4)).^2)]);
A(1:N,1:N) = Abar;

%% Y
AF = A(1:3*N-2,1:2*N-1)*F_ext;
Y = xi_nu(1:3*N-2)*abs(AF);
disp(['Y = ',num2str(i2f(Y))])

%% Z1
u_ext2 = [u;zeros(2*N-2,1)];
[~,DF_ext2] = F_DF_PM(u_ext2,alpha,beta,g);
ADF = A(:,1:3*N-2) * DF_ext2(:,1:2*N-1);
B = spdiags(ones(2*N-1,1),0,4*N-3,2*N-1) - ADF;
Z1_finite = max( ( xi_nu * abs(B) ) ./ xi_nu(1:2*N-1) );
disp(['Z1_finite = ',num2str(i2f(Z1_finite))])

e_ext = [e;zeros(N-1,1)];
Z1_tail = xi_nu(1:2*N-1) * abs(e_ext-convo(w,2*u,2*N-1)) ...
          + 1/(pi*N)^2 * xi_nu(1:2*N-1) * abs(convo(w,alpha*e-2*beta*u,2*N-1));
disp(['Z1_tail = ',num2str(i2f(Z1_tail))])
      
Z1 = max(Z1_finite,Z1_tail);
disp(['Z1 = ',num2str(i2f(Z1))])
   
%% Z2
normA_finite = max( ( xi_nu(1:2*N-1) * abs(A(1:2*N-1,1:N)) ) ./ xi_nu(1:N) );
normA_tail = 1/(pi*N)^2  * xi_nu(1:N) * abs(w);
normA = max(normA_finite,normA_tail);

normALap_finite = max( ( xi_nu(1:2*N-1) * abs(A(1:2*N-1,1:N)*diag(ipi^2*(0:N-1).^2)) ) ./ xi_nu(1:N) );
normALap_tail = xi_nu(1:N) * abs(w);
normALap = max(normALap_finite,normALap_tail);

Z2 = normALap*2 + normA*2*beta;
disp(['Z2 = ',num2str(i2f(Z2))])


%% Cheking that we have a contraction

Delta = (1-Z1)^2-2*Y*Z2;
if i2f(Z1)<1 && i2f(Delta,'inf')>0
    rmin = i2f( (1-Z1-sqrt(Delta)) / Z2 );
    rmax = i2f( (1-Z1) / Z2, 'inf' );
    disp(['[r_min,r_max] = ',num2str(rmin),',',num2str(rmax)])
    fprintf('\n')
else
    rmin = NaN;
    rmax = NaN;
    fprintf('The validation failed :(\n')
end
