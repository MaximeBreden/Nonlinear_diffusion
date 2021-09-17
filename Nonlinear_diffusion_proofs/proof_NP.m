function [rmin,rmax,Abar,w,sigma]=proof_NP(u,alpha,beta,gamma,g,nu,Abar,w,sigma)

%The last 3 arguments are optionals, and are computed inside this function
%if not given as input.

if exist('intval','file') && isintval(u(1))
    fprintf('\nRigorous validation with interval arithmetic\n')
    ipi = intval('pi');
else
    fprintf('\nPrevalidation without interval arithmetic\n')
    ipi = pi;
end

%% Init
N = length(u)/2;
u1 = u(1:N);
u2 = u(N+1:2*N);
e = zeros(N,1);
e(1) = 1;
I = eye(N);

xi_nu = [1, 2*nu.^(1:4*N-4)];

u_ext = zeros(2*(2*N-1),1);
if exist('intval','file') && isintval(u(1))
    u_ext = intval(u_ext);
end
u_ext(1:N) = u1;
u_ext(2*N-1+(1:N)) = u2;
[F_ext,DF_ext] = F_DF_NP(u_ext,alpha,beta,gamma,g);

if nargin<7
    Mu2 = convomat(u2);
    w = (I-Mu2)\e;
    sigma = convo(w,gamma*e+u1);
    Abar = inv(DF_ext);
    Abar = [Abar(1:N,1:N), Abar(1:N,2*N-1+(1:N));
            Abar(2*N-1+(1:N),1:N), Abar(2*N-1+(1:N),2*N-1+(1:N))];
end

A = [convomat(w,4*N-3), convomat(sigma,4*N-3) * diag([0,-1./(ipi^2*(1:(4*N-4)).^2)]);
     zeros(4*N-3,4*N-3), diag([0,-1./(ipi^2*(1:(4*N-4)).^2)])];
A(1:N,1:N) = Abar(1:N,1:N);
A(1:N,4*N-3+(1:N)) = Abar(1:N,N+1:2*N);
A(4*N-3+(1:N),1:N) = Abar(N+1:2*N,1:N);
A(4*N-3+(1:N),4*N-3+(1:N)) = Abar(N+1:2*N,N+1:2*N);

%% Y
AF = A([1:3*N-2,4*N-3+(1:3*N-2)], [1:2*N-1,4*N-3+(1:2*N-1)]) * F_ext;
Y = [xi_nu(1:3*N-2), xi_nu(1:3*N-2)] * abs(AF);
disp(['Y = ',num2str(i2f(Y))])

%% Z1
u_ext2 = zeros(2*(3*N-2),1);
if exist('intval','file') && isintval(u(1))
    u_ext2 = intval(u_ext2);
end
u_ext2(1:N) = u1;
u_ext2(3*N-2+(1:N)) = u2;
[~,DF_ext2] = F_DF_NP(u_ext2,alpha,beta,gamma,g);
ADF = A(:, [1:3*N-2, 4*N-3+(1:3*N-2)]) * DF_ext2(:, [1:2*N-1, 3*N-2+(1:2*N-1)]);
I_ext = spdiags(ones(2*N-1,1),0,4*N-3,2*N-1);
Zeros_ext = zeros(4*N-3,2*N-1);
B = [I_ext, Zeros_ext; Zeros_ext, I_ext] - ADF;
Z1_finite = max( ( [xi_nu,xi_nu] * abs(B) ) ./ [xi_nu(1:2*N-1),xi_nu(1:2*N-1)] );
disp(['Z1_finite = ',num2str(i2f(Z1_finite))])

e_ext = [e;zeros(N-1,1)];
sigma_ext = [sigma;zeros(N-1,1)];
Z1_tail = max( xi_nu(1:2*N-1) * abs(e_ext-convo(w,e-u2,2*N-1)) ...
               + 1/(ipi*N)^2 * (1 + xi_nu(1:N) * abs(sigma)) * xi_nu(1:N) * abs(alpha*e-2*beta*u1), ...
               xi_nu(1:2*N-1) * abs(sigma_ext-convo(w,gamma*e+u1,2*N-1)) );
disp(['Z1_tail = ',num2str(i2f(Z1_tail))])
      
Z1 = max(Z1_finite,Z1_tail);
disp(['Z1 = ',num2str(i2f(Z1))])
   
%% Z2
normAblock_finite = [max( ( xi_nu(1:2*N-1) * abs(A(1:2*N-1,1:N)) ) ./ xi_nu(1:N) ),...
                     max( ( xi_nu(1:2*N-1) * abs(A(1:2*N-1,4*N-3+(1:N))) ) ./ xi_nu(1:N) );
                     max( ( xi_nu(1:N) * abs(A(4*N-3+(1:N),1:N)) ) ./ xi_nu(1:N) ),...
                     max( ( xi_nu(1:N) * abs(A(4*N-3+(1:N),4*N-3+(1:N))) ) ./ xi_nu(1:N) )];
normAblock_tail = [xi_nu(1:N) * abs(w), 1/(ipi*N)^2 * xi_nu(1:N) * abs(sigma);
                   0, 1/(ipi*N)^2];
                
normAblock = max(normAblock_finite, normAblock_tail);

Z2 = max( normAblock(1,1)+normAblock(2,1), (normAblock(1,2)+normAblock(2,2))*(2*beta) );
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




