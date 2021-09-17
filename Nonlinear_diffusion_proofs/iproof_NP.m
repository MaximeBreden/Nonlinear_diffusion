function [rmin,rmax]=iproof_NP(X,alpha,beta,gamma,g,nu)

fprintf('\nRigorous validation with interval arithmetic\n')

%% Init
N = length(X)/2;
u = X(1:N);
v = X(N+1:2*N);
e = zeros(N,1);
e(1) = 1;
I = eye(N);

Mv = convomat(v);
w = (I-Mv)\e;
sigma = convo(w,gamma*e+u);
xi_nu = [1, 2*nu.^(1:4*N-4)];

X_ext = zeros(2*(2*N-1),1);
X_ext(1:N) = u;
X_ext(2*N-1+(1:N)) = v;
[F_ext,DF_ext] = F_DF_NP(X_ext,alpha,beta,gamma,g);
Abar = inv(DF_ext);
Abar = [Abar(1:N,1:N), Abar(1:N,2*N-1+(1:N));
        Abar(2*N-1+(1:N),1:N), Abar(2*N-1+(1:N),2*N-1+(1:N))];

A = [convomat(w,4*N-3), convomat(sigma,4*N-3) * diag([0,-1./(pi^2*(1:(4*N-4)).^2)]);
     zeros(4*N-3,4*N-3), diag([0,-1./(pi^2*(1:(4*N-4)).^2)])];
A(1:N,1:N) = Abar(1:N,1:N);
A(1:N,4*N-3+(1:N)) = Abar(1:N,N+1:2*N);
A(4*N-3+(1:N),1:N) = Abar(N+1:2*N,1:N);
A(4*N-3+(1:N),4*N-3+(1:N)) = Abar(N+1:2*N,N+1:2*N);

%% Y
AF = A([1:3*N-2,4*N-3+(1:3*N-2)], [1:2*N-1,4*N-3+(1:2*N-1)]) * F_ext;
Y = [xi_nu(1:3*N-2), xi_nu(1:3*N-2)] * abs(AF);
disp(['Y = ',num2str(Y)])

%% Z1
X_ext2 = zeros(2*(3*N-2),1);
X_ext2(1:N) = u;
X_ext2(3*N-2+(1:N)) = v;
[~,DF_ext2] = F_DF_NP(X_ext2,alpha,beta,gamma,g);
ADF = A(:, [1:3*N-2, 4*N-3+(1:3*N-2)]) * DF_ext2(:, [1:2*N-1, 3*N-2+(1:2*N-1)]);
I_ext = spdiags(ones(2*N-1,1),0,4*N-3,2*N-1);
Zeros_ext = zeros(4*N-3,2*N-1);
B = [I_ext, Zeros_ext; Zeros_ext, I_ext] - ADF;
Z1_finite = max( ( [xi_nu,xi_nu] * abs(B) ) ./ [xi_nu(1:2*N-1),xi_nu(1:2*N-1)] );
disp(['Z1_finite = ',num2str(Z1_finite)])

e_ext = [e;zeros(N-1,1)];
sigma_ext = [sigma;zeros(N-1,1)];
Z1_tail = max( xi_nu(1:2*N-1) * abs(e_ext-convo(w,e-v,2*N-1)) ...
               + 1/(pi*N)^2 * (1 + xi_nu(1:N) * abs(sigma)) * xi_nu(1:N) * abs(alpha*e-2*beta*u), ...
               xi_nu(1:2*N-1) * abs(sigma_ext-convo(w,gamma*e+u,2*N-1)) );
disp(['Z1_tail = ',num2str(Z1_tail)])
      
Z1 = max(Z1_finite,Z1_tail);
disp(['Z1 = ',num2str(Z1)])
   
%% Z2
normAblock_finite = [max( ( xi_nu(1:2*N-1) * abs(A(1:2*N-1,1:N)) ) ./ xi_nu(1:N) ),...
                     max( ( xi_nu(1:2*N-1) * abs(A(1:2*N-1,4*N-3+(1:N))) ) ./ xi_nu(1:N) );
                     max( ( xi_nu(1:N) * abs(A(4*N-3+(1:N),1:N)) ) ./ xi_nu(1:N) ),...
                     max( ( xi_nu(1:N) * abs(A(4*N-3+(1:N),4*N-3+(1:N))) ) ./ xi_nu(1:N) )];
normAblock_tail = [xi_nu(1:N) * abs(w), 1/(pi*N)^2 * xi_nu(1:N) * abs(sigma);
                   0, 1/(pi*N)^2];
                
normAblock = max(normAblock_finite, normAblock_tail);

Z2 = max( normAblock(1,1)+normAblock(2,1), (normAblock(1,2)+normAblock(2,2))*(2*beta) );
disp(['Z2 = ',num2str(Z2)])


%% Cheking that we have a contraction

Delta = (1-Z1)^2-2*Y*Z2;
if Z1<1 && Delta>0
    rmin = (1-Z1-sqrt(Delta)) / Z2;
    rmax = (1-Z1) / Z2;
    disp(['[r_min,r_max] = ',num2str(rmin),',',num2str(rmax)])
    fprintf('\n')
else
    rmin = NaN;
    rmax = NaN;
    fprintf('The validation failed :(\n')
end




