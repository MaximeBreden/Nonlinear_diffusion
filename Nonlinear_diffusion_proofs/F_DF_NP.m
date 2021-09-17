function [F,DF] = F_DF_NP(u,alpha,beta,gamma,g)

if exist('intval','file') && isintval(u(1))
    ipi = intval('pi');
else
    ipi = pi;
end

N = length(u)/2;
u1 = u(1:N);
u2 = u(N+1:2*N);

Mu1 = convomat(u1);
u1u1 = Mu1*u1;
u1u2 = Mu1*u2;
Mu2 = convomat(u2);
Lap = -ipi^2*(0:N-1)'.^2;
I = eye(N);

if length(g)<N
    g=[g;zeros(N-length(g),1)];
end

F = [u1 - gamma*u2 - u1u2;
     Lap.*u2 + alpha*u1 - beta*u1u1 + g];
 
DF = [I-Mu2, -(gamma*I+Mu1);
      alpha*I-2*beta*Mu1, diag(Lap)];