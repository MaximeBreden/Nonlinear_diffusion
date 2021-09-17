function [F,DF] = F_DF_PM(u,alpha,beta,g)

if exist('intval','file') && isintval(u(1))
    ipi = intval('pi');
else
    ipi = pi;
end

N = length(u);
Mu = convomat(u);
u2 = Mu*u;
Lap = -ipi^2*(0:N-1)'.^2;

if length(g)<N
    g=[g;zeros(N-length(g),1)];
end

F = Lap.*u2 + alpha*u - beta*u2 +g;

DF = 2*diag(Lap)*Mu + alpha*eye(N) -2*beta*Mu;
