function [F,DF] = F_DF_SKT(u,para)

if exist('intval','file') && isintval(u(1))
    ipi = intval('pi');
else
    ipi = pi;
end
    
N = length(u)/2;
u1 = u(1:N);
u2 = u(N+1:2*N);

d1 = para.d1;
d2 = para.d2;
d11 = para.d11;
d12 = para.d12;
d21 = para.d21;
d22 = para.d22;
r1 = para.r1;
a1 = para.a1;
b1 = para.b1;
r2 = para.r2;
a2 = para.a2;
b2 = para.b2;

Mu1 = convomat(u1);
Mu2 = convomat(u2);
u1u1 = Mu1*u1;
u2u2 = Mu2*u2;
u1u2 = Mu1*u2;
Lap = -ipi^2*(0:N-1)'.^2;
MLap = spdiags(Lap,0,N,N); %MLap = diag(Lap);
I = eye(N);

Phi1 = d1*u1 + d11*u1u1 + d12*u1u2;
Phi2 = d2*u2 + d21*u1u2 + d22*u2u2;
R1 = r1*u1 - a1*u1u1 - b1*u1u2;
R2 = r2*u2 - b2*u1u2 - a2*u2u2;

F = [Lap.*Phi1 + R1;
     Lap.*Phi2 + R2];
 
if nargout>1 %if we also want DF 
    DPhi11 = d1*I + 2*d11*Mu1 + d12*Mu2;
    DPhi12 = d12*Mu1;
    DPhi21 = d21*Mu2;
    DPhi22 = d2*I + d21*Mu1 + 2*d22*Mu2;

    DR11 = r1*I - 2*a1*Mu1 - b1*Mu2;
    DR12 = -b1*Mu1;
    DR21 = -b2*Mu2;
    DR22 = r2*I - b2*Mu1 - 2*a2*Mu2;

    DF = [MLap*DPhi11 + DR11, MLap*DPhi12 + DR12;
          MLap*DPhi21 + DR21, MLap*DPhi22 + DR22];
end