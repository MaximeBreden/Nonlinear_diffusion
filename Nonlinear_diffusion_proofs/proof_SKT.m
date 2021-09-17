function [rmin,rmax,Abar,w11,w12,w21,w22]=proof_SKT(u,para,nu,Abar,w11,w12,w21,w22)

%The last 5 arguments are optionals, and are computed inside this function
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
xi_nu = [1, 2*nu.^(1:4*N-4)];

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


DPhi11 = d1*e + 2*d11*u1 + d12*u2;
DPhi12 = d12*u1;
DPhi21 = d21*u2;
DPhi22 = d2*e + d21*u1 + 2*d22*u2;

u_ext = zeros(2*(2*N-1),1);
if exist('intval','file') && isintval(u(1))
    u_ext = intval(u_ext);
end
u_ext(1:N) = u1;
u_ext(2*N-1+(1:N)) = u2;

if nargin<4
    MDPhi = [convomat(DPhi11), convomat(DPhi12);
             convomat(DPhi21), convomat(DPhi22)];

    w_dot1 = MDPhi \ [e;zeros(N,1)];
    w11 = w_dot1(1:N);
    w21 = w_dot1(N+1:2*N);

    w_dot2 = MDPhi \ [zeros(N,1);e];
    w12 = w_dot2(1:N);
    w22 = w_dot2(N+1:2*N);
   
    [F_ext,DF_ext] = F_DF_SKT(u_ext,para);
    Abar = inv(DF_ext);
    Abar = [Abar(1:N,1:N), Abar(1:N,2*N-1+(1:N));
            Abar(2*N-1+(1:N),1:N), Abar(2*N-1+(1:N),2*N-1+(1:N))];
else
    F_ext = F_DF_SKT(u_ext,para);
end

%InvLap = diag([0,-1./(ipi^2*(1:(4*N-4)).^2)]);
InvLap = spdiags([0,-1./(ipi^2*(1:(4*N-4)).^2)]',0,4*N-3,4*N-3);
A = [convomat(w11,4*N-3)*InvLap, convomat(w12,4*N-3)*InvLap;
     convomat(w21,4*N-3)*InvLap, convomat(w22,4*N-3)*InvLap];

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
[~,DF_ext2] = F_DF_SKT(u_ext2,para);
ADF = A(:, [1:3*N-2, 4*N-3+(1:3*N-2)]) * DF_ext2(:, [1:2*N-1, 3*N-2+(1:2*N-1)]);
I_ext = spdiags(ones(2*N-1,1),0,4*N-3,2*N-1);
Zeros_ext = zeros(4*N-3,2*N-1);
B = [I_ext, Zeros_ext; Zeros_ext, I_ext] - ADF;
Z1_finite = max( ( [xi_nu,xi_nu] * abs(B) ) ./ [xi_nu(1:2*N-1),xi_nu(1:2*N-1)] );
disp(['Z1_finite = ',num2str(i2f(Z1_finite))])

norm_w = [ xi_nu(1:N)*abs(w11), xi_nu(1:N)*abs(w12);
           xi_nu(1:N)*abs(w21), xi_nu(1:N)*abs(w22); ];
       
DR11 = r1*e - 2*a1*u1 - b1*u2;
DR12 = -b1*u1;
DR21 = -b2*u2;
DR22 = r2*e - b2*u1 - 2*a2*u2;

norm_DR = [ xi_nu(1:N)*abs(DR11), xi_nu(1:N)*abs(DR12);
           xi_nu(1:N)*abs(DR21), xi_nu(1:N)*abs(DR22); ];
       
norm_w_DR = norm_w * norm_DR;

e_ext = [e;zeros(N-1,1)];
Z1_tail = max( xi_nu(1:2*N-1) * abs(e_ext-convo(w11,DPhi11,2*N-1)-convo(w12,DPhi21,2*N-1)) ...
               + xi_nu(1:2*N-1) * abs(convo(w21,DPhi11,2*N-1)+convo(w22,DPhi21,2*N-1)) ...
               + 1/(ipi*N)^2 * ( norm_w_DR(1,1) + norm_w_DR(2,1) ), ...
               xi_nu(1:2*N-1) * abs(convo(w11,DPhi12,2*N-1)+convo(w12,DPhi22,2*N-1)) ...
               + xi_nu(1:2*N-1) * abs(e_ext-convo(w21,DPhi12,2*N-1)-convo(w22,DPhi22,2*N-1)) ...
               + 1/(ipi*N)^2 * ( norm_w_DR(1,2) + norm_w_DR(2,2) ) );
disp(['Z1_tail = ',num2str(i2f(Z1_tail))])
           
Z1 = max(Z1_finite,Z1_tail);
disp(['Z1 = ',num2str(i2f(Z1))])
   
%% Z2
normAblock_finite = [max( ( xi_nu(1:2*N-1) * abs(A(1:2*N-1,1:N)) ) ./ xi_nu(1:N) ),...
                     max( ( xi_nu(1:2*N-1) * abs(A(1:2*N-1,4*N-3+(1:N))) ) ./ xi_nu(1:N) );
                     max( ( xi_nu(1:2*N-1) * abs(A(4*N-3+(1:2*N-1),1:N)) ) ./ xi_nu(1:N) ),...
                     max( ( xi_nu(1:2*N-1) * abs(A(4*N-3+(1:2*N-1),4*N-3+(1:N))) ) ./ xi_nu(1:N) )];
normAblock_tail = 1/(ipi*N)^2 * norm_w;
normAblock = max(normAblock_finite, normAblock_tail);

Lap = -ipi^2*diag((0:N-1).^2);
normALapblock_finite = [max( ( xi_nu(1:2*N-1) * abs(A(1:2*N-1,1:N)*Lap) ) ./ xi_nu(1:N) ),...
                     max( ( xi_nu(1:2*N-1) * abs(A(1:2*N-1,4*N-3+(1:N))*Lap) ) ./ xi_nu(1:N) );
                     max( ( xi_nu(1:2*N-1) * abs(A(4*N-3+(1:2*N-1),1:N)*Lap) ) ./ xi_nu(1:N) ),...
                     max( ( xi_nu(1:2*N-1) * abs(A(4*N-3+(1:2*N-1),4*N-3+(1:N))*Lap) ) ./ xi_nu(1:N) )];
normALapblock_tail = norm_w;
normALapblock = max(normALapblock_finite, normALapblock_tail);

Z2_11 = 2*d11*(normALapblock(1,1)+normALapblock(2,1)) + 2*a1*(normAblock(1,1)+normAblock(2,1));
Z2_22 = 2*d22*(normALapblock(1,2)+normALapblock(2,2)) + 2*a2*(normAblock(1,2)+normAblock(2,2));
Z2_12 = d12*(normALapblock(1,1)+normALapblock(2,1)) + b1*(normAblock(1,1)+normAblock(2,1)) ...
        + d21*(normALapblock(1,2)+normALapblock(2,2)) + b2*(normAblock(1,2)+normAblock(2,2));
Z2 = max([Z2_11, Z2_22, Z2_12]);
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

