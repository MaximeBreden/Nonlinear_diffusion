clear variables
close all
clc

%parameters of the model
alpha = 1;
beta = 1;
gamma = 0.1;

%number of Fourier modes
N = 50;

%weight for the proof (the nu in \ell^1\nu)
nu = 1.1;
if exist('intval','file')
    inu = intval('1.1');
end

%forcing term g
g = zeros(N,1);
g(1:5) = [1/2;3/2;1;-1/2;3];

%domain (only affects the plots)
a = 0;
b = 1;
% discretization for the plots
nb = 500;

plot_cos(g,a,b,nb)
xlabel('$x$', 'Interpreter', 'latex')
title('g')
axis tight

%loading a precomputed solution
if gamma == 3
    load('dataNP3.mat', 'u')
    if exist('intval','file')
        igamma = intval('3');
    end
elseif gamma == 0.1
    load('dataNP01.mat', 'u')
    if exist('intval','file')
        igamma = intval('0.1');
    end
else
    disp("No precomputed solution for this value of gamma, but you might try to find one by playing around with Newton's method")
    return
end
Ndata = length(u)/2;
if Ndata<N
    u = [u(1:Ndata); zeros(N-Ndata,1); u(Ndata+1:2*Ndata); zeros(N-Ndata,1)];
else
    u = [u(1:N); u(Ndata+1:Ndata+N)];
end

%parameters for Newton's method
it_max = 20;
tol = 10^-12;

%refinement of the numerical solution using Newton's method
fprintf("\nRefinement of the numerical solution using Newton's method\n")
it = 0;
[F,DF] = F_DF_NP(u,alpha,beta,gamma,g);
err = norm(F,1)
while err>tol && it<it_max && err<10^10
    u = u -DF\F;
    [F,DF] = F_DF_NP(u,alpha,beta,gamma,g);
    err = norm(F,1)
    it = it + 1;
end

u1 = u(1:N);
u2 = u(N+1:2*N);
figure
semilogy(abs(u1),'b')
hold on
semilogy(abs(u2),'r')
title('Decay of the Fourier coefficients')

fig = 3;
plot_cos(u1,a,b,nb,'b-',fig)
hold on
plot_cos(u2,a,b,nb,'r--',fig)
xlabel('$x$', 'Interpreter', 'latex')
legend('u1','u2')
set(gca,'FontSize',15)
axis tight
ylim([0 inf])

%"prevalidation" (without interval arithmetic)
nu = 1.1;
[rmin,rmax,Abar,w,sigma]=proof_NP(u,alpha,beta,gamma,g,nu);

%rigorous proof (with interval arithmetic)
if exist('intval','file') && not(isnan(rmin))
        ialpha = intval(alpha);
        ibeta = intval(beta);
        ig = intval(g);
        inu = intval('1.1');
        iu = intval(u);
        iAbar = intval(Abar);
        iw = intval(w);
        isigma = intval(sigma);
        [irmin,irmax]=proof_NP(iu,ialpha,ibeta,igamma,ig,inu,iAbar,iw,isigma);
else
    fprintf("You need Intlab in order to run the rigorous proof\n")
end
