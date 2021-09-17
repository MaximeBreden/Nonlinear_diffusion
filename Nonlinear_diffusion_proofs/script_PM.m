clear variables
close all
clc

%parameters of the model
alpha = 1;
beta = 1;

%number of Fourier modes
N = 20;

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
set(gca,'FontSize',15)
title('g')
axis tight

%loading a precomputed solution
load('dataPM.mat', 'u')
Ndata = length(u);
if Ndata<N
    u = [u; zeros(N-Ndata,1)];
else
    u = u(1:N);
end

%parameters for Newton's method
it_max = 20;
tol = 10^-12;

%refinement of the numerical solution using Newton's method
fprintf("\nRefinement of the numerical solution using Newton's method\n")
it = 0;
[F,DF] = F_DF_PM(u,alpha,beta,g);
err = norm(F,1)
while err>tol && it<it_max
    u = u -DF\F;
    [F,DF] = F_DF_PM(u,alpha,beta,g);
    err = norm(F,1)
    it = it + 1;
end

figure
semilogy(abs(u))
title('Decay of the Fourier coefficients')

plot_cos(u,a,b,nb)
xlabel('$x$', 'Interpreter', 'latex')
title('u')
set(gca,'FontSize',15)
axis tight

%"prevalidation" (without interval arithmetic)
[rmin,rmax,Abar,w]=proof_PM(u,alpha,beta,g,nu);

%rigorous proof (with interval arithmetic)
if exist('intval','file') && not(isnan(rmin))
        ialpha = intval(alpha);
        ibeta = intval(beta);
        ig = intval(g);
        iu = intval(u);
        iAbar = intval(Abar);
        iw = intval(w);
        [irmin,irmax]=proof_PM(iu,ialpha,ibeta,ig,inu,iAbar,iw);
else
    fprintf("You need Intlab in order to run the rigorous proof\n")
end


