clear variables
close all
clc

%number of Fourier modes
N = 100;

%weight for the proof (the nu in \ell^1\nu)
nu = 1.01;
if exist('intval','file')
    inu = intval('1.01');
end

%ind is used to select the parameters (and the associated solutions we are 
%interested in. ind should be taken in {1,..,4}. The four different cases 
%correspond to the four cases studied in the paper. See just below for a 
%brief description of each case. Depending on the case you choose, you may 
%have to change N and/or nu for the proof to be successful. Appropriate
%values for each cases are given in the paper.
ind = 4;

%parameters of the model and precomputed solutions
switch ind
    
    case 1 %"usual" triangular + weak competition case, from Mimura & al.
        d1=0.005;d2=0.005;d12=3;d21=0;d11=0;d22=0;r1=5;r2=2;a1=3;a2=3;b1=1;b2=1;
        load('dataSKT1.mat', 'Tab_u')
        
    case 2 %weak competition but with d12=d21 large, from Breden & al.
        d1=0.005;d2=0.005;d12=100;d21=100;d11=0;d22=0;r1=15/2;r2=16/7;a1=4;a2=2;b1=6;b2=1;
        load('dataSKT2.mat', 'Tab_u')
        
    case 3 %triangular but outside of the usual weak/strong range, from Breden & al.
        d1=0.05;d2=0.05;d12=3;d21=0;d11=0;d22=0;r1=15;r2=5;a1=1;a2=1;b1=0.5;b2=3;
        load('dataSKT3.mat', 'Tab_u')
        
    case 4 %With self-diffusion, and negative standard diffusion d1 and d2, inspired from Breden & al.
        d1=-0.007;d2=-0.007;d12=3;d21=0.002;d11=0.05;d22=0.05;r1=5;r2=2;a1=3;a2=3;b1=1;b2=1;
        load('dataSKT4.mat', 'Tab_u')
        
    otherwise
        disp("Invalid value of 'ind'")
        return
        
end

para.d1 = d1;
para.d2 = d2;
para.d11 = d11;
para.d12 = d12;
para.d21 = d21;
para.d22 = d22;
para.r1 = r1;
para.a1 = a1;
para.b1 = b1;
para.r2 = r2;
para.a2 = a2;
para.b2 = b2;

if exist('intval','file')
    ipara = para;
    %We have to take care of the parameter values which are not exactly
    %representable with floatting point numbers.
    switch ind
        case 1
            ipara.d1=intval('0.005');ipara.d2=intval('0.005');
        case 2
            ipara.d1=intval('0.005');ipara.d2=intval('0.005');ipara.r2=16/intval(7);
        case 3 
            ipara.d1=intval('0.05');ipara.d2=intval('0.05');
        case 4
            ipara.d1=intval('-0.007');ipara.d2=intval('-0.007');ipara.d21=intval('0.002');ipara.d11=intval('0.05');ipara.d22=intval('0.05');
    end
end


Ndata = size(Tab_u,1)/2;
L = size(Tab_u,2);
if Ndata<N
    Tab_u = [Tab_u(1:Ndata,:); zeros(N-Ndata,L); Tab_u(Ndata+1:2*Ndata,:); zeros(N-Ndata,L)];
else
    Tab_u = [Tab_u(1:N,:); Tab_u(Ndata+1:Ndata+N,:)];
end

%parameters for Newton's method
it_max = 10;
tol = 10^-12;

for l=1:L
    u = Tab_u(:,l);
    fprintf("\nSolution number %d:\n",l)
    
    %refinement of the numerical solution using Newton's method
    fprintf("\nRefinement of the numerical solution using Newton's method\n")
    it = 0;
    [F,DF] = F_DF_SKT(u,para);
    err = norm(F,1)
    while err>tol && it<it_max && err<10^10
        u = u -DF\F;
        [F,DF] = F_DF_SKT(u,para);
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

    %domain (only affects the plots)
    a = 0;
    b = 1;
    % discretization for the plots
    nb = 500;
    fig = 2*l;

    plot_cos(u1,a,b,nb,'b-',fig)
    hold on
    plot_cos(u2,a,b,nb,'r--',fig)
    xlabel('$x$', 'Interpreter', 'latex')
    legend('u1','u2')
    set(gca,'FontSize',15)
    axis tight

    %"prevalidation" (without interval arithmetic)
    [rmin,rmax,Abar,w11,w12,w21,w22]=proof_SKT(u,para,nu);

    %rigorous proof (with interval arithmetic)
    if exist('intval','file') && not(isnan(rmin))
        iu = intval(u);
        iAbar = intval(Abar);
        iw11 = intval(w11);
        iw12 = intval(w12);
        iw21 = intval(w21);
        iw22 = intval(w22);
        [irmin,irmax]=proof_SKT(iu,ipara,inu,iAbar,iw11,iw12,iw21,iw22);
    end
end
