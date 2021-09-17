function plot_cos(f,a,b,nb,col,fig)
%Plots the function represented by the cos coefficients in f, on [a,b]
%discretized with step size nb.
%Normalization: f(t)=f_0 + 2\sum_{k\geq 2} f_k cos(k*pi*(t-a)/(b-a)).

if nargin<6
    figure
else
    figure(fig)
end
if nargin<5
    col='b';
end
if nargin<4
    nb=500;
end
if nargin<2
    a=0;
    b=pi;
end

k=0:length(f)-1;
t=(a:(b-a)/nb:b)';

M=cos((t-a)/(b-a)*k*pi);
f(2:end)=2*f(2:end);
ft=M*f;

plot(t,ft,col,'Linewidth',2)

