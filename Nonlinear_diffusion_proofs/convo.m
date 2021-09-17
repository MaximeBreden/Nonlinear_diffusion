function s=convo(a1,a2,m_output)
%Compute the convolution product of two vectors of cosine coeffs

%The assumed normalizations are f=f_0+2*\sum_{k\geq 1} f_k cos(kx)

%The two vector can be of different lengths.

%The third input prescribes the length of the output and is optional. It
%can be used to enforce zero-padding. By default, there is only
%zero-padding to ensure that the two vectors have the same length, and the
%output is of that length. This part is probably not optimal memory-wise.

m1=length(a1);
m2=length(a2);
M=max(m1,m2);
if nargin==3
    M=max(M,m_output);
end
if exist('intval','file') && (isintval(a1(1)) || isintval(a2(1)))
    a1=[a1;intval(zeros(M-m1,1))];
    a2=[a2;intval(zeros(M-m2,1))];
else
    a1=[a1;zeros(M-m1,1)];
    a2=[a2;zeros(M-m2,1)];
end

s=convomat(a1,M)*a2;

if nargin==3 && M>m_output
    s=s(1:m_output);
end
