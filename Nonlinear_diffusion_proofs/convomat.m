function M=convomat(u,m_output)
%Compute the matrix associated to the convolution product of a vector of 
%cosine coeffs with a vector of cosine coeffs. That is, u is assumed to be
%a vector of cosine coeffs, and M is such that, for any vector v of cosine
%coeffs, u*v=Mv.

%The assumed normalizations is f=f_0+2*\sum_{k\geq 1} f_k cos(kx) 

%The third input is optional, and can be used to enforce the size of M. By
%default the size is the one of u.

m=length(u);

if nargin==2
    if exist('intval','file') && isintval(u(1))
        u=[u;intval(zeros(m_output-m,1))];
    else
        u=[u;zeros(m_output-m,1)];
    end
end

A=hankel(u);

C=toeplitz(u,u);%!! If u is complex, toeplitz(u) gives toeplitz(u,u'), which is not what we want here !!
C(:,1)=0;

M=A+C;

if nargin==2 && m>m_output
    M=M(1:m_output,1:m_output);
end
