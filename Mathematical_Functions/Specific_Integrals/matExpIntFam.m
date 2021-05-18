function [H,Q,M,W]=matExpIntFam(T,A,B,Qc)
%%MATEXPINTFAM This function explicitly solves a number of real matrix
%       integral problems that arise when solving for optimal linear
%       regulators. The integrals are:
%       H=int_0^T expm(A*s)*B ds
%       Q=int_0^T expm(A'*s)*Qc*expm(A*s) ds
%       M=int_0^T expm(A'*s)*Qc*H(s) ds
%       W=int_0^T H(s)'*Qc*H(s) ds
%       where s is the scalar parameter of integration, expm is the matrix
%       exponential function int refers to a definite integral, and the
%       function
%       H(s) is defined to be
%       H(s)=(int_0^sexpm(A*tau) dtau)*B
%
%INPUTS: T The real, positive, scalar upper bound of the integrals.
%        A A real nXn matrix.
%        B A real nXp matrix. If this parameter is omitted or an empty
%          matrix is passed, then it is assumed that p=0, in which case H,
%          M, and W will be empty matrices.
%       Qc A symmetric positive (semi)definite nXn matrix. If omitted or an
%       empty matrix is passed, then it is assumed that Qc=eye(n,n);
%
%OUTPUTS: H,Q,M,W The values of the desired matrix integrals, which are
%         defined above
%
%Expressions for obtaining the solutions to the integrals by evaluating a
%single matrix exponential are given in [1]. The type of linear regulator
%problem addressed by these integrals is given in [2].
%
%REFERENCES:
%[1] C. F. Van Loan, "Computing Integrals Involving the Matrix
%    Exponential," IEEE Transactions on Automatic Control, vol. AC-23, no.
%    3, pp. 395-404, Jun. 1967.
%[2] E. S. Armstrong and A. K. Caglayan, "An algorithm for the weighting
%    matrices in the sampled-data optimal linear regulator problem,"
%    National Aeronautics and Space Administration, Langley Research
%    Center, Hampton, VA, Tech. Rep. NASA TN D-8372, Dec. 1976.
%
%September 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(A,1);

if(nargin<4||isempty(Qc))
    Qc=eye(n,n); 
end

if(isempty(B)||nargin<2)
    p=0;
    B=zeros(n,p);
else
    p=size(B,2);
end

C=[-A.',        eye(n,n), zeros(n,n+p);
    zeros(n,n), -A.',     Qc,   zeros(n,p);
    zeros(n,2*n),         A,    B;
    zeros(p,3*n+p)];

CExp=expm(C*T);

F3=CExp((2*n+1):(3*n),(2*n+1):(3*n));
G2=CExp((n+1):(2*n),(2*n+1):(3*n));
G3=CExp((2*n+1):(3*n),(3*n+1):(3*n+p));
H2=CExp((n+1):(2*n),(3*n+1):(3*n+p));
K1=CExp(1:n,(3*n+1):(3*n+p));

H=G3;
Q=F3.'*G2;
M=F3.'*H2;

temp=B.'*F3.'*K1;
W=temp+temp.';

end

%LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.
