function F=matExpInt(T,A,B)
%%MATEXPINT Evaluate the integral from 0 to T of expm(A*t)*B dt, where expm
%           is a matrix exponential and A and B are matrices. This integral
%           arises when discretizing a continuous-time linear dynamic
%           system with a control input.
%
%INPUTS: T The period of the integral, presumably positive.
%        A An NXN matrix.
%        B An NXM matrix. If omitted, the identity matrix is substituted.
%
%The algorithm is based on the technique in [1], which solves for a large
%variety of such solutions. Simplifications to avoid computing unneeded
%solutions were undertaken.
%
%REFERENCES:
%[1] C. F. Van Loan, "Computing Integrals Involving the Matrix
%    Exponential," IEEE Transactions on Automatic Control, vol. AC-23, no.
%    3, pp. 395-404, Jun. 1967.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(A,1);

if(nargin<3)
    B=eye(N);
end
M=size(B,2);

Z1=zeros(N,M);
Z2=zeros(M,M);

%Build the augmented matrix.
C=[A,B;
   Z1',Z2];

CT=expm(C*T);

F=CT(1:N,(N+1):(N+M));
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
