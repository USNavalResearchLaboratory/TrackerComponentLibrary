function [a,F,f]=solveAutoregressiveSys(u,L,useTLS)
%%SOLVEAUTOREGRESSIVESYS This solves an autoregressive linear system or
%           order L for the parameters given the order of the system. This
%           is solving the system of equations
%           [u(L+1)]  [u(1), ..., u(L-1), u(L)  ] [a(1)]
%           |u(L+2)|  |u(2), ..., u(L),   u(L+1)| |a(2)|
%           |u(L+3)|= |u(3), ..., u(L+1), u(L+2)|*|a(3)|
%           |...   |  |...   ...   ...    ...   | |... |
%           [u(N)  ]  [u(N-L)..., u(N-2), u(N-1)] [a(L)]
%           Given u. This is equivalent to the difference equation
%           u(N)=a(1)*u(N-1)+a(2)*u(N-2)+...+a(L)*u(N-L)
%           This function provides the least-squares solution.
%           
%INPUTS: u An NXnumSnapshots matrix of numSnapshots of independent sets of
%          data. This function provides a least squares solution over all
%          snapshots. The data can be real or complex.
%        L The order of the model, L<=fix(N/2). If this parameter is
%          omitted or an empty matrix is passed, then L=fix(N/2) is used.
%   useTLS This function just solves a linear system. The system can be
%          solved either using the pseudoinverse function pinv or using the
%          totalLeastSquares function to solve the total least squares
%          problem. If this prameter is true, then totalLeastSquares is
%          used. The default if this parameter is omitted or an empty
%          matrix is passed is true.
%
%OUTPUTS: a The LX1 set of parameters of the system.
%         F The matrix multipled by the a vector in the linear system
%           described above.
%         f The subset of the u that is on the left-side of the linear
%           system described above.
%
%Autoregressive systems arise in a number of filtering problems, such as in
%Prony's direction estimation in [2], and are discussed in Chapter 1.5 of
%[1].
%
%This function just solves the linear system using the pinv or
%totalLeastSquares functions. For multiple snapshots, the extra rows from
%the new data are stacked.
%
%EXAMPLE:
% L=5;
% N=12;
% a=((L:-1:1)/10)';%The true parameters
% u=zeros(N,1);
% u(1:L)=rand(L,1);%The inital values to start
% %Fill in the rest of u via the recursion described above.
% for k=(L+1):N
%     u(k)=sum(u((k-L):(k-1)).*a);
% end
% aEst=solveAutoregressiveSys(u,L)
%One will see that aEst is the same as a, within finite precision limits.
%
%REFERENCES:
%[1] S. Haykin, Adaptive Filter Theory, 4th ed. Upper Saddle River, NJ:
%    Prentice Hall, 2002.
%[2] G. M. Pitstick, J. R. Cruz, and R. J. Mulholland, "A novel
%    interpretation of Prony's method," Proceedings of the IEEE, vol. 76,
%    no. 8, pp. 1052-1053, Aug. 1988.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(u,1);
numSnapshots=size(u,2);

if(nargin<2||isempty(L))
    L=fix(N/2);
end

if(nargin<3||isempty(useTLS))
   useTLS=true; 
end

numRows=numSnapshots*(N-L);

F=zeros(numRows,L);
f=zeros(numRows,1);
curRow=1;
for curSnapshot=1:numSnapshots
    for k=1:(N-L)
        F(curRow,:)=u((k:(L+k-1)),curSnapshot);
        f(curRow)=u(L+k,curSnapshot);
        curRow=curRow+1;
    end
end

if(useTLS)
    a=totalLeastSquares(F,f);
else
    a=lsqminnorm(F,f);
end


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
