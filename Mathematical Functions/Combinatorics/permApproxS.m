function retVal=permApproxS(A,numIter)
%PERMAPPROXS Compute a stochastic approximation to the permanent of a real
%            matrix.
%
%INPUTS: A An mXn matrix. If m<=n, then the standard matrix permanent bound
%          is found. If m>n, then the permanent bound of A' is found to be
%          consistent with the permanents of A and A' being equal in square
%          matrices. Empty matrices have a permanent of one by definition.
%  numIter The number of iterations of the stochastic permanent
%          approximation algorithm to use. If this parameter is omitted or
%          an empty matrix is passed, then the default value of n*n*m is
%          used. 
%
%OUTPUTS: val An approxmation of the matrix permanent of A.
%
%The algorithm is the real algorithm described in [1]. The algorithm
%utilizes random sampling, so the results will not be equal when the
%algorithm is run multiple times.
%
%REFERENCES:
%[1] A. Barvinok, "Polynomial time algorithms to approximate permanents
%    and mixed discriminants within a simply exponential factor," Random
%    Structures and Algorithms, vol. 14, no. 1, pp. 29-61, Jan. 1999.
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Empty matrices have a permanent of 1 by definition.
    if(isempty(A))
        retVal=1;
        return; 
    elseif(numel(A)==1) 
       val=A; 
       return; 
    end 

    m=size(A,1);
    n=size(A,2);
    
    if(m>n)
        A=A';
        temp=m;
        m=n;
        n=temp;
    end
    
    if(nargin<2||isempty(numIter))
        numIter=n*n*m;
    end
    
    curComb=0:(m-1);
    
    retVal=0;
    while(~isempty(curComb))
        retVal=retVal+permApproxStochSquare(A(:,curComb+1),numIter);
        curComb=getNextCombo(curComb,n);
    end
end

function val=permApproxStochSquare(A,numIter)
    N=size(A,1);

    alpha=zeros(numIter,1);
    for curIter=1:numIter
        U=randn(N,N);
        B=U.*sqrt(A);
        alpha(curIter)=det(B)^2;
    end

    val=mean(alpha);
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
