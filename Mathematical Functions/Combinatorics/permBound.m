function retVal=permBound(A)
%%PERMBOUND  Calculate an upper bound to the matrix permanent. This bound
%            might be useful when approximating target-measurement
%            association probabilities using matrix permanents.
%
%INPUTS: A An mXn matrix. If m<=n, then the standard matrix permanent bound
%          is found. If m>n, then the permanent bound of A' is found to be
%          consistent with the permanents of A and A' being equal in square
%          matrices. Empty matrices have a permanent of one by definition.
%
%OUTPUTS: val An approxmation of the matrix permanent of A.
%
%This function implements a modified version of the the matrix permanent
%approximation e4 described in [1]. In the case of a zero column sum, the
%e4 permanent bound by column is bad, but the same approximation by row is
%good. In general, since the bound is an upper bound regardless of whether
%it is taken by row or by column, it should be evaluated both by row and by
%column and the smaller value should be taken. That is what is done here.
%
%The above reference also discusses "conditioning techniques" for the
%matrix under question before finding the bound on the permanent. No
%conditioning is done here.
%
%When given an mXn rectangular matrix with n>m (more rows than columns),
%the matrix permanent is the sum of all permanents of all possible square
%matrices formed by choosing m columns. Thus, the matrix permanent
%upper bound here for rectangular matrices is obtained through a similar
%sum. For the case of m>n, the permanent bound of the transposed matrix is
%returned. The permanent of an empty matrix is defined to be one.
%
%REFERENCES:
%[1] J. K. Uhlmann, "Matrix permanent inequalities for approximating joint
%    assignment matrices in tracking systems," Journal of the Franklin
%    Institute, vol. 341, no. 7, pp. 569-593, Nov. 2004.
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
    
    curComb=0:(m-1);
    
    retVal=0;
    while(~isempty(curComb))
        ACur=A(:,curComb+1);
        
        bound1=permApproxSquare(ACur);
        bound2=permApproxSquare(ACur');

        retVal=retVal+min(bound1,bound2);
        
        curComb=getNextCombo(curComb,n);
    end
end

function retVal=permApproxSquare(A)
n=size(A,1);

%Alllocate space
c=sum(A,2);%Holds the column sums.

retVal=1;
for curRow=1:n
    retVal=retVal*c(curRow)*min(sum(A(curRow,:)./c'),1);
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
