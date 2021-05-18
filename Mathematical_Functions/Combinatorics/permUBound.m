function retVal=permUBound(A,algorithm)
%%PERMUBOUND Calculate an upper bound to the matrix permanent. This bound
%           might be useful when approximating target-measurement
%           association probabilities using matrix permanents. The bound
%           tends to be a poor approximation of the matrix permanent.
%           However, a ratio of bounds is used in tracking, and the ratio
%           can be comparably accurate.
%
%INPUTS: A An mXn matrix. If m<=n, then the standard matrix permanent bound
%          is found. If m>n, then the permanent bound of A' is found to be
%          consistent with the permanents of A and A' being equal in square
%          matrices. Empty matrices have a permanent of one by definition.
% algorithm An optional selecting the algorithm to use. Possible values
%          are:
%          0 (The default if omitted or an empty matrix is passed) Use a
%            modified version of the matrix approximation e4 described in
%            [1].
%          1 Use the Hadamard-type bound that is given in Lemma 4.1a in
%            [2].
%          2 Use the Brégman-Minc bound that is given in Lemma 4.1b in [2].
%            This bound is only valid if A has binary elements.
%
%OUTPUTS: val An approxmation of the matrix permanent of A.
%
%In the case of a zero column sum, the e4 permanent bound in [1] by column
%is very loose, but the same approximation by row is good. In general,
%since the bound is an upper bound regardless of whether it is taken by row
%or by column, it should be evaluated both by row and by column and the
%smaller value should be taken. That is what is done in the implementation
%of algorithm 0.
%
%In [1] "conditioning techniques" for the matrix under question before
%finding the bound on the permanent are discussed. No conditioning is done
%here.
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
%[2] B. Roos, "New permanent approximation inequalities via identities,"
%    arXiv, 21 Feb. 2018. [Online]. Available:
%    https://arxiv.org/abs/1612.03702
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    algorithm=0;
end

%Empty matrices have a permanent of 1 by definition.
if(isempty(A))
    retVal=1;
    return; 
elseif(numel(A)==1) 
    retVal=A; 
    return; 
end 

switch(algorithm)
    case 0%The bound of [1].
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
    case 1%The bound of Lemma 4.1a of [2].
        N=size(A,1);
        n=size(A,2);

        %If it is not a thin matrix, then use the transpose.
        if(n>N)
            A=A';
            temp=N;
            N=n;
            n=temp;
        end
        
        retVal=exp(gammaln(N+1)-gammaln(N-n+1)+sum((1/2)*log((1/N)*sum(A.^2,1))));
    case 2%The bound of Lemma 4.1b of [2] --only for binary matrices.
        N=size(A,1);
        n=size(A,2);

        %If it is not a thin matrix, then use the transpose.
        if(n>N)
            A=A';
            temp=N;
            N=n;
            n=temp;
        end
        
        NzTilde=sum(A,1);

        zetaLn=(1./NzTilde).*gammaln(NzTilde+1);
        zetaLn(NzTilde==0)=0;
        
        retVal=exp(gammaln(N+1)-gammaln(N-n+1)-(n/N)*gammaln(N+1)+sum(zetaLn));
    otherwise
        error('Unknown bound selected.')
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
