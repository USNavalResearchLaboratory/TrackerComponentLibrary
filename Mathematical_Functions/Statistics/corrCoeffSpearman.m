function [C,R]=corrCoeffSpearman(x1,x2,assumeUnique)
%%CORRCOEFFSPEARMAN Compute the Spearman rank correlation coefficient 
%                   matrix of a number of variables across columns of a
%                   matrix or of two vectors. Introduced in [1], the
%                   Spearman rank correlation coefficient between samples
%                   (x_i, y_i) of two random variables x and y is a
%                   correlation measure of the variables that is useful
%                   even if the variables in question have very
%                   non-Gaussian distributions.
%
%INPUTS: x1,x2 If x2 is omitted or an empty matrix is passed, then x1 is a
%              real NXxDim matrix of N samples of xDim variables whose
%              pairwise rank correlations are desired. If both x1 and x2
%              are given, then they are both length N vectors of separate
%              variables, so xDim=2.
% assumeUnique The computation of the rank correlation coefficients is
%              faster if it can be assumed that all of the elements in each
%              of the xDim vectors are unique. If they are known to be
%              unique, then one can set this to true. The default if
%              omitted or an empty matrix is passed is false.
%
%OUTPUTS: C An xDimXxDim Spearman rank correlation coefficient matrix.
%           C(i,j) is the rank correlation coefficient between variables i
%           and j. The diagonal is always all ones. The off-diagonal
%           elements range from -1 to 1.
%         R The NXxDim matrix of the ranks of each of the input values.
%
%Spearman's rank correlation coefficient is discussed in Chapter 1.14 of
%[2]. Equation 3.8 of Chapter 3.8 of [2] is implemented here to obtain the
%correlation coefficients. When there might be duplicates, the function
%rankOrder is used to obtain the rankings.
%
%When there are no duplicates, the Spearman rank correlation coefficient is
%equivalent to replacing the data in x1 with its rank and then computing
%the Pearson correlation coefficient of that (e.g. the corrcoef function in
%Matlab).
%
%EXAMPLE 1:
%This is the example in Chapter 1.14 of [1].
% x1=[7,4,3,10,6,2,9,8,1,5];
% x2=[5,7,3,10,1,9,6,2,8,4];
% C=corrCoeffSpearman(x1,x2)
%One will find that C is a 2X2 matrix with 1's on the main diagonal and
%-0.103030303030303 on the off-diagonals.
%
%EXAMPLE 2:'
%In this example, there is a repeated value, meaning that the ranks are
%repeated.
% xMat=[ 4,   90;
%        6, -120;
%       12,   82;
%        9,   64;
%       12,   pi];
% C=corrCoeffSpearman(xMat) 
%One find the off-diagonal elements of C to be -0.205195670417031.
%
%REFERENCES:
%[1] C. Spearman, "The proof and measure of association between two
%    things," The American Journal of Psychology, vol. 15, no. 1, pp.
%    72-101, Jan. 1904.
%[2] M. G. Kendall, Rank Correlation Methods, 3rd ed. LiverPool: Charles
%    Birch and Sons Ltd., 1962.
%
%August 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(assumeUnique))
    assumeUnique=false;
end

if(nargin<2||isempty(x2))
    xDim=size(x1,2);
    n=size(x1,1);
    nConst=n*(n^2-1)/6;

    %Replace x1 with the rank orders.
    if(assumeUnique)
        [~,R]=sort(x1,'ascend');
        for i=1:xDim
            R(:,i)=inversePermutation(R(:,i));
        end

        s=zeros(1,xDim);
    else
        [R,s]=rankOrder(x1,false,false);
        s=s/12;
    end

    %Allocate space and set the diagonal elements to 1.
    C=eye(xDim,xDim);
    for i=1:(xDim-1)
        for j=(i+1):xDim
            d2Sum=sum((R(:,i)-R(:,j)).^2);
            
            %The correlation coefficient in Equation 3.8 in Chapter 3.8 of
            %[2].
            rho=(nConst-d2Sum-s(i)-s(j))/sqrt((nConst-2*s(i))*(nConst-2*s(j)));
            
            C(i,j)=rho;
            C(j,i)=rho;
        end
    end
else
    x1=x1(:);
    x2=x2(:);
    
    n=length(x1);
    nConst=n*(n^2-1)/6;
    
    R=zeros(n,2);
    
    if(assumeUnique)
        [~,r]=sort(x1,'ascend');
        R(:,1)=inversePermutation(r);
        [~,r]=sort(x2,'ascend');
        R(:,2)=inversePermutation(r);
        
        s1=0;
        s2=0;
    else    
        [r,s1]=rankOrder(x1,false,false);
        R(:,1)=r;
        [r,s2]=rankOrder(x2,false,false);
        R(:,2)=r;
    end
    s1=s1/12;
    s2=s2/12;
    
    d2Sum=sum((R(:,1)-R(:,2)).^2);
    
    %The correlation coefficient in Equation 3.8 in Chapter 3.8 of [2].
    rho=(nConst-d2Sum-s1-s2)/sqrt((nConst-2*s1)*(nConst-2*s2));

    C=[1,  rho;
       rho,1];
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
