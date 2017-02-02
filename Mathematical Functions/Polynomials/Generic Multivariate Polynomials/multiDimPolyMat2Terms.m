function termMat=multiDimPolyMat2Terms(coeffs,ordering,padToLength)
%%MULTIDIMPOLYMAT2TERMS Given a multivariate polynomial represented as a
%         hypermatrix of its coefficients (including cross terms), obtain a
%         2D matrix representing the same polynomial, where only the
%         nonzero terms are kept.
%
%INPUTS: coeffs A hypermatrix of the coefficients for the multivariate
%               polynomial. These are arranged such that
%               coeffs(a1,a2,a3...an) corresponds to the coefficient of an
%               x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1) term.  Thus, the
%               number of indices coeffs takes is equal to the
%               dimensionality of x (not counting singleton dimensions at
%               the end of coeffs). Note that this ordering is the reverse
%               that used in the 1D polyval function that is built into
%               Matlab. The number of elements for each index in coeffs is
%               the maximum order of that dimension +1.
%      ordering Specifies the type of ordering to use for the terms in
%               termMat. Possible values are
%               0 (The default if omitted or an empty matrix is passed) Use
%                 lexicographic order in the powers a1,a2,a3...an taking an
%                 to be the most significant power. This means that an is
%                 taken to be the most significant value, a(n-1) the next =
%                 most and so on. For three indices, this means that
%                 (1,0,0)<(2,0,0)<(0,1,0)<(1,2,0), for example.
%                 Lexicographic order means that the LOWEST value in
%                 a1,a2,a3...an comes first.
%               1 Use graded lexicographic ordering. In this instance, the
%                 coefficients sets are in order of increasing values of
%                 (a1+a2+...+an), and ties are broken with lexicographic
%                 order with an being the most significant element. This si
%                 also known as degree negative lexicographic ordering.
%               2 Use lexicographic ordering where a1 is the most
%                 significant element. This option provides the fastest
%                 conversion from a coefficient hypermatrix to a matrix of
%                 terms.
%               3 Use graded lexicographic ordering where a1 is the most
%                 significant element
%   padToLength If this parameter is omitted, then termMat is an
%               (n+1)XnumTerms matrix. However, if padToLength is provided,
%               then termMat is a max(n+1,padToLength+1)Xnumterms matrix,
%               with the extra terms all zero.
%
%OUTPUTS: termMat An (n+1)XnumTerms matrix (or with more rows, depending on
%                 padToLength) such that termMat(:,i)=[c,a1,a2,...,an] is
%                 the a term from coeffs where c is the value of
%                 coeffs(a1,a2,a3...an). Terms with c=0 are not included in
%                 termMat.
%
%EXAMPLE:
%Consider the multivariate polynomial 14+3*x1^2-18*x2+12*x1*x3-3*x2*x3:
% coeffs=zeros(3,2,2);
% coeffs(0+1,0+1,0+1)=14;
% coeffs(2+1,0+1,0+1)=3;
% coeffs(0+1,1+1,0+1)=-18;
% coeffs(1+1,0+1,1+1)=12;
% coeffs(0+1,1+1,1+1)=-3;
% %Lexicographic ordering  with the last element being the most signfiicant
% %is
% termMat=multiDimPolyMat2Terms(coeffs,2)
% %where one will get
% termMat =
%     14     3   -18    12    -3
%      0     2     0     1     0
%      0     0     1     0     1
%      0     0     0     1     1
% %Graded lexicographic ordering with the last element being the most
% %significant is
% termMat=multiDimPolyMat2Terms(coeffs,3)
% %where one will get
% termMat =
%     14   -18     3    12    -3
%      0     0     2     1     0
%      0     1     0     0     1
%      0     0     0     1     1
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(ordering))
    ordering=0;
end

if(nargin<3||isempty(padToLength))
    padToLength=0;
end

if(ordering>1)
    MSESwap=false;
    ordering=ordering-2; 
else
    MSESwap=true;
end

%Find will give us the indices of the nonzero coefficients in lexicographic
%order.
idx=find(coeffs(:));
dimVec=size(coeffs);
numDims=length(dimVec);
numTerms=length(idx);

%Allocate space
termMat=zeros(max(numDims+1,padToLength+1),numTerms);

%First, construct termMat in lexicographic order with the first index being
%the LEAST significant.
for i=1:numTerms
    curIdx=idx(i);
    indices=index2NDim(dimVec,curIdx);
    termMat(1:(numDims+1),i)=[coeffs(curIdx);indices-1];
end

%If one desires a lexicographic ordering where the first index is the MOST
%significant, then the data must be resorted. There are most likely more
%efficient ways to do this, but repeated sorting is the simplest.
if(MSESwap==true)
    gtCompareFunc=@(a,b)(lexOrderGT(a(2:end),b(2:end),true)>0);
    termMat=heapSort(termMat,true,gtCompareFunc);
end

switch(ordering)
    case 0
        return;
    case 1%A graded ordering is desired.
        %Sort by the order (grade)  using a stable sorting algorithm so
        %that the secondary lexicographic or reverse lexicographic ordering
        %remains.
        orders=sum(termMat(2:end,:),1);
        %We use bubble sort, because it is a stable sorting algorithm.
        %Thus, the (reverse) lexicographic order of the data will be
        %preserved to break ties. bubbleSort sorts in ascending order.
        [~,idxList]=bubbleSort(orders);
        termMat=termMat(:,idxList);
    otherwise
        error('Unknown ordering option specified.')
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
