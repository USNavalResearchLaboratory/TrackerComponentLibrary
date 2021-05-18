function coeffsRet=shrinkMultiDimPoly2Fit(coeffs)
%%SHRINKMULTIDIMPOLY2FIT Given a hypermatrix of coefficients of a
%                  multivariate polynomial, many of which might be zero,
%                  return the smallest hypermatrix possible that can
%                  represent the polynomial. This function effectively gets
%                  rid of high-order terms that are all zero.
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
%
%OUTPUTS: coeffsRet The same multivariate polynomial as represented by
%               coeffs, but shrunk as small as possible.
%
%As an example, consider a 2D polynomial:
%440-288*x2+60*x2^2-4*x2^3-110*x1+72*x1*x2-15*x1*x2^2+x1*x2^3
% coeffs=zeros(5,5);
% coeffs(0+1,0+1)=440;
% coeffs(0+1,1+1)=-288;
% coeffs(0+1,2+1)=60;
% coeffs(0+1,3+1)=-4;
% coeffs(1+1,0+1)=-110;
% coeffs(1+1,1+1)=72;
% coeffs(1+1,2+1)=-15;
% coeffs(1+1,3+1)=1;
% coeffsRet=shrinkMultiDimPoly2Fit(coeffs)
%One will find that coeffsRet is only 2X4 in size, versus the 5X5 for the
%original matrix. Both represent the same 2D polynomial.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

dimSizeList=size(coeffs);
numIdx=length(dimSizeList);

newDimSizeList=zeros(numIdx,1);

%Go through to see which if any of the coefficients including other terms
%is nonzero.
idxVec(1:numIdx)={':'};
for curIdx=1:numIdx
    %idxVec(1:(curIdx-1))={':'};
    %idxVec((curIdx+1):end)={':'};
    
    maxIdx=1;
    for i=1:dimSizeList(curIdx)
        idxVec{curIdx}=i;
        %This type of element selection is faster than loops in Matlab, but
        %actually copying the elements into els just to use the any
        %function is very inefficient.
        els=coeffs(idxVec{:});
        if(any(els(:))~=0)
            maxIdx=i;
        end
        
    end
    idxVec{curIdx}=':';
    newDimSizeList(curIdx)=maxIdx;
end

%Allocate space for the reduced-size array to return and make it the
%correct size.
totalNewEls=prod(newDimSizeList);
coeffsRet=zeros(totalNewEls,1);
coeffsRet=reshape(coeffsRet,newDimSizeList(:)');
%Now, we must copy all of the elements in coeffs that are being kept into
%the proper positions in coeffsRet.

%This time, we need to use a loop and perform the copy one element at a
%time. The copying is performed by going through the linear indices of the
%destination matrix, translating them into the indices for the source
%matrix, one at a time.
for curEl=1:totalNewEls
    coeffsRet(curEl)=coeffs(nDim2Index(dimSizeList,index2NDim(newDimSizeList,curEl)));
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
