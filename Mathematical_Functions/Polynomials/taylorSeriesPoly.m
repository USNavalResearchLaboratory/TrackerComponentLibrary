function polyMat=taylorSeriesPoly(a,multiVarDerivs,order)
%%TAYLORSERIESPOLY Obtain the multivariate polynomial corresponding to a
%               Taylor series expansion of a scalar function that takes a
%               multivariate parameter.
%
%INPUTS: a The nX1 vector about which the Taylor series expansion is
%          performed. For a Maclauren series, set this to zero.
% multiVarDerivs A function handle or multidimensional matrix providing the
%          necessary derivatives for the Taylore series expansion. If this
%          is a function, it takes a numDimX1 vector where each element
%          corresponds to the number of derivatives taken for each element
%          (from 0 to order). If this is a matrix, then
%          multiVarDerivs(a1,a2,a3,..an) should correspond to a derivative
%          with respect to x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1).
%    order The order of the Taylor series aproximation. If this
%          parameter is omitted or an empty matrix is passed, a default
%          value of 2 is used.
%
%OUTPUTS: polyMat A hypermatrix of the coefficients for the multivariate
%               polynomial representign the Taylor series expansion. These
%               are arranged such that polyMat(a1,a2,a3...an) corresponds
%               to the coefficient of an
%               x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1) term.  Thus, the
%               number of indices polyMat takes is equal to the
%               dimensionality of x (not counting singleton dimensions at
%               the end of polyMat). Note that this ordering is the reverse
%               that used in the 1D polyval function that is built into
%               Matlab. The number of elements for each index in polyMat is
%               the order+1.
%
%A (truncated) Taylor series is an approximation to a continuous,
%differentiable function in the neighborhood of the point a based on the
%derivatives of the function. Expressions for multivariate Taylor series
%expansions above order 2 are seldom given explicitly. However, one can
%find expressions for general multivariate Taylor series expansions in
%Chapter 2.4 of [1].
%
%This function builds the polynomial from the Taylor series expansion just
%by using the given derivatives and the rules of algebra. Most of the loops
%deal with enumerating the different multivariate derivatives.
%
%EXAMPLE:
%In this example, we would like to consider the error in Taylor series
%polynomial approximations of the bivaraite function exp(x1^2*x2) in the
%neighborhood of the point x1=2, x2=3. We will consider polynomial
%approximations up to the fourth order. To do this, we need a function that
%will provide derivatives up to the fourth order. This is:
% 
% function val=testDerivs(deg,x1,x2)
%     if(deg(1)==0&&deg(2)==0)%Zeroth-order
%         val=exp(x1^2*x2);
%     elseif(deg(1)==1&&deg(2)==0)%First-order
%         val=2*exp(x1^2*x2)*x1*x2;
%     elseif(deg(1)==0&&deg(2)==1)
%         val=exp(x1^2*x2)*x1^2;
%     elseif(deg(1)==2&&deg(2)==0)%Second-order
%         val=2*exp(x1^2*x2)*x2+4*exp(x1^2*x2)*x1^2*x2^2;
%     elseif(deg(1)==1&&deg(2)==1)
%         val=2*exp(x1^2*x2)*x1+2*exp(x1^2*x2)*x1^3*x2;
%     elseif(deg(1)==0&&deg(2)==2)
%         val=exp(x1^2*x2)*x1^4;
%     elseif(deg(1)==3&&deg(2)==0)%Third-order
%         val=12*exp(x1^2*x2)*x1*x2^2+8*exp(x1^2*x2)*x1^3*x2^3;
%     elseif(deg(1)==2&&deg(2)==1)
%         val=2*exp(x1^2*x2)+10*exp(x1^2*x2)*x1^2*x2+4*exp(x1^2*x2)*x1^4*x2^2;
%     elseif(deg(1)==1&&deg(2)==2)
%         val=4*exp(x1^2*x2)*x1^3+2*exp(x1^2*x2)*x1^5*x2;
%     elseif(deg(1)==0&&deg(2)==3)
%         val=exp(x1^2*x2)*x1^6;
%     elseif(deg(1)==4&&deg(2)==0)%Fourth-order
%         val=12*exp(x1^2*x2)*x2^2+48*exp(x1^2*x2)*x1^2*x2^3+16*exp(x1^2*x2)*x1^4*x2^4;
%     elseif(deg(1)==3&&deg(2)==1)
%         val=24*exp(x1^2*x2)*x1*x2+36*exp(x1^2*x2)*x1^3*x2^2+8*exp(x1^2*x2)*x1^5*x2^3;
%     elseif(deg(1)==2&&deg(2)==2)
%         val=12*exp(x1^2*x2)*x1^2+18*exp(x1^2*x2)*x1^4*x2+4*exp(x1^2*x2)*x1^6*x2^2;
%     elseif(deg(1)==1&&deg(2)==3)
%         val=6*exp(x1^2*x2)*x1^5+2*exp(x1^2*x2)*x1^7*x2;
%     elseif(deg(1)==0&&deg(2)==4)
%         val=exp(x1^2*x2)*x1^8;
%     else
%         error('Bad degrees given')
%     end
% end
%
% %Given the above function, we create the function handle to pass the
% multiVarDerivs function as
% a=zeros(2,1);
% a(1)=2;
% a(2)=3;
% multiVarDerivs=@(deg)testDerivs(deg,a(1),a(2));
% %Next, we  create the interpolating polynomials from first to fourth
% %order as
% polyMat1=taylorSeriesPoly(a,multiVarDerivs,1);
% polyMat2=taylorSeriesPoly(a,multiVarDerivs,2);
% polyMat3=taylorSeriesPoly(a,multiVarDerivs,3);
% polyMat4=taylorSeriesPoly(a,multiVarDerivs,4);
%
% %We would like to demonstrate that near x1=2, x2=3, the approximation
% %becomes better with increasing polynomial order. Thus, we look at the
% %magnitudes of the differences between the estimates and the true value
% %at the point x1=2.01, x2=3.01
% x=zeros(2,1);
% x(1)=2.01;
% x(2)=3.01;
% 
% abs(polyValMultiDim(polyMat1,x)-testDerivs([0;0],x(1),x(2)))
% abs(polyValMultiDim(polyMat2,x)-testDerivs([0;0],x(1),x(2)))
% abs(polyValMultiDim(polyMat3,x)-testDerivs([0;0],x(1),x(2)))
% abs(polyValMultiDim(polyMat4,x)-testDerivs([0;0],x(1),x(2)))
%
%The values returned from the above example should be about 2.3329e+03,
%135.7034, 6.2048, and 0.2363 demonstrating how higher polynomial orders
%improve estimation accuracy.
%
%REFERENCES:
%[1] K. Königsberger, Analysis 2, 5th ed. Berlin: Springer, 2004, (In
%    German).
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(order))
    order=2;
end

if(isa(multiVarDerivs,'function_handle'))
    derivFun=multiVarDerivs;
else%If a matrix of values is passed.
    matDims=size(multiVarDerivs);
    derivFun=@(numDerivVec)multiVarDerivs(nDim2Index(matDims,numDerivVec+1));
end

numDim=size(a,1);

%Allocate space for the multivariate polynomial. This requires a matrix
%involving numDim indices. 
polyMat=zeros(repmat(order+1,1,numDim));

%Set the zeroth-order term.
polyMat(1)=derivFun(zeros(numDim,1));

for curOrder=1:order
    %For each order, we go through all combinations of which set of
    %variables are differentiated (we choose curOrder variables to
    %differentiate). We effectively form loops where the outermost (i1)
    %goes from 1 to numDim, the second (i2) from i1 to numDim, the third
    %(i3) from i2 to numDim and so on. Thus, we eliminate repeated
    %derivatives.
    
    %The initial set of derivatives.
    vars2Diff=ones(curOrder,1);
    curIdx=curOrder;
    isAscending=false;
    while(curIdx>0)
        if(curIdx==curOrder)
            if(curIdx==1)
                minK=1;
            else
                minK=vars2Diff(curIdx-1);
            end
            
            for k=minK:numDim
                vars2Diff(curIdx)=k;
                
                %Now, we have to build the multivariate polynomial
                %involving the (x(i)-a(i)) terms.
                multiVarPoly=1;
                for curTerm=1:curOrder
                    polyCur=zeros(repmat(1+1,1,numDim));
                    
                    %The zeroth-order term
                    polyCur(1)=-a(vars2Diff(curTerm));
                    idxVec=ones(numDim,1);
                    idxVec(vars2Diff(curTerm))=2;
                    curEl=nDim2Index(2,idxVec);
                    polyCur(curEl)=1;%The linear term.
                    
                    %Multiply the term.
                    multiVarPoly=convn(multiVarPoly,polyCur);
                end
                
                %The value for this set of partial derivatives. 
                %orderVec now contains a list of variables with which we
                %should take derivatives. We need to translate that into
                %the number of derivatives per variable.
                numDerivVec=zeros(numDim,1);
                for curTerm=1:curOrder
                    idx=vars2Diff(curTerm);
                    numDerivVec(idx)=numDerivVec(idx)+1;
                end
                
                derivVal=derivFun(numDerivVec);
                factorialTerm=1./prod(factorial(numDerivVec));
                multiVarPoly=factorialTerm*derivVal*multiVarPoly;
                
                %Add the current polynomial terms to the running sum.
                polyMat=polySumMultiDim(multiVarPoly,polyMat);
            end
            isAscending=false;
            curIdx=curIdx-1;
        elseif(isAscending==false)
            vars2Diff(curIdx)=vars2Diff(curIdx)+1;
            if(vars2Diff(curIdx)<=numDim)
                %Ascend in index
                isAscending=true;
                curIdx=curIdx+1;
            else%Keep descending in index
                curIdx=curIdx-1;
            end
        else%We are ascending in index
            vars2Diff(curIdx)=vars2Diff(curIdx-1);
            curIdx=curIdx+1;
        end
    end
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
