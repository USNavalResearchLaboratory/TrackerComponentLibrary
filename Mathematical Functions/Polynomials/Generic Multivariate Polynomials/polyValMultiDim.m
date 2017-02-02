function val=polyValMultiDim(coeffs,x)
%%POLYVALMULTIDIM Evaluate a multidimensional polynomial at one or more
%                 points given a hypermatrix of its coefficients
%                 (including cross terms).
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
%             x The xDimXnumPoints set of points at which the multivariate
%               polynomial should be evaluated. Note that if xDim is larger
%               than the number of indices that coeffs takes, it is assumed
%               that the polynomial does not depend on the extra elements
%               of x. This deals with issues regarding trailing singleton
%               dimensions.
%
%OUTPUTS: val A 1XnumPoints vector of the polynomial evaluated at the given
%               points.
%
%In 1D, one typically uses some form of Horner's method to efficiently
%evaluate a multivariate polynomial. However, devising an optimal Horner's
%method for multiple variables is difficult as there is no one unique way
%to factor things; different orderings can lead to algorithms with
%different complexities. Moreover such multivariate Horner's methods often
%have a bit of overhead. Thus, this function just directly evaluates every
%term in the sum, using something of a recusion so that the computation of
%the powers of the multivariate x is more efficient.
%
%As an example, consider the multivariate polynomial
%x1^3+12*x1+24*x2-36*x3-3*x1*x2*x3+16*x3^3+6*x1^2*x3-18*x1*x2+64
%This can be evaluated at the point [12;24;36] as
% coeffs=zeros(4,2,4);
% coeffs(3+1,0+1,0+1)=1;
% coeffs(1+1,0+1,0+1)=12;
% coeffs(0+1,1+1,0+1)=24;
% coeffs(0+1,0+1,1+1)=-36;
% coeffs(1+1,1+1,1+1)=-3;
% coeffs(0+1,0+1,3+1)=16;
% coeffs(1+2,0+1,1+1)=6;
% coeffs(1+1,1+1,1+0)=-18;
% coeffs(0+1,0+1,0+1)=64;
% x=[12;24;36];
% val=polyValMultiDim(coeffs,x)
%One will find that the result is 742528.
%
%Note that when evaluating a 1D polynomial the order of the coefficients is
%reverse that of the polyval function in Matlab. For example, to evaluate 
%2*x^2-3*x+4 at the point x=-8, one uses
% coeffs=[4;-3;2];
% x=-8;
% val=polyValMultiDim(coeffs,x)
%to get a result of 156.
%
%Note that the convn function can be used for multivariate polynomial
%multiplication.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The maximum degrees of the polynomials are size(coeffs)-1.
highestDegsP1=size(coeffs);
numDim=length(highestDegsP1);

%Matlab suppresses trailing singleton dimensions except when something is
%1D. Thus, we shall determine the actual dimensionality.
if(numDim==2)
    if(highestDegsP1(2)==1&&highestDegsP1(1)==1)
        %If all dimensions are singleton, then that just represents a
        %constant value and does not depend on the input.
        val=repmat(coeffs(1),[1,size(x,2)]);
        return;
    elseif(highestDegsP1(2)==1)
        %The last dimension is singleton and thus does not contribute to
        %the polynomial value.
        numDim=numDim-1;
    end
end

%Trailing singleton dimensions do not play a role in the polynomial value.
%Other singleton dimensions also do not play a role, but it is harder to
%omit them.
numDim=min(size(x,1),numDim);
numX=size(x,2);

degList=zeros(numDim,1);
cumProds=ones(numDim,numX);
curDim=numDim;
isBacktracking=false;
val=zeros(1,numX);
curCoeffIdx=1;

%The loop uses degList to go through all of the powers of all of the
%products of the x terms.
while(curDim<=numDim)
    if(curDim==1)
        %If at a leaf node, add the terms together
        val=val+cumProds(1,:)*coeffs(curCoeffIdx);%Add the x^0 term.
        curCoeffIdx=curCoeffIdx+1;
        
        %Add the higher order terms.
        for curOrder=1:(highestDegsP1(1)-1)
            cumProds(1,:)=cumProds(1,:).*x(1,:);            
            val=val+cumProds(1,:)*coeffs(curCoeffIdx);
            curCoeffIdx=curCoeffIdx+1;
        end
        
        %Backtrack
        isBacktracking=true;
        curDim=curDim+1;
    elseif(isBacktracking)
        degList(curDim)=degList(curDim)+1;
        if(degList(curDim)<highestDegsP1(curDim))
            cumProds(curDim,:)=cumProds(curDim,:).*x(curDim,:);
            isBacktracking=false;
            curDim=curDim-1;
            cumProds(curDim,:)=cumProds(curDim+1,:);
        else
            curDim=curDim+1;%Backtrack another level.
        end
    else%Not backtracking, just entering a new level
        degList(curDim)=0;
        %Continue downward.
        curDim=curDim-1;
        cumProds(curDim,:)=cumProds(curDim+1,:);
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
