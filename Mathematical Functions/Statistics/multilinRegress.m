function [coeffs,sigma2]=multilinRegress(zSamp,y,w,hasIntercept)
%%MULTILINREGRESS Perform multilinear regression. This means, given a set
%         of multidimensional points of the form [x1;x2;x3;...;xn],
%         determine the line/plane/hyperplane best fitting the points
%         such that
%         xn=coeffs(1)*x(2)+coeffs(2)*x(3)+....coeffs(n-1)*x(n-1)+coeffs(n)
%         if there is an intercept or
%         xn=coeffs(1)*x(2)+coeffs(2)*x(3)+....coeffs(n)*x(n).
%         if it is known there is a zero intercept.
%         The "best" fitting line/plane/hyperplane is the one that
%         minimizes
%         sum_{i=1}^N w(i)*(dot(coeffs(1:(n-1)),zSamp(1:(n-1),i))+coeffs(n)-zSamp(n,i))^2
%         if there is an intercept and y is not provided or that minimizes
%         sum_{i=1}^N w(i)*(dot(coeffs(1:(n-1)),zSamp(:,i))+coeffs(n)-y(i))^2
%         if y is provided and there is a y intercept. If there is no y
%         intercept, then the coeffs(n) term goes away in both instances.
%         The w(i) terms are all positive weights. For unweighted fusion,
%         w(i)=1.
%
%INPUTS: There are two parameterizations for zSamp and y. If y is an empty
%        matrix, then zSamp is a numDimsXnumPoints set of numPoints data
%        samples numDims>=2, as given above. On the other hand, if y is
%        provided, then y corresponds to zSamp(n,:) in the previous
%        parameterization (and can be a row or column vector) and zSamp thus
%        lacks the final row that it has in the initial parameterization.
%        Thus, zSamp will be a (numDims-1)XnumPoints matrix when y is
%        given.
%      w If weighted linear least squares optimization is supposed to be
%        performed, this is a length numPoints vector of values >0
%        containing the weight of each sample. If the regression is
%        performed with un unweighted cost function, this parameter can be
%        omitted or an empty matrix can be passed.
% hasIntercept A boolean value indicating whether the fit has a nonzero y
%        intercept. The default if omitted or an empty matrix is passed is
%        true.
%
%OUTPUTS: coeffs The numDimsX1 set of weights for the linear equation
%                described above. If hasIntercept=true, then the last
%                element is the additive constant.
%         sigma2 The residual mean square.  This is SSE/(numPoints-numDims)
%                where SSE is the vector of ordinary residuals e'*diag(w)*e
%                where e=(y(:)-x'*coeffs), where x is
%                [zSamp;ones(1,numPoints)] if an intercept is desired and y
%                is given, or x=zSamp if an intercept is not desired and y
%                is given. Similar expressions for when y is not given
%                hold, whereby y is replaced by the last row of zSamp. This
%                is the quantity being minimized in the optimization.
%
%Formulae for this type of regression are given in Appendix A.7.1 and A.7.3
%of [1] for 2D and 3D measurements without the weights. This function just
%continues the pattern implementing the results for an arbitrary number of
%dimensions, for weights and a modification allowing for there to be a
%known zero y intercept is also added. The derivation is just basic
%differential calculus. Note that this function will fail if the plane
%being estimated does not depend on the independent parameter. In 2D, this
%would correspond to a vertical line (e.g. x=4).
%
%EXAMPLE 1:
%Here is an example with a nonzero y intercept using both input formats.
% coeffs=[4;6;2];
% numPts=20;
% xVals=linspace(-4,4,numPts);
% yVals=linspace(-4,4,numPts);
% [x,y]=ndgrid(xVals,yVals);
% xyPts=[x(:).';y(:).'];
% zPts=sum(bsxfun(@times,coeffs(1:2),xyPts),1)+coeffs(3);
% zPts=zPts+0.1*randn(1,numPts^2); %Add noise.
% zSamp=[xyPts;zPts];
% [coeffsEst,sigma2]=multilinRegress(zSamp)
% [coeffsEstAlt,sigma2Alt]=multilinRegress(xyPts,zPts)
%One finds that coeffsEst is the same as coeffsEstAlt and is close to
%coeffs. A perfect fit (wihtin finite precision limitaitions) is achieved
%if there is zero noise.
%
%EXAMPLE 2:
%Here is an example when it is known that the y intercept is zero.
% coeffs=[4;6];
% numPts=20;
% xVals=linspace(-4,4,numPts);
% yVals=linspace(-4,4,numPts);
% [x,y]=ndgrid(xVals,yVals);
% xyPts=[x(:).';y(:).'];
% zPts=sum(bsxfun(@times,coeffs(1:2),xyPts),1);
% zPts=zPts+0.1*randn(1,numPts^2); %Add noise.
% zSamp=[xyPts;zPts];
% [coeffsEst,sigma2]=multilinRegress(zSamp,[],[],false)
% [coeffsEstAlt,sigma2Alt]=multilinRegress(xyPts,zPts,[],false)
%Again, the results are comparable to those in problem 1 with the two input
%formulations equal and close to coeffs.
%
%REFERENCES:
%[1] P. J. Schneider and D. H. Eberly, Geometric Tools for Computer
%    Graphics. Amsterdam: Morgan Kaufmann Publishers, 2003.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDims=size(zSamp,1);
numPoints=size(zSamp,2);

if(nargin<4||isempty(hasIntercept))
    hasIntercept=true; 
end

if(nargin<3)
   w=[]; 
end

if(nargin>1&&~isempty(y))
    if(hasIntercept)
        x=[zSamp;ones(1,numPoints)];
        p=numDims+1;
    else
        x=zSamp;
        p=numDims;
    end
else
    if(hasIntercept)
        x=[zSamp(1:(numDims-1),:);ones(1,numPoints)];
        p=numDims;
    else
        x=zSamp(1:(numDims-1),:);
        p=numDims-1;
    end
    
    y=zSamp(numDims,:);
    y=y(:);
end

if(isempty(w))
    A=x*x';
    b=x*y(:);
else
    %In Matlab, x*eye(size(x,2))*x'-x*x' and
    %bsxfun(@times,x,ones(numPoints,1).')*x'-x*x' are not numerically zero,
    %which means that just using w with weights all equal to 1 will change
    %the last few bits of the result.

    w=w(:);
    A=bsxfun(@times,w.',x)*x';
    b=x*(w.*y(:));
end

coeffs=A\b;

if(nargout>1)
    if(nargin>1&&~isempty(y))
        e=y(:)-x'*coeffs;
    else
        e=vec(zSamp(numDims,:))-x'*coeffs;
    end

    if(isempty(w))
        %The residual sum of squares.
        sigma2=e'*e/(numPoints-p);
    else%The weighted residual sum of squares.
        sigma2=(w.*e)'*e/(numPoints-p);
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
