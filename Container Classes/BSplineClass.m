classdef BSplineClass < matlab.mixin.Copyable
%BSPLINECLASS This class simplifies the interpolation of real or complex
%             functions with multivariate inputs and scalar or multivariate
%             outputs along with derivatives thereof.
%
%The class calls BSplinePolyFitMultiDim and then stores the interpolation
%surface for the BSplineInterpValMultiDim functions. However, in addition
%to that, the class can also store additional derivative interpolation
%polynomials that one might want to interpolate separately or at the same
%time. This can simplify the interpolation of a function as well as its
%gradient and Hessian all at once.
%    
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

properties
%The properties are class array holding arrays or matrices. The number of
%dimensions of the cell arrays equals the number of dimensions of the input
%parameters to the interpolants. For 2D, {i,j} returns the interpolant
%parameter of the (i-1)the derivative in the first dimensions and the
%(j-1)th derivative in the second dimension.
    
    aSpline
    tSpline
    tLength
    k
end
    
methods
function newObj=BSplineClass(tau,tauLengths,y,k,maxDerivs,numCentralKnots,useSparse)
%BSPLINECLASS The constructor method. Given a set of real vector (or
%         scalar) points as well as real or complex scalar function values,
%         this method initializes the class to be realy to perform
%         b-spline interpolation to a desired polynomial degree within the
%         sampled region of the function. This also precomputes
%         interpolation coefficients for the desired number of derivatives
%         in each dimension.
%
%INPUTS: tau A numPointsXnumDims matrix of the coordinates of the points in
%            each dimension in ascending order. The full grid of points
%            in all dimensions could be found using
%            [tau1,tau2,etc.]=ndgrid(tau(1:tauLengths(1),1),tau(1:tauLengths(2),2),...)
%            and the dimensions can be put together into the full set of
%            numDimsXtotalNumPoints points as
%            tauTotal=[tau1(:).';tau2(:).';...]. Note that using meshgrid
%            in 2D/3D will put the points in the wrong order; one has to
%            use the ordering of ndgrid. These values must be real.
% tauLengths A numDimsX1 or a 1XnumDims vector this lists the number of
%            points in each dimension of tau. The total number of sample
%            points is totalNumPoints=prod(tauLengths).
%          y This is the totalNumPointsXnumSets full sets of function
%            values taken at all of the points implied by the vectors in
%            tau for numSets different interpolation problems. The ordering
%            of the points is the same as the ordering of tauTotal, which
%            is described above. These values can be real or complex.
%          k A numDimsX1 or 1XnumDims set of the order of the b-splines in
%            each dimension. The value k-1 is the polynomial order of the
%            approximation being performed. If the order is the same in all
%            dimensions, then a single scalar can be passed. If this
%            parameter is omitted or an empty matrix is passed, then the
%            default of maxDerivs+2 is used.
%  maxDerivs A numDimsX1 or 1XnumDims set of the maximum number of
%            derivatives for which one wishes tocommpute interpolation
%            coefficients in each dimension. For example, in 2D,
%            maxDerivs=[2;2] means that interpolation coefficients for
%            derivatives (0,0), (1,0), (0,1), (0,2), (2,0), (1,2), (2,1)
%            and (2,2) are computed, where (a,b) indicated a derivatives in
%            the first dimension and b derivatives in the second dimension.
%            The default if this parameter is omitted or an empty matrix is
%            passed is [0;0], indicating no derivatives.
% numCentralKnots This optional parameter is used if one passes a higher
%            density of points in tau than is needed, because one wants to
%            perform least squares interpolation. This is the numDimsX1 or
%            1XnumDims list of the number of knots in the central region of
%            knots in each dimension. It varies from 0 to n-k. The higher
%            the number, the more knots are used total. A total of
%            numCentralKnots+2*k will be used by this class.
%  useSparse The linear system being solved for the coefficients can be
%            quite large, but is not necessarily dense. If useSparse is
%            true, then the matrix for the system will be allocated as a
%            sparse matrix. This can slow down the algorithm, but it can
%            also keep Matlab from running out of memory on large systems.
%            The default if this parameter is omitted or an empty matrix is
%            passed is true.
%
%OUTPUTS: newObj The initialized BSPline class object.
%
%This function format sin the input and calls BSplinePolyFitMultiDim. It
%then store
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    numDims=size(tau,2);

    if(nargin<7||isempty(useSparse))
        useSparse=true;
    end
    
    if(nargin<5||isempty(maxDerivs))
        maxDerivs=zeros(numDims,1); 
    end
    
    if(isscalar(maxDerivs))
        maxDerivs=maxDerivs*ones(numDims,1);
    end
    maxDerivs=maxDerivs(:).';
    
    if(nargin<4||isempty(k))
       k=maxDerivs+2; 
    end
    
    if(isscalar(k))
        k=repmat(k,[numDims,1]); 
    end

    if(nargin<6||isempty(numCentralKnots))
        numCentralKnots=tauLengths(:)-k(:);
    end
    
    if(any(k-1<maxDerivs))
       error('The maximum number of derivatives desired is more than supported by k.') 
    end
    
    %Allocate space for the polynomial and all desired derivatives.
    newObj.aSpline=cell(maxDerivs+1);
    newObj.tSpline=cell(maxDerivs+1);
    newObj.tLength=cell(maxDerivs+1);
    newObj.k=cell(maxDerivs+1);
    
    %Determine the coefficients for the polynomial itself.
    [a,t,tLength]=BSplinePolyFitMultiDim(tau,tauLengths,y,k,numCentralKnots,useSparse);

    newObj.aSpline{1,1}=a;
    newObj.tSpline{1,1}=t;
    newObj.tLength{1,1}=tLength;
    newObj.k{1,1}=k;
    
    newObj.precomputeDerivCoeffs(maxDerivs);
end

function vals=interpVals(theObj,points,derivList)
%INTERPVALS Interpolate values of the function value and/or its derivatives
%           at the points in points for specified derivatives. Multiple
%           sets of values can be evaluted at once.
%           
%INPUTS: theObj The implicitly passed BSplineClass object.
%       points The numDimXnumPoints set of points where one wishes to
%               perform interpolation. These values must be real.
%     derivList The numDimXnumDerivs list of derivatives that should be
%               interpolated. For examples in 3D, [0;0;0] is the function
%               value and [1;2;0] takes the 1st derivative with respect to
%               the first component and the second derivative with respect
%               to the second component. If this parameter is omitted or an
%               empty matrix is passed, then the all zero vector is used.
%
%OUTPUTS: vals A numSetsXnumDerivsXnumPoints set of interpolated values for
%              all of the points for each of the sets in the class and for
%              all of the specified derivatives.
%
%This function evaluates all of the polynomials using
%BSplineInterpValMultiDim. One cannot request more derivatives than were
%set when the class was created.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    numDims=length(theObj.tLength{1});
    numA=size(theObj.aSpline{1});
    if(length(numA)==numDims)
        numSets=1;
    else
        numSets=numA(end);
    end
    
    if(nargin<3||isempty(derivList))
        derivList=zeros(numDims,1);
    end
    
    numDerivVals=size(derivList,2);
    numPoints=size(points,2);
    
    vals=zeros(numSets,numDerivVals,numPoints);

    %Add 1 to simplify indexation.
    derivList=derivList+1;
    for curDerivIdx=1:numDerivVals
        derivIdx=derivList(:,curDerivIdx);
        
        tCur=theObj.tSpline{derivIdx(1),derivIdx(2)};
        tLengthCur=theObj.tLength{derivIdx(1),derivIdx(2)};
        aCur=theObj.aSpline{derivIdx(1),derivIdx(2)};
        kCur=theObj.k{derivIdx(1),derivIdx(2)};
        
        vals(:,curDerivIdx,:)=reshape(BSplineInterpValMultiDim(points,tCur,tLengthCur,aCur,kCur),[numSets,1,numPoints]);
    end
end

function vals=interpAllVals(theObj,points)
%INTERPVALS Interpolate values of the function value and/or its derivatives
%           at the points in x for all instances if precomputed
%           derivatives. Multiple sets of values can be evaluted at once.
%           
%INPUTS: theObj The implicitely passed BSplineClass object.
%             x The numDimXnumPoints set of points where one wishes to
%               perform interpolation. These values must be real.
%
%OUTPUTS: vals A numSetsXtotalNumDerivsXnumPoints set of interpolated
%              values for all of the points for each of the sets in the
%              class and for all of the derivatives where the coefficients
%              in the class have been precomputed.
%
%This function uses the BSplineInterpValMultiDim fucntion on all of the
%precomputed interpolation polynomials.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    numDims=length(theObj.tLength{1});

    numPoints=size(points,2);
    numA=size(theObj.aSpline{1});
    if(length(numA)==numDims)
        numSets=1;
    else
        numSets=numA(end);
    end
    
    totalNumDerivs=numel(theObj.aSpline);
    
    vals=zeros(numSets,totalNumDerivs,numPoints);
    for curDeriv=1:totalNumDerivs
        tCur=theObj.tSpline{curDerivIdx};
        tLengthCur=theObj.tLength{curDerivIdx};
        aCur=theObj.aSpline{curDerivIdx};
        kCur=theObj.k{curDerivIdx};

        vals(:,curDeriv,:)=reshape(BSplineInterpValMultiDim(points,tCur,tLengthCur,aCur,kCur),[numSets,1,numPoints]);
    end
end
end

methods(Access=private)
function precomputeDerivCoeffs(theObj,maxDerivs)
%PRECOMPUTEDERIVCOEFFS Compute interpolation polynomals for all
%             multivaraite derivatives of the function up to a given order.
%
%INPUTS: theObj The implicitely passed BSplineClass object.
%            default of maxDerivs+2 is used.
%  maxDerivs A numDimsX1 or 1XnumDims set of the maximum number of
%            derivatives for which one wishes to precompute interpolation
%            coefficients in each dimension. For example, in 2D,
%            maxDerivs=[2;2] means that interpolation coefficients for
%            derivatives (0,0), (1,0), (0,1), (0,2), (2,0), (1,2), (2,1)
%            and (2,2) are computed, where (a,b) indicated a derivatives in
%            the first dimension and b derivatives in the second dimension.
%
%OUTPUTS: None. The class itself is changed.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.  
    
    numDims=length(theObj.tLength{1});

    if(~isempty(theObj.aSpline{end}))
       %The test for whether the size of the arrays implies that the values
       %should be there is only performed if the last element is empty,
       %which implies that the class was just allocated and the values have
       %not been filled in.
        curMaxVals=size(theObj.aSpline)-1;
        
        %If the values have already been computed.
        if(all(curMaxVals>=maxDerivs))
            return; 
        end
    end
    
    a0=theObj.aSpline{1};
    tLength0=theObj.tLength{1};
    t0=theObj.tSpline{1};
    k0=theObj.k{1};
    
    if(any(k0(:)-1<maxDerivs(:)))
       error('The maximum number of derivatives desired is more than supported by k.') 
    end
    
    totalDerivEls=prod(maxDerivs+1);
    
    curDeriv=zeros(numDims,1);
    for curIdx=2:totalDerivEls
        curDeriv=getNextTuple(curDeriv,maxDerivs);
        
        %The flipud is because the last element in the tuple is the least
        %signaificant, but the ordering of elements in aSpline from Matlab
        %makes the first element the least significant.
        [aD,tD,tLengthD,kD]=BSplineInterpDerivMultiDim(t0,tLength0,a0,k0,flipud(curDeriv));
        theObj.aSpline{curIdx}=aD;
        theObj.tSpline{curIdx}=tD;
        theObj.tLength{curIdx}=tLengthD;
        theObj.k{curIdx}=kD;
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
