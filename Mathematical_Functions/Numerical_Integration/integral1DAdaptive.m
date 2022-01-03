function [intEst,totalError,exitCode]=integral1DAdaptive(f,bounds,n,algorithm,c1,RelTol,AbsTol,maxSearchReg)
%%INTEGRAL1DADAPTIVE Perform adaptive numerical integration over a function
%                 that takes a scalar input and returns a scalar value. The
%                 integration is based on the use of embedded cubature
%                 points, which can use weighting functions, if necessary.
%                 The integral function built into Matlab does not support
%                 the use of high-order embedded cubature points. This
%                 function evaluates the integral from bounds(1) to
%                 bounds(2) of w(x)*f(x) dx where the bounds are finite and
%                 w(x) is a weighting function selected by the algorithm
%                 option (default w(x)=1).
%
%INPUTS: f The handle to the function. f(x) can take a vector of points and
%          return the scalar value of the function at all of the points.
%          The output should be the same size as the input.
%   bounds A 2X1 or 1X2 vector such that bounds(1) is the lower bound of
%          integration and bounds(2) is the upper bound of integration. it
%          is not required that bounds(2)>bounds(1). Integration bounds
%          must be finite.
%        n An optional input specifying the number of points to use in the
%          base set of cubature points to use for integration (not the
%          extended set). The total order of integration should be 2*n-1.
%          This parameter is not used if algorithm=-1. The default if this
%          parameter is omitted or an empty matrix is passed is 8.
% algorithm This optionally specifies the type of weighting function used
%          for the integral. This corresponds to the value in the
%          extQuadraturePoints1D function except option 10 is not
%          available and methods with unbounded integration regions cannot
%          be used. Additionally, algorithm=-1 means that n is not used
%          and the following input, c1, is used to directly provide the
%          cubature points and weights needed. See the comments to the
%          function extQuadraturePoints1D for options. The default if this
%          parameter is omitted or an empty matrix is passed is 2, which
%          corresponds to a weighting function of 1 (no weighting).
%       c1 This input is only needed if the algorithm input requires an
%          additional parameter in the extQuadraturePoints1D function, or
%          if algorithm =-1. If algorithm=-1, this is a structure
%          containing the following members:
%           wG The nX1 set of weights associated with the first n entires
%              in xiK that would provide integration of an n order 2*n-1
%              polynomial accuracy.
%          xiK The Gauss-Kronrod extension of an n-point integration rule.
%              The first n points are those for an order 2*n-1 rule using
%              weights wG and the next n+1 points are different. This with
%              the weights wK have a polynomial accuracy of 2*n+1.
%           wG The weights of xiK. These are all different from the
%              weights in wG.
%          bounds The 2X1 or 1X2 vector of bounds for the integration
%              using the points in xiK
%   RelTol The maximum relative error tolerance allowed. If omitted or an
%          empty matrix is passed, the default value of 1e-8 is used.
%   AbsTol The absolute error tolerance allowed. If omitted or an empty
%          matrix is passed, the default value of 1e-11 is used.
% maxSearchReg The algorithm works by splitting space into increasingly
%          small regions. This optional parameter is the maximum
%          number of regions allowed. If omitted or an empty matrix
%          is passed, the default value of 1000 is used.
%
%OUTPUTS: intEst The estimate of the integral.
%     totalError The estimate of the error of the integral. This is
%                usually an underestimate. It was used to determine
%                convergence.
%       exitCode A value indicating the status of the algorithm on
%                termination. Possible values are
%                0 Termination occurred due to the absolute or relative
%                  error tolerances being fulfilled.
%                1 Termination occurred due to the search region having
%                  been split into the maximum number of subregions.
%                2 A singularity of a NaN was encountered during a step;
%                  the returned values of intEst and totalError are from
%                  the previous step and do not meet the error criterion.
%                3 A singularity or NaN was encountered. The values of
%                  intEst and totalError are most likely invalid.
%
%Note that singularities at the ends of the interval are often not a
%problem as the ends of the interval are typically not evaluated unless
%finite precision errors begin to dominate near the singularity.
%
%The algorithm utilizes embedded cubature points. The function evaluated at
%a subset of points, multiplied by associated weights and summed provides
%an approximation of a certain polynomial order. Using the entire set of
%points and a different set of weights provides an approximation of a
%higher polynomial order. The magnitude of the differences in values is an
%estimate of the estimation error. If the region has been broken into
%multiple parts, then the estimation error is approximated as the sum of
%the errors of the parts.
%
%The algorithm starts by breaking the integration region into two parts. If
%the total error does not satisfy either the relative or the absolute error
%tolerances, then the region with the highest error is split in two. That
%continues until either an error tolerance is satisfied or the maximum
%number of splits occurs.
%
%EXAMPLE 1:
%This is a very simple integral. The correct answer is 2.
% bounds=[0;pi];
% f=@(x)sin(x);
% [intEst,totalError,exitCode]=integral1DAdaptive(f,bounds)
%One will see that the intEst is 2.
%
%EXAMPLE 2:
%This is a difficult integral, because it is infinitely oscillatory near
%the origin. However, it can be solved exactly. The exact solution is
%-pi/2+cos(1)+sinInt(1), which is -0.084410950559573886889031770373595 to
%extended precision. Here, we compute it to an _estimated_ seven digits of
%accuracy. A very high order cubature formula is used.
% bounds=[0;1];
% f=@(x)cos(1./x);
% AbsTol=1e-7;
% RelTol=1e-7;
% n=120;
% [intEst,totalError,exitCode]=integral1DAdaptive(f,bounds,n,[],[],AbsTol,RelTol)
%One finds a solution who true accuracy is still good, but is close to
%double the estimated accuracy. This is a very difficult integral. If one
%were to use the integral function built into Matlab (as of Matlab2016b):
% intEst=integral(f,0,1,'AbsTol',AbsTol,'RelTol',RelTol)
%one would simply get a NaN as the returned value.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(bounds))
   bounds=[-1;1]; 
end

if(nargin<3||isempty(n))
    n=8;
end

if(nargin<4||isempty(algorithm))
    algorithm=2; 
end

if(nargin<5||isempty(c1))
    c1=[]; 
end

if(nargin<6||isempty(RelTol))
    RelTol=1e-8; 
end

if(nargin<7||isempty(AbsTol))
    AbsTol=1e-11;
end

if(nargin<8||isempty(maxSearchReg))
   maxSearchReg=1000; 
end

if(algorithm~=-1)
    %Get quadrature points for integrating over -1 to 1.
    [~,wG,xiK,wK,boundsOrig]=extQuadraturePoints1D(n,algorithm,c1);
else
    wG=c1.wG;
    xiK=c1.xiK;
    wK=c1.wK;
    boundsOrig=c1.bounds;
end

%This binary heap stores the regions that are to be searched in order of
%the maximum error of the region.
searchRegions=BinaryHeap(maxSearchReg,true);

%First, we integrate the two halves of the search region and push the
%results onto the heap.

%The total error is used to determine termination contstraints.
totalError=0;
%The integral estimate could be computed at the end if one only cared about
%absolute error bounds. For relative error bounds, we must keep a running
%approximation of the integral value.
intEst=0;
center=sum(bounds)/2;
boundsDesired=[bounds(1);center];
[val,err]=integrate1DWithBounds(f,boundsDesired,boundsOrig,wG,xiK,wK);
if(~isfinite(val)||~isfinite(err))
    totalError=Inf;
    intEst=val;
    exitCode=2;%Encountered a singularity or a NaN.
    return;
end
totalError=totalError+err;
intEst=intEst+val;

regionParams.regionValue=val;
regionParams.bounds=boundsDesired;
searchRegions.insert(err,regionParams);

boundsDesired=[center;bounds(2)];
[val,err]=integrate1DWithBounds(f,boundsDesired,boundsOrig,wG,xiK,wK);
if(~isfinite(val)||~isfinite(err))
    totalError=Inf;
    intEst=val;
    exitCode=2;%Encountered a singularity or a NaN.
    return;
end
totalError=totalError+err;
intEst=intEst+val;

regionParams.regionValue=val;
regionParams.bounds=boundsDesired;
searchRegions.insert(err,regionParams);

if(~isfinite(totalError)||~isfinite(intEst))
    exitCode=3;%Other overflow error
    return;
end

%Next, we check whether the convergence criterion is fulfilled and if not,
%start dividing the search region.
while(1)
    %Determine whether either of the error bounds has been met.
    if(totalError<AbsTol||totalError<RelTol*abs(intEst))
        exitCode=0;
        break;
    end

    %If the heap of search regions is full, then terminate the algorithm.
    if(searchRegions.heapSize()>=maxSearchReg)
        exitCode=1;
        return;
    end

    %Remove the region with the highest error and split it.
    topRegion=searchRegions.deleteTop();
    
    curVal=topRegion.value.regionValue;
    curBounds=topRegion.value.bounds;
    curErr=topRegion.key;
    
    %Remove the effects of the last integration step
    intEstPrev=intEst;
    totalErrorPrev=totalError;
    intEst=intEst-curVal;
    totalError=totalError-curErr;
    
    %Now add the integrals over the two smaller regions.
    center=sum(curBounds)/2;
    boundsDesired=[curBounds(1);center];
    [val,err]=integrate1DWithBounds(f,boundsDesired,boundsOrig,wG,xiK,wK);
    if(~isfinite(val)||~isfinite(err))
        intEst=intEstPrev;
        totalError=totalErrorPrev;
        exitCode=2;%Encountered a singularity or a NaN.
        return;
    end
    
    totalError=totalError+err;
    intEst=intEst+val;
    
    regionParams.regionValue=val;
    regionParams.bounds=boundsDesired;
    searchRegions.insert(err,regionParams);

    boundsDesired=[center;curBounds(2)];
    [val,err]=integrate1DWithBounds(f,boundsDesired,boundsOrig,wG,xiK,wK);
    if(~isfinite(val)||~isfinite(err))
        intEst=intEstPrev;
        totalError=totalErrorPrev;
        exitCode=2;%Encountered a singularity or a NaN.
        return;
    end
    
    totalError=totalError+err;
    intEst=intEst+val;

    regionParams.regionValue=val;
    regionParams.bounds=boundsDesired;
    searchRegions.insert(err,regionParams);
    
    if(~isfinite(totalError)||~isfinite(intEst))
        intEst=intEstPrev;
        totalError=totalErrorPrev;
        exitCode=3;%Other overflow error
        return;
    end
end
end

function [valHigh,err]=integrate1DWithBounds(f,boundsDesired,boundsOrig,wG,xiK,wK)
%%INTEGRATE1DWITHBOUNDS This function uses the embedded cubature points to
%               perform 1-dimensional integration over the function f in
%               the desired finite region with upper and lower bounds in
%               boundsDesired. The finite bounds for which the cubature
%               points were desired are given in boundsOrig. This function
%               is simply implemented using a change of variables
%               substitution.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

%Extract upper and lower desired integration bounds as well as the bounds
%given by the points xiG and xiK, the original integration bounds.
lbOrig=boundsOrig(1);
ubOrig=boundsOrig(2);
lbDes=boundsDesired(1);
ubDes=boundsDesired(2);

deltaOrig=ubOrig-lbOrig;
deltaDes=ubDes-lbDes;
deltaRatio=(deltaDes/deltaOrig);

numLowerApprox=length(wG);
points=lbDes+deltaRatio*(xiK-lbOrig);

vals=f(points);

valsG=vals(1:numLowerApprox);
valLo=sum(wG(:).*valsG(:))*deltaRatio;
valHigh=sum(wK(:).*vals(:))*deltaRatio;

err=abs(valHigh-valLo);
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
