function [xi,w]=spherSurfPoints2SpherPoints(xiSurf,wSurf,order,method)
%%SPHERSURFPOINTS2SPHERPOINTS Given cubature points of a certain order for
%               integration over the surface of a unit (hyper-)sphere,
%               obtain cubature points for integration over the unit sphere
%               (turn cubature points for a surface integral into those for
%               a volumetric integral). 
%
%INPUTS: xiSurf, wSurf A numDimXnumPoints set of cubture points and the
%               corresponding numPointsX1 set of cubature weights for
%               integration over the surface of the unit sphere.
%         order The positive order of the spherical surface cubature points
%               passed. This will also be the order of spherical points
%               produced by the algorithm.
%        method This specified the choice of the radial rule that is used
%               for transforming the cubature points for the spherical
%               shell into cubature points for the sphere. Possible values
%               are:
%               0 (The default if omitted or an empty matrix is passed and
%                 numDim is odd). Use half of the points for an integral
%                 over -1 to 1 with weighting function w(x)=x^(2*c1) on
%                 (-1,1), c1>=0 and c1 is an integer. This is akin to the
%                 technique used in S14 14-1 in [1], pg. 292.
%               1 Use cubature points obtained directly over 0 to 1 with
%                 weighting function w(x)=x^c1 on 0-1.
%               2 (The default if omitted or an empty matrix is passed and
%                 numDim is even) Use half the cubature points for an
%                 integral over -1 to 1 with weighting function 1. This
%                 required an adjustment to the weights as compared to
%                 the approach used in algorithm S14-1 in [1]. 
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%The conversion is based on the development of Chapter 2.8 of [1], which is
%used in developing formula S3 14-1 in [1], pg. 292.
%
%Given spherical surface points of a given order, and 1D points of the same
%or higher order for integrating over the region 0-1 of x^(numDim-1)*f(x),
%then a product formula can be used to obtain the spherical points. Points
%for integrating over a 1D region from 0-1 with weighting function
%w(x)=1 can be obtained using a formula for an order that is an extra
%numDim-1 higher than would be normally needed and by multiplying the
%weights by the radial values used in the product rule as used in S14-1 in
%[1]. To make them over the region 0-1, half are discarded, which is valid
%as the points are symmetric about the origin. This is the approach of
%algorithm 2. Algorithm 0 uses a formula for integrating over the
%region -1 to 1 of x^(numDim-1)*f(x) and half of the points can be
%similarly discarded to make the algorithm valid for an integral from 0 to
%1. Algorithm 1 directly uses points designed for a region of 0 to 1.
%
%As an example, consider obtaining fourteenth-order spherical points in 3D:
% [xiSurf,wSurf]=fourteenthOrderSpherSurfCubPoints(3);
% [xi,w]=spherSurfPoints2SpherPoints(xiSurf,wSurf,14);
% %As a test for the points, consider the integral over the sphere of z^14:
% sum(bsxfun(@times,sum(xi(1,:).^14,1),w'))
%One will find that the result is the same as the exact analytic solution
%of 4*pi/255.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(xiSurf,1);
deltaStep=length(wSurf);

if(nargin<4||isempty(method))
    if(mod(numDim,2)==1)
        method=0; 
    else
        method=2;
    end
end

switch(method)
    case 0%Use half the orthogonal polynomial points for w(x)=x^(2*c1) for
          %integration between -1 and 1.
        if(mod(numDim,2)==0)
           error('This formula can only be used when the number of dimensions is odd');
        end

        [r,A]=quadraturePoints1D(ceil((order+1)/2),8,numDim-1);
        %We discard the negative half of the points as described in the
        %development to Formula S3 14-1 in [1], pg. 292. That changes the
        %points from integrating from -1 to 1 to integrating from 0 to 1.
        %This only works, because of the symmetry of the points about 0.
        sel=(r>=0);
        r=r(sel);
        A=A(sel);
        num1DPoints=length(r);
    case 1%Use orthogonal polynomials directly for w(x)=x^c1 on 0-1.
        num1DPoints=ceil((order+1)/2);
        [r,A]=quadraturePoints1D(num1DPoints,7,numDim-1);
    case 2%Half a set of Legendre polynomials is to be used.
        linPointOrder=order+numDim-1;
        [r,A]=GaussLegendrePoints1D(ceil((linPointOrder+1)/2));
        sel=(r>=0);
        r=r(sel);
        A=A(sel);
        num1DPoints=length(r);
    otherwise
        error('Unknown method specified')
end

%Allocate space
numPoints=num1DPoints*deltaStep;
xi=zeros(numDim,numPoints);
w=zeros(numPoints,1);

curStart=1;
if(method==0||method==1)
    for i=1:num1DPoints
        xi(:,curStart:(curStart+deltaStep-1))=r(i)*xiSurf;
        w(curStart:(curStart+deltaStep-1))=A(i)*wSurf;
        curStart=curStart+deltaStep;
    end
else%Method=2
    for i=1:num1DPoints
        xi(:,curStart:(curStart+deltaStep-1))=r(i)*xiSurf;
        w(curStart:(curStart+deltaStep-1))=r(i)^(numDim-1)*A(i)*wSurf;
        curStart=curStart+deltaStep;
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
