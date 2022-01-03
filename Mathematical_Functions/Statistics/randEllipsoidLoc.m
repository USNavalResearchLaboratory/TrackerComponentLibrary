function [CartPoint,latLong]=randEllipsoidLoc(N,a,f)
%%RANDELLIPSOIDLOC Generate a random point uniformly distributed on the
%              surface of a reference ellipsoid, such as the WGS-84
%              reference ellipsoid. The point is provided in Cartesian
%              coordinates as well as in terms of latitude and longitude.
%
%INPUTS: N The number of random samples on the reference ellipsoid to draw.
%        a The semi-major axis of the reference ellipsoid. If this argument
%          is omitted, the value in Constants.WGS84SemiMajorAxis is used.
%        f The flattening factor of the reference ellipsoid. If this
%          argument is omitted, the value in Constants.WGS84Flattening is
%          used.
%
%OUTPUTS: CartPoint A 3XN matrix of random Cartesian [x;y;z] points on the
%                   ellipsoid with the given semi-major axis and flattening
%                   factor.
%           latLong A 2XN matrix of ellipsoidal latitude and longitude
%                   corresponding to points in CartPoint.
%
%The problem of generating a uniformly distributed sample on a reference
%ellipsoid is number 28 in Chapter 3.4.1 of [1], where the algorithm to do
%so is given in the back of the book. That algorithm is implemented here,
%where the general ellipsoid formulation has been modified to support the
%typical parameterization in terms of a semimajor axis and a flattening
%factor. The algorithm is a rejection sampling method.
%
%Large flattening factor values can make the algorithm slow.
%
%REFERENCES:
%[1] D. Knuth, The Art of Computer Programming: Seminumerical Algorithms,
%    3rd ed. Reading, MA: Addison-Wesley, 1998, vol. 2.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    f=Constants.WGS84Flattening;
end

if(nargin<2)
    a=Constants.WGS84SemiMajorAxis;
end

%The equation for a reference ellipsoid is
%(x^2+y^2)/a^2+z^2/(a^2*(1-f)^2)=1.
%The algorithm describe in Knuth's book wants things paramterized in terms
%of c(1)*x^2+c(2)*y^2+c(3)*z^2=1. Thus, the coefficients are
c(1,1)=1/a^2;
c(2,1)=1/a^2;
c(3,1)=1/(a^2*(1-f)^2);

%The c's have to be sorted. This means that the ordering of the components
%in x, y, and z will be wrong, so indices to undo the sorting when saving
%the final result need to be found.
[c,sortIdx]=sort(c,'descend');
unsortIdx=1:3;
unsortIdx=unsortIdx(sortIdx);

n=3;%The number of dimensions.

%Allocate space for the return variables.
CartPoint=zeros(3,N);

for curSamp=1:N
    while(1)
        %First, generate a random point on the unit sphere
        y=randDirVec(3);
        rho=sqrt(sum(c.*y.^2));

        %Find K
        if(n*c(n)>=c(1))
            K=sqrt(c(n)^(n-1));
        else
            K=sqrt(((n+1)/(c(1)+c(n)))^(n+1)*(c(1)*c(n)/n)^n);
        end

        %Generate a uniform random variable U
        U=rand(1);

        if(rho^(n+1)*U<K*sqrt(sum(c.^2.*y.^2)))
            CartPoint(unsortIdx,curSamp)=y/rho;
            break;
        end
    end
end

if(nargout>1)
    point=Cart2Ellipse(CartPoint,[],a,f);
    latLong=point(1:2,:);
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
