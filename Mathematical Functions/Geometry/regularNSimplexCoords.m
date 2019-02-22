function points=regularNSimplexCoords(numDim,method)
%%REGULARNSIMPLEXCOORDS Find the coordinates of a regular n-dimensional
%               simplex. That is, a hypertetrahedron in n-dimensions.
%               All sides in the simplex have equal length, the centroid of
%               the simplex is zero, and the distance of all of the points
%               in the simplex from the origin is 1.
%
%INPUTS: numDim The number of dimensions of the simplex. numDim>=1.
%        method The simplex is not unique (it can be rotated). If the
%               method parameter is specified, it selects which method of
%               finding a simplex is used. Possible values are
%               0 (The default if omitted or an empty matrix is passed. Use
%                  the method of [1], which utilizes cosines, is used.
%               1 The method of [1] utilizing square roots instead of
%                 cosines is used.
%
%OUTPUTS: points numDimXnumPoints set of points that form a
%                numDim-simplex. numPoints is always equal to numDim+1.
%
%The two methods given are from Chapter 9 of [1]. The methods do not
%actually provide points that are a unit distance from the origin, though
%they are all equidistant from the origin. Thus, this scales the points to
%be a unit distance form the origin.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(method))
    method=0;
end

switch(method)
    case 0%The method in [1] utilizing cosines
        points=zeros(numDim,numDim+1);

        r=1:fix(numDim/2);

        for k=0:numDim
            points(2*r-1,k+1)=cos(2*r*k*pi/(numDim+1));
            points(2*r,k+1)=sin(2*r*k*pi/(numDim+1));
            %For an odd number of dimensions
            if(mod(numDim,2)~=0)
                points(numDim,k+1)=(-1)^k/sqrt(2);
            end
        end
        points=points/sqrt(numDim/2);
    case 1%The method in [1] utilizing square roots.
        points=zeros(numDim,numDim+1);
        for k=1:(numDim+1)
            i=1:k;
            points(i,k)=-sqrt((numDim+1)./((numDim-i+2).*(numDim-i+1)));
            points(k,k)=sqrt((numDim+1)*(numDim-k+1)/(numDim-k+2));
            %Components for i>k are zero.
        end
        points=points/sqrt(numDim);
    otherwise
       error('Invalid method specified') 
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
