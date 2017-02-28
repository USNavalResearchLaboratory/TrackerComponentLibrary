function [xi, w]=secondOrderCubPoints(numDim,w0)
%SECONDDORDERCUBPOINTS Generate second order cubature points for 
%           integration involving a multidimensional Gaussian probability
%           density function (PDF).
%
%INPUTS:    numDim  An integer specifying the dimensionality of the points
%                   to be generated.
%               w0 An optional input parameter. This is the value of w(1),
%                  the coefficient associated with a cubature point placed
%                  at the origin. It is required that 1>w0>0. If omitted or
%                  an empty matrix is passed, a default value of 1/3 is
%                  used.
%
%OUTPUTS:   xi   A numDim X numCubaturePoints matrix containing the
%                cubature points. (Each "point" is a vector)
%           w    A numCubaturePoints X 1 vector of the weights
%                 associated with the cubature points.
%
%The algorithm implemented is for the spherical simplex sigma points from
%Appendix III of [1]. It is modified because all of the quadrature points,
%as described in the paper, need to be multiplied by
%(1/sqrt((1/w0-1)/(numDim+1))) to assure that the sample covariance matrix
%of the points has a unit diagonal.
%
%REFERENCES:
%[1] S. J. Julier and J. K. Uhlmann, "Unscented filtering and nonlinear
%    estimation," Proceedings of the IEEE, vol. 92, no. 3, pp. 401-422,
%    Mar. 2004.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<2||isempty(w0))
       w0=1/3; 
    end

    numPoints=numDim+2;
    w=zeros(numPoints,1);
    xi=zeros(numDim,numPoints);
    w(1)=w0;
    w(2:end)=(1-w0)/(numDim+1);
    
    i=0;
    xi(1,i+1)=0;
    i=1;
    xi(1,i+1)=-1/sqrt(2*w(1));
    i=2;
    xi(1,i+1)=1/sqrt(2*w(1));

    for j=2:numDim
        for i=1:j
            xi(j,i+1)=-1/sqrt(j*(j+1)*w(1));
        end
        
        i=j+1;
        xi(j,i+1)=j/sqrt(j*(j+1)*w(1));
    end

    xi=xi/sqrt((1/w0-1)/(numDim+1));
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
