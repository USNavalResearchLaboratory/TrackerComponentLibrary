function [xi,w]=secondOrderNDimCubPoints(numDim,algorithm)
%%SECONDORDERNDIMCUBPOINTS Generate first-order cubature points for
%               integration over a numDim-dimensional cube with bounds in
%               coordinates of (-1,1).
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated. numDim>1.
%     algorithm An optional parameter specifying the algorithm to be used
%               to generate the points. Possible values are:
%               0 (The default if omitted or an empty matrix is passed)
%                 Formula Cn 2-1 in [1], pg 229, numDim+1 points.
%               1 Formula Cn 2-2 in [1], pg. 230, 2*numDim+1 points.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points (Each "point" is a vector).
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    algorithm=0; 
end

switch(algorithm)
    case 0%Cn 2-1 in [1], pg 229, numDim+1 points.
        V=2^numDim;
        
        w=(1/(numDim+1))*ones(numDim+1,1)*V;
        xi=zeros(numDim,numDim+1);
        i=0:numDim;
        for k=1:fix(numDim/2)
            xi(2*k-1,i+1)=sqrt(2/3)*cos(2*i*k*pi/(numDim+1));
            xi(2*k,i+1)=sqrt(2/3)*sin(2*i*k*pi/(numDim+1));
        end
        
        if(mod(numDim,2)~=0)
            xi(numDim,i+1)=(-1).^i/sqrt(3);
        end
    case 1%Cn 2-2 in [1], pg. 230, 2*numDim+1 points.
        r=sqrt(3)/6;
        V=2^numDim;
        B1=V;
        B2=-r*V;
        B3=r*V;
    
        xi=[2*r*ones(numDim,1),genAllMultisetPermutations([1;r*ones(numDim-1,1)]),genAllMultisetPermutations([-1;r*ones(numDim-1,1)])];
        w=[B1;B2*ones(numDim,1);B3*ones(numDim,1)];
    otherwise
        error('Unknown algorithm specified');
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
