function [xi,w]=thirdOrderNDimShellCubPoints(numDim,K1,K2,algorithm)
%%THIRDORDERNDIMCHELLCUBPOINTS Generate third-order cubature points for
%               integration over a numDim-dimensional cubic shell. The
%               shell is defined as the region K1<=sum(abs(x))<=K2.
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated. numDim>1.
%         K1,K2 The upper and lower bounds defining the shell 0<K1<K2<Inf.
%     algorithm An optional parameter specifying the algorithm to be
%               used to generate the points. Possible values are:
%               0 (The default if omitted or an empty matrix is passed)
%                 Formula C_n^{shell} 3-1 in [1], pg. 266, 2*numDim points.
%               1 Formula C_n^{shell} 3-3 in [1], pg. 267, 2^(numDim+1)
%                 points.
%               2 Formula C_n^{shell} 3-4 in [1], pg. 267, numDim*2^numDim
%                 points.
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
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(algorithm))
    algorithm=0; 
end

n=numDim;

switch(algorithm)
    case 0%C_n^{shell} 3-1 in [1], pg. 266, 2*numDim points.
        r=sqrt(n*(K2^(n+2)-K1^(n+2))/(3*(K2^n-K1^n)));
        xi=fullSymPerms([r;zeros(n-1,1)]);

        V=2^n*(K2^n-K1^n);
        B=(1/(2*n))*V;

        w=B*ones(2*n,1);
    case 1%C_n^{shell} 3-3 in [1], pg. 267, 2^(numDim+1) points.
        xi=[PMCombos(K1*ones(n,1)),PMCombos(K2*ones(n,1))];
        B1=(2*K2^(n+2)-3*K2^2*K1^n+K1^(n+2))/(3*(K2^2-K1^2));
        B2=(2*K1^(n+2)-3*K1^2*K2^n+K2^(n+2))/(3*(K2^2-K1^2));
        
        w=[B1*ones(2^n,1);B2*ones(2^n,1)];
    case 2%C_n^{shell} 3-4 in [1], pg. 267, numDim*2^numDim points.
        u=sqrt(n*(K2^(n+2)-K1^(n+2))/((n+2)*(K2^n-K1^n)));
        v=sqrt(1/3)*u;
        
        xi=fullSymPerms([u;v*ones(n-1,1)]);
        
        V=2^n*(K2^n-K1^n);
        B1=V/(n*2^n);
        w=B1*ones(n*2^n,1);
    otherwise
        error('Unknown algorithm specified.')
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
