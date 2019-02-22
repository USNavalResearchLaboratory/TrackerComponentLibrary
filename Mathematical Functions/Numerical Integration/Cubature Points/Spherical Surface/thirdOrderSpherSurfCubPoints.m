function [xi,w]=thirdOrderSpherSurfCubPoints(numDim,algorithm)
%%THIRDORDERSPHERSURFCUBPOINTS Obtain third-order cubature points for
%                   integration over the surface of the unit sphere
%                   (weighting function is just 1).
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated. numDim>=1.
%     algorithm An optional parameter specifying the algorithm to be
%               used to generate the points. Possible values are
%               0 (The default if omitted or an empty matrix is passed)
%                 Use algorithm Un 3-1 on page 294 of [1]. This requires
%                 2*numDim points.
%               1 Use algorithm Un 3-2 on page 294 of [1]. This requires
%                 2^numDim points.
%               2 Use the degree 3 algorithm of [2] requiring 2*(numDim+1)
%                 points. The algorithm is also in [3], which is easier to
%                 get than [2].
%               3 Algorithm U3 3-1 of [1], requiring 12 points and that
%                 numDim=3.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%REFERENCES:
%[1] A. H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%[2] I. P. Mysovskikh, The approximation of multiple integrals by using
%    interpolatory cubature formulae, in Quantitative Approximation, R. A.
%    DeVore and K. Scherer, eds., Academic Press, New York, 1980, pp.
%    217-243.
%[3] A. Genz and J. Monahan, "Stochastic integration rules for infinite
%    regions," SIAM Journal on Scientific Computing, vol. 19, no. 2, pp.
%    426-439, Mar. 1998.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    algorithm=0;
end

switch(algorithm)
    case 0%Stroud's method, 3-1, 2*n points.
        absUm=2*pi^(numDim/2)/gamma(numDim/2);
        w=repmat(absUm/(2*numDim),[2*numDim,1]);
        xi=[eye(numDim,numDim),-eye(numDim,numDim)];
    case 1%Strouds method, 3-2 2^n points
        absUm=2*pi^(numDim/2)/gamma(numDim/2);
        w=repmat(absUm*2^(-numDim),[2^numDim,1]);
        
        r=1/sqrt(numDim);
        xi=PMCombos(repmat(r,[numDim,1]));
    case 2%Mysovskikh's method
        absUm=2*pi^(numDim/2)/gamma(numDim/2);
        v=regularNSimplexCoords(numDim);
        xi=[v,-v];
        w=repmat(absUm/(2*(numDim+1)),[2*(numDim+1),1]);
    case 3%Algorithm U3 3-1 of [1], requiring 12 points
        if(numDim~=3)
           error('This formula requires that numDim=3') 
        end
        V=2*pi^(3/2)/gamma(3/2);
        
        r=1/sqrt(2);
        xi=fullSymPerms([r;r;0]);
        w=(V/12)*ones(12,1);
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
