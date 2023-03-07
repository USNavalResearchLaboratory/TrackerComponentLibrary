function [xi,w]=arbOrderGaussCubPoints(numDim,n,beta,algorithm,randomize)
%%ARBORDERGAUSSCUBPOINTS Generate cubature points for integrating with a
%               weighting function of a standard 0-I multivariate Gaussian
%               PDF times |x|^beta. That is,
%               w(x)=1/(2*pi)^(numDim/2)*norm(x)^beta*exp(-x'*x/2).
%               For beta=0, the weighting function is just the standard
%               normal PDF. For 1D points, one should just use the
%               quadraturePoints1D function.
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated. numDim>1.
%             n A positive integer such that 2n-1 is the highest degree to
%               which the cubature  points are accurate.
%          beta An integer specifying the exponent of the norm(x) term in
%               the weighting function (the thing times the normal PDF in
%               the weighting function). beta>-numDim. If omitted or an
%               empty matrix is passed, beta=0 is used.
%     algorithm An optional parameter specifying which algorithm to use to
%               generate arbitrary order Gaussian cubature points. Possible
%               values are:
%               0 (The default if omitted or an empty matrix is passed) Use
%                 the algorithm described in Chapter 2.6 of [1]. The
%                 specialization for Gaussian problems starts near Equation
%                 2.6-16. The transformation to a Gaussian PDF is performed
%                 by scaling the points and weights.
%               1 Combine the functions arbOrderSpherSurfCubPoints and
%                 spherSurfPoints2GaussPoints to obtain arbitrary order
%                 points.
%               2 Use the tensor product rule of linCubPoints2MultiDim
%                 applied to 1D points. This option only supports beta=0. 
%               3 Use Smolyak's algorithm via linCubPoints2MultiDim
%                 applied to 1D points. This option only supports beta=0.
%     randomize If this parameter is true, then the points will be
%               multiplied by a random orthonormal rotation matrix. This
%               does not change the moments up to the order of the points.
%               This randomization is done in [2] and [3] to lessen various
%               effects that arise when using points in the same
%               orientation repeatedly in tracking. The default if this
%               parameter is omitted or an empty matrix is passed is false.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%Algorithm 0 is as described in Chapter 2.6 of [1]. The specialization
%for Gaussian problems starts near Equation 2.6-16. The transformation to
%a Gaussian PDF is performed by scaling the points and weights.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%[2] O. Straka, D. Duník, and M. Simandl, "Randomized unscented Kalman
%    filter in tracking," in Proceedings of the 15th International
%    Conference on Information Fusion, Singapore, 9-12 Jul. 2012, pp.
%    503-510.
%[3] J. Duník, O. Straka, and M. Simandl, "The development of a randomised
%    unscented Kalman filter," in Proceedings of the 18th World Congress,
%    The International Federation of Automatic Control, Milan, Italy, 28
%    Aug. - 2 Sep. 2011, pp. 8-13.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

if(nargin<5||isempty(randomize))
    randomize=false;
end

if(nargin<4||isempty(algorithm))
   algorithm=0; 
end

if(nargin<3||isempty(beta))
   beta=0; 
end

if(numDim==1)
    error('numDim must be >1.') 
end

switch(algorithm)
    case 0
        [r,Ar]=quadraturePoints1D(n,9,numDim-1+beta);
        numRVals=length(r);

        %Allocate space for the second set of points.
        numGegenbauerPoints=n;
        y=zeros(numGegenbauerPoints,numDim-1);
        Ak=zeros(numGegenbauerPoints,numDim-1);

        %First, find all of the yk and Ak values
        for k=1:(numDim-1)
            [yCur,AkCur]=quadraturePoints1D(n,3,(k-1)/2);

            y(:,k)=yCur(:);
            Ak(:,k)=AkCur;
        end

        %Allocate space
        numPoints=numRVals*numGegenbauerPoints;
        xi=zeros(numDim,numPoints);
        w=zeros(numPoints,1);

        dims=numGegenbauerPoints*ones(numDim-1,1);
        curCubPoint=1;
        %Go through all of the ranges
        for curR=1:numRVals
            curIdx=1;
            i=index2NDim(dims,curIdx);

            %Now, we go through all the tuples
            while(~isempty(i))
                xi(numDim,curCubPoint)=r(curR)*y(i(numDim-1),numDim-1);
                wProd=Ak(i(numDim-1),numDim-1);

                prodRecur=1;
                for k=(numDim-1):-1:2
                    prodRecur=prodRecur*sqrt(1-y(i(k),k)^2);
                    xi(k,curCubPoint)=r(curR)*prodRecur*y(i(k-1),k-1);
                    wProd=wProd*Ak(i(k-1),k-1);
                end
                prodRecur=prodRecur*sqrt(1-y(i(1),1)^2);
                xi(1,curCubPoint)=r(curR)*prodRecur;
                w(curCubPoint)=wProd*Ar(curR);

                curCubPoint=curCubPoint+1;
                curIdx=curIdx+1;
                i=index2NDim(dims,curIdx);
            end
        end

        %We now have cubature points and weights for integrating with a
        %weighting function of w(x)=|x|^beta*exp(-x'*x). We will transform
        %these to integrating over the normal 0-I distribution times
        %|x|^beta.
        w=pi^(-numDim/2)*sqrt(2)^(beta)*w;
        xi=sqrt(2)*xi;
    case 1
    %Obtain arbitrary order spherical surface points and turn them into
    %points for integrating over a Gaussian distribution using the
    %spherSurfPoints2GaussPoints function.
        [xiSurf,wSurf]=arbOrderSpherSurfCubPoints(numDim,n);
        [xi,w]=spherSurfPoints2GaussPoints(xiSurf,wSurf,2*n-1,beta);
    case 2
        if(beta~=0)
           error('This algorithmic option is not supported with beta~=0') 
        end
        [xi,w]=linCubPoints2MultiDim(numDim,n,0);
    case 3
        if(beta~=0)
           error('This algorithmic option is not supported with beta~=0') 
        end
        [xi,w]=linCubPoints2MultiDim(numDim,n,1);
    otherwise
        error('Unknown algorithm specified')
end

if(randomize)
    R=randOrthoMat(numDim);

    numPoints=length(w);
    for curPoint=1:numPoints
        xi(:,curPoint)=R*xi(:,curPoint);
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
