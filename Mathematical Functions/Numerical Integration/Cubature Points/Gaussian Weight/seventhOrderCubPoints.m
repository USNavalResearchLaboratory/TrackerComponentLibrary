function [xi,w]=seventhOrderCubPoints(numDim,algorithm,randomize)
%%SEVENTHORDERCUBPOINTS Generate seventh order cubature points for
%                       integration involving a multidimensional Gaussian
%                       probability density function (PDF). Note that if
%                       numDim=8, then 32 of the points will have zero
%                       weights (and are thus eliminated) and if numDim>8,
%                       then some points will have negative weights.
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated; numDim>=1.
%     algorithm An optional parameter specifying the algorithm to be used
%               to generate the points. Possible values are
%               0 (The default if omitted or an empty matrix is passed and
%                 numDim~=2 and numDim~=1) Use  algorithm E_n^{r^2} 7-3 on
%                 page 319 of [1], requiring 2*(2^numDim+2*numDim^2) points.
%                 Note, however, that the formula in the book contains a
%                 typo. Using [2] and [3], one can get the correct
%                 formulae, which are summarized in [4].
%               1 E_n^{r^2} 7-1 on page 318 of [1], requires
%                 2^numDim+2*numDim^2+1 points, numDim=3,4,6,7.
%               2 (The default if numDim==2) E_2^{r^2} 7-1, page 324 of
%                 [1], 12 points, numDim=2.
%               3 E_2^{r^2} 7-2, page 324 of [1], 16 points, numDim=2. This
%                 includes a correction as a scale factor of (4/3) was
%                 missing from the squared expressions for r and s.
%               4 E_3^{r^2} 7-1, page 327 of [1], 27 points, numDim=3,
%                 first variant using the upper signs.
%               5 E_3^{r^2} 7-1, page 327 of [1], 27 points, numDim=3,
%                 second variant using the lower signs.
%               6 E_3^{r^2} 7-2, page 328 of [1], 33 points, numDim=3,
%                 first variant using the upper signs.
%               7 E_3^{r^2} 7-2, page 328 of [1], 33 points, numDim=3,
%                 second variant using the lower signs.
%               8 E_4^{r^2} 7-1, page 329 of [1], 49 points, numDim=4,
%                 first variant using the upper signs, corrected by
%                 multiplying s,t, and r by a factor of sqrt(4/5).
%               9 (The default if numDim==1) Call quadraturePoints1D(4), 4
%                  points, numDim=1.
%     randomize If this parameter is true, then the points will be
%               multiplied by a random orthonormal rotation matrix. This
%               does not change the moments up to the order of the points.
%               This randomization is done in [5] and [6] to lessen various
%               effects that arise when using points in the same
%               orientation repeatedly in tracking. The default if this
%               parameter is omitted or an empty matrix is passed is false.
%               This parameter is ignored if numDim==1.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%%For more details on how to use these points, see the comments in the
%function fifthOrderCubPoints.m.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%[2] A. H. Stroud, "Some seventh degree integration formulas for symmetric
%    regions," SIAM Journal on Numerical Analysis, vol. 4, no. 1, pp.
%    37-44, Mar. 1967.
%[3] A. H. Stroud, "Some seventh degree integration formulas for the
%    surface of an n-sphere," Numerische Mathematik, vol. 11, no. 3, pp.
%    273-276, Mar. 1968.
%[4] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%[5] O. Straka, D. Duník, and M. Simandl, "Randomized unscented Kalman
%    filter in tracking," in Proceedings of the 15th International
%    Conference on Information Fusion, Singapore, 9-12 Jul. 2012, pp.
%    503-510.
%[6] J. Duník, O. Straka, and M. Simandl, "The development of a randomised
%    unscented Kalman filter," in Proceedings of the 18th World Congress,
%    The International Federation of Automatic Control, Milan, Italy, 28
%    Aug. - 2 Sep. 2011, pp. 8-13.
%
%July 2012 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    n=numDim;
    
    if(nargin<2||isempty(algorithm))
        if(n>2)
            algorithm=0;
        elseif(n==2)
            algorithm=2;
        else%Assume %n==1
            algorithm=9;
        end
    end
    
    if(nargin<3||isempty(randomize))
        randomize=false;
    end
    
    switch(algorithm)
        case 0%Algorithm E_n^{r^2} 7-1 on page 318 of [1], requiring
              %2^numDim+2*numDim^2+1 points.
            if(n<3)
               error('The selected algorithm requires numDim>=3') 
            end
            
            %First, we will generate points for the n-sphere.
            [u,wu]=seventhOrderSpherSurfCubPoints(n,0);
            VSphere=2*pi^(numDim/2)/gamma(numDim/2);
            wu=wu/VSphere;
            numUPoints=size(u,2);

            %Finally, all of this must be transformed as described in Stroud's
            %book to make it work with the desired Gaussian weighting function.
            r1=sqrt((n+2-sqrt(2*(n+2)))/2);
            r2=sqrt((n+2+sqrt(2*(n+2)))/2);

            %The denominators were wrong in the paper if we want them to sum to
            %1. They are correct here.
            A(1)=(n+2+sqrt(2*(n+2)))/(2*(n+2));
            A(2)=(n+2-sqrt(2*(n+2)))/(2*(n+2));

            w=zeros(2*numUPoints,1);
            xi=zeros(n,2*numUPoints);

            w(1:numUPoints)=wu*A(1);
            w(numUPoints+(1:numUPoints))=wu*A(2);

            xi(:,1:numUPoints)=r1*u;
            xi(:,numUPoints+(1:numUPoints))=r2*u;

            %To make it approximate the standard normal distribution.
            xi=sqrt(2)*xi;

            if(numDim==8)
                %Get rid of the zero weights when numDim==8.
                sel=(w~=0);
                w=w(sel);
                xi=xi(:,sel);
            end
        case 1%E_n^{r^2} 7-1 on page 318 of [1], 2^numDim+2*numDim^2+1 
              %points, numDim=3,4,6,7.
            if(n~=3&&n~=4&&n~=6&&n~=7)
               error('The selected algorithm only supports numDim=3, 4, 6,or 7.')
            end
            
            r=sqrt((3*(8-n)-(n-2)*sqrt(3*(8-n)))/(2*(5-n)));
            s=sqrt((3*n-2*sqrt(3*(8-n)))/(2*(3*n-8)));
            t=sqrt((6+sqrt(3*(8-n)))/2);

            B=(8-n)/(8*r^6);
            C=1/(2^(n+3)*s^6);
            D=1/(16*t^6);
            A=1-2*n*B-2^n*C-2*n*(n-1)*D;
            
            %The sqrt(2) makes it for the normal 0-I distribution.
            xi=sqrt(2)*[zeros(n,1),fullSymPerms([r;zeros(numDim-1,1)]),PMCombos(s*ones(n,1)),fullSymPerms([t;t;zeros(numDim-2,1)])];
            w=[A;B*ones(2*n,1);C*ones(2^n,1);D*ones(2*n*(n-1),1)];
        case 2%E_2^{r^2} 7-1, page 324 of [1], 12 points, numDim=2.
            if(n~=2)
               error('The selected algorithm only supports numDim=2.')
            end
            
            r=sqrt(3);
            s=sqrt((9-3*sqrt(5))/8);
            t=sqrt((9+3*sqrt(5))/8);
            A=1/36;
            B=(5+2*sqrt(5))/45;
            C=(5-2*sqrt(5))/45;
            
            %The sqrt(2) makes it for the normal 0-I distribution.
            xi=sqrt(2)*[fullSymPerms([r;0]),PMCombos([s;s]),PMCombos([t;t])];
            w=[A*ones(4,1);B*ones(4,1);C*ones(4,1)];
        case 3%E_2^{r^2} 7-2, page 324 of [1], 16 points, numDim=2. This
              %includes a correction as a scale factor of (4/3) was missing
              %from the squared expressions for r and s.
            if(n~=2)
               error('The selected algorithm only supports numDim=2.')
            end
            
            r=sqrt((4/3)*(3+sqrt(6))/2);
            s=sqrt((4/3)*(3-sqrt(6))/2);
            
            A=(5-2*sqrt(6))/48;
            B=(5+2*sqrt(6))/48;
            C=1/48;
            
            %The sqrt(2) makes it for the normal 0-I distribution.
            xi=sqrt(2)*[fullSymPerms([r;0]),fullSymPerms([s;0]),fullSymPerms([r;s])];
            w=[A*ones(4,1);B*ones(4,1);C*ones(8,1)];
        case 4%E_3^{r^2} 7-1, page 327 of [1], 27 points, numDim=3, first
              %variant using the upper signs.
            if(n~=3)
               error('The selected algorithm only supports numDim=3.')
            end
            
            r=sqrt((15+sqrt(15))/4);
            s=sqrt((6-sqrt(15))/2);
            t=sqrt((9+2*sqrt(15))/2);
            A=(720+8*sqrt(15))/2205;
            B=(270-46*sqrt(15))/15435;
            C=(162+41*sqrt(15))/6174;
            D=(783-202*sqrt(15))/24696;
            
            %The sqrt(2) makes it for the normal 0-I distribution.
            xi=sqrt(2)*[zeros(3,1),fullSymPerms([r;0;0]),fullSymPerms([s;s;0]),PMCombos([t;t;t])];
            w=[A;B*ones(6,1);C*ones(12,1);D*ones(8,1)];
        case 5%E_3^{r^2} 7-1, page 327 of [1], 27 points, numDim=3, second
              %variant using the lower signs.         
            if(n~=3)
               error('The selected algorithm only supports numDim=3.')
            end
            
            r=sqrt((15-sqrt(15))/4);
            s=sqrt((6+sqrt(15))/2);
            t=sqrt((9-2*sqrt(15))/2);
            A=(720-8*sqrt(15))/2205;
            B=(270+46*sqrt(15))/15435;
            C=(162-41*sqrt(15))/6174;
            D=(783+202*sqrt(15))/24696;
            
            %The sqrt(2) makes it for the normal 0-I distribution.
            xi=sqrt(2)*[zeros(3,1),fullSymPerms([r;0;0]),fullSymPerms([s;s;0]),PMCombos([t;t;t])];
            w=[A;B*ones(6,1);C*ones(12,1);D*ones(8,1)];  
        case 6%E_3^{r^2} 7-2, page 328 of [1], 33 points, numDim=3, first
              %variant using the upper signs.
            if(n~=3)
               error('The selected algorithm only supports numDim=3.')
            end
 
            r=sqrt((25+15*sqrt(2)+5*sqrt(5)+3*sqrt(10))/4);
            s=sqrt((25+15*sqrt(2)-5*sqrt(5)-3*sqrt(10))/4);
            t=sqrt((3-sqrt(2))/2);
            u=sqrt((9-3*sqrt(2)-3*sqrt(5)+sqrt(10))/4);
            v=sqrt((9-3*sqrt(2)+3*sqrt(5)-sqrt(10))/4);

            A=(80+8*sqrt(2))/245;
            B=(395-279*sqrt(2))/13720;
            C=(45+29*sqrt(2))/2744;

            %The sqrt(2) makes it for the normal 0-I distribution.
            xi=sqrt(2)*[zeros(3,1),PMCombos([r;s;0]),PMCombos([0;r;s]),PMCombos([s;0;r]),PMCombos([u;v;0]),PMCombos([0;u;v]),PMCombos([v;0;u]),PMCombos([t;t;t])];
            w=[A;B*ones(12,1);C*ones(20,1)];
        case 7%E_3^{r^2} 7-2, page 328 of [1], 33 points, numDim=3, second
              %variant using the lower signs.
            if(n~=3)
               error('The selected algorithm only supports numDim=3.')
            end
            r=sqrt((25-15*sqrt(2)+5*sqrt(5)-3*sqrt(10))/4);
            s=sqrt((25-15*sqrt(2)-5*sqrt(5)+3*sqrt(10))/4);
            t=sqrt((3+sqrt(2))/2);
            u=sqrt((9+3*sqrt(2)-3*sqrt(5)-sqrt(10))/4);
            v=sqrt((9+3*sqrt(2)+3*sqrt(5)+sqrt(10))/4);

            A=(80-8*sqrt(2))/245;
            B=(395+279*sqrt(2))/13720;
            C=(45-29*sqrt(2))/2744;

            %The sqrt(2) makes it for the normal 0-I distribution.
            xi=sqrt(2)*[zeros(3,1),PMCombos([r;s;0]),PMCombos([0;r;s]),PMCombos([s;0;r]),PMCombos([u;v;0]),PMCombos([0;u;v]),PMCombos([v;0;u]),PMCombos([t;t;t])];
            w=[A;B*ones(12,1);C*ones(20,1)];
        case 8%E_4^{r^2} 7-1, page 329 of [1], 49 points, numDim=4, first
              %variant using the upper signs, corrected by multiplying s,t,
              %and r by a factor of sqrt(4/5).
            if(n~=4)
               error('The selected algorithm only supports numDim=4.')
            end
            s=sqrt((3-sqrt(3))/2);
            t=(3+sqrt(3));
            r=2*s;
            
            %Correction from paper.
            s=sqrt(4/5)*s;
            t=sqrt(4/5)*t;
            r=sqrt(4/5)*r;
            
            A=1/4;
            B=(9+5*sqrt(3))/576;
            C=(9-5*sqrt(3))/576;
            
            %The sqrt(2) makes it for the normal 0-I distribution.
            xi=sqrt(2)*[zeros(4,1),fullSymPerms([r;0;0;0]),PMCombos([s;s;s;s]),fullSymPerms([t;t;0;0])];
            w=[A;B*ones(24,1);C*ones(24,1)];
        case 9%Use quadraturePoints1D(4) for 1D points
            if(n~=1)
                error('This algorithm requires numDim=1.')
            end
            [xi,w]=quadraturePoints1D(4);
        otherwise
            error('Unknown algorithm selected')
    end
    
    if(numDim>1&&randomize)
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
