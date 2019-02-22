function [xi,w]=fifthOrderSpherCubPoints(numDim,algorithm,alpha)
%%FIFTHORDERSPHERCUBPOINTS Generate fifth-order cubature points for
%               integration over a unit spherical (or hyperspherical)
%               region. The weighting function can either be 1 or
%               sum(x.^2)^(-alpha/2) for alpha>-numDim, depending on the
%               chosen algorithm.
%
%INPUTS:  numDim An integer specifying the dimensionality of the points
%                to be generated.
%      algorithm A value indicating which algorithm should be used.
%                Possible values are
%                0 (The default if omitted or an empty matrix is passed)
%                  Formula Sn 5-2 in [1], pg. 269, 2*numDim^2+1 points.
%                1 Formula Sn 5-3 in [1], pg. 270, 2^numDim+2*numDim
%                  points, alpha=0, numDim>=4 for all points to be internal
%                  to the unit sphere
%                2 Formula Sn 5-4 in [1], pg. 270, 2^(n+1)-1 points,
%                  numDim>=2 for all points inside the unit sphere.
%                3 Formula Sn 5-5 in [1], pg. 270, numDim*2^numDim+1
%                  points, alpha=0, numDim>=2 for all points within the
%                  unit sphere.
%                4 Formula Sn 5-6 in [1], pg. 271, 2^numDim*(numDim+1)
%                  points, alpha=0, numDim>=2 for all points within the
%                  unit sphere.
%                5 Formula S2 5-1 in [1], pg. 279, 7 points, numDim=2, with
%                  corrections for a missing +/- term in front of the t in
%                  the third set of points and a missing V in the
%                  expression for B.
%                6 Formula S2 5-2 in [1], pg. 280, 9 points, numDim=2,
%                  alpha=0.
%                7 Formula S3 5-1 in [1], pg. 289, 13 points, numDim=3.
%                8 Formula S3 5-2 in [1], pg. 290, 21 points, numDim=3,
%                  alpha=0.
%          alpha A parameter specifying the exponent of the weighting
%                function over the sphere. The weighting function is 
%                sum(x.^2)^(-alpha/2). If omitted or an empty matrix is
%                passed, alpha=0 is used. Not all algorithms support
%                nonzero weighting functions.
%
%OUTPUTS:   xi      A numDim X numCubaturePoints matrix containing the
%                   cubature points. (Each "point" is a vector)
%           w       A numCubaturePoints X 1 vector of the weights
%                   associated with the cubature points.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

if(nargin<2||isempty(algorithm))
   algorithm=0; 
end

if(nargin<3||isempty(alpha))
   alpha=0; 
end

switch(algorithm)
    case 0%Sn 5-2 in [1], pg. 269, 2*numDim^2+1 points.
        V=2/(numDim+alpha)*(pi^(numDim/2)/gamma(numDim/2));
        r=sqrt(3*(numDim+alpha+2)/((numDim+2)*(numDim+alpha+4)));
        
        B2=V*(numDim+2)*(numDim+alpha)*(numDim+alpha+4)/(36*numDim*(numDim+alpha+2)^2);
        B1=V*(4-numDim)*(numDim+2)*(numDim+alpha)*(numDim+alpha+4)/(18*numDim*(numDim+alpha+2)^2);
        B0=V-2*numDim*B1-2*numDim*(numDim-1)*B2;
        
        xi=[zeros(numDim,1),fullSymPerms([r;zeros(numDim-1,1)]),fullSymPerms([r;r;zeros(numDim-2,1)])];
        w=[B0;B1*ones(2*numDim,1);B2*ones(2*numDim*(numDim-1),1)];
    case 1%Sn 5-3 in [1], pg. 270, 2^numDim+2*numDim points, alpha=0,
          %numDim>=4 for all pionts to be internal to the unit sphere.
        if(alpha~=0)
            error('This algorithm requires that alpha=0')
        end
        
        r=sqrt((numDim+4-sqrt(2*(numDim+4)))/(numDim+4));
        s=sqrt((numDim*(numDim+4)+2*sqrt(2*(numDim+4)))/((numDim^2+2*numDim-4)*(numDim+4)));
        V=(2/numDim)*(pi^(numDim/2)/gamma(numDim/2));
        
        B1=V/((numDim+2)*(numDim+4)*r^4);
        B2=V/(2^numDim*(numDim+2)*(numDim+4)*s^4);
        
        xi=[fullSymPerms([r;zeros(numDim-1,1)]),PMCombos(s*ones(numDim,1))];
        w=[B1*ones(2*numDim,1);B2*ones(2^numDim,1)];
    case 2%Sn 5-4 in [1], pg. 270, 2^(n+1)-1 points, numDim>=2 for all
          %points inside the unit sphere.
        s=sqrt((numDim+alpha+2)/((numDim+2)*(numDim+alpha+4)));
        V=2/(numDim+alpha)*(pi^(numDim/2)/gamma(numDim/2));
        
        numPoints=2^(numDim+1)-1;
        xi=zeros(numDim,numPoints);
        w=zeros(numPoints,1);
        
        vec=s*ones(numDim,1);
        BPowSum=0;
        curStart=1;
        for k=1:numDim
            r=sqrt((k+2)*(numDim+alpha+2)/((numDim+2)*(numDim+alpha+4)));
            B=V*2^(k-numDim)*(numDim+2)*(numDim+alpha)*(numDim+alpha+4)/(numDim*(k+1)*(k+2)*(numDim+alpha+2)^2);
            
            vec(k)=r;
            xiCur=PMCombos(vec(k:end));
            num2Add=size(xiCur,2);
            xi(k:end,curStart:(curStart+num2Add-1))=xiCur;
            w(curStart:(curStart+num2Add-1))=B;
            curStart=curStart+num2Add;
            
            BPowSum=BPowSum+2^(numDim-k+1)*B;
        end
        
        if(alpha==0)
            w(end)=4*V/(numDim+2)^2;
        else
            w(end)=V-BPowSum;
        end
    case 3%Sn 5-5 in [1], pg. 270, numDim*2^numDim+1 points, alpha=0,
          %numDim>=2 for all points within the unit sphere.
        if(alpha~=0)
            error('This algorithm requires that alpha=0')
        end
        V=(2/numDim)*(pi^(numDim/2)/gamma(numDim/2));
        
        r=sqrt((numDim+2+(numDim-1)*sqrt(2*(numDim+2)))/(numDim*(numDim+4)));
        s=sqrt((numDim+2-sqrt(2*(numDim+2)))/(numDim*(numDim+4)));
        B0=4*V/(numDim+2)^2;
        B1=V*(numDim+4)/(2^numDim*(numDim+2)^2);
        
        xi=[zeros(numDim,1),fullSymPerms([r;s*ones(numDim-1,1)])];
        w=[B0;B1*ones(numDim*2^numDim,1)];
    case 4%Sn 5-6 in [1], pg. 271, 2^numDim*(numDim+1) points, alpha=0,
          %numDim>=2 for all points within the unit sphere.
        if(alpha~=0)
            error('This algorithm requires that alpha=0')
        end
        V=(2/numDim)*(pi^(numDim/2)/gamma(numDim/2));
        
        n=numDim;
        r=sqrt((n*(n+4)+2*sqrt(n+4)+(n-1)*sqrt(2*(n+1)*(n+2)*(n+4)))/(n*(n+2)*(n+4)));
        s=sqrt((n*(n+4)+2*sqrt(n+4)-sqrt(2*(n+1)*(n+2)*(n+4)))/(n*(n+2)*(n+4)));
        t=sqrt((n+4-2*sqrt(n+4))/((n+2)*(n+4)));
        B=V/(2^n*(n+1));
        
        xi=[fullSymPerms([r;s*ones(numDim-1,1)]),PMCombos(t*ones(numDim,1))];
        w=B*ones(2^n*(n+1),1);
    case 5%S2 5-1 in [1], pg. 279, 7 points, numDim=2, with a correction
          %for a missing +/- term in front of the t in the third set of
          %points and a missing V in the expression for B.
        if(numDim~=2)
            error('This algorithm requires numDim=2');
        end
        
        V=2/(numDim+alpha)*(pi^(numDim/2)/gamma(numDim/2));
        
        r=sqrt((alpha+4)/(alpha+6));
        s=sqrt((alpha+4)/(4*(alpha+6)));
        t=sqrt(3*(alpha+4)/(4*(alpha+6)));
        A=4*V/(alpha+4)^2;
        B=V*(alpha+2)*(alpha+6)/(6*(alpha+4)^2);
        
        xi=[[0;0],PMCombos([r;0]),PMCombos([s;t])];
        w=[A;B*ones(6,1)];
    case 6%S2 5-2 in [1], pg. 280, 9 points, numDim=2, alpha=0.
        if(alpha~=0)
            error('This algorithm requires that alpha=0')
        end
        V=(2/numDim)*(pi^(numDim/2)/gamma(numDim/2));
        
        r=1/sqrt(2);
        A=V/6;
        B=V/24;
        
        xi=[[0;0],fullSymPerms([r;0]),PMCombos([r;r])];
        w=[A*ones(5,1);B*ones(4,1)];
    case 7%S3 5-1 in [1], pg. 289, 13 points, numDim=3.
        if(numDim~=3)
            error('This algorithm requires numDim=3');
        end
        
        V=2/(numDim+alpha)*(pi^(numDim/2)/gamma(numDim/2));
        
        r=sqrt((alpha+5)*(5+sqrt(5))/(10*(alpha+7)));
        s=sqrt((alpha+5)*(5-sqrt(5))/(10*(alpha+7)));
        B0=V*4/(alpha+5)^2;
        B1=V*(alpha+3)*(alpha+7)/(12*(alpha+5)^2);
        
        xi=[[0;0;0],PMCombos([r;s;0]),PMCombos([0;r;s]),PMCombos([s;0;r])];
        w=[B0;B1*ones(12,1)];
    case 8%S3 5-2 in [1], pg. 290, 21 points, numDim=3, alpha=0.
        if(numDim~=3)
            error('This algorithm requires numDim=3');
        end
        if(alpha~=0)
            error('This algorithm requires that alpha=0')
        end
        V=(2/numDim)*(pi^(numDim/2)/gamma(numDim/2));
        
        r=sqrt((15+5*sqrt(5))/42);
        s=sqrt((15-5*sqrt(5))/42);
        t=sqrt(5/21);
        B0=4*V/25;
        B1=21*V/500;
        
        xi=[[0;0;0],PMCombos([r;s;0]),PMCombos([0;r;s]),PMCombos([s;0;r]),PMCombos([t;t;t])];
        w=[B0;B1*ones(20,1)];
    case 9%S4 5-1 in [1], pg. 292, 25 points, numDim=4, alpha=0.
        if(numDim~=4)
            error('This algorithm requires numDim=4');
        end
        if(alpha~=0)
            error('This algorithm requires that alpha=0')
        end
        V=(2/numDim)*(pi^(numDim/2)/gamma(numDim/2));
        r=sqrt(3/8);
        
        B1=V/27;
        B0=V/9;
        
        xi=[[0;0;0;0],fullSymPerms([r;r;0;0])];
        w=[B0;B1*ones(24,1)];
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
