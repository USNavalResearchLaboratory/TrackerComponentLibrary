function [xi,w]=eleventhOrderSpherCubPoints(algorithm)
%%ELEVENTHORDERSPHERCUBPOINTS Generate second-order cubature points for
%               integration over a unit spherical (or hyperspherical)
%               region. The weighting function is 1. Currently only
%               integration in 2D is supported.
%
%INPUTS:algorithm A value indicating which algorithm should be used.
%                Possible values are
%                0 (The default if omitted or an empty matrix is passed)
%                  Formula S2 11-1 in [1], pg. 285, 28 points, numDim=2.
%                1 Formula S2 11-4 in [1], pg.286, 32 points, numDim=2.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<1||isempty(algorithm))
   algorithm=0; 
end

switch(algorithm)
    case 0%S2 11-1 in [1], pg. 285, 28 points, numDim=2.
        numDim=2;
        V=2/(numDim)*(pi^(numDim/2)/gamma(numDim/2));
        
        r1=sqrt((31-sqrt(601))/60);
        r2=sqrt(3/5);
        r3=sqrt((31+sqrt(601))/60);
        r4=sqrt((10-sqrt(10))/20);
        r5=sqrt((10+sqrt(10))/20);
        
        s4=sqrt((10-sqrt(10))/60);
        s5=sqrt((10+sqrt(10))/60);
        
        B1=V*(857*sqrt(601)+12707)/(20736*sqrt(601));
        B2=V*125/3456;
        B3=V*(857*sqrt(601)-12707)/(20736*sqrt(601));
        B4=V*(340+25*sqrt(10))/10368;
        B5=V*(340-25*sqrt(10))/10368;
        
        xi=[fullSymPerms([r1;0]),fullSymPerms([r2;0]),fullSymPerms([r3;0]),fullSymPerms([r4;s4]),fullSymPerms([r5;s5])];
        w=[B1*ones(4,1);B2*ones(4,1);B3*ones(4,1);B4*ones(8,1);B5*ones(8,1)];
    case 1%S2 11-4 in [1], pg.286, 32 points, numDim=2.
        numDim=2;
        V=2/(numDim)*(pi^(numDim/2)/gamma(numDim/2));
        
        r1=sqrt((5-sqrt(15))/10);
        r2=1/sqrt(2);
        r3=sqrt((5+sqrt(15))/10);
        u1=sqrt((5+sqrt(45-10*sqrt(15)))/20);
        v1=sqrt((5-sqrt(45-10*sqrt(15)))/20);
        u2=sqrt((5+sqrt(15)+2*sqrt(-150+40*sqrt(15)))/20);
        v2=sqrt((5+sqrt(15)-2*sqrt(-150+40*sqrt(15)))/20);
        t=sqrt((5-sqrt(15))/20);

        B1=V*5/144;
        B2=V*(34-5*sqrt(15))/396;
        B3=V*(4805-620*sqrt(15))/103824;
        C1=V*(10+5*sqrt(15))/792;
        C2=V*(2405+620*sqrt(15))/207648;
        
        xi=[fullSymPerms([r1;0]),fullSymPerms([r2;0]),fullSymPerms([r3;0]),fullSymPerms([u1;v1]),fullSymPerms([u2;v2]),fullSymPerms([t;t])];
        w=[B1*ones(4,1);B2*ones(4,1);B3*ones(4,1);C1*ones(8,1);C2*ones(8,1);B1*ones(4,1)];
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
