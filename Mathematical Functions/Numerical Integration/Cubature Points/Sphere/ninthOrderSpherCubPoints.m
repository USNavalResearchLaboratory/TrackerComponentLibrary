function [xi,w]=ninthOrderSpherCubPoints(algorithm)
%%NINTHORDERSPHERCUBPOINTS Generate second-order cubature points for
%               integration over a unit spherical (or hyperspherical)
%               region. The weighting function is 1. Currently only
%               integration in 2D is supported.
%
%INPUTS:algorithm A value indicating which algorithm should be used.
%                Possible values are
%                0 (The default if omitted or an empty matrix is passed)
%                  Formula S2 9-1 in [1], pg. 281, 19 points, numDim=2.
%                1 Formula S2 9-3 in [1], pg.284, 21 points, numDim=2.
%                2 Formula S2 9-5 in [1], pg. 285, 28 points, numDim=2.
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
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.


if(nargin<1||isempty(algorithm))
   algorithm=0; 
end

switch(algorithm)
    case 0%S2 9-1 in [1], pg. 281, 19 points, numDim=2.
        numDim=2;
        V=2/(numDim)*(pi^(numDim/2)/gamma(numDim/2));
        
        rho1=sqrt((96-4*sqrt(111))/155);
        rho2=sqrt((96+4*sqrt(111))/155);
        B0=251*V/2304;
        B1=V*(110297+5713*sqrt(111))/2045952;
        B2=V*(110297-5713*sqrt(111))/2045952;
        C=125*V/3072;

        k=1:6;
        r1=rho1*cos(k*pi/3);
        r2=rho2*cos(k*pi/3);
        s1=rho1*sin(k*pi/3);
        s2=rho2*sin(k*pi/3);
        u=sqrt(4/5)*cos((2*k-1)*pi/6);
        v=sqrt(4/5)*sin((2*k-1)*pi/6);
        
        xi=[[0;0],[r1;s1],[r2;s2],[u;v]];
        w=[B0;B1*ones(6,1);B2*ones(6,1);C*ones(6,1)];
    case 1%S2 9-3 in [1], pg.284, 21 points, numDim=2.
        numDim=2;
        V=2/(numDim)*(pi^(numDim/2)/gamma(numDim/2));
    
        rho1=sqrt((6-sqrt(6))/10);
        rho2=sqrt((6+sqrt(6))/10);
        B1=V*(16+sqrt(6))/360;
        B2=V*(16-sqrt(6))/360;
        B0=V/9;
        
        k=1:10;
        r1=rho1*cos(k*pi/5);
        r2=rho2*cos(k*pi/5);
        s1=rho1*sin(k*pi/5);
        s2=rho2*sin(k*pi/5);
        
        xi=[[0;0],[r1;s1],[r2;s2]];
        w=[B0;B1*ones(10,1);B2*ones(10,1)];
    case 2%S2 9-5 in [1], pg. 285, 28 points, numDim=2.
        numDim=2;
        V=2/(numDim)*(pi^(numDim/2)/gamma(numDim/2));
        
        r=sqrt((5+sqrt(15))/10);
        u1=sqrt((5-sqrt(15))/10)*cos(pi/8);
        v1=sqrt((5-sqrt(15))/10)*sin(pi/8);
        u2=sqrt(1/2)*cos(pi/8);
        v2=sqrt(1/2)*sin(pi/8);
        u3=sqrt((5+sqrt(15)+sqrt(-700+185*sqrt(15)))/20);
        v3=sqrt((5+sqrt(15)-sqrt(-700+185*sqrt(15)))/20);
        
        B1=V*(12060-1440*sqrt(15))/254088;
        B2=5*V/144;
        B3=V/18;
        B4=V*(5585+1440*sqrt(15))/508176;
        
        xi=[fullSymPerms([r;0]),fullSymPerms([u1;v1]),fullSymPerms([u2;v2]),fullSymPerms([u3;v3])];
        w=[B1*ones(4,1);B2*ones(8,1);B3*ones(8,1);B4*ones(8,1)];
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
