function [xi,w]=seventhOrderSpherCubPoints(numDim,algorithm)
%%SEVENTHORDERSPHERCUBPOINTS Generate second-order cubature points for
%               integration over a unit spherical (or hyperspherical)
%               region. The weighting function is 1.
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated.
%     algorithm A value indicating which algorithm should be used.
%               Possible values are
%               0 (The default if omitted or an empty matrix is passed)
%                 Formula Sn 7-2 in [1], pg. 271, 2^(numDim+1)+4*numDim^2
%                 points, numDim>=3.
%               1 Formula S2 7-1 in [1], pg. 281, 12 points, numDim=2,
%                 with a correction of a coefficient for B2 from 4 to 41.
%               2 Formula S2 7-2 in [1], pg. 281, 16 points, numDim=2.
%               3 Formula S3 7-2 in [1], pg. 291, 32 points, numDim=3.
%               4 Formula S3 7-3 in [1], pg. 291, 33 points, numDim=3.
%               5 Formula S4 7-2 in [1], pg. 292, 72 points, numDim=4.
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
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

if(nargin<2||isempty(algorithm))
   algorithm=0; 
end

switch(algorithm)
    case 0%Sn 7-2 in [1], pg. 271, 2^(numDim+1)+4*numDim^2 points,
          %numDim>=3.
        if(numDim<3)
           error('This algorithm requires numDim>=3') 
        end
        %Get points for the surface of the sphere.
        [u,B]=seventhOrderSpherSurfCubPoints(numDim,0);
        r=zeros(2,1);
        A=zeros(2,1);
        
        n=numDim;
        
        r(1)=sqrt(((n+2)*(n+4)-2*sqrt(2*(n+2)*(n+4)))/((n+4)*(n+6)));
        r(2)=sqrt(((n+2)*(n+4)+2*sqrt(2*(n+2)*(n+4)))/((n+4)*(n+6)));
        A(1)=(2*(n+2)^2-(n-2)*sqrt(2*(n+2)*(n+4)))/(4*n*(n+2)^2);
        A(2)=(2*(n+2)^2+(n-2)*sqrt(2*(n+2)*(n+4)))/(4*n*(n+2)^2);

        xi=[r(1)*u,r(2)*u];
        w=[A(1)*B;A(2)*B];
    case 1%S2 7-1 in [1], pg. 281, 12 points, numDim=2, with a correction
          %of a coefficient for B2 from 4 to 41.
        if(numDim~=2)
           error('This algorithm requires numDim=2') 
        end
        V=(2/numDim)*(pi^(numDim/2)/gamma(numDim/2));
        r=sqrt(3/4);
        s=sqrt((27-3*sqrt(29))/104);
        t=sqrt((27+3*sqrt(29))/104);
        B1=V*2/27;
        B2=V*(551+41*sqrt(29))/6264;
        B3=V*(551-41*sqrt(29))/6264;
        
        xi=[fullSymPerms([r;0]),PMCombos([s;s]),PMCombos([t;t])];
        w=[B1*ones(4,1);B2*ones(4,1);B3*ones(4,1)];
    case 2%S2 7-2 in [1], pg. 281, 16 points, numDim=2.
        if(numDim~=2)
           error('This algorithm requires numDim=2') 
        end
        V=(2/numDim)*(pi^(numDim/2)/gamma(numDim/2));
        
        r1=sqrt((3-sqrt(3))/6);
        r2=sqrt((3+sqrt(3))/6);
        
        w=V/16*ones(16,1);
        
        xi=zeros(2,16);
        k=1:8;
        
        xi(1,1:8)=r1*cos((2*k-1)*pi/8);
        xi(2,1:8)=r1*cos((2*k-1)*pi/8);
        xi(1,9:16)=r2*cos((2*k-1)*pi/8);
        xi(2,9:16)=r2*cos((2*k-1)*pi/8);
    case 3%S3 7-2 in [1], pg. 291, 32 points, numDim=3.
        if(numDim~=3)
           error('This algorithm requires numDim=3') 
        end
        V=(2/numDim)*(pi^(numDim/2)/gamma(numDim/2));
        
        r=sqrt((1715-7*sqrt(17770))/2817);
        s=sqrt((1715+7*sqrt(17770))/2817);
        t=sqrt(7/18);
        u=sqrt(7/27);
        B1=V*(2965*sqrt(17770)+227816)/(72030*sqrt(17770));
        B2=V*(2965*sqrt(17770)-227816)/(72030*sqrt(17770));
        B3=324*V/12005;
        B4=2187*V/96040;
        
        xi=[fullSymPerms([r;0;0]),fullSymPerms([s;0;0]),fullSymPerms([t;t;0]),PMCombos([u;u;u])];
        w=[B1*ones(6,1);B2*ones(6,1);B3*ones(12,1);B4*ones(8,1)];
    case 4%S3 7-3 in [1], pg. 291, 33 points, numDim=3.
        if(numDim~=3)
           error('This algorithm requires numDim=3') 
        end
        V=(2/numDim)*(pi^(numDim/2)/gamma(numDim/2));
        
        r=sqrt((5+sqrt(5))/18);
        s=sqrt((5-sqrt(5))/18);
        u=sqrt((3-sqrt(5))/6);
        v=sqrt((3+sqrt(5))/6);
        t=1/sqrt(3);
        B0=16*V/175;
        B1=81*V/1400;
        B2=3*V/280;
        
        xi=[[0;0;0],PMCombos([r;s;0]),PMCombos([0;r;s]),PMCombos([s;0;r]),PMCombos([u;v;0]),PMCombos([0;u;v]),PMCombos([v;0;u]),PMCombos([t;t;t])];
        w=[B0;B1*ones(12,1);B2*ones(20,1)];
    case 5%S4 7-2 in [1], pg. 292, 72 points, numDim=4.
        if(numDim~=4)
           error('This algorithm requires numDim=4') 
        end
        V=(2/numDim)*(pi^(numDim/2)/gamma(numDim/2));
        
        r=sqrt((39-3*sqrt(41))/64);
        s=sqrt((39+3*sqrt(41))/64);
        t=1/sqrt(2);
        u=1/2;
        
        B1=V*(33*sqrt(41)+109)/(1440*sqrt(41));
        B2=V*(33*sqrt(41)-109)/(1440*sqrt(41));
        B3=V/240;
        B4=V/60;
        
        xi=[fullSymPerms([r;0;0;0]),fullSymPerms([s;0;0;0]),fullSymPerms([t;t;0;0]),fullSymPerms([u;u;u;0])];
        w=[B1*ones(8,1);B2*ones(8,1);B3*ones(24,1);B4*ones(32,1)];
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
