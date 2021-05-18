function [xi,w]=fifthOrderCrossPolyCubPoints(numDim,algorithm)
%%FIFTHORDERCROSSPOLYCUBPOINTS Generate third order cubature points for
%               integration over the region sum(abs(x))<=1 with weighting
%               function w(x)=1. In two dimensions, this is a square, in
%               three an octohedron, in n a cross polytope (or
%               n-dimensional octahedron).
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated.
%     algorithm A value indicating which algorithm should be used.
%               Possible values are
%               0 (The default if omitted or an empty matrix is passed)
%                 Formula G_n 5-2 in [1], pg. 305, 4*numDim^2-4*numDim+1
%                 points, 2<=numDim<=7.
%               1 G_n 5-3 in [1], pg. 305, 2^numDim+2*numDim points,
%                 numDim=3,4,5, with a correction that the first term with
%                 (r;0;0...0) be put through fullSymPerms.
%               2 G_n 5-4 in [1], pg. 306, 2^{numDim+1)-1 points, n<4.
%               3 G_3 5-1 in [1], pg. 306, 13 points, numDim=3.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%              points. (Each "point" is a vector)
%            w A numCubaturePoints X 1 vector of the weights associated
%              with the cubature points.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

if(nargin<2||isempty(algorithm))
   algorithm=0; 
end

switch(algorithm)
    case 0%G_n 5-2 in [1], pg. 305, 4*numDim^2-4*numDim+1 points,
          %2<=numDim<=7.
        if(numDim<2||numDim>7)
           error('This formula requires 2<=numDim<=7.') 
        end
        V=exp(numDim*log(2)-gammaln(numDim+1));
        n=numDim;
        
        r=sqrt((n+5-sqrt((n+5)*(7-n)))/((n+3)*(n+4)));
        s=sqrt((n+5+sqrt((n+5)*(7-n)))/((n+3)*(n+4)));
        
        B0=V*(n^2+5*n+10)/((n+1)*(n+2)*(n+5));
        B1=V*(n+3)*(n+4)/(4*(n-1)*(n+1)*(n+2)*(n+5));
        
        xi=[zeros(numDim,1),fullSymPerms([r;s;zeros(numDim-2,1)])];
        w=[B0;B1*ones(4*numDim^2-4*numDim,1)];
    case 1%G_n 5-3 in [1], pg. 305, 2^numDim+2*numDim points, numDim=3,4,5,
          %with a correction that the first term with (r;0;0...0) must be
          %put through fullSymPerms.
        if(numDim~=3 &&numDim~=4&&numDim~=5)
           error('This formula requires numDim=3, 4, or 5.') 
        end
        V=exp(numDim*log(2)-gammaln(numDim+1));
        n=numDim;
        
        r=sqrt((5*(n+3)*(n+4)+sqrt(5*(n+3)*(n+4)*(n^2+5*n+10)))/((n+3)*(n+4)*(2*n+5)));
        s=sqrt((2*n*(n+3)*(n+4)-2*sqrt(5*(n+3)*(n+4)*(n^2+5*n+10)))/((n-2)*(n+3)*(n+4)*(n^2+4*n+5)));
        
        B1=10*V/((n+1)*(n+2)*(n+3)*(n+4)*r^4);
        B2=2^(2-n)*V/((n+1)*(n+2)*(n+3)*(n+4)*s^4);
        
        xi=[fullSymPerms([r;zeros(n-1,1)]),PMCombos(s*ones(numDim,1))];
        w=[B1*ones(2*n,1);B2*ones(2^n,1)];
    case 2%G_n 5-4 in [1], pg. 306, 2^{numDim+1)-1 points, n<4.
        if(numDim>=4)
           error('This formula requires numDim<4') 
        end
        V=exp(numDim*log(2)-gammaln(numDim+1));
        n=numDim;
        
        s=sqrt(2/((n+3)*(n+4)));
        B0=V*(n^2+5*n+10)/((n+1)*(n+2)*(n+5));
        
        numPoints=2^(n+1)-1;
        xi=zeros(numDim,numPoints);
        w=zeros(numPoints,1);
        w(1)=B0;
        
        curStart=2;
        vec=s*ones(numDim,1);
        for k=1:numDim
            r=sqrt(2*(k+5)/((n+3)*(n+4)));
            B=V*2^(k-n-1)*5*(n+3)*(n+4)/((n+1)*(n+2)*(k+4)*(k+5));

            vec(k)=r;
            xiCur=PMCombos(vec(k:end));
            num2Add=size(xiCur,2);
            xi(k:end,curStart:(curStart+num2Add-1))=xiCur;
            w(curStart:(curStart+num2Add-1))=B;

            curStart=curStart+num2Add;
        end
    case 3%G_3 5-1 in [1], pg. 306, 13 points, numDim=3.
        if(numDim~=3)
           error('This formula requires numDim=3') 
        end
        V=exp(numDim*log(2)-gammaln(numDim+1));
        
        r=sqrt((4+2*sqrt(2))/21);
        s=sqrt((4-2*sqrt(2))/21);
        B0=V*68/320;
        B1=V*21/320;
        
        xi=[[0;0;0],PMCombos([r;s;0]),PMCombos([0;r;s]),PMCombos([s;0;r])];
        w=[B0;B1*ones(12,1)];
    otherwise
        error('Unknown algorithm specified');   
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
