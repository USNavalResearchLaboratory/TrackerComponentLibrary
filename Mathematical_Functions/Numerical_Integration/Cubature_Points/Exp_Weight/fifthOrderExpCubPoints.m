function [xi,w]=fifthOrderExpCubPoints(numDim,algorithm)
%%FIFTHORDEREXPCUBPOINTS Generate fifth-order cubature points for
%               integration over real space involving the weighting
%               function w(x)=exp(-sqrt(sum(x.*x))).
%
%INPUTS:  numDim An integer specifying the dimensionality of the points
%                to be generated.
%      algorithm A value indicating which algorithm should be used.
%                Possible values are
%                0 (The default if omitted or an empty matrix is passed)
%                  Formula E_n^r 5-1 in [1], pg. 329, 2*numDim^2+1 points.
%                1 Formula E_n^r 5-3 in [1], pg. 330, 2^(numDim+1)-1
%                  points.
%                2 Formula E_n^r 5-4 in [1], pg. 330, numDum*2^numDum+1
%                  points.
%                3 Formula E_2^r 5-1 in [1], pg. 331, 7 points, numDim=2.
%                4 Formula E_3^r 5-1 in [1], pg. 334, 13 points, numDim=3.
%                5 Formula E_3^r 5-2 in [1], pg. 334, 15 points, numDim=3.
%                6 Formula E_4^r 5-1 in [1], pg. 335, 25 points, numDim=4.
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
    case 0%E_n^r 5-1 in [1], pg. 329, 2*numDim^2+1 points.
        V=2*pi^(numDim/2)*exp(gammaln(numDim)-gammaln(numDim/2));
        
        r=sqrt((numDim+2)*(numDim+3));
        s=sqrt((numDim+2)*(numDim+3)/2);
        
        A=V*2*(2*numDim+3)/((numDim+2)*(numDim+3));
        B=V*(4-numDim)*(numDim+1)/(2*(numDim+2)^2*(numDim+3));
        C=V*(numDim+1)/((numDim+2)^2*(numDim+3));
        
        xi=[zeros(numDim,1),fullSymPerms([r;zeros(numDim-1,1)]),fullSymPerms([s;s;zeros(numDim-2,1)])];
        w=[A;B*ones(2*numDim,1);C*ones(2*numDim*(numDim-1),1)];
    case 1%E_n^r 5-3 in [1], pg. 330, 2^(numDim+1)-1 points.
        V=2*pi^(numDim/2)*exp(gammaln(numDim)-gammaln(numDim/2));
        n=numDim;
        
        s=sqrt(n+3);
        vec=s*ones(numDim,1);
        
        numPoints=2^(n+1)-1;
        xi=zeros(numDim,numPoints);
        w=zeros(numPoints,1);
        curStart=1;
        
        BPowSum=0;
        for k=1:numDim
            B=V*(2^(k-n)*(n+1))/((k+1)*(k+2)*(n+3));
            r=sqrt((k+2)*(n+3));
            
            BPowSum=BPowSum+2^(n-k+1)*B;
            vec(k)=r;
            xiCur=PMCombos(vec(k:end));
            num2Add=length(xiCur);
            xi(k:end,curStart:(curStart+num2Add-1))=xiCur;
            w(curStart:(curStart+num2Add-1))=B;
            
            curStart=curStart+num2Add;
        end
        
        B0=V-BPowSum;
        w(end)=B0;
    case 2%E_n^r 5-4 in [1], pg. 330, numDum*2^numDum+1 points.
        V=2*pi^(numDim/2)*exp(gammaln(numDim)-gammaln(numDim/2));
        
        n=numDim;
        r=sqrt(((n+2)*(n+3)+(n-1)*(n+3)*sqrt(2*(n+2)))/n);
        s=sqrt(((n+2)*(n+3)-(n+3)*sqrt(2*(n+2)))/n);
        A=V*(4*n+6)/((n+2)*(n+3));
        B=V*2^(-n)*(n+1)/((n+2)*(n+3));
        
        xi=[zeros(numDim,1),fullSymPerms([r;s*ones(numDim-1,1)])];
        w=[A;B*ones(numDim*2^numDim,1)];
    case 3%E_2^r 5-1 in [1], pg. 331, 7 points, numDim=2.
        if(numDim~=2)
            error('This algorithm requires numDim=2.')
        end
        V=2*pi^(numDim/2)*exp(gammaln(numDim)-gammaln(numDim/2));
        
        r=2*sqrt(5);
        s=sqrt(5);
        t=sqrt(15);
        A=7*V/10;
        B=V/20;
        
        xi=[[0;0],PMCombos([r;0]),PMCombos([s;t])];
        w=[A;B*ones(6,1)];
    case 4%E_3^r 5-1 in [1], pg. 334, 13 points, numDim=3.
        if(numDim~=3)
            error('This algorithm requires numDim=3.')
        end
        V=2*pi^(numDim/2)*exp(gammaln(numDim)-gammaln(numDim/2));
        
        r=sqrt(15+3*sqrt(5));
        s=sqrt(15-3*sqrt(5));
        A=3*V/5;
        B=V/30;
        
        xi=[[0;0;0],PMCombos([r;s;0]),PMCombos([0;r;s]),PMCombos([s;0;r])];
        w=[A;B*ones(12,1)];
    case 5%E_3^r 5-2 in [1], pg. 334, 15 points, numDim=3.
        if(numDim~=3)
            error('This algorithm requires numDim=3.')
        end
        V=2*pi^(numDim/2)*exp(gammaln(numDim)-gammaln(numDim/2));
        
        r=sqrt(30);
        s=sqrt(10);
        A=3*V/5;
        B=2*V/75;
        C=3*V/100;
        
        xi=[[0;0;0],fullSymPerms([r;0;0]),PMCombos([s;s;s])];
        w=[A;B*ones(6,1);C*ones(8,1)];
    case 6%E_4^r 5-1 in [1], pg. 335, 25 points, numDim=4.
        if(numDim~=4)
            error('This algorithm requires numDim=4.')
        end
        V=2*pi^(numDim/2)*exp(gammaln(numDim)-gammaln(numDim/2));
        
        r=sqrt(21);
        A=11*V/21;
        B=5*V/252;
        
        xi=[[0;0;0;0],fullSymPerms([r;r;0;0])];
        w=[A;B*ones(24,1)];
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
