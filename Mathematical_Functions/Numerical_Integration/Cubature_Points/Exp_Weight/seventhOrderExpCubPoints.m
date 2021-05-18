function [xi,w]=seventhOrderExpCubPoints(numDim,algorithm)
%%SEVENTHORDEREXPCUBPOINTS Generate seventh-order cubature points for
%               integration over real space involving the weighting
%               function w(x)=exp(-sqrt(sum(x.*x))).
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated.
%     algorithm A value indicating which algorithm should be used. Possible
%               values are:
%               0 (The default if omitted or an empty matrix is passed)
%                 Formula E_n^r 7-1 in [1], pg. 331, 2^numDim+2*numDim^2+1
%                 points, 3<=numDim<=7.
%               1 Formula E_2^r 7-1 in [1], pg. 332, 12 points, numDim=2.
%               2 Formula E_3^r 7-1 in [1], pg. 334, 27 points, numDim=3,
%                 with a correction for the denominator of the D term.
%               3 Formula E_3^r 7-2 in [1], pg. 335, 33 ppints, numDim=3.
%               4 Formula E_4^r 7-1 in [1], pg. 335, 49 points, numDim=4.
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
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

if(nargin<2||isempty(algorithm))
   algorithm=0; 
end

switch(algorithm)
    case 0%E_n^r 7-1 in [1], pg. 331, 2^numDim+2*numDim^2+1 points,
          %3<=numDim<=7.
        if(numDim<3||numDim>7)
            error('This algorithm requires 3<=numDim<=7.')
        end
        V=2*pi^(numDim/2)*exp(gammaln(numDim)-gammaln(numDim/2));
        
        n=numDim;
        r=sqrt((n+5)*(3*(n+3)*(8-n)-(n-2)*sqrt(3*(n+3)*(2*n+7)*(8-n)))/(-2*n^2+6*n+11));
        s=sqrt((n+5)*(3*n*(n+3)-2*sqrt(3*(n+3)*(2*n+7)*(8-n)))/(3*n^2+5*n-56));
        t=sqrt((n+5)*(6*(n+3)+sqrt(3*(n+3)*(2*n+7)*(8-n)))/(2*n-5));
        B=V*(8-n)*(n+1)*(n+3)*(n+5)/r^6;
        C=V*(n+1)*(n+3)*(n+5)/(2^n*s^6);
        D=V*(n+1)*(n+3)*(n+5)/(2*t^6);
        A=V-2*n*B-2^n*C-2*n*(n-1)*D;
        
        xi=[zeros(numDim,1),fullSymPerms([r;zeros(numDim-1,1)]),PMCombos(s*ones(numDim,1)),fullSymPerms([t;t;zeros(numDim-2,1)])];
        w=[A;B*ones(2*numDim,1);C*ones(2^numDim,1);D*ones(2*numDim*(numDim-1),1)];
    case 1%E_2^r 7-1 in [1], pg. 332, 12 points, numDim=2.
        if(numDim~=2)
            error('This algorithm requires numDim=2.')
        end
        V=2*pi^(numDim/2)*exp(gammaln(numDim)-gammaln(numDim/2));
        
        r=sqrt(42);
        s=sqrt((6615+21*sqrt(74255))/454);
        t=sqrt((6615-21*sqrt(74255))/454);
        A=V*5/588;
        B=V*(5272105-18733*sqrt(74255))/43661940;
        C=V*(5272105+18733*sqrt(74255))/43661940;
        
        xi=[fullSymPerms([r;0]),PMCombos([s;s]),PMCombos([t;t])];
        w=[A*ones(4,1);B*ones(4,1);C*ones(4,1)];
    case 2%E_3^r 7-1 in [1], pg. 334, 27 points, numDim=3, with a
          %correction for the denominator of the D term.
        if(numDim~=3)
            error('This algorithm requires numDim=3.')
        end
        V=2*pi^(numDim/2)*exp(gammaln(numDim)-gammaln(numDim/2));
        
        r=sqrt((720-24*sqrt(130))/11);
        s=sqrt(288+24*sqrt(130));
        t=sqrt((-216+24*sqrt(130))/7);
        
        A=V*(5175-13*sqrt(130))/8820;
        B=V*(3870+283*sqrt(130))/493920;
        C=V*(3204-281*sqrt(130))/197568;
        D=V*(4239+373*sqrt(130))/197568;
        
        xi=[[0;0;0],fullSymPerms([r;0;0]),fullSymPerms([s;s;0]),PMCombos([t;t;t])];
        w=[A;B*ones(6,1);C*ones(12,1);D*ones(8,1)];
    case 3%E_3^r 7-2 in [1], pg. 335, 33 ppints, numDim=3.
        if(numDim~=3)
            error('This algorithm requires numDim=3.')
        end
        V=2*pi^(numDim/2)*exp(gammaln(numDim)-gammaln(numDim/2));
        
        r=sqrt(-50-10*sqrt(5)+10*sqrt(39)+2*sqrt(195));
        s=sqrt(-50+10*sqrt(5)+10*sqrt(39)-2*sqrt(195));
        t=sqrt(36+4*sqrt(39));
        u=sqrt(54-18*sqrt(5)+6*sqrt(39)-2*sqrt(195));
        v=sqrt(54+18*sqrt(5)+6*sqrt(39)+2*sqrt(195));
        A=V*(1725-26*sqrt(39))/2940;
        B=V*(1065+171*sqrt(39))/54880;
        C=V*(297-47*sqrt(39))/32928;
        
        xi=[[0;0;0],PMCombos([r;s;0]),PMCombos([0;r;s]),PMCombos([s;0;r]),PMCombos([u;v;0]),PMCombos([0;u;v]),PMCombos([v;0;u]),PMCombos([t;t;t])];
        w=[A;B*ones(12,1);C*ones(20,1)];
    case 4%E_4^r 7-1 in [1], pg. 335, 49 points, numDim=4.
        if(numDim~=4)
            error('This algorithm requires numDim=4.')
        end
        V=2*pi^(numDim/2)*exp(gammaln(numDim)-gammaln(numDim/2));
        
        s=sqrt(63-9*sqrt(35)); 
        r=2*s;
        t=sqrt(126+18*sqrt(35));
        
        A=V*53/108;
        B=V*(385+65*sqrt(35))/36288;
        C=V*(385-65*sqrt(35))/36288;
        
        xi=[[0;0;0;0],fullSymPerms([r;0;0;0]),PMCombos([s;s;s;s]),fullSymPerms([t;t;0;0])];
        w=[A;B*ones(24,1);C*ones(24,1)];
        
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
