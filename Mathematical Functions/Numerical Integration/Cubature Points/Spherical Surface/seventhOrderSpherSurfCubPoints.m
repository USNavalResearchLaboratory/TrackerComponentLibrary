function [xi,w]=seventhOrderSpherSurfCubPoints(numDim,algorithm)
%%SEVENTHORDERSPHERSURFCUBPOINTS Generate seventh-order cubature points for
%               integration over the surface of a unit (hyper-)sphere
%               (weighting function is just 1).
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated, numDim>=2.
%     algorithm A value indicating which algorithm should be used.
%               Possible values are
%               0 (The default if omitted or an empty matrix is passed)
%                 Use formula I from [1], 2^numDim+2*numDim^2 points,
%                 numDim>=3. For numDim>8, negative weights are present.
%               1 Formula II from [1],
%                 2^numDim+(4/3)*numDim*(numDim-1)*(numDim-2)+2*numDim
%                 points, numDim>=4.
%               2 Formula III from [1],
%                 2^numDim+(4/3)*numDim*(numDim-1)*(numDim-2)+2*numDim*(numDim-1)
%                 points, numDim>=4.
%               3 Formula from [2], 2^n+2*n^2 points, numDim>=3
%               4 Un 7-1 from [3], 2^n+2*n^2 points, with the correction
%                 from Table I of [4], numDim>=3.
%               5 Un 7-2 from [3], 2^n*(n+1) points
%               6 U3 7-1 in [3], 24 points, numDim=3
%               7 U3 7-2 in [3], 26 points, numDim=3
%               8 U4 7-1 in [3], 48 points, numDim=4
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%REFERENCES:
%[1] A. H. Stroud, "Some seventh degree integration formulas for the
%    surface of an n-sphere," Numerische Mathematik, vol. 11, no. 3, pp.
%    273-276, Mar. 1968.
%[2] A. H. Stroud, "Some seventh degree integration formulas for symmetric
%    regions," SIAM Journal on Numerical Analysis, vol. 4, no. 1, pp.
%    37-44, Mar. 1967.
%[3] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%[4] R. Cools, "An encyclopedia of cubature formulas," Journal of
%    Complexity, vol. 19, no. 3, pp. 445-453, Jun. 2003.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
   algorithm=0; 
end

switch(algorithm)
    case 0%Formula I from [1], 2^n+2*n^2 points, numDim>=3
        if(numDim<3)
           error('numDim must be >=3 to use this formula.') 
        end
        k1=1;
        k2=2;
        k3=numDim;
        
        I1=2*pi^(numDim/2)/gamma(numDim/2);
        
        A1=(8-numDim)/(numDim*(numDim+2)*(numDim+4))*I1;
        A2=4/(numDim*(numDim+2)*(numDim+4))*I1;
        A3=2^(-numDim)*numDim^3/(numDim*(numDim+2)*(numDim+4))*I1;
        
        xi1=fullSymPerms([repmat(1/sqrt(k1),[k1,1]);zeros(numDim-k1,1)]);
        xi2=fullSymPerms([repmat(1/sqrt(k2),[k2,1]);zeros(numDim-k2,1)]);
        xi3=fullSymPerms([repmat(1/sqrt(k3),[k3,1]);zeros(numDim-k3,1)]);
        
        xi=[xi1,xi2,xi3];
        w=[repmat(A1,[size(xi1,2),1]);repmat(A2,[size(xi2,2),1]);repmat(A3,[size(xi3,2),1])];
    case 1%Formula II from [1], 2^n+(4/3)*n*(n-1)*(n-2)+2*n points,
          %numDim>=4
        if(numDim<4)
           error('numDim must be >=4 to use this formula.') 
        end
        
        k1=1;
        k2=3;
        k3=numDim;
        
        I1=2*pi^(numDim/2)/gamma(numDim/2);
        
        A1=(14-numDim)/(2*numDim*(numDim+2)*(numDim+4))*I1;
        A2=27/(4*numDim*(numDim+2)*(numDim+4)*(numDim-3))*I1;
        A3=numDim^3*(numDim-5)/((2^numDim)*numDim*(numDim+2)*(numDim+4)*(numDim-3))*I1;
        
        xi1=fullSymPerms([repmat(1/sqrt(k1),[k1,1]);zeros(numDim-k1,1)]);
        xi2=fullSymPerms([repmat(1/sqrt(k2),[k2,1]);zeros(numDim-k2,1)]);
        xi3=fullSymPerms([repmat(1/sqrt(k3),[k3,1]);zeros(numDim-k3,1)]);
        
        xi=[xi1,xi2,xi3];
        w=[repmat(A1,[size(xi1,2),1]);repmat(A2,[size(xi2,2),1]);repmat(A3,[size(xi3,2),1])];
    case 2%Formula III from [1], 2^n+(4/3)*n*(n-1)*(n-2)+2*n*(n-1) points,
          %numDim>=4
        if(numDim<4)
           error('numDim must be >=4 to use this formula.') 
        end
        
        k1=2;
        k2=3;
        k3=numDim;
        
        I1=2*pi^(numDim/2)/gamma(numDim/2);
        
        A1=4*(14-numDim)/(numDim*(numDim+2)*(numDim+4)*(numDim-2))*I1;
        A2=27*(numDim-8)/(2*numDim*(numDim+2)*(numDim+4)*(numDim-2)*(numDim-3))*I1;
        A3=numDim^3*(numDim^2-9*numDim+38)/((2^numDim)*numDim*(numDim+2)*(numDim+4)*(numDim-2)*(numDim-3))*I1;
        
        xi1=fullSymPerms([repmat(1/sqrt(k1),[k1,1]);zeros(numDim-k1,1)]);
        xi2=fullSymPerms([repmat(1/sqrt(k2),[k2,1]);zeros(numDim-k2,1)]);
        xi3=fullSymPerms([repmat(1/sqrt(k3),[k3,1]);zeros(numDim-k3,1)]);
        
        xi=[xi1,xi2,xi3];
        w=[repmat(A1,[size(xi1,2),1]);repmat(A2,[size(xi2,2),1]);repmat(A3,[size(xi3,2),1])];
    case 3%The formula from [2], 2^n+2*n^2 points, numDim>=3
        if(numDim<3)
           error('numDim must be >=3 to use this formula.') 
        end
        
        V=2*pi^(numDim/2)/gamma(numDim/2);
        
        A1=(8-numDim)/(numDim*(numDim+2)*(numDim+4))*V;
        A2=2^(-numDim)*numDim^3/(numDim*(numDim+2)*(numDim+4))*V;
        A3=4/(numDim*(numDim+2)*(numDim+4))*V;
        
        k1=1;
        k2=numDim;
        k3=2;
        
        xi1=fullSymPerms([repmat(1/sqrt(k1),[k1,1]);zeros(numDim-k1,1)]);
        xi2=fullSymPerms([repmat(1/sqrt(k2),[k2,1]);zeros(numDim-k2,1)]);
        xi3=fullSymPerms([repmat(1/sqrt(k3),[k3,1]);zeros(numDim-k3,1)]);
        
        xi=[xi1,xi2,xi3];
        w=[repmat(A1,[size(xi1,2),1]);repmat(A2,[size(xi2,2),1]);repmat(A3,[size(xi3,2),1])];
    case 4%Un 7-1 from [3], 2^n+2*n^2 points
        n=numDim;
        V=2*pi^(n/2)/gamma(n/2);
        
        r=1;
        s=sqrt(1/n);
        t=sqrt(1/2);
        
        xi=[fullSymPerms([r;zeros(n-1,1)]),PMCombos(s*ones(n,1)),fullSymPerms([t;t;zeros(n-2,1)])];

        B=((8-n)/(n*(n+2)*(n+4)))*V;
        C=((2^(-n)*n^3)/(n*(n+2)*(n+4)))*V;
        D=(4/(n*(n+2)*(n+4)))*V;
        
        w=[B*ones(2*n,1);
           C*ones(2^n,1);
           D*ones(2*n*(n-1),1)];
    case 5%Un 7-2 from [3], 2^n*(n+1) points

        r=sqrt(1/numDim);
        s=sqrt(5/(numDim+4));
        t=sqrt(1/(numDim+4));
        
        V=2*pi^(numDim/2)/gamma(numDim/2);
        
        A=-numDim^2/(2^(numDim+3)*(numDim+2))*V;
        B=(numDim+4)^2/(2^(numDim+3)*numDim*(numDim+2))*V;
        
        xi=[PMCombos(repmat(r,[numDim,1])),fullSymPerms([s;t*ones(numDim-1,1)])];
        w=[A*ones(2^numDim,1);B*ones(numDim*2^numDim,1)];
    case 6%U3 7-1 in [3], 24 points, numDim=3
        if(numDim~=3)
           error('numDim must be 3 to use this formula.') 
        end
        
        rstVals=roots([1;0;-1;0;1/5;0;-1/105]);
        rstVals=sort(rstVals(rstVals>0),'descend');
        r=rstVals(1);
        s=rstVals(2);
        t=rstVals(3);
        
        u=[r;-r;s;-s;t;-t];
        v=[s;t;t;r;r;s];
        w=[t;s;r;t;s;r];
        
        xi=zeros(3,24);
        for curPoint=1:6
           xi(:,curPoint)=[u(curPoint);v(curPoint);w(curPoint)];
           xi(:,curPoint+6)=[u(curPoint);-v(curPoint);-w(curPoint)];
           xi(:,curPoint+12)=[u(curPoint);w(curPoint);-v(curPoint)];
           xi(:,curPoint+18)=[u(curPoint);-w(curPoint);v(curPoint)];
        end

        V=2*pi^(numDim/2)/gamma(numDim/2);
        B=V/24;
        w=ones(24,1)*B;
    case 7%U3 7-2 in [3], 26 points, numDim=3
        if(numDim~=3)
           error('numDim must be 3 to use this formula.') 
        end
        
        r=1;
        s=1/sqrt(2);
        t=1/sqrt(3);
        
        V=2*pi^(numDim/2)/gamma(numDim/2);
        
        xi=zeros(3,26);
        w=zeros(26,1);
        xi(:,1:6)=fullSymPerms([r;0;0]);
        w(1:6)=(40/840)*V;
        xi(:,7:18)=fullSymPerms([s;s;0]);
        w(7:18)=(32/840)*V;
        xi(:,19:26)=PMCombos([t;t;t]);
        w(19:26)=(27/840)*V;
    case 8%U4 7-1 in [3], 48 points, numDim=4
        if(numDim~=4)
           error('numDim must be 4 to use this formula.') 
        end
        
        V=2*pi^(numDim/2)/gamma(numDim/2);
        r=1;
        s=1/2;
        t=1/sqrt(2);
        B=V/48;
        
        xi=[fullSymPerms([r;0;0;0]),PMCombos([s;s;s;s]),fullSymPerms([t;t;0;0])];
        w=B*ones(48,1);
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
