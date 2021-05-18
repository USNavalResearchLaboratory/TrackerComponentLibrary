function [xi,w]=fifthOrderSimplexCubPoints(numDim,algorithm)
%%FIFTHORDERSIMPLEXCUBPOINTS Generate fifth-order cubature points for
%               integration over the simplex of points containing all
%               points such that sum(x)<=1 and all x>=0, where x is a
%               numDimX1 vector.
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               to be generated.
%     algorithm A value indicating which algorithm should be used.
%               Possible values are
%               0 (The default if omitted or an empty matrix is passed,
%                 unless numDim=2 or 3) Formula T_n 5-2 in [1], pg. 312,
%                 factorial(numDim+5)/(factorial(5)*factorial(numDim))
%                 points, numDim>=4, where the denominator of B1 has been
%                 corrected.
%               1 (The default if numDim=2) Formula T_2 5-1 in [1], pg.
%                 314, 7 points, numDim=2.
%               2 (The default if numDim=3) Formula T_3 5-1 in [1], pg.
%                 315, 15 points, numDim=3.
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
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

if(nargin<2||isempty(algorithm))
    if(numDim>=4)
        algorithm=0;
    elseif(numDim==2)
        algorithm=1;
    else
       algorithm=2; 
    end
end

n=numDim;
V=1/factorial(n);

switch(algorithm)
    case 0%T_n 5-2 in [1], pg. 312,
          %factorial(numDim+5)/(factorial(5)*factorial(numDim)) points,
          %numDim>=4, where the denominator of B1 has been corrected.
        if(numDim<4)
           error('This algorithm requires that numDim>=4.') 
        end
          
        r=1/5;
        s=4/5;
        u=2/5;
        v=3/5;
        
        nFactRat=exp(gammaln(n+1)-gammaln(n+5+1));
        
        B1=V*(12*n^4-82*n^3+477*n^2-1277*n+1440)*nFactRat/12;
        B2=V*5^2*(-3*n^3+19*n^2-96*n+170)*nFactRat/12;
        B3=V*5^2*(-n^3+13*n^2-47*n+65)*nFactRat/6;
        B4=V*5^3*(n^2-6*n+20)*nFactRat/3;
        B5=V*5^3*(n^2-11*n+20)*nFactRat/4;
        B6=V*5^4*(5-n)*nFactRat/2;
        B7=V*5^5*nFactRat;
        
        numPoints=fix(exp(gammaln(n+5+1)-gammaln(5+1)-gammaln(n+1)));
        xi=zeros(n,numPoints);
        w=zeros(numPoints,1);
        curStart=1;
        
        xiExtended=genAllMultisetPermutations([zeros(n,1);1]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=B1;
        curStart=curStart+num2Add;
        
        xiExtended=genAllMultisetPermutations([zeros(n-1,1);r;s]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=B2;
        curStart=curStart+num2Add;
        
        xiExtended=genAllMultisetPermutations([zeros(n-1,1);u;v]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=B3;
        curStart=curStart+num2Add;
        
        xiExtended=genAllMultisetPermutations([zeros(n-2,1);r;r;v]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=B4;
        curStart=curStart+num2Add;
        
        xiExtended=genAllMultisetPermutations([zeros(n-2,1);r;u;u]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=B5;
        curStart=curStart+num2Add;
        
        xiExtended=genAllMultisetPermutations([zeros(n-3,1);r;r;r;u]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=B6;
        curStart=curStart+num2Add;
        
        xiExtended=genAllMultisetPermutations([zeros(n-4,1);r;r;r;r;r]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=B7;
    case 1%T_2 5-1 in [1], pg. 314, 7 points, numDim=2.
        if(numDim~=2)
           error('This algorithm requires that numDim=2.') 
        end
        
        t=1/3;
        r=(6-sqrt(15))/21;
        u=(6+sqrt(15))/21;
        s=(9+2*sqrt(15))/21;
        v=(9-2*sqrt(15))/21;
        
        A=V*9/40;
        B=V*(155-sqrt(15))/1200;
        C=V*(155+sqrt(15))/1200;
        
        xi=zeros(2,7);
        w=zeros(7,1);
        
        xi(:,1)=[t;t];
        w(1)=A;
        
        xi(:,2)=[r;r];
        xi(:,3)=[r;s];
        xi(:,4)=[s;r];
        w(2:4)=B;
        
        xi(:,5)=[u;u];
        xi(:,6)=[u;v];
        xi(:,7)=[v;u];
        w(5:7)=C;
    case 2%T_3 5-1 in [1], pg. 315, 15 points, numDim=3.
        if(numDim~=3)
           error('This algorithm requires that numDim=3.') 
        end
        
        r=1/4;
        s1=(7-sqrt(15))/34;
        s2=(7+sqrt(15))/34;
        u=(10-2*sqrt(15))/40;
        t1=(13+3*sqrt(15))/34;
        t2=(13-3*sqrt(15))/34;
        v=(10+2*sqrt(15))/40;
        
        A=V*16/135;
        B1=V*(2665+14*sqrt(15))/37800;
        B2=V*(2665-14*sqrt(15))/37800;
        C=V*20/378;
        
        xi=zeros(3,15);
        w=zeros(15,1);
        
        xi(:,1)=[r;r;r];
        w(1)=A;
        
        xi(:,2)=[s1;s1;s1];
        xi(:,3)=[s1;s1;t1];
        xi(:,4)=[s1;t1;s1];
        xi(:,5)=[t1;s1;s1];
        w(2:5)=B1;
        
        xi(:,6)=[s2;s2;s2];
        xi(:,7)=[s2;s2;t2];
        xi(:,8)=[s2;t2;s2];
        xi(:,9)=[t2;s2;s2];
        w(6:9)=B2;
        
        xi(:,10)=[u;u;v];
        xi(:,11)=[u;v;u];
        xi(:,12)=[v;u;u];
        xi(:,13)=[v;v;u];
        xi(:,14)=[v;u;v];
        xi(:,15)=[u;v;v];
        w(10:15)=C;
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
