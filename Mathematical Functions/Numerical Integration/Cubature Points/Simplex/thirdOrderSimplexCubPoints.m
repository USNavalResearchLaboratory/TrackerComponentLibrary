function [xi,w]=thirdOrderSimplexCubPoints(numDim,algorithm)
%%THIRDORDERSIMPLEXCUBPOINTS Generate third-order cubature points for
%               integration over the simplex of points containing all
%               points such that sum(x)<=1 and all x>=0, where x is a
%               numDimX1 vector.
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated.
%     algorithm A value indicating which algorithm should be used. Possible
%               values are:
%               0 (The default if omitted or an empty matrix is passed)
%                 Formula T_n 3-1 in [1], pg. 308, numDim+2 points.
%               1 Formula T_n 3-2 in [1], pg. 308, 2*numDim+2 points.
%               2 Formula T_n 3-3 in [1], pg. 308, 2*n+3 points.
%               3 Formula T_n 3-4 in [1], pg. 308,
%                 (numDim+1)*(numDim+2)/2 points, numDim<7 with a
%                 correction to the denominator of B.
%               4 FormulaT_n 3-5 in [1], pg. 309, (numDim+1)*(numDim+4)/2
%                 points, numDim>=3.
%               5 Formula T_n 3-7 in [1], pg. 309,
%                 (numDim^3+5*numDim+12)/6 points.
%               6 Formula T_n 3-8 in [1], pg. 310,
%                 (numDim^3+11*numDim+12)/6 points.
%               7 Formula T_n 3-9 in [1], pg. 310,
%                 (numDim+1)*(numDim+2)*(numDim+3)/6 points.
%               8 Formula T_n 3-10 in [1], pg. 310,
%                 (numDim^3-numDim+3)/3 points, numDim>=3 and numDim~=5.
%               9 Formula T_n 3-11 in [1], pg. 311,
%                 (numDim^3+2*numDim+3)/3 points, numDim>=3 and numDim~=5.
%              10 Formula T_2 3-1 in [1], pg. 314, 6 points, numDim=2.
%              11 Formula T_5 3-1 in [1], pg. 315, 16 points, numDim=5.
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
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

if(nargin<2||isempty(algorithm))
   algorithm=0; 
end

n=numDim;
V=1/factorial(n);

switch(algorithm)
    case 0%T_n 3-1 in [1], pg. 308, numDim+2 points.
        r=1/(n+1);
        s=1/(n+3);
        t=3/(n+3);
        
        B=-V*(n+1)^2/(4*(n+2));
        C=V*(n+3)^2/(4*(n+1)*(n+2));
        
        numPoints=n+2;
        xi=zeros(n,numPoints);
        w=zeros(numPoints,1);
        
        xi(:,1)=r*ones(numDim,1);
        w(1)=B;
        
        xiExtended=genAllMultisetPermutations([s*ones(n,1);t]);
        xi(:,2:end)=xiExtended(1:n,:);
        w(2:end)=C;
    case 1%T_n 3-2 in [1], pg. 308, 2*numDim+2 points.
        r=(2*n+5-sqrt(4*n+13))/(2*(n+1)*(n+3));
        s=1-n*r;
        
        B=V*(1-sqrt(4*n+13))/(2*(n+1)*(n+2)*(n+3));
        C=V*(2*n^2+10*n+11+sqrt(4*n+13))/(2*(n+1)*(n+2)*(n+3));
        
        numPoints=2*n+2;
        xi=zeros(n,numPoints);
        w=zeros(numPoints,1);
        
        xiExtended=genAllMultisetPermutations([zeros(n,1);1]);
        xi(:,1:(n+1))=xiExtended(1:n,:);
        w(1:(n+1))=B;
        
        xiExtended=genAllMultisetPermutations([r*ones(n,1);s]);
        xi(:,(n+2):end)=xiExtended(1:n,:);
        w((n+2):end)=C;
    case 2%T_n 3-3 in [1], pg. 308, 2*n+3 points.
        r=1/(n+1);
        s=1/n;
        A=V*(3-n)*(n+1)^2/((n+2)*(n+3));
        B=V*3/((n+1)*(n+2)*(n+3));
        C=V*n^3/((n+1)*(n+2)*(n+3));
        
        numPoints=2*n+3;
        xi=zeros(n,numPoints);
        w=zeros(numPoints,1);
        xi(:,1)=r*ones(n,1);
        w(1)=A;
        xi(:,2)=0;
        xi(:,3:(n+2))=eye(n,n);
        w(2:(n+2))=B;
        
        xiExtended=genAllMultisetPermutations([s*ones(n,1);0]);
        xi(:,(n+3):end)=xiExtended(1:n,:);
        w((n+3):end)=C;
    case 3%T_n 3-4 in [1], pg. 308, (numDim+1)*(numDim+2)/2 points,
          %numDim<7 with a correction to the denominator of B.
        if(numDim>=7)
            error('This algorithm requires that numDim<7')
        end
        
        rVals=roots([2*(n-2)*(n+1)*(n+3);-(5*n^2+5*n-18);4*n;-1]);
        %We just choose the first real root.
        rVals=rVals(imag(rVals)==0);
        r=rVals(1);
        t=1/2;
        s=1-n*r;
        
        B=V*(n-2)/((n+1)*(n+2)*(1-2*n*r^2-2*(1-n*r)^2));
        C=(2/n)*(V/(n+1)-B);

        numPoints=(n+1)*(n+2)/2;
        xi=zeros(n,numPoints);
        w=zeros(numPoints,1);
        
        xiExtended=genAllMultisetPermutations([r*ones(n,1);s]);
        xi(:,1:(n+1))=xiExtended(1:n,:);
        w(1:(n+1))=B;
        
        xiExtended=genAllMultisetPermutations([zeros(n-1,1);t;t]);
        xi(:,(n+2):end)=xiExtended(1:n,:);
        w((n+2):end)=C;
    case 4%T_n 3-5 in [1], pg. 309, (numDim+1)*(numDim+4)/2 points,
          %numDim>=3.
        if(numDim<3)
            error('This algorithm requires that numDim>=3')
        end
        r=1/2;
        s=1/n;
        B=V*(6-n)/((n+1)*(n+2)*(n+3));
        C=V*8*(n-3)/((n-2)*(n+1)*(n+2)*(n+3));
        D=V*n^3/((n-2)*(n+1)*(n+2)*(n+3));
        
        numPoints=(n+1)*(n+4)/2;
        xi=zeros(n,numPoints);
        w=zeros(numPoints,1);
        
        xi(:,1)=0;
        xi(:,2:(n+1))=eye(n,n);
        w(1:(n+1))=B;
        curStart=n+2;
        
        xiExtended=genAllMultisetPermutations([zeros(n-1,1);r;r]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=C;
        curStart=curStart+num2Add;
        
        xiExtended=genAllMultisetPermutations([s*ones(n,1);0]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=D;
    case 5%T_n 3-7 in [1], pg. 309, (numDim^3+5*numDim+12)/6 points.
        r=1/(n+1);
        s=1/3;
        A=V*(n+1)^2*(n-3)/((n-2)*(n+2)*(n+3));
        B=V*(9-n)/(2*(n+1)*(n+2)*(n+3));
        C=V*27/((n-2)*(n+1)*(n+2)*(n+3));
        
        numPoints=(n^3+5*n+12)/6;
        xi=zeros(n,numPoints);
        w=zeros(numPoints,1);
        
        xi(:,1)=r;
        w(1)=A;
        
        xi(:,2)=0;
        xi(:,3:(n+2))=eye(n,n);
        w(2:(n+2))=B;
        
        xiExtended=genAllMultisetPermutations([zeros(n-2,1);s;s;s]);
        xi(:,(n+3):end)=xiExtended(1:n,:);
        w((n+3):end)=C;
    case 6%T_n 3-8 in [1], pg. 310, (numDim^3+11*numDim+12)/6 points.
        r=1/n;
        s=1/3;
        A=V*(-n^2+11*n-12)/(2*(n-1)*(n+1)*(n+2)*(n+3));
        B=V*n^3/((n-1)*(n+1)*(n+2)*(n+3));
        C=V*27/((n-1)*(n+1)*(n+2)*(n+3));
        
        numPoints=(n^3+11*n+12)/6;
        xi=zeros(n,numPoints);
        w=zeros(numPoints,1);
        curStart=1;
        
        xiExtended=genAllMultisetPermutations([zeros(n,1);1]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=A;
        curStart=curStart+num2Add;
        
        xiExtended=genAllMultisetPermutations([r*ones(n,1);0]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=B;
        curStart=curStart+num2Add;
        
        xiExtended=genAllMultisetPermutations([zeros(n-2,1);s;s;s]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=C;
    case 7%T_n 3-9 in [1], pg. 310, (numDim+1)*(numDim+2)*(numDim+3)/6
          %points.
        r=1/3;
        s=2/3;
        B=V*(n^2-4*n+6)/((n+1)*(n+2)*(n+3));
        C=V*(27-9*n)/(2*(n+1)*(n+2)*(n+3));
        D=V*27/((n+1)*(n+2)*(n+3));
          
        numPoints=(n+1)*(n+2)*(n+3)/6;
        xi=zeros(n,numPoints);
        w=zeros(numPoints,1);
        curStart=1;
       
        xiExtended=genAllMultisetPermutations([zeros(n,1);1]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=B;
        curStart=curStart+num2Add;
        
        xiExtended=genAllMultisetPermutations([zeros(n-1,1);r;s]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=C;
        curStart=curStart+num2Add;
        
        xiExtended=genAllMultisetPermutations([zeros(n-2,1);r;r;r]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=D;
    case 8%T_n 3-10 in [1], pg. 310, (numDim^3-numDim+3)/3 points,
          %numDim>=3 and numDim~=5.
        if(numDim<3||numDim==5)
         error('This algorithm requires that numDim>=3 and numDim~=5.') 
        end
        r=1/(n+1);
        s=1/3;
        t=1/(n-2);
        
        A=V*(3-n)*(n-12)*(n+1)^2/(3*(n-2)*(n+2)*(n+3));
        B=V*54*(3*n-11)/((n-5)*(n-2)*(n-1)*(n+1)*(n+2)*(n+3));
        C=V*2*(n-2)*(n-2)*(n-9)/((n-5)*(n-1)*(n+1)*(n+2)*(n+3));
        
        numPoints=(n^3-n+3)/3;
        xi=zeros(n,numPoints);
        w=zeros(numPoints,1);
        
        xi(:,1)=r;
        w(1)=A;
        curStart=2;
        
        xiExtended=genAllMultisetPermutations([zeros(n-2,1);s;s;s]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=B;
        curStart=curStart+num2Add;
        
        xiExtended=genAllMultisetPermutations([t*ones(n-2,1);0;0;0]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=C;
    case 9%T_n 3-11 in [1], pg. 311, (numDim^3+2*numDim+3)/3 points,
          %numDim>=3 and numDim~=5.
        if(numDim<3||numDim==5)
            error('This algorithm requires that numDim>=3 and numDim~=5.') 
        end
        
        s=1/3;
        t=1/(n-2);
        A=V*(12-n)/(2*(n+1)*(n+2)*(n+3));
        B=V*27*(n-7)/((n-5)*(n-1)*(n+1)*(n+2)*(n+3));
        C=V*6*(n-2)*(n-2)/((n-5)*(n-1)*(n+1)*(n+2)*(n+3));
        
        numPoints=(n^3+2*n+3)/3;
        xi=zeros(n,numPoints);
        w=zeros(numPoints,1);
        curStart=1;
        
        xiExtended=genAllMultisetPermutations([zeros(n,1);1]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=A;
        curStart=curStart+num2Add;
          
        xiExtended=genAllMultisetPermutations([zeros(n-2,1);s;s;s]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=B;
        curStart=curStart+num2Add;
        
        xiExtended=genAllMultisetPermutations([t*ones(n-2,1);0;0;0]);
        num2Add=size(xiExtended,2);
        xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
        w(curStart:(curStart+num2Add-1))=C;
    case 10%T_2 3-1 in [1], pg. 314, 6 points, numDim=2.
        if(numDim~=2)
            error('This algorithm requires that numDim=2.') 
        end 
        r=1/2;
        s=0;
        u=1/6;
        v=4/6;
        
        B=V/30;
        C=V*9/30;
        
        xi=zeros(2,6);
        w=zeros(6,1);
        
        xi(:,1)=[r;r];
        xi(:,2)=[r;s];
        xi(:,3)=[s;r];
        w(1:3)=B;
        
        xi(:,4)=[u;u];
        xi(:,5)=[u;v];
        xi(:,6)=[v;u];
        w(4:6)=C;
    case 11%T_5 3-1 in [1], pg. 315, 16 points, numDim=5.
        if(numDim~=5)
            error('This algorithm requires that numDim=5.') 
        end
        
        r=1/6;
        s=1/2;
        
        xi=zeros(5,16);
        w=zeros(16,1);
        
        xi(:,1)=[r;r;r;r;r];
        w(1)=V*27/42;
        
        xiExtended=genAllMultisetPermutations([0;0;0;0;s;s]);
        xi(:,2:16)=xiExtended(1:n,:);
        w(2:16)=V/42;
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
