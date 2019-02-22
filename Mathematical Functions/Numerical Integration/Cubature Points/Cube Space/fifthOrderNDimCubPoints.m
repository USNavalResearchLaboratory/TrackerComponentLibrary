function [xi,w]=fifthOrderNDimCubPoints(numDim,algorithm)
%%FIFTHORDERNDIMCUBPOINTS Generate fifth-order cubature points for
%               integration over a numDim-dimensional cube with bounds in
%               coordinates of (-1,1).
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated. numDim>1.
%     algorithm An optional parameter specifying the algorithm to be used
%               to generate the points. Possible values are:
%               0  (The default if omitted or an empty matrix is passed)
%                  Formula Cn 5-2 in [1], pg. 231, 2*numDum^2+1 points.
%               1  Formula Cn 5-3 in [1], pg. 232, 3*numDim^2+3*numDim+1
%                  points.
%               2  Formula Cn 5-4 in [1], pg. 233 2^numDim+2*numDim points.
%               3  Formula Cn 5-5 in [1], pg. 233, 2^numDim+2*numDim+1
%                  points.
%               4  Formula Cn 5-6 in [1], pg. 233, 2^(numDim+1)-1 points.
%               5  Formula Cn 5-7 in [1], pg. 234, numDim*2^numDim+1
%                  points.
%               6  Formula Cn 5-8 in [1], pg. 235, 2^numDim*(numDim+1)
%                  points.
%               7  Formula Cn 5-9 in [1], pg. 234, 3^numDim points.
%               8  Formula C2 5-1 in [1], pg. 246, 7 points, numDim=2.
%               9  Formula C2 5-2 in [1], pg. 247, 7 points, numDim=2.
%               10 Formula C2 5-3 in [1], pg. 248, 8 points, numDim=2.
%               11 Formula C2 5-4 in [1], pg. 249, 9 points, numDim=2.
%               12 Formula C2 5-5 in [1], pg. 250, 13 points, numDim=2.
%               13 Formula C2 5-6 in [1], pg. 251, 13 points, numDim=2.
%               14 Formula C3 5-1 in [1], pg. 262, 13 points, numDim=3.
%                  This is actually the same as the first solution in [4].
%               15 Formula C3 5-2 in [1], pg. 263, 14 points, numDim=3.
%               16 Formula C3 5-3 in [1], pg. 263, 21 points, numDim=3.
%               17 Formula C3 5-4 in [1], pg. 263, 23 points, numDim=3.
%               18 Formula C3 5-5 in [1], pg. 263, 25 points, numDim=3.
%               19 Formula C3 5-6 in [1], pg. 264, 27 points, numDim=3.
%               20 Formula C3 5-7 in [1], pg. 264, 29 points, numDim=3.
%               21 Formula C3 5-8 in [1], pg. 264, 42 points, numDim=3,
%                  with corrections taken from Sadowsky's original paper
%                  [2].
%               22 Formula C4 5-1 in [1], pg. 266, 31 points, numDim=4,
%                  with corrections from the orignal paper by Stroud and
%                  Goit [3].
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%[2] M. Sadowsky, "A formula for approximate computation of a triple
%    integral," The American Mathematical Monthly, vol. 47, no. 8, pp.
%    539-543, Oct. 1940.
%[3] A. H. Stroud and E. H. Goit Jr., "Some extensions of integration
%    formulas," SIAM Journal on Numerical Analysis, vol. 5, no. 2, pp.
%    243-251, Jun. 1968.
%[4] J. W. Peterson, "Analytical formulae for two of A. H. Stroud's
%    quadrature rules," arXiv, 28 Sep. 2009. [Online]. Available:
%    http: //arxiv.org/abs/0909.5106
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    algorithm=0; 
end

switch(algorithm)
    case 0%Cn 5-2 in [1], pg. 231, 2*numDum^2+1 points.
        V=2^numDim;
        
        r=sqrt(3/5);
        B0=V*(25*numDim^2-115*numDim+162)/162;
        B1=V*(70-25*numDim)/162;
        B2=V*25/324;
        
        xi=[zeros(numDim,1),fullSymPerms([r;zeros(numDim-1,1)]),fullSymPerms([r;r;zeros(numDim-2,1)])];
        w=[B0;B1*ones(2*numDim,1);B2*ones(2*numDim*(numDim-1),1)];
    case 1%Cn 5-3 in [1], pg. 232, 3*numDim^2+3*numDim+1 points.
        r=sqrt(7/15);
        s=sqrt((7+sqrt(24))/15);
        t=sqrt((7-sqrt(24))/15);
        V=2^numDim;
        B0=V*(5*numDim^2-15*numDim+14)/14;
        B1=V*25/168;
        B2=-V*25*(numDim-2)/168;
        B3=V*5/48;
        B4=-V*5*(numDim-2)/48;
        
        numPoints=3*numDim^2+3*numDim+1;
        xi=zeros(numDim,numPoints);
        w=zeros(numPoints,1);
        
        %xi(:,1) is all zeros.
        w(1)=B0;
        curStart=2;
        
        xiCur=genAllMultisetPermutations([r;r;zeros(numDim-2,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStart:(curStart+num2Add-1))=xiCur;
        w(curStart:(curStart+num2Add-1))=B1;
        curStart=curStart+num2Add;
        
        xiCur=genAllMultisetPermutations([-r;-r;zeros(numDim-2,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStart:(curStart+num2Add-1))=xiCur;
        w(curStart:(curStart+num2Add-1))=B1;
        curStart=curStart+num2Add;
        
        xiCur=fullSymPerms([r;zeros(numDim-1,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStart:(curStart+num2Add-1))=xiCur;
        w(curStart:(curStart+num2Add-1))=B2;
        curStart=curStart+num2Add;
        
        xiCur=genAllMultisetPermutations([s;-t;zeros(numDim-2,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStart:(curStart+num2Add-1))=xiCur;
        w(curStart:(curStart+num2Add-1))=B3;
        curStart=curStart+num2Add;
        
        xiCur=genAllMultisetPermutations([-s;t;zeros(numDim-2,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStart:(curStart+num2Add-1))=xiCur;
        w(curStart:(curStart+num2Add-1))=B3;
        curStart=curStart+num2Add;
    
        xiCur=fullSymPerms([s;zeros(numDim-1,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStart:(curStart+num2Add-1))=xiCur;
        w(curStart:(curStart+num2Add-1))=B4;
        curStart=curStart+num2Add;
        
        xiCur=fullSymPerms([t;zeros(numDim-1,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStart:(curStart+num2Add-1))=xiCur;
        w(curStart:(curStart+num2Add-1))=B4;
    case 2%Cn 5-4 in [1], pg. 233 2^numDim+2*numDim points.
        V=2^numDim;
        r=sqrt((5*numDim+4)/30);
        s=sqrt((5*numDim+4)/(15*numDim-12));
        B1=40*V/(5*numDim+4)^2;
        B2=2^(-numDim)*((5*numDim-4)/(5*numDim+4))^2*V;
        
        xi=[fullSymPerms([r;zeros(numDim-1,1)]),PMCombos(s*ones(numDim,1))];
        w=[B1*ones(2*numDim,1);B2*ones(2^numDim,1)];
    case 3%Cn 5-5 in [1], pg. 233, 2^numDim+2*numDim+1 points.
        V=2^numDim;
        
        r=sqrt(2/5);
        B0=V*(8-5*numDim)/9;
        B1=V*5/18;
        B2=V*1/(9*2^numDim);
        
        xi=[zeros(numDim,1),fullSymPerms([r;zeros(numDim-1,1)]),PMCombos(ones(numDim,1))];
        w=[B0;B1*ones(2*numDim,1);B2*ones(2^numDim,1)];
    case 4%Cn 5-6 in [1], pg. 233, 2^(numDim+1)-1 points.
        V=2^numDim;
        
        numPoints=2^(numDim+1)-1;
        
        xi=zeros(numDim,numPoints);
        w=zeros(numDim,1);
        
        s=1/sqrt(3);
        
        %The first point is all zeros.
        w(1)=(4/(5*numDim+4))*V;
        
        curStart=2;
        vec=s*ones(numDim,1);
        for k=1:numDim
            vec(k)=sqrt((5*k+4)/15);
            xiCur=PMCombos(vec(k:end));
            numCur=size(xiCur,2);
            
            xi(k:end,curStart:(curStart+numCur-1))=xiCur;
            w(curStart:(curStart+numCur-1))=5*2^(k-numDim+1)*V/((5*k-1)*(5*k+4));
            
            curStart=curStart+numCur;
        end
    case 5%Cn 5-7 in [1], pg. 234, numDim*2^numDim+1 points.
        V=2^numDim;
        B0=V*4/(5*numDim+4);
        B1=V*5*2^(-numDim)/(5*numDim+4);
        
        r=sqrt((5*numDim+4+2*(numDim-1)*sqrt(5*numDim+4))/(15*numDim));
        s=sqrt((5*numDim+4-2*sqrt(5*numDim+4))/(15*numDim));
        
        xi=[zeros(numDim,1),fullSymPerms([r;s*ones(numDim-1,1)])];
        w=[B0;B1*ones(numDim*2^numDim,1)];
    case 6%Cn 5-8 in [1], pg. 235, 2^numDim*(numDim+1) points.
        V=2^numDim;
        r=sqrt((5*numDim-2*sqrt(5)+2*(numDim-1)*sqrt(5*numDim+5))/(15*numDim));
        s=sqrt((5*numDim-2*sqrt(5)-2*sqrt(5*numDim+5))/(15*numDim));
        t=sqrt((5+2*sqrt(5))/15);
        
        w=V/(2^numDim*(numDim+1))*ones(2^numDim*(numDim+1),1);
        xi=[fullSymPerms([r;s*ones(numDim-1,1)]),PMCombos(t*ones(numDim,1))];
    case 7%Cn 5-9 in [1], pg. 234, 3^numDim points.
        r=[-sqrt(3/5);0;sqrt(3/5)];
        A=[5/9;8/9;5/9];
        
        xi=zeros(numDim,3^numDim);
        w=zeros(3^numDim,1);
        
        dims=3*ones(numDim,1);
        for curIdx=1:(3^numDim)
            idxVec=index2NDim(dims,curIdx);
            
            xi(:,curIdx)=r(idxVec);
            w(curIdx)=prod(A(idxVec));
        end
    case 8%C2 5-1 in [1], pg. 246, 7 points, numDim=2.
        if(numDim~=2)
           error('The selected algorithm requires that numDim=2') 
        end
        V=2^numDim;
        
        r=sqrt(3/5);
        s=sqrt(1/3);
        t=sqrt(14/15);
        
        xi=[PMCombos([r;s]),PMCombos([0,t]),[0;0]];
        w=[5/36*ones(4,1);5/63*ones(2,1);2/7]*V;
    case 9%C2 5-2 in [1], pg. 247, 7 points, numDim=2.
        if(numDim~=2)
           error('The selected algorithm requires that numDim=2') 
        end
        V=2^numDim;
        
        r=sqrt(7/15);
        s=sqrt((7+sqrt(24))/15);
        t=sqrt((7-sqrt(24))/15);
        
        xi=[[0;0],[r;r],[-r;-r],[s;-t],[-s;t],[t;-s],[-t;s]];
        w=[2/7;25/168;25/168;5/48;5/48;5/48;5/48]*V;
    case 10%C2 5-3 in [1], pg. 248, 8 points, numDim=2.
        if(numDim~=2)
           error('The selected algorithm requires that numDim=2') 
        end
        V=2^numDim;
        
        r=sqrt(7/15);
        s=sqrt(7/9);
        
        xi=[fullSymPerms([r;0]),PMCombos([s;s])];
        w=[10/49*ones(4,1);9/196*ones(4,1)]*V;
    case 11%C2 5-4 in [1], pg. 249, 9 points, numDim=2.
        if(numDim~=2)
           error('The selected algorithm requires that numDim=2') 
        end
        V=2^numDim;
        
        r=sqrt(3/5);
        xi=[[0;0],fullSymPerms([r;0]),PMCombos([r;r])];
        w=[16/81;10/81*ones(4,1);25/324*ones(4,1)]*V;
    case 12%C2 5-5 in [1], pg. 250, 13 points, numDim=2.
        if(numDim~=2)
           error('The selected algorithm requires that numDim=2') 
        end
        V=2^numDim;
        
        xi=[[0;0],PMCombos([1;1]),fullSymPerms([1;0]),fullSymPerms([0.5;0])];
        w=[-28/45;1/36*ones(4,1);1/45*ones(4,1);16/45*ones(4,1)]*V;
    case 13%C2 5-6 in [1], pg. 251, 13 points, numDim=2.
        if(numDim~=2)
           error('The selected algorithm requires that numDim=2') 
        end
        V=2^numDim;
        
        xi=[[0;0],fullSymPerms([1;0]),PMCombos([1;1]), PMCombos([0.5;0.5])];
        w=[2/45;2/45*ones(4,1);1/60*ones(4,1);8/45*ones(4,1)]*V;
    case 14%C3 5-1 in [1], pg. 262, 13 points, numDim=3, which is the same
           %as the first formula in [4].
        if(numDim~=3)
           error('The selected algorithm requires that numDim=3') 
        end
        
        t=sqrt(71440+6802*sqrt(19));
        
        eta=0;
        lambda=sqrt((1919-148*sqrt(19)+4*t)/3285);
        xi=-sqrt((1121+74*sqrt(19)-2*t)/3285);
        mu=sqrt((1121+74*sqrt(19)+2*t)/3285);
        gammaVal=sqrt((1919-148*sqrt(19)-4*t)/3285);
        
        A=32/19;
        B=133225/(260072-1520*sqrt(19)+(133-37*sqrt(19))*t);
        C=133225/(260072-1520*sqrt(19)-(133-37*sqrt(19))*t);
        
        xi=[[eta;eta;eta],[lambda;xi;xi],-[lambda;xi;xi],[xi;lambda;xi],-[xi;lambda;xi],[xi;xi;lambda],-[xi;xi;lambda],[mu;mu;gammaVal],-[mu;mu;gammaVal],[mu;gammaVal;mu],-[mu;gammaVal;mu],[gammaVal;mu;mu],-[gammaVal;mu;mu]];
        
        w=[A;B;B;B;B;B;B;C;C;C;C;C;C];
    case 15%C3 5-2 in [1], pg. 263, 14 points, numDim=3.
        if(numDim~=3)
           error('The selected algorithm requires that numDim=3') 
        end
        V=2^numDim;
        
        r=sqrt(19/30);
        s=sqrt(19/33);
        B1=V*40/361;
        B2=V*121/2888;
        
        xi=[fullSymPerms([r;0;0]),PMCombos([s;s;s])];
        w=[B1*ones(6,1);B2*ones(8,1)];
    case 16%C3 5-3 in [1], pg. 263, 21 points, numDim=3.
        if(numDim~=3)
           error('The selected algorithm requires that numDim=3') 
        end
        V=2^numDim;
        
        r=1/2;
        
        xi=[[0;0;0],fullSymPerms([r;0;0]),fullSymPerms([1;0;0]),PMCombos([1;1;1])];
        w=[-62/45;16/45*ones(6,1);1/45*ones(6,1);1/72*ones(8,1)]*V;
    case 17%C3 5-4 in [1], pg. 263, 23 points, numDim=3.
        if(numDim~=3)
           error('The selected algorithm requires that numDim=3') 
        end
        V=2^numDim;
        
        r=1/2;
        
        xi=[[0;0;0],fullSymPerms([1;0;0]),PMCombos([r;r;r]),PMCombos([1;1;1])];
        w=[-2/45;2/45*ones(6,1);4/45*ones(8,1);1/120*ones(8,1)]*V;
    case 18%C3 5-5 in [1], pg. 263, 25 points, numDim=3.
        if(numDim~=3)
           error('The selected algorithm requires that numDim=3') 
        end
        V=2^numDim;
        
        r=1/2;
        
        xi=[[0;0;0],fullSymPerms([r;0;0]),fullSymPerms([1;0;0]),fullSymPerms([1;1;0])];
        w=[-19/45;16/45*ones(6,1);-1/30*ones(6,1);1/36*ones(12,1)]*V;
    case 19%C3 5-6 in [1], pg. 264, 27 points, numDim=3.
        if(numDim~=3)
           error('The selected algorithm requires that numDim=3') 
        end
        V=2^numDim;
        
        r=1/2;
        
        xi=[[0;0;0],fullSymPerms([r;0;0]),fullSymPerms([1;1;0]),PMCombos([1;1;1])];
        w=[-4/3;16/45*ones(6,1);1/90*ones(12,1);1/120*ones(8,1)]*V;
    case 20%C3 5-7 in [1], pg. 264, 29 points, numDim=3.
        if(numDim~=3)
           error('The selected algorithm requires that numDim=3') 
        end
        V=2^numDim;
        
        r=1/2;
        
        xi=[[0;0;0],fullSymPerms([1;1;0]),PMCombos([r;r;r]),PMCombos([1;1;1])];
        w=[2/45;1/45*ones(12,1);4/45*ones(8,1);-1/360*ones(8,1)]*V;
    case 21%C3 5-8 in [1], pg. 264, 42 points, numDim=3, with corrections
           %taken from Sadowsky's original paper.
        if(numDim~=3)
           error('The selected algorithm requires that numDim=3') 
        end

        r=sqrt(5/8);
        xi=[fullSymPerms([1;0;0]),fullSymPerms([1;1;0]),fullSymPerms([1;r;r])];
        w=[91*ones(6,1);-40*ones(12,1);16*ones(24,1)]*4/225;
    case 22%C4 5-1 in [1], pg. 266, 31 points, numDim=4, with corrections
           %from the orignal paper by Stroud and Goit.
        if(numDim~=4)
           error('The selected algorithm requires that numDim=4') 
        end
        V=2^numDim;
        
        q=sqrt(8/5);
        r=sqrt(13/30);
        s=sqrt(7/6);
        t=sqrt(14/15);
        u=sqrt(2/3);
        v=sqrt(3/5);
        w=sqrt(1/3);
        
        xi=[[0;0;0;0],PMCombos([0;0;0;q]),PMCombos([0;0;r;s]),PMCombos([0;t;u;0]),PMCombos([0;t;0;u]),PMCombos([v;w;u;0]),PMCombos([v;w;0;u])];
        w=[91/546;-55/1092*ones(2,1);5/91*ones(4,1);5/252*ones(8,1);5/144*ones(16,1)]*V;
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
