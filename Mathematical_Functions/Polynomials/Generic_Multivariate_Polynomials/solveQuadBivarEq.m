function sols=solveQuadBivarEq(coeffsA,coeffsB,AbsTol)
%%SOLVEQUADBIVAREQ Find the zeros of a set of two quadratic bivariate
%                  equations. This assumes that the systems are well
%                  defined so there exists a finite number of solutions.
%
%INPUTS: coeffsA, coeffsB 3X3 matrices of the coefficients for the
%               multivariate polynomials whose zeros are desired. These are
%               arranged such that coeffsA(a1,a2) corresponds to the
%               coefficient of an x1^(a1-1)*x2^(a2-1) term in the first
%               polynomial. Entries corresponding to degrees higher than
%               two are ignored. For example coeffsA(3+1,3+1) is ignored if
%               a matrix larger than 3X3 is passed.
%        AbsTol A tolerance value for determining whether the sum of the
%               magnitudes of the two equations is zero. If omitted or an
%               empty matrix is passed, the default of 1e-7 is used. There
%               is a loss of precision when repeated roots arise.
%
%OUTPUTS: sols A 2XnumSol set of the solutions. sols(1,:) are the x
%              solutions; sols(2,:) are the corresponding y solutions.
%              There can be up to four solutions.
%
%The basic idea behind this function is discussed in Appendix A of [1]. The
%problem is solved by eliminating one variable resulting in a quartic
%equation that can be easily solved. Given the solution to the quartic, the
%results can be substituted back into the quadratic equations. However, the
%quadratics produce extra solutions. Normally, one would just choose the
%solutions that agrees with both quadratic equations. However, that is
%problematic if y has repeated roots. Thus, we save all solutions below
%AbsTol in agreement, but we can only return four solutions, because
%Bézout's number is four. Thus, if more than four solutions are saved, we
%take the lowest-cost solution and then sequentially take the others such
%that they maximize the minimum distance to previous solutions (This
%heuristic helps avoid near-duplicate solutions). In the end, however,
%finite precision issues can occasionally cause problems when repeated
%roots are present.
%
%EXAMPLE 1:
%Consider the system:
%0=-x1^2+2*x1*x2+x2^2+5*x1-3*x2-4
%0=x1^2+2*x1*x2+x2^2-1
% coeffsA=zeros(2+1,2+1);
% coeffsA(0+1,0+1)=-4;
% coeffsA(2+1,0+1)=-1;
% coeffsA(1+1,1+1)=2;
% coeffsA(0+1,2+1)=1;
% coeffsA(1+1,0+1)=5;
% coeffsA(0+1,1+1)=-3;
% 
% coeffsB=zeros(2+1,2+1);
% coeffsB(0+1,0+1)=-1;
% coeffsB(2+1,0+1)=1;
% coeffsB(1+1,1+1)=2;
% coeffsB(0+1,2+1)=1;
% sols=solveQuadBivarEq(coeffsA,coeffsB)
%One obtains the solutions
%(4,-5), (1,0), (3,-2), (0,-1).
%
%EXAMPLE 2:
%In this instance, we consider where repeated solutions are present.
%0=(1/4)*x^2+y^2-1
%0=1-2*x+x^2+y^2-1
% coeffsA=zeros(2+1,2+1);
% coeffsA(0+1,0+1)=-1;
% coeffsA(2+1,0+1)=1/4;
% coeffsA(0+1,2+1)=1;
% 
% coeffsB=zeros(2+1,2+1);
% coeffsB(1+1,0+1)=-2;
% coeffsB(2+1,0+1)=1;
% coeffsB(0+1,2+1)=1;
% sols=solveQuadBivarEq(coeffsA,coeffsB)
%The solutions are (2,0), (2,0), (2/3, 2*sqrt(2)/3), (2/3, -2*sqrt(2)/3)
%and in this instance, the loss of precision due to repeated roots is not
%that serious.
%
%EXAMPLE 3:
%This is the same as example 2, but switching x and y
% coeffsA=zeros(2+1,2+1);
% coeffsA(0+1,0+1)=-1;
% coeffsA(0+1,2+1)=1/4;
% coeffsA(2+1,0+1)=1;
% 
% coeffsB=zeros(2+1,2+1);
% coeffsB(0+1,1+1)=-2;
% coeffsB(0+1,2+1)=1;
% coeffsB(2+1,0+1)=1;
% sols=solveQuadBivarEq(coeffsA,coeffsB)
%Here, one can see that there is a slight loss of precision and the results
%veer slightly complex. Of course, if one knows that only real solutions
%should exist, one can simple discard the imaginary component.
%
%REFERENCES:
%[1] D. F. Crouse, "General multivariate polynomial target localization and
%    initial estimation," Journal of Advances in Infromation Fusion, vol.
%    13, no. 1, pp. 68-91, Jun. 2018.
%
%April 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(AbsTol))
    AbsTol=1e-7;%Absolute tolerance for declaring equations zero.
end

%Normalize coefficients to reduce the influence of scaling issues.
coeffsA=coeffsA/max(abs(coeffsA(:)));
coeffsB=coeffsB/max(abs(coeffsB(:)));

a1=coeffsA(0+1,2+1);
a2=coeffsA(1+1,1+1);
a3=coeffsA(0+1,1+1);
a4=coeffsA(2+1,0+1);
a5=coeffsA(1+1,0+1);
a6=coeffsA(0+1,0+1);

b1=coeffsB(0+1,2+1);
b2=coeffsB(1+1,1+1);
b3=coeffsB(0+1,1+1);
b4=coeffsB(2+1,0+1);
b5=coeffsB(1+1,0+1);
b6=coeffsB(0+1,0+1);

c0=a6^2*b4^2+b6*(a5^2*b4-a4*a5*b5+a4^2*b6)+a6*(-a5*b4*b5+a4*(b5^2-2*b4*b6));
c1=(a5^2*b3*b4+a6*b4*(2*a3*b4-a2*b5)+2*a4^2*b3*b6-a5*(a6*b2*b4+a4*b3*b5+a3*b4*b5+a4*b2*b6-2*a2*b4*b6)+a4*(-2*a6*b3*b4+2*a6*b2*b5+a3*b5^2-2*a3*b4*b6-a2*b5*b6));
c2=(b4*(a5^2*b1+(a3^2+2*a1*a6)*b4-a5*(a3*b2-2*a2*b3+a1*b5)-a2*(a6*b2+a3*b5)+a2^2*b6)+a4^2*(b3^2+2*b1*b6)-a4*(-a6*b2^2+a5*b2*b3+2*a6*b1*b4+2*a3*b3*b4+a5*b1*b5-2*a3*b2*b5+a2*b3*b5-a1*b5^2+a2*b2*b6+2*a1*b4*b6));
c3=(2*a4^2*b1*b3-a4*(a5*b1*b2-a3*b2^2+a2*b2*b3+2*a3*b1*b4+2*a1*b3*b4+a2*b1*b5-2*a1*b2*b5)+b4*(-a1*a5*b2+a2^2*b3+2*a1*a3*b4+a2*(2*a5*b1-a3*b2-a1*b5)));
c4=(a4^2*b1^2+b4*(a2^2*b1-a1*a2*b2+a1^2*b4)+a4*(-a2*b1*b2+a1*(b2^2-2*b1*b4)));

y=roots([c4;c3;c2;c1;c0]);

sols=zeros(2,2*4);
costVals=zeros(2*4,1);
numSol=0;
for curY=1:4
    yCur=y(curY);
    %There are two equations that we can try to solve. Solve the one that
    %has the highest weight on the quadratic term.
    if(abs(a4)>abs(b4))
        x1=(-a5-a2*yCur+sqrt((a5+a2*yCur)^2-4*a4*(a6+a3*yCur+a1*yCur^2)))/(2*a4);
        x2=(-a5-a2*yCur-sqrt((a5+a2*yCur)^2-4*a4*(a6+a3*yCur+a1*yCur^2)))/(2*a4);
    else
        x1=(-b5-b2*yCur+sqrt((b5+b2*yCur)^2-4*b4*(b6+b3*yCur+b1*yCur^2)))/(2*b4);
        x2=(-b5-b2*yCur-sqrt((b5+b2*yCur)^2-4*b4*(b6+b3*yCur+b1*yCur^2)))/(2*b4);
    end
    
    curSolVal=[x1;yCur];
    costVal=abs(polyValMultiDim(coeffsA,curSolVal))+abs(polyValMultiDim(coeffsB,curSolVal));
    if(costVal<AbsTol)
        numSol=numSol+1;
        sols(:,numSol)=curSolVal;
        costVals(numSol)=costVal;
    end

    curSolVal=[x2;yCur];
    costVal=abs(polyValMultiDim(coeffsA,curSolVal))+abs(polyValMultiDim(coeffsB,curSolVal));
    if(costVal<AbsTol)
        numSol=numSol+1;
        sols(:,numSol)=curSolVal;
        costVals(numSol)=costVal;
    end
end

sols=sols(:,1:numSol);

%If sols contains more than four solutions (the maximum allowed), then
%some of them are false. We will only keep the one that has the lowest
%cost and add others the maximize the minimum distance from the other
%solutions.
if(numSol>4)
    distMat=Inf*ones(numSol,numSol);

    for curSol1=1:(numSol-1)
       for curSol2=(curSol1+1):numSol
           diff=sols(:,curSol1)-sols(:,curSol2);
           distMat(curSol1,curSol2)=diff'*diff;
           distMat(curSol2,curSol1)=distMat(curSol1,curSol2);
       end
    end

    solsKeptIdx=zeros(4,1);
    solsKeptBool=false(numSol,1);

    [~,minIdx]=min(costVals);
    %Keep the minimum cost solution.
    solsKeptIdx(1)=minIdx;
    solsKeptBool(minIdx)=true;

    %The next solutions added maximize the minimum distance from the
    %solutions already added.
    numSolFound=1;
    for k=2:4
        maxDist=-Inf;
        maxIdx=0;

        for curSol=1:numSol
            if(solsKeptBool(curSol)==true)
                continue;
            end
            %For this solution, find the previously chosen solution to
            %which it is closest. 
            minDist=Inf;
            minIdx=0;
            for prevSol=1:numSolFound
                prevSolIdx=solsKeptIdx(prevSol);
                if(distMat(curSol,prevSolIdx)<minDist)
                    minDist=distMat(curSol,prevSolIdx);
                    minIdx=curSol;
                end
            end

            if(minDist>maxDist)
               maxDist=minDist;
               maxIdx=minIdx;
            end
        end

        solsKeptIdx(k)=maxIdx;
        solsKeptBool(maxIdx)=true;
        numSolFound=numSolFound+1;
    end
    sols=sols(:,solsKeptIdx);
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
