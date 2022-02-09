function [rootVals,exitCode]=polyRootsMultiDim(polyCoeffMats,maxDegIncreases,useMotzkinNull)
%%POLYROOTSMULTIDIM Find the roots of a system of simultaneous multivariate
%                   polynomials. Only the affine roots are found (which are
%                   generally the only ones desired), not the roots at
%                   infinity. The algorithm used transforms the
%                   multivariate root finding problem into a generalized
%                   eigenvalue problem involving the nullspace of a
%                   Macaulay matrix. Due to finite precision errors and the
%                   combinational complexity of the algorithm, it is
%                   generally best if the maximum number of variables is
%                   <=3 and the maximum degree of the polynomials is <=3.
%                   Sparser polynomials are much easier than dense ones
%                   (where all possible monomial coefficients are nonzero).
%                   When dealing with dense polynomials, it is generally
%                   best if the maximum number of variables is 2, where the
%                   maximum degree can get as high as sixth order without
%                   notable problems.
%
%INPUTS: polyCoeffMats An nX1 or 1Xn cell array containing n hypermatrices
%               of polynomial coefficients for a multivariate polynomial. A
%               hypermatrix of the coefficients for a multivariate
%               polynomial contains elements arranged such that
%               coeffs(a1,a2,a3...an) corresponds to the coefficient of an
%               x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1) term.  Thus, the
%               number of indices coeffs takes is equal to the
%               dimensionality of x (not counting singleton dimensions at
%               the end of coeffs). Note that this ordering is the reverse
%               that used in the 1D polyval function that is built into
%               Matlab. The number of elements for each index in coeffs is
%               the maximum order of that dimension +1.
% maxDegIncreases An optional parameter specifying the maximum number
%               of degree increases of the Motzkin matrix to use when
%               trying to solve the multivariate polynomial system. If
%               maxDegIncreases is too small, then the polynomial root
%               solver will fail. This parameter is provided as an option
%               to help identify when solving large systems would be very
%               slow. The default if this parameter is omitted or an empty
%               matrix is passed is 10*n. If one wishes there to be no
%               limit, then one can set maxDegIncreases to Inf and degree
%               increases will be allowed until Matlab runs out of memory.
% useMotzkinNull A key step in the algorithm is the determination of
%               nullspaces. By defauly, the nullspace function is used.
%               However, if desired, the Motzkin algorithm described in [1]
%               can be used. The Motzkin algorithm is generally not as
%               numerically stable, but it can be useful if one wishes to
%               step through the code and compare the results to the values
%               in [1].
%
%OUTPUTS: rootVals An nXnumSol matrix of the numSol zeros of the polynomial
%                  system found, or an empty matrix if exitCode does not
%                  equal zero.
%         exitCode A parameter indicating how the algorithm terminated.
%                  Possible values:
%                  0 The algorithm was a success; the roots were found.
%                  1 The maximum number of degree increases elapsed.
%                  2 A finite precision error occured either causing the
%                    nullity of the Macaulay matrix to change after
%                    stabilization, or to decrease when the degree was
%                    increased.
%
%This function implements Algorithm 3 of [1], an affine null space based
%root finding method for multivariate poylnomials. The polynomials in
%polyCoeffMats are converted into term matrices. The degree negative
%lexicographic ordering in the paper, for a fixed degree, coincides with
%the ordering of the terms provided by unrankTComposition with
%firstElMostSig=true. Thus, ranking and unranking of compositions is used
%to keep track of the monomials in the columns of the Macaulay matrix.
%
%The algorithm requires that a "shift function" be chosen (g(x) in [1]).
%The shift function was arbitrarily chosen to be the sum_{i=1}^n i*xi,
%where xi is the ith unknown variable. That is, the sum of all first degree
%terms multiplied by their position.
%
%The use of the Motzkin algorithm for finding the nullspace computation
%(Sections 3.2.1, 3.2.2 and 6.2.5 of [1]) allows one to avoid computation
%of degree-augmented Macaulay matrices, and provides a basis for the
%nullspace havign the nice sparse structure as used in [1]. However,
%Section 6.2.3 demonstrated that the choice of basis for the nullspace does
%not change the final solution. Thus, the option to use the (possibly more
%stable) null command (which uses singular value decomposition) is offered,
%though it does require keeping track of the entire Macaulay matrix.
%
%Though the Macaulay matrix has a sparse quasi-Toeplitz structure, the
%sparsity was not taken advantage of in the implementation of this
%algorithm.
%
%Roots at infinity are not "true" solutions to the problem; they only arise
%when the system of polynomials has been homogenized (all monomial terms in
%each equation made to have the same total degree) through the introduction
%of an additional variable. Section 2.3 of [1] discusses what they are.
%This algorithm does not find them.
%
%EXAMPLE 1:
%An example is the third-order polynomial system on page 94 of [1]. The
%equations are
%0=x1^2+5*x1*x2+4*x2*x3-10
%0=x2^3+3*x1^2*x2-12
%0=x3^3+4*x1*x2*x3-8
%Thus, the solution is 
% p=zeros(4,4,4);%First Equation
% p(2+1,0+1,0+1)=1;
% p(1+1,1+1,0+1)=5;
% p(0+1,1+1,1+1)=4;
% p(0+1,0+1,0+1)=-10;
% 
% q=zeros(4,4,4);%Second Equation
% q(0+1,3+1,0+1)=1;
% q(2+1,1+1,0+1)=3;
% q(0+1,0+1,0+1)=-12;
% 
% k=zeros(4,4,4);%Third Equation
% k(0+1,0+1,3+1)=1;
% k(1+1,1+1,1+1)=4;
% k(0+1,0+1,0+1)=-8;
% 
% polyCoeffMats=cell(3,1);
% polyCoeffMats{1}=p;
% polyCoeffMats{2}=q;
% polyCoeffMats{3}=k;
% rootVals=polyRootsMultiDim(polyCoeffMats)
%Where 18 solutions are found. Two are real; the rest complex.
%
%EXAMPLE 2:
%This is the system of equations used as an example in Section 2.1. of [1].
%0=-x1^2+2*x1*x2+x2^2+5*x1-3*x2-4
%0=x1^2+2*x1*x2+x2^2-1
%Thus, the solution is 
% p=zeros(2+1,2+1);
% p(0+1,0+1)=-4;
% p(2+1,0+1)=-1;
% p(1+1,1+1)=2;
% p(0+1,2+1)=1;
% p(1+1,0+1)=5;
% p(0+1,1+1)=-3;
% 
% q=zeros(3+1,3+1);
% q(0+1,0+1)=-1;
% q(2+1,0+1)=1;
% q(1+1,1+1)=2;
% q(0+1,2+1)=1;
% polyCoeffMats=cell(2,1);
% polyCoeffMats{1}=p;
% polyCoeffMats{2}=q;
% rootVals=polyRootsMultiDim(polyCoeffMats)
%In this instance, all four roots of the system are real. The solutions are
%(4,-5), (1,0), (3,-2), (0,-1).
%
%REFERENCES:
%[1] P. Dreesen, "Back to the roots: Polynomial system solving using linear
%    algebra," Ph.D. dissertation, Katholieke Universiteit Leuven, Leuven,
%    Flanders, Belgium, Sep. 2013.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The number polynomials and the number of variables.
n=length(polyCoeffMats);

if(nargin<2||isempty(maxDegIncreases))
    maxDegIncreases=10*n;
end

if(nargin<3||isempty(useMotzkinNull))
    useMotzkinNull=false;
end

%Determine the number of variables. Not all polynomials have to have all of
%the same variables, so 
termMats=cell(n,1);

%Holds the degrees of each polynomial.
d=zeros(n,1);

%Convert the polynomial hypermatrices to term matrices. The order of the
%monomials in the term matrix does not matter here (though the ordering of
%the monomials when added to the Motzkin matrix is important).
for curPoly=1:n
    termMat=multiDimPolyMat2Terms(polyCoeffMats{curPoly},2,n);
    
    %Normalize the equation so that the highest coefficient has magnitude
    %one. This helps reduce finite precision problems.
    maxVal=max(abs(termMat(1,:)));
    termMat(1,:)=termMat(1,:)/maxVal;

    termMats{curPoly}=termMat;
    
    %The total degree of the current polynomial.
    d(curPoly)=max(sum(termMat(2:end,:),1));
end

%The maximum total degree of all of the polynomials.
d0=max(d);

%This keeps track of the number of columns of the Macaulay matrix  before
%the first column containing a certain degree monomial term (the
%number of monomial terms before the first monomial term of the given
%degree. The array numBeforeDeg is extended each time another degree block
%is added to the Macaulay matrix. numBeforeDeg(i+1) is the number of items
%before degree i, starting with degree 0, where numBeforeDeg(0+1)=0.
numBeforeDeg=getNumElsBeforeDeg(n,d0);

%Construct the initial Macaulay matrix of degree d0, which is the highest
%total degree of all of the polynomials going in.
[N,numBeforeDeg]=buildInitialMacaulayMat(termMats,d,d0,numBeforeDeg);
p=size(N,1);
q=size(N,2);

%Compute the inital basis for the nullspace
if(useMotzkinNull)
    Z=MotzkinMatrix(N);
else
    Z=nullspace(N);
end
N=[];%N is no longer needed.

%We keep track of the nullity to determine when the size of the nullspace
%has stabilized. Only after it has stabilized do we consider whether we can
%identify where any solutions at infinity are (so as to eliminate them).
nullity=size(Z,2);

dG=[];
dCur=d0;

nullityStabilized=false;
for curDegIncrease=1:maxDegIncreases
    %Get the new rows that would have to be augmented.
    [NRows,numBeforeDeg]=newRows4MacaulayMatrix(p,numBeforeDeg,termMats,d,dCur);
    qNew=size(NRows,2);
    
    %Compute the augmented nullspace for the Macaulay matrix using
    %Motzkin's algorithm (Algorithm 4). This used the block method of
    %Section 6.2.5 to expand the nullspace.
    N1=NRows(:,1:q);
    N2=NRows(:,(q+1):qNew);
    if(useMotzkinNull)
        XY=MotzkinMatrix([N1*Z,N2]);
    else
        XY=nullspace([N1*Z,N2]);
    end
    numZPrev=size(Z,2);
    X=XY(1:numZPrev,:);
    Y=XY((numZPrev+1):end,:);
    Z=[Z*X;Y];

    nullityNew=size(Z,2);

    dCur=dCur+1;
    p=size(N,1);
    q=qNew;
    
    %If we have found the degree for which the nullity no longer increases.
    if(nullityNew==nullity)
        nullityStabilized=true;
        
        %Check whether we are at dG -the degree where we can detect the
        %affine basis set (and thus know that there are ma solutions).
        [isAtdG,ma,degOfGap]=checkFordG(Z,numBeforeDeg,dCur);
        
        if(isAtdG)
            dG=dCur;
            break;
        end
    elseif(nullityStabilized==true||nullityNew<nullity)
        %The nullity has stabilized and then become unstabilized or the
        %nullity decreased. Neither should happen and both are indicative
        %of finite precision problems.
        rootVals=[];
        exitCode=2;
        return;
    else
        nullity=nullityNew;
    end
end

if(isempty(dG))
    rootVals=[];
    exitCode=1;%The maximum number of iterations elapsed.
    return;
end

%S1 will select all monomials up to one less than the degree of the gap. We
%will use all monomials in the nullspace up to the degree of the gap.
k=numBeforeDeg(degOfGap+2);
S1=eye(numBeforeDeg(degOfGap+1),k);

%Perform column compression on Z to get W11 using Theorem 6.9
[~,~,Q]=svd(Z(1:k,:));
W=Z*Q;
W11=W(1:k,1:ma);

Sg=constructSg(n,degOfGap,numBeforeDeg);

%The generalized eigenvalue problem to solve is S1*W11*V11*D=Sg*W11*V11
%However, eig cannot solve that as it is a rectangular problem. We use the
%method of Section 6.2.2 to transform it into a square eigenvalue problem.
[V11,~]=eig(lsqminnorm(S1*W11,Sg*W11));

%Using Corollary 6.11 to extract the actual solutions.
Ka1=W11*V11;
Ka1=bsxfun(@rdivide,Ka1,Ka1(1,:));
rootVals=Ka1(2:(n+1),:);

exitCode=0;%Successful return

end

function Sg=constructSg(n,maxDeg,numBeforeDeg)
%This function constructs Sg where g is the shift function, chosen to be
%the sum of all first order monomial terms weighted by their index.
%That is x1+2*x2+3*x3+...+n*xn
%These values are defined through Proposition 6.3. Our chose of sg was
%arbitrary.

termMats{1}=[(1:n);eye(n)];

Sg=buildInitialMacaulayMat(termMats,1,maxDeg,numBeforeDeg);
end

function [dGHasBeenReached,ma,degOfGap]=checkFordG(Z,numBeforeDeg,maxDeg)
%This implements Corollary 6.12 to determine whether we have reached a
%point where we can identify the affine basis set.

    curRank=1;%The rank of the degree-0 block.
    
    for curDeg=1:maxDeg
        selRows=(1:numBeforeDeg(curDeg+2));
        newRank=rank(Z(selRows,:));
        if(newRank==curRank)
            dGHasBeenReached=true;
            ma=curRank;%The number of affine bases.
            degOfGap=curDeg;
            return;
        end
        curRank=newRank;
    end
    dGHasBeenReached=false;%No gap detected.
    ma=[];
    degOfGap=[];
end

function W=MotzkinRow(b,epsVal)
%Chapter 3.2.1. A subroutine in the Motzkin method for finding the null
%space of a matrix. This implements the row-subroutine used by the
%MotzkinMatrix function.

n=length(b);
W=zeros(n,n-1);

%Explicitly zero values that are too small in magnitude.
b(abs(b)<epsVal)=0;

%Step 1: Find the rightmost pivot. This is the rightmost nonzero value.
ip=find(b,1,'last');

%If there is no nonzero element.
if(isempty(ip))
    W=eye(n,n);
    return;
end

%Step 2: Rescale b
b=b/b(ip);

%Step 3
e=zeros(n,1);
for i=(ip+1):n
    e(i)=1;
    W(:,i-1)=e;
    e(i)=0;
end

%Step 4
for i=(ip-1):-1:1
    W(i,i)=1;
    W(ip,i)=-b(i);
end
end

function H=MotzkinMatrix(A)
%Chapter 3.2.2 to find the Motzkin (canonical) nullspace of a matrix.
%An example with A=[0,2,1,0,2,0;1,0,3,2,0,1;4,0,3,2,1,0] in the paper.

m=size(A,1);

%A value for declaring elements of b zero for purposes of the Motzkin row
%algorithm.
epsVal=max(size(A))*eps(norm(A));

b=A(1,:);

H=MotzkinRow(b,epsVal);
for i=2:m
   b=A(i,:)*H;
   H=H*MotzkinRow(b,epsVal);
end
end

function [MRows,numBeforeDeg]=newRows4MacaulayMatrix(p,numBeforeDeg,termMats,d,dCur)
%To enlarge the Macaulay matrix to one degree higher, new rows must be
%appended. This computes those rows. The Macaulay matrix is defined in
%Section 5.1.2. This function only computes the rows of a particular new
%degree. When increasing the degree of the Macaulay matrix, it is not
%necessary to respect the same row ordering as in a Macaulay matrix built
%to a particular degree from scratch.

s=length(termMats);%The number of polynomials
n=size(termMats{1},1)-1;%The number of variables

%The new maximum degree if these rows are appended.
d0=dCur+1;

%Update offsets with the current degree.
numCompositionsInDeg=binomial(dCur+1+n-1,n-1);
numBeforeDeg(end+1)=numBeforeDeg(end)+numCompositionsInDeg;

[pNew,qNew]=MacaulayMatrixSize(d0,d);

%Allocate space for the new rows.
MRows=zeros(pNew-p,qNew);

curRow=1;
for i=1:s%Going through all the polynomials
    termMat=termMats{i};%The current polynomial being manipulated.
    numTerms=size(termMat,2);

    degree=d0-d(i);
    
    %For each value n, all possible monomials are given by all
    %compositions of the degree into n parts (as n is the number of
    %variables). Of course, these compositions are defined such that 
    %the minimum value in each element is 1 and we want it to be zero,
    %so we have to adjust for that.

    numMonomials=binomial(degree+n-1,n-1);
    for j=0:(numMonomials-1)
        curMonomial=unrankTComposition(j,n,degree+n,true)-1;

        %Add the polynomial multiplied by the current monomial.
        for curTerm=1:numTerms
            monomial2Add=curMonomial+termMat(2:end,curTerm);

            %The degree of the term to add.
            deg=sum(monomial2Add);
            offset=numBeforeDeg(deg+1);

            idx=rankTComposition(monomial2Add+1,true)+offset+1;
            MRows(curRow,idx)=termMat(1,curTerm);
        end
        curRow=curRow+1;
    end
end
end

function numBeforeDeg=getNumElsBeforeDeg(n,dMax)
%A composition of n into d parts is the number of ways of splitting the
%integer n into the sum of d nonzero values. Compositions can also be
%offset to split n into d values, which are >=0. The lexicographic rank of
%a composition (starting from 0) is the offset within graded lexicographic
%ordering of a value with a particular fixed degree. We are going to have
%to keep track of these offsets to be able to index vectors.

numBeforeDeg=zeros(dMax+1,1);
numBeforeDeg(0+1)=0;
numBeforeDeg(1+1)=1;
for degree=1:dMax
    %The number of compositions of degree into n parts, where the parts can
    %be zero. We are including d0 as the function checkFordG will need it
    %even though this function only needs it up to degree d0-1.
    numCompositionsInDeg=binomial(degree+n-1,n-1);
    
    numBeforeDeg(degree+2)=numBeforeDeg(degree+1)+numCompositionsInDeg;
end

end

function [M,numBeforeDeg]=buildInitialMacaulayMat(termMats,d,d0,numBeforeDeg)
%%BUILDINITIASLMACAULAYMAT Construct the minimum-sized Macaulay matrix for
%             the system given. The degree of the matrix is the maximum
%             total degree of all of the polynomials. The Macaulay matrix
%             is defined in Section 5.1.2. This function is also used as a
%             convenience function to construct Sg.

s=length(termMats);%The number of polynomials
n=size(termMats{1},1)-1;%The number of variables

%Initially a Macaulay matrix of degree d0 is constructed. The number of
%rows p and columns q of the matrix are given 
[p,q]=MacaulayMatrixSize(d0,d);

%Allocate space for the d0th order Macaulay matrix.
M=zeros(p,q);

%Now, fill in the Macaulay matrix of order d0. Here, each polynomial is
%multiplied by all monomials of degree <=d0-d(i) and added. This includes
%the zeroth-degree, which contains
curRow=1;
for i=1:s%Going through all the polynomials
    termMat=termMats{i};%The current polynomial being manipulated.
    
    %First, add the zeroth-order term (the polynomial itself.
    numTerms=size(termMat,2);
    for curTerm=1:numTerms
        %The degree of the term to add.
        deg=sum(termMat(2:end,curTerm));
        offset=numBeforeDeg(deg+1);
        
        idx=rankTComposition(termMat(2:end,curTerm)+1,true)+offset+1;
        M(curRow,idx)=termMat(1,curTerm);
    end
    curRow=curRow+1;
    
    %Now, add all other terms associated with this polynomial.
    for degree=1:(d0-d(i))%All the orders for the current polynomial.
        %For each order, all possible monomials are given by all
        %compositions of degree into n parts (as n is the number of
        %variables). Of course, these compositions are defined such that 
        %the minimum value in each element is 1 and we want it to be zero,
        %so we have to adjust for that
        
        numMonomials=binomial(degree+n-1,n-1);
        for j=0:(numMonomials-1)
            curMonomial=unrankTComposition(j,n,degree+n,true)-1;
            
            %Add the polynomial multiplied by the current monomial.
            for curTerm=1:numTerms
                monomial2Add=curMonomial+termMat(2:end,curTerm);
                
                %The degree of the term to add.
                deg=sum(monomial2Add);
                offset=numBeforeDeg(deg+1);
                
                idx=rankTComposition(monomial2Add+1,true)+offset+1;
                M(curRow,idx)=termMat(1,curTerm);
            end
            curRow=curRow+1;
        end
    end
end
end

function [p,q]=MacaulayMatrixSize(d0,d)
%%%MACAULAYMATRIXSIZE The size of a d0th order Macaulay matrix when the
%                   polynomials have degree d(i). This implements Lemma 5.8
%                   in Section 5.2.1 of [1].

%The number of polynomials.
n=length(d);

q=binomial(n+d0,d0);

p=0;
for i=1:n
   p=p+binomial(n+d0-d(i),d0-d(i)); 
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
