function [xi,w]=GenzKeisterPoints(numDim,m,algorithm,epsVal,randomize)
%%GENZKEISTERPOINTS Generate multiple dimensional cubature points and
%             weights for integration involving a multidimensional Gaussian
%             probability density function (PDF) using the method of Genz
%             and Keister in [1].
%
%INPUTS: numDim A positive integer specifying the dimensionality of the
%               points to be generated.
%             m An integer related to the number of generator points used.
%               This must be >0 and <=17 for algorithm 0 and <=15 for
%               algorithm 1 unless algorithm provides a different set of
%               parameters from which higher order generators can be found.
%     algorithm An optional parameter indicating how the cubature points
%               are generated. Possible values are
%               0 (The default if omitted or an empty matrix is passed)
%                 Use the Q(25)[1+2+6+10+16] generator that is given in
%                 the first column of Table 3.4 in [1].
%               1 Use the Q(25)[1+2+8+20] generator that is given in the
%                 second column of Table 3.4 in [1].
%              nu This is a vector of increment values nu as in [1] so that
%                 a new generator can be computed. The generator using
%                 these values is computed based off the first two
%                 generator points being fixed at 0 and sqrt(3). More
%                 details on this input are given in the comments below.
%       epsVal An optional parameter specifying how small the absolute
%              value of a weight must be before the weight and associated
%              cubature point are omitted. The generators given in the
%              paper were designed so that many weights are numerically
%              zero, meaning that cubature formula need far fewer points.
%              If this parameter is omitted or an empty matrix is passed,
%              the default value of eps(1) is used.
%    randomize If this parameter is true, then the points will be
%              multiplied by a random orthonormal rotation matrix. This
%              does not change the moments up to the order of the points.
%              This randomization is done in [3] and [4] to lessen various
%              effects that arise when using points in the same
%              orientation repeatedly in tracking. The default if this
%              parameter is omitted or an empty matrix is passed is false.
%              
%OUTPUTS: xi A d X numCubaturePoints matrix containing the cubature points.
%            (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%The algorithm of [1], which is an extension of that of [2] is used. For
%more details on how to use these points when integrating over a Gaussian
%PDF, see the comments in the function fifthOrderCubPoints.m.
%
%The algorithm is based on initially determining a set of generating
%coefficients lambda, where the first two are fixed at 1 and sqrt(3).
%The coefficients are determined by recursively expanding the number of
%coefficients in intervals. The generator Q_P obtained the coefficients by
%expanding them in intervals of nu=[3;5;8]. The generator \hat{Q}_P
%expanded them in intervals of nu=[4;10]; The total number of points for an
%interval vector nu is 2+sum(nu) Though the option is given to generate
%other coefficients, which would be done using the 
%getGenzKeisterGenerators subroutine in this function, this is generally a
%bad idea. First, not all intervals produce valid coefficients, so an error
%can occur when a bad combination is attempted. More importantly, however,
%is the fact that extended precision is usually needed to compute the
%generators to a reasonable accuracy. The tabulated values used for
%algorithm=0 and algorithm=1 were computed using extended precision by the
%authors of [1] and are taken from the paper. Though the subroutine
%getGenzKeisterGenerators can recreate the results, there is a severe loss
%of precision with the higher-order lamba terms.
%
%REFERENCES:
%[1] A. Genz and B. D. Keister, "Fully symmetric interpolatory rules for
%    multiple integrals over infinite regions with Gaussian weight,"
%    Journal of Computational and Applied mathematics, vol. 71, no. 2, pp.
%    299-309, Jul. 1996.
%[2] A. Genz, "Fully symmetric interpolatory rules for multiple integrals,"
%    SIAM Journal on Numerical Analysis, vol. 23, no. 6, pp. 1273-1283,
%    Dec. 1986.
%[3] O. Straka, D. Duník, and M. Simandl, "Randomized unscented Kalman
%    filter in tracking," in Proceedings of the 15th International
%    Conference on Information Fusion, Singapore, 9-12 Jul. 2012, pp.
%    503-510.
%[4] J. Duník, O. Straka, and M. Simandl, "The development of a randomised
%    unscented Kalman filter," in Proceedings of the 18th World Congress,
%    The International Federation of Automatic Control, Milan, Italy, 28
%    Aug. - 2 Sep. 2011, pp. 8-13.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
   algorithm=0; 
end

if(nargin<4||isempty(epsVal))
    epsVal=eps(1);
end

if(nargin<5||isempty(randomize))
    randomize=false;
end

n=numDim;

%First, we shall determine the maximum number of points that might be
%required so that we can allocate memory for xi and w. The formula is given
%in the text after Equation 1.
maxNumPoints=0;
for curCard=0:m
    pRecur=getNextMPartition(n+curCard,n);
    while(~isempty(pRecur))%While we haven't passed the final partition.
        %We need to determine the multiplicities of all of the unique
        %elements in pRecur that are not zero. These are the i-terms in the
        %text.
        p=pRecur(pRecur~=1);
        numReps=diff(find(diff([Inf;sort(p(:),'descend');-Inf])));
        %The inner diff([Inf;E;-Inf]) marks where repetitions end. The
        %outer find finds the indices of the nonzero elements. The
        %differences between those indices is the number of elements in
        %each repetitions. This is efficient as it avoids loops in Matlab,
        %but it is an inefficient method of determining the multiplicity of
        %the unique elements in E.
        
        iCard=sum(numReps);
        
        %factorial(x)=gamma(x+1). The logarithms of the gamma functions are
        %used for the factorials to avoid overflow in the intermediate
        %results.
        maxNumPoints=fix(maxNumPoints+(2^iCard)*exp(gammaln(n+1)-sum(gammaln(numReps+1))-gammaln(n-iCard+1)));

        pRecur=getNextMPartition(pRecur);
    end
end

%Allocate space for the points and weights
xi=zeros(n,maxNumPoints);
w=zeros(maxNumPoints,1);

switch(algorithm)
    case 0%Generator for Q_P
        lambda=[0;
                1.7320508075688773;
                4.1849560176727319;
                0.74109534999454084;
                2.8612795760570581;
                6.3633944943363700;
                1.2304236340273060;
                5.1870160399136561;
                2.5960831150492022;
                3.2053337944991945;
                9.0169397898903025;
                0.24899229757996061;
                7.9807717985905609;
                2.2336260616769417;
                7.1221067008046167;
                3.6353185190372782;
                5.6981777684881096;
                4.7364330859522971];
        %The function computeAValues could be used to find the values of a
        %in Equation 1. However, the values below were computed using
        %extended precision given the 16-bit values of lambda in the paper.
        a=[1;
           1;
           0;
           6;
           -48.378475125832451;
           0;
           0;
           0;
           34020;
           -986064.53173677489;
           0;
           0;
           0;
           0;
           0;
           1.2912054173706603e12;
           -1.1268664521456168e14;
           2.9248520348796280e15];
    case 1 %Generator for \hat{Q}_P
        lambda=[0;
                1.7320508075688773;
                4.9791465117195582;
                0.84628809835102170;
                3.7355715460409573;
                2.6840395601585692;
                9.0508037980317400;
                0.47371420996884380;
                8.0130130598043254;
                1.2435457006528093;
                7.1482776511870860;
                2.2210157242456798;
                6.3725842092196923;
                3.1782891110545301;
                5.6545621267720157;
                4.3394221426603945]; 
       a=[1;
          1;
         0;
         6;
         -93.0486211834777976;
         504.496566347049718;
         0;
         0;
         0;
         0;
         1.93536000000001776e6;
         -2.38763644847775079e8;
         4.62442819320708296e9;
         8.22485150843440875e9;
         -4.90446886942675039e12;
         7.61797098142559229e13];
    otherwise%Create the generator from the sequence of nu values.
        lambda=getGenzKeisterGenerators(algorithm);
        a=computeAValues(lambda);
end

%In the unnumbered equation for computing the weights, before Equation 1 in
%[1], which is the same as Equation 2.4 in [2], there is a term involving
%the ratio of a to a product of differences of lambda's. Those terms are
%the same for all pairs of (p(i),k(i)+p(i)). Thus, here, we compute all of
%the values ahead of time and put them into a table that can be addressed
%as innerTerms4Weights(p(i)+1,k(i)+p(i)+1), where the plus one terms are
%because Matlab addresses things from 1 instead of from 0. This table os
%passed to the function to compute the weights so as to avoid needed
%recomputation of the same terms.
innerTerms4Weights = zeros(m+1,m+1);
innerTerms4Weights(1,1) = a(1);
lambda2=lambda.*lambda;
for pSubi = 1:(m+1)
    prodVal=1;
    for pkSum=2:(m+1)
        if(pkSum>pSubi)
            prodVal=prodVal*(lambda2(pSubi) - lambda2(pkSum));
        else
            prodVal=prodVal*(lambda2(pSubi) - lambda2(pkSum-1));
        end
        
        if(pkSum>=pSubi)
            innerTerms4Weights(pSubi,pkSum)=a(pkSum)/prodVal;
        end
    end
end

numPoints=0;
%The loops below implement the sums in the equation for Q^{(m,n)}(f) on
%page 300.

%The outer loop is over the cardinality of the partitions p
for curCard=0:m
    %The second loop goes through all partitions of curCard+1 elements into
    %n parts. However, we want the partitions to be from 0, not 1. This
    %means we must increase the set from which we draw by n. Then pRecur-1
    %is the desired permutation set p. The partitions are produced in
    %colexicographic order, which means that the ordering of the magnitudes
    %of p(1) to p(n) match that desired in the paper.
    pRecur=getNextMPartition(n+curCard,n);

    while(~isempty(pRecur))%While we haven't passed the final partition.
        p=pRecur-1;%Get the actual partition with elements starting from zero.
    
        %The weight for all of the values involved in this partition.
        wp=computeW(m,p,innerTerms4Weights);
        %Now, we must go through all possible permutations of the elements
        %in p. Because elements are repeated, these are multiset
        %permutations.
        pPerm=genAllMultisetPermutations(p);

        numPerm=size(pPerm,2);
        for i=1:numPerm
            curPerm=pPerm(:,i);
            %For the current multiset permutation of the lambdas, get all
            %possible +/- combinations of the lambdas.
            xiNew=PMCombos(lambda(curPerm+1));
            numVals=size(xiNew,2);
            
            %Add the points.
            xi(:,(numPoints+1):(numPoints+numVals))=xiNew;
            w((numPoints+1):(numPoints+numVals))=wp;
            numPoints=numPoints+numVals;
        end

        pRecur=getNextMPartition(pRecur);
    end
end

%Get rid of points with weights that are numerically zero.
sel=abs(w)>epsVal;
w=w(sel);
xi=xi(:,sel);

if(randomize)
    R=randOrthoMat(numDim);

    numPoints=length(w);
    for curPoint=1:numPoints
        xi(:,curPoint)=R*xi(:,curPoint);
    end
end
end

function w = computeW(m,p,innerTerms4Weights)
%%COMPUTEW Compute the weight w for a given partition using the unnumbered
%          equation before Equation 1 in [1], which is Equation 2.4 in [2].
%          The efficient method of going through the values of k without
%          explicitly computing k is based on correspondece with Dr. Genz
%          of [1], and [2].

n=length(p);
    
%K is defined before Equation 1 in [1].
K=sum(p~=0);

%The vector k is never generated. Rather, the kpSum, the sum of k and p is
%tracked. If we start with k=0, then that means that the cardinality of k
%can increase by m-sum(p). The cardinalityLeft value is how much
%cardinality we have left for
cardinalityLeft=m-sum(p);

%If k is the maximum cardinality, then the maximum value of k(i)+p(i) is
%such that the sum of the cardinalities equals m. The array kpSum will hold
%the values of k(i)+p(i) for all i. If k=0, the initial implicit value of
%k, then k+p=p.
kpSum=p;

%prodSums accumulates the values of the products over i AND accumulates the
%sum of the products over all k.
prodSums=zeros(1,n+1); 
while(1)
    for i=1:n
        prodSums(1)=1;
        if(cardinalityLeft>=0)
            prodSums(i+1)=prodSums(i+1)+innerTerms4Weights(p(i)+1,kpSum(i)+1)*prodSums(i);
            prodSums(i)=0;
            kpSum(i)=kpSum(i)+1;
            cardinalityLeft=cardinalityLeft-1;
            break;
        end
        
        cardinalityLeft=cardinalityLeft+kpSum(i)-p(i);
        kpSum(i)=p(i);
    end

    if(i==n&&kpSum(n)==p(n))
        break;
    end
end
    
w=2^(-K)*prodSums(n+1);
end


function a=computeAValues(lambda)
%%COMPUTEAVALUES This computes all of the a values from Equation 1 in [1]
%           given the generator set. Note that finite precision errors will
%           add up quickly.
%
%Expanding the product inside of the integral inquestion, one can see that
%there are terms involving moments of the standard normal N(0,1)
%distribution times products involving lambda terms. The nth moment of the
%standard normal distribution is doubleFactorial(n-1) if n is even and is
%0 if n is odd. Thus, we can evaluate the function term by term. This is
%done here by expanding the polynomial in x via convolutions and then
%evaluating the values of the moments for appropriate terms.

numVals=length(lambda);
a=zeros(numVals,1);%Space for the return values.

%The maximum degree of the polynomial is 2*numVals. However, since we know
%that all odd moments are zero, we only need to save the even moments
%(including the zeroth-order moment).
x2Moments=zeros(numVals,1);
%Precompute all of the necessary moments
x2Moments(0+1)=1;%The zeroth-order moment is 1.
for n=1:numVals
    x2Moments(n+1)=doubleFactorial(2*n-1);
end

%Polynomials only have squared lambda terms in them.
lambda2=lambda.*lambda;

%Because the polynomial is only in terms of x^2 terms, there will be no odd
%terms. Thus, the polynomial is build from second moments (we do not save
%the even moments). For simplicity, we use the opposite ordering of the
%polyval function such that the lowest element in the polynomial vector is
%the zero-th order element. 

p=zeros(numVals+1,1);%To hold the polynomial terms.

p(1:2)=[-lambda2(1);1];
a(1)=1;%The zeroth-order moment.
a(2)=sum(p(1:2).*x2Moments(1:2));
for curMult=2:(numVals-1)
    %Add in the other product term.
    p(1:(curMult+1))=conv(p(1:curMult),[-lambda2(curMult);1]);
    %Evaluate the integrals.
    a(curMult+1)=sum(p(1:(curMult+1)).*x2Moments(1:(curMult+1)));
end

end

function lambda=getGenzKeisterGenerators(nuList)
%%GETGENZKEISTERGENERATORS Obtain the generators (the lambda terms).
%
%Use the method of Section 2 to find the generators (lambdas) for
%Genz-Keister's algorithm for a given number of generators.
%nu must be at least 2.

numStages=length(nuList);
maxOrder=2+sum(nuList)*2;

%The maximum degree of the polynomial in x in the second to last equation
%on page 301 is 2*maxOrder. However, since we know that all odd moments are
%zero, we only need to save the even moments. We include the zeroth-order
%moment, because the polynomial will have a zeroth-order term.
x2Moments=zeros(maxOrder,1);
%Precompute all of the necessary moments
x2Moments(0+1)=1;%The zeroth-order moment is 1.
for n=1:maxOrder
    x2Moments(n+1)=doubleFactorial(2*n-1);
end

lambda=zeros(2+sum(nuList),1);
lambda(1)=1;%lamda_0
lambda(2)=sqrt(3);%lambda_1
mu=2;

%The accumulated polynomial consisting of known values of lambda. This is
%the product from 0 to mu in the un-numbered equation in the middle of
%page 3.
polyAccum=zeros(3,1);
polyAccum(3)=1;%x^(2*2) term.
polyAccum(2)=-3;%x^2 term
polyAccum(1)=0;%x^0 term

for curAdd=1:numStages
    nu=nuList(curAdd);%The number of points to add.
    
    %We need to keep track of the coefficients of the s term. These are
    %determined below by the powers of x multiplying the s terms, because
    %once the integral in the equation in the middle of page 3 is
    %evaluated, all of the x terms go away and the solution is linear in s.
    %The constant term is the negative of the one that has no s terms.
    SCoeffMat=zeros(nu,nu);
    constTerm=zeros(nu,1);

    for k=0:(nu-1)        
        %First, we evaluate the constant term. This is x^(2*nu) multiplied
        %by the powers in curPolyX, taking into account the offset by the
        %x^(2*k) term.
        accumLen=length(polyAccum);
        sumVal=0;
        for curVal=1:accumLen
            sumVal=sumVal+x2Moments(curVal+nu+k)*polyAccum(curVal);
        end
        constTerm(k+1)=-sumVal;
        
        %The accumulated polynomial in x is multiplied by the polynomial in
        %the unknowns s. We will perform the multiplication for each s term
        %and put the result into the SCoeffMat.
        for curS=1:nu
            sumVal=0;
            for curVal=1:accumLen
                sumVal=sumVal+x2Moments(curVal+k+curS-1)*polyAccum(curVal);
            end
            SCoeffMat(k+1,curS)=sumVal;
        end
    end
    %Now, solve the linear equation to get the s terms.
    s=SCoeffMat\constTerm;
    
    %Put the s terms into an array so that the roots of the s equation can
    %be found. This means inserting zeros
    sTotal=zeros(2*nu+1,1);
    sTotal(1)=1;
    for i=1:nu
       sTotal(2*nu-(2*(i-1))+1)=s(i);
    end
    
    %The new lambda values are the roots of the s terms.
    lambdaNew=roots(sTotal);
    if(any(imag(lambdaNew)~=0))
        error('The sequence of nu chosen to determine the generators lambda led to an imaginary generator.')
    end
    
    %The roots should appear in pairs of opposite sign. We will just keep
    %the positive ones. However, that is easier said than done due to
    %finite precision errors in the computation of the roots (i.e. one
    %cannot just use the unique function after taking the absolute value).
    %Thus, we will choose closest pairs and average them.
    lambdaNew=abs(lambdaNew);
    lambdaNewChosen=zeros(nu,1);
    unassignedVals=true(2*nu,1);
    numChosen=0;
    for curChosen=1:nu
        %Find the next unassigned value and assign it.
        temp=find(unassignedVals==true);
        idx2Assign=temp(1);
        unassignedVals(idx2Assign)=false;
        
        %Now, find the next closest point. This is an inefficient way to do
        %it, but it lets us avoid using a loop in matlab.
        [~,pairIdx]=min(abs(lambdaNew-lambdaNew(idx2Assign))./(unassignedVals));
        numChosen=numChosen+1;
        lambdaNewChosen(numChosen)=(lambdaNew(idx2Assign)+lambdaNew(pairIdx))/2;
        unassignedVals(pairIdx)=false;
    end
    lambda((mu+1):(mu+nu))=lambdaNewChosen;
    
    %Now, add in the s terms to the accumulated polynomial.
    polyAccum=conv(polyAccum,[s;1]);
    
    %Record the total number of terms present.
    mu=mu+nu;
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
