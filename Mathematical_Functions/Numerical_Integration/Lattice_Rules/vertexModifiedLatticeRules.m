function [xi,w,varMin]=vertexModifiedLatticeRules(numDim,method,tabIdx,lowerBounds,upperBounds)
%%VERTEXMODIFIEDLATTICERULES Obtain lattice points for multidimensional
%               numerical integration based on Niederreiter and Sloan
%               vertex modified lattice rules for non periodic integrands.
%               The standard points are meant for integrating a function
%               from 0-1 in all dimensions. Inputs can change the
%               integration region to other finite or infinite bounds.
%
%INPUTS: numDim The number of dimensions of the lattice points that are to
%               be generated. numDim>=2. Valid values depend on the method
%               chosen.
%       methods The method of obtaining the lattice points that should be
%               used. Possible values are
%               0 Use tabulated points taken from Appendix B of [1]. These
%                 are available for numDim=2 to 12. The next input tabIdx
%                 is a value from 1 to 8 specifying the formula in the
%                 table to use. Higher indices produce lattices with more
%                 points. The number of points is 2^numDim+N-1. For
%                 tabIdx=1:2, N is 79, 157, 313, 619, 1249, 2503, 5003, and
%                 10007.
%               1 Determine the optimal vertex modified lattice by
%                 minimizing the vertex variance of Equation 8.13 in [1].
%                 This uses Equation 8.10 for the vertices and the formulae
%                 of Chapter 8.3 and Apendix B. In this instance, tabIdx is
%                 a prime number that affects the number of points. The
%                 Matlab function primes can be used to find prime numbers.
%                 This option allows the use of untabulated numbers of
%                 dimensions and points.
%        tabIdx A parameter whose value depends on the method selected. See
%               the comments to the method input for more information.
%lowerBounds,upperBounds If these p[arameters are provided (and are not
%               empty matrices), then the lattice points will be
%               transformed for integration over regions other than just
%               +/-1 in all dimensions. These are numDimX1 or 1XnumDim
%               vectors where lowerBounds are the lower integration bounds
%               and upperBounds are the upper integration bounds. If
%               transformed, the following substitutions are used:
%               Considering only the current dimensions of xi where a is
%               the lower bound and b is the upper bound:
%               * If a and b are finite:
%                 xi(curDim,:)=(b-a)*xi(curDim,:)+a, w=(b-a)*w
%               * If a=-Inf and b is finite:
%                 xi(curDim,:)=-1./xi(curDim,:)+b+1; w=w./xi(curDim,:).^2;
%                 --Note this may cause problems if any elements of xi are
%                 zero.
%               * If a is finite and b=Inf:
%                 xi(curDim,:)=1./(1-xi(curDim,:))+a-1;
%                 w=w./(1-xi(curDim,:)).^2;
%                 --Note this may cause problems if any elements of xi are
%                 1.
%               * If a=-Inf and b=Inf:
%                 xi(curDim,:)=1./(1-xi(curDim,:))-1./xi(curDim,:);
%                 w=w.*(1./(1-xi(curDim,:)).^2+1./xi(curDim,:).^2);
%                 --Note this may cause problems if any elements of xi are
%                 zero or 1.
%               * Same as any above, but the lower bound is greater than
%                 the upper bound:
%                 The bounds are flipped, but the same substitution is
%                 made, just w is multiplied by -1.
%
%OUTPUTS: xi A numDimXnumPoints matrix of the lattice points generated
%            using the selected method.
%          w A numPointsX1 set of weights associated with the lattice
%            points.
%     varMin The vertex variance associated with the points. This is from
%            Equation 8.13 in Chapter 8.3 of [1].
%
%EXAMPLE:
%Consider the examples in 4D of Section 6.6 of [1]. We start with the
%problem of integrating over the 4-dimensional cube of the function
%exp(x1*x2*x3*x4).
% %The function
% f=@(x)exp(prod(x,1));
% %Integrated over +/-1 in all dimensions, for 4 dimensions has an "exact"
% %solution of
% exactSol=1.0693976088597706235;
% %(Which comes from evaluating a hypergeometric function)
% %First, consider using tabulated formulae.
% numDim=4;
% tabIdx=2;
% [xi,w]=vertexModifiedLatticeRules(numDim,0,tabIdx);
% est1=sum(bsxfun(@times,f(xi),w'));
% error1=abs(exactSol-est1)
% %Next deriving the optimal formula for the same number of points.
% N=length(w)-2^numDim+1;%A prime number.
% [xi,w,l]=vertexModifiedLatticeRules(numDim,1,N);
% est2=sum(bsxfun(@times,f(xi),w'));
% error2=abs(exactSol-est2)
%One might notice that the "optimal" solution is not identical and that
%it is marginally worse in performance. The difference in the cost
%function was only by a few digits. The tabulated values are
%"recommended" values, not necessarily those found numerically optimal
%using double precision.
%
%EXAMPLE 2:
%Here, we consider integration over a region that is not the unit cube. We
%consider the function sum(abs(x-1/2)) in five dimensions with bounds as
%follows:
% f=@(x)sum(abs(x-1/2),1);
% lowerBounds=[-2;-1;0;4;3];
% upperBounds=[2;1;1;6;8];
% %The exact solution against which we will compare is easily found:
% exactSol=915;
% numDim=5;
% method=0;
% tabIdx=1;
% [xi,w]=vertexModifiedLatticeRules(numDim,method,tabIdx,lowerBounds,upperBounds);
% est1=sum(bsxfun(@times,f(xi),w'));
% error1=abs(exactSol-est1)
% %Increasing the number of points used, we get
% tabIdx=4;
% [xi,w]=vertexModifiedLatticeRules(numDim,method,tabIdx,lowerBounds,upperBounds);
% est2=sum(bsxfun(@times,f(xi),w'));
% error2=abs(exactSol-est2)
%And one observed that error2 is less than error1, as one would expect when
%increasing the number of points used.
%
%REFERENCES:
%[1] I. H. Sloan and S. Joe, Lattice Methods for Multiple integration.
%    Oxford, United Kingdom: Clarendon Press, 1994.
%[2] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

s=numDim;

switch(method)
    case 1%If the points should be found via optimization.

        N=tabIdx;
        if(~isprime(N))
            error('The number of points must be prime.');
        end

        [xi,w,varMin]=findOptRule(s,tabIdx);
    case 0%IF tabulated points and weights should be used.


    switch(s)
        case 2
            paramTable=[79,     15;
                        157,    28;
                        313,    25;
                        619,    74;
                        1249,   585;
                        2503,   340;
                        5003,   2318;
                        10007,  4722];
        case 3
            paramTable=[79,     23;
                        157,    22;
                        313,    51;
                        619,    239;
                        1249,   507;
                        2503,   521;
                        5003,   1086;
                        10007,  4073];
        case 4
            paramTable=[79,     29;
                        157,    22;
                        313,    51;
                        619,    173;
                        1249,   338;
                        2503,   183;
                        5003,   1033;
                        10007,  1481];
        case 5
            paramTable=[79,     15;
                        157,    22;
                        313,    38;
                        619,    132;
                        1249,   377;
                        2503,   814;
                        5003,   2036;
                        10007,  1481];
        case 6
            paramTable=[79,     15;
                        157,    22;
                        313,    59;
                        619,    132;
                        1249,   143;
                        2503,   347;
                        5003,   708;
                        10007,  357];
        case 7
            paramTable=[79,     6;
                        157,    42;
                        313,    59;
                        619,    53;
                        1249,   131;
                        2503,   530;
                        5003,   687;
                        10007,  1855];
        case 8
            paramTable=[79,     15;
                        157,    42;
                        313,    59;
                        619,    53;
                        1249,   235;
                        2503,   223;
                        5003,   687;
                        10007,  1855];
        case 9
            paramTable=[79,     28;
                        157,    36;
                        313,    59;
                        619,    53;
                        1249,   235;
                        2503,   223;
                        5003,   687;
                        10007,  2521];
        case 10
            paramTable=[79,     28;
                        157,    36;
                        313,    59;
                        619,    53;
                        1249,   235;
                        2503,   372;
                        5003,   1292;
                        10007,  4474];
        case 11
            paramTable=[79,     28;
                        157,    36;
                        313,    62;
                        619,    40;
                        1249,   235;
                        2503,   372;
                        5003,   1180;
                        10007,  4474];
        case 12
            paramTable=[79,     28;
                        157,    66;
                        313,    62;
                        619,    40;
                        1249,   328;
                        2503,   140;
                        5003,   201;
                        10007,  3708];
        otherwise
            error('Untabulated number of dimensions given')
    end

    N=paramTable(tabIdx,1);
    l=paramTable(tabIdx,2);

    i=(0:(s-1))';
    z=mod(l.^i,N);

    numPoints=2^s+N-1;

    %Allocate space
    w=zeros(numPoints,1);
    xi=zeros(s,numPoints);

    %First, fill in the N-1 points.
    for curPoint=1:(N-1)
        w(curPoint)=1/N;
        xi(:,curPoint)=mod(curPoint/N*z,1);
    end

    %Next, fill in the vertex points.
    maxIdxVals=ones(s,1);
    for i=1:(2^s)
        curPoint=curPoint+1;
        %The values of the indices of the sums in Equation 8.8 in [1].
        iVals=unrankTuple(i-1,maxIdxVals);
        xi(:,curPoint)=iVals;

        %Construct the current w using Equation 8.10 in [1].
        sumVal=0;

        for j=1:(N-1)
            lCur=1;
            for curlIdx=1:s
                if(iVals(curlIdx)==1)
                    a=1-xi(curlIdx,j);
                else
                    a=xi(curlIdx,j);
                end
                lCur=lCur*a;
            end
            sumVal=sumVal+lCur;
        end
        w(curPoint)=1/2^s-sumVal/N;
    end

    varMin=sum(w(N:end).^2);
end


%If lower and upper bounds are provided.
if(nargin>3&&~isempty(lowerBounds))
    for curDim=1:numDim
        a=lowerBounds(curDim);
        b=upperBounds(curDim);
        
        if(a<b)
            flipSign=1;
        else%Make a<b and flip the sign of the result.
            temp=a;
            a=b;
            b=temp;
            flipSign=-1;
        end
        %Now, a<b.

        if(isfinite(a)&&isfinite(b))%Both bounds finite
            xi(curDim,:)=a+(b-a)*xi(curDim,:);
            w=flipSign*w*(b-a);
        elseif(~isfinite(a)&&isfinite(b))%a=-Inf
            xi(curDim,:)=-1./xi(curDim,:)+b+1;
            w=flipSign*w./xi(curDim,:).^2;
        elseif(isfinite(a)&&~isfinite(b))%b=Inf
            xi(curDim,:)=1./(1-xi(curDim,:))+a-1;
            w=flipSign*w./(1-xi(curDim,:)).^2;
        else%a=-Inf and b=Inf
            xi(curDim,:)=1./(1-xi(curDim,:))-1./xi(curDim,:);
            w=flipSign*w.*(1./(1-xi(curDim,:)).^2+1./xi(curDim,:).^2);
        end
    end
end
end

function [xiMin,wMin,varMin]=findOptRule(s,N)
%%FINDOPTRULE Find the optimum vertex modified lattice rule. This minimizes
%             Equation 8.13 of Chapter 8.3 of [1] via brute force. As
%             noted in Appendix B, only values of l from 1:N/2 need be
%             searched.
%
%INPUTS:      s The number of dimensions of the lattice points that are to
%               be generated. numDim>=2.
%             N A prime number that affects how many lattice points are
%               generalted. In all 2^s+N-1 lattice points will be
%               generated.
%
%OUTPUTS: xiMin A numDimXnumPoints matrix of the optimal lattice points.
%          wMin A numPointsX1 set of weights associated with the optimal 
%               lattice points.
%        varMin The vertex variance associated with the points. This is
%               from Equation 8.13 in Chapter 8.3 of [1]. 
%
%REFERENCES:
%[1] I. H. Sloan and S. Joe, Lattice Methods for Multiple integration.
%    Oxford, United Kingdom: Clarendon Press, 1994.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
numPoints=2^s+N-1;

%Allocate space
w=zeros(numPoints,1);
xi=zeros(s,numPoints);

%The first N-1 weights do not depeond on the optimization.
w(1:(N-1))=1/N;

varMin=Inf;
for l=1:fix(N/2)
    i=(0:(s-1))';
    z=mod(l.^i,N);
    
    %First, fill in the N-1 points.
    for curPoint=1:(N-1)
        xi(:,curPoint)=mod(curPoint/N*z,1);
    end
    
    %Next, fill in the vertex points and weights.
    maxIdxVals=ones(s,1);
    for i=1:(2^s)
        curPoint=curPoint+1;
        %The values of the indices of the sums in Equation 8.8 in [1].
        iVals=unrankTuple(i-1,maxIdxVals);
        xi(:,curPoint)=iVals;

        %Construct the current w using Equation 8.10 in [1].
        sumVal=0;

        for j=1:(N-1)
            lCur=1;
            for curlIdx=1:s
                if(iVals(curlIdx)==1)
                    a=1-xi(curlIdx,j);
                else
                    a=xi(curlIdx,j);
                end
                lCur=lCur*a;
            end
            sumVal=sumVal+lCur;
        end
        w(curPoint)=1/2^s-sumVal/N;
    end
    
    %Finally, compute the vertes variance.
    varCur=sum(w(N:end).^2);
    
    if(varCur<varMin)
        varMin=varCur;
        wMin=w;
        xiMin=xi;
        lMin=l;
    end
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
