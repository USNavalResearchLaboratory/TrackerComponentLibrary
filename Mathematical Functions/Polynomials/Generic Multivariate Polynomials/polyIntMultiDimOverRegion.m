function intVal=polyIntMultiDimOverRegion(coeffs,region,numVars,c1,c2)
%POLYINTMULTIDIMOVERREGION Perform a definite integral over a multivariate
%               polynomial (including cross terms) over a certain (real)
%               region with a specified (real) weighting function. These
%               are common regions for which explicit solutions are
%               available.
%
%INPUTS: coeffs A hypermatrix of the coefficients for the multivariate
%               polynomial. These are arranged such that
%               coeffs(a1,a2,a3...an) corresponds to the coefficient of an
%               x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1) term.  Thus, the
%               number of indices coeffs takes is equal to the
%               dimensionality of x (not counting singleton dimensions at
%               the end of coeffs). Note that this ordering is the reverse
%               that used in the 1D polyval function that is built into
%               Matlab. The number of elements for each index in coeffs is
%               the maximum order of that dimension +1.
%        region A parameter specifying the region over which the integral
%               should be taken. When a weighting function is present, that
%               means that the integral is of the weighting function times
%               the polynomial in coeffs. Possible values (where x(i) is
%               the ith variable -ith index in coeffs) are:
%               0 The n-dimensional cube -1<=x(i)<=1
%               1 The n-dimensional cubic shell c1<=|x(i)<=c2, 0<c1<c2<Inf
%               2 The n-dimensional sphere sum(x.^2)<=1 with weighting
%                 function sum(x.^2)^(c1/2) where c1>-numVars.
%               3 The n-dimensional shell c1<=sum(x.^2)<=1, 0<c1<1
%               4 The surface of the n-dimensional unit sphere sum(x.^2)=1
%               5 The n-dimensional cross polytope (n-dimensional
%                 octahedron) sum(abs(x))<=1.
%               6 The n-dimensional simplex sum(x)<=1, x(i)>=0
%               7 The entire n-dimensional space weighted by exp(-x'*x)
%               8 The entire n-dimensional space weighted by exp(-sqrt(x'x))
%       numVars A parameter specifying the number of variables present.
%               This is here to deal with singleton dimensions in coeffs.
%               When coeffs has more than two indices, Matlab automatically
%               suppresses those trailing ones having size 1. Also, given a
%               1D value, Matlab automatically adds a trailing singleton
%               dimension of 1. Thus, one cannot always determine the
%               number of variables over which integration should be
%               performed from the size of coeffs alone.
%         c1,c2 Constants that are used with regions 1,2, and 3 and are
%               otherwise ignored.
%
%OUTPUTS: intVal The value of the integral over the selected region.
%
%This evaluated the integral term-by term, using the formulae in Chapter 7
%of [1].
%
%As an example, consider the multivariate polynomial:
%14+3*x1^2-18*x2+12*x1*x3-3*x2*x3-7*x1^2*x2^2*x3^2+x3^4
% numVars=3;
% coeffs=zeros(3,3,5);
% coeffs(0+1,0+1,0+1)=14;
% coeffs(2+1,0+1,0+1)=3;
% coeffs(0+1,1+1,0+1)=-18;
% coeffs(1+1,0+1,1+1)=12;
% coeffs(0+1,1+1,1+1)=-3;
% coeffs(2+1,2+1,2+1)=-7;
% coeffs(0+1,0+1,4+1)=1;
% %Integral over the unit hypercube:
% intVal=polyIntMultiDimOverRegion(coeffs,0,numVars)
% %The result is intVal=16136/135.
% %Integral over the cubical shell from 1 to 10.
% intVal=polyIntMultiDimOverRegion(coeffs,1,numVars,1,10)
% %The result is intval=-(10285810968/5).
% %Integral over the weighted unit sphere with c1=1.
% intVal=polyIntMultiDimOverRegion(coeffs,2,numVars,1)
% %The result is intVal=737*pi/50.
% %Integral over the n-dimensional shell with c1=0.75.
% intVal=polyIntMultiDimOverRegion(coeffs,3,numVars,0.75)
% %The result is intVal=710568101*pi/61931520.
% %Integral over the surface of the n-dimensional unit sphere 
% intVal=polyIntMultiDimOverRegion(coeffs,4,numVars)
% %The result is intVal=908*pi/15.
% %Integral over the n-dimensional cross polytope
% intVal=polyIntMultiDimOverRegion(coeffs,5,numVars)
% %The result is intVal=108317/5670.
% %Integral over the n-dimensional simplex
% intVal=polyIntMultiDimOverRegion(coeffs,6,numVars)
% %The result is intVal=77699/45360.
% %Integral over the entire n-dimensional space weighted by exp(-x'*x)
% intVal=polyIntMultiDimOverRegion(coeffs,7,numVars)
% %The result is intVal=(123*pi^(3/2))/8.
% %Integral over the entire n-dimensional space weighted by exp(-sqrt(x'x))
% intVal=polyIntMultiDimOverRegion(coeffs,8,numVars)
% %The result is intVal=-9968*pi.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(coeffs))
    intVal=[];
    return;
end

numDims=size(coeffs);
numIdx=length(numDims);

switch(region)
    case 0%The n-dimensional cube -1<=x(i)<=1
        V=2^numVars;
        I=@(alpha)cubeMonomialInt(alpha);
    case 1%The n-dimensional cubic shell c1<=|x(i)<=c2, 0<c1<c2<Inf
        V=2^numVars*(c2^numVars-c1^numVars);
        I=@(alpha)cubeShellMonomialInt(alpha,c1,c2);
    case 2%The n-dimensional sphere sum(x.^2)<=1 with weighting function
          %sum(x.^2)^(c1/2) where c1>-numVars.
        V=2/(numVars+c1)*pi^(numVars/2)/gamma(numVars/2);
        I=@(alpha)spherMonomialInt(alpha,c1);
    case 3%The n-dimensional shell c1<=sum(x.^2)<=1, 0<c1<1
        V=2*pi^(numVars/2)/(numVars*gamma(numVars/2))*(1-c1^numVars);
        I=@(alpha)spherShellMonomialInt(alpha,c1);
    case 4%The surface of the n-dimensional unit sphere sum(x.^2)=1
        V=2*pi^(numVars/2)/gamma(numVars/2);
        I=@(alpha)spherSurfMonomialInt(alpha);
    case 5%The n-dimensional cross polytope (n-dimensional octahedron)
          %sum(abs(x))<=1
        V=exp(numVars*log(2)-gammaln(numVars+1));
        I=@(alpha)crossPolyMonomialInt(alpha);
    case 6%The n-dimensional simplex sum(x)<=1, x(i)>=0
        V=1/factorial(numVars);
        I=@(alpha)simplexMonomialInt(alpha);
    case 7%The entire n-dimensional space weighted by exp(-x'*x)
        V=pi^(numVars/2);
        I=@(alpha)exp2MonomialIntegral(alpha);
    case 8%The entire n-dimensional space weighted by exp(-sqrt(x'x))
        V=2*pi^(numVars/2)*exp(gammaln(numVars)-gammaln(numVars/2));
        I=@(alpha)expMonomialIntegral(alpha);
    otherwise
        error('Unknown region specified')
end

numEls=numel(coeffs);
alpha=zeros(numVars,1);

%The integral of each term shall be evaluated and then summed.
intVal=coeffs(1)*V;%The zeroth-order term.
for curEl=2:numEls
    idxVec=index2NDim(numDims,curEl);
    alpha(1:numIdx)=idxVec-1;
    intVal=intVal+coeffs(curEl)*I(alpha);
end
end

function val=cubeMonomialInt(alpha)
    if(any(mod(alpha,2)~=0))
        val=0;
    else
        numVars=length(alpha);
        val=2^numVars/prod(alpha+1);
    end
end

function val=cubeShellMonomialInt(alpha,c1,c2)
    if(any(mod(alpha,2)~=0))
        val=0;
    else
        numVars=length(alpha);
        alphaSum=sum(alpha);
        val=2^numVars*(c2^(numVars+alphaSum)-c1^(numVars+alphaSum))/prod(alpha+1);
    end
end

function val=spherMonomialInt(alpha,c1)
    if(any(mod(alpha,2)~=0))
        val=0;
    else
        numVars=length(alpha);
        alphaSum=sum(alpha);
        val=2/(numVars+c1+alphaSum)*exp(sum(gammaln((alpha+1)/2))-gammaln((numVars+alphaSum)/2));
    end
end

function val=spherShellMonomialInt(alpha,c1)
    if(any(mod(alpha,2)~=0))
        val=0;
    else
        numVars=length(alpha);
        alphaSum=sum(alpha);
        val=2/(numVars+alphaSum)*exp(sum(gammaln((alpha+1)/2))-gammaln((numVars+alphaSum)/2))*(1-c1^(numVars+alphaSum));
    end
end

function val=spherSurfMonomialInt(alpha)
    if(any(mod(alpha,2)~=0))
        val=0;
    else
        numVars=length(alpha);
        val=2*exp(sum(gammaln((alpha+1)/2))-gammaln((numVars+sum(alpha))/2));
    end
end

function val=crossPolyMonomialInt(alpha)
    if(any(mod(alpha,2)~=0))
        val=0;
    else
        numVars=length(alpha);
        val=exp(numVars*log(2)+sum(gammaln(alpha+1))-gammaln(numVars+sum(alpha)+1));
    end
end

function val=simplexMonomialInt(alpha)
    numVars=length(alpha);
    val=exp(sum(gammaln(alpha+1))-gammaln(numVars+sum(alpha)+1));
end

function val=exp2MonomialIntegral(alpha)
    if(any(mod(alpha,2)~=0))
        val=0;
    else
        val=prod(gamma((alpha+1)/2));
    end
end

function val=expMonomialIntegral(alpha)
    if(any(mod(alpha,2)~=0))
        val=0;
    else
        numVars=length(alpha);
        alphaSum=sum(alpha);
        val=2*exp(gammaln(numVars+alphaSum)+sum(gammaln((alpha+1)/2))-gammaln((numVars+alphaSum)/2));
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
