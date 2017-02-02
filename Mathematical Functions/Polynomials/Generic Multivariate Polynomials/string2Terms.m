function termMat=string2Terms(theString,numVars)
%%STRING2TERMS Given a string representing a real multidimensional
%           polynomial equation, convert it into a matrix of terms of the
%           polynomial.
%
%INPUTS: theString The string representing the multivariate polynomial with
%           real coefficients. All coefficients come before variables.
%           Variables are all x followed by a number, such as x3. All
%           numbers for variables must be >0. No parenthesis should be
%           used. This formatting is less permissive than in string2Terms.
%           For example, 
%           '60*x2^2-4*x2^3+72*x1*x2+12' is a valid string.
%   numVars This is an optional parameter indicating the number of
%           spaces to leave for variables in termMat. This must be equal to
%           or larger than the highest number of the variable in the
%           string. This is used to make termMat have a larger number of
%           indices for exponents (rows) than the maximum exponent number.
%           If omitted or an empty matrix is passed, only enough spaces for
%           the variables in theString is used. numVars less than the
%           required minimum amount of space means that numVars will be
%           ignored.
%
%OUTPUTS: termMat An (n+1)XnumTerms matrix such that
%           termMat(:,i)=[c,a1,a2,...,an] is a monomial term where c is the
%           value of of the monomial coefficient and the monomial is
%           x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1).
%           n is the highest numbered used to index a variable in
%           theString, or is numVars, whichever is larger. Thus, one
%           should not use needlessly high indices in the string.
%
%EXAMPLE 1:
%  theString='x1*x3-2*x2*x3-6*x1-x2-4*x4+12';
%  termMat=string2Terms(theString)
%Here termMat is 5X1, having space for coefficients of each of the four
%variables.
%
%EXAMPLE 2:
%  theString='x2*x3-2*x2*x3-x2-4*x4+12';
%  termMat=string2Terms(theString)
%Here, termMat is again 5X1, even though there are three variables (x1 is
%gone). This is because the coefficient number corresponds to a position.
%
%EXAMPLE 3
%  theString='x2*x3-2*x2*x3-x2-4*x4+12';
%  termMat=string2Terms(theString,5)
%This is the same as the second example, except there is an extra empty
%row.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(numVars))
   numVars=0; 
end

[monomials,signs]=strsplit(theString,{'+','-'});
%If the first element in the string was a sign like +/-
if(isempty(monomials{1}))
   monomials(1)=[]; 
end
numMonomials=length(monomials);

%Add an initial + if no leading +/- was provided.
if(numMonomials>length(signs))
    signs=['+',signs];
end

signs=cell2mat(signs);

%We will convert the signs to +/-1
signVals=zeros(numMonomials,1);
sel=(signs=='+');
signVals(sel)=1;
signVals(~sel)=-1;

%First, scan through and split the monomials into parts around the
%multiplication and extract the coefficients.
coeffs=zeros(numMonomials,1);
maxXIdx=0;
for curMonomial=1:numMonomials
    splitMonomial=strsplit(monomials{curMonomial},'*');
    
    %Extract the coefficient.
    val=sscanf(splitMonomial{1},'%f');
    
    %If the coefficient is just +/-one.
    if(isempty(val))
        coeffs(curMonomial)=signVals(curMonomial);
    else
        coeffs(curMonomial)=signVals(curMonomial)*val;
        splitMonomial(1)=[];
    end
    
    %Next, get the index of each term an its order
    numSplitMonomials=length(splitMonomial);
    
    if(numSplitMonomials==0)
        degList=[];
        varIdxList=[];
    else
        degList=zeros(numSplitMonomials,1);
        varIdxList=zeros(numSplitMonomials,1);
    end
    
    for curSplitMonomial=1:numSplitMonomials
        scanVals=sscanf(splitMonomial{curSplitMonomial},'x%i^%i');
        if(length(scanVals)==1)
            varIdxList(curSplitMonomial)=scanVals;
            degList(curSplitMonomial)=1;
        else
            varIdxList(curSplitMonomial)=scanVals(1);
            degList(curSplitMonomial)=scanVals(2);
        end
    end
    
    maxXIdx=max([varIdxList;maxXIdx]);
    
    %Save the variable indices and degrees. 
    monomials{curMonomial}=[varIdxList,degList];
end

maxXIdx=max(maxXIdx,numVars);

%Knowing the maximum index of x, we can allocate a term matrix.
termMat=zeros(maxXIdx+1,numMonomials);
termMat(1,:)=coeffs';
for curMonomial=1:numMonomials
    monomialVals=monomials{curMonomial};
    numTerms=size(monomialVals,1);
    
    if(numTerms>0)
        varIdxList=monomialVals(:,1);
        degList=monomialVals(:,2);
        
        for curTerm=1:numTerms
            termMat(1+varIdxList(curTerm),curMonomial)=degList(curTerm);
        end
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
