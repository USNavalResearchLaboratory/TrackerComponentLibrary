function [termMats,alreadyHomogenized]=homogenizeMultiDimPolys(termMats,makeProjTransform)
%%HOMOGENIZEMULTIDIMPOLYS Given a system of multivariate polynomials,
%              produce a homogenized system. The degree of a monomial term
%              c*x1^a1*x2^a2*...xn^an for a1...an>=0 is the sum of all
%              a1...an. Most multivariate polynomials will contain terms of
%              different powers. This function adds an extra variable so
%              that all of the terms have the same power. That is for each
%              polynomial f(x1,x2,...,xn), the homogenized polynomial is
%              x0^d*f(x1/x0,x2/x0,...,xn/x0), where d is the highest degree
%              of any monomial term in the polynomial.
%
%INPUTS: termMats Either a single term matrix of a cell array containing
%           multiple term matrices. Each term matrix is An (n+1)XnumTerms
%           matrix such that termMat(:,i)=[c,a1,a2,...,an] is a monomial
%           term where c is the value of of the monomial coefficient and
%           the monomial is x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1). The
%           number of variables and the maximum degree need not be the same
%           for each term matrix. Differing numbers of variables implies
%           that the missing coefficients all have zero exponents.
% makeProjTransform An optional  boolean variable indicating whether the
%           final output system should be a projective transformation of
%           the input system. In such an instance, a cell array must be
%           returned as termMats is expanded by one equation, which is a
%           linear function consisting of random complex numbers times each
%           of the variables. As noted in Section 8 of [1], for a system of
%           n simultaneous polynomial equations in n unknowns, the
%           homogenized system/ projective transform system has a number of
%           solutions exactly equal to the B�zout bound and has no
%           solutions at infinity. If omitted or an empty matrix is passed,
%           the default is false.
%
%OUTPUTS: termMats The homogenized equations. If a sinle term matrix is
%            passed, then this is a matrix. otherwise, this is a cell array
%            of term matrices like the input. All of the term matrices now
%            contain the same number of variables (with zeros for unused
%            variables) with one more variable than the equation with the
%            most variables on the input. This new variable is the
%            homogenizing variable.
%  alreadyHomogenized A boolean value indicating that the input system was
%            already homogenized. In such an instance, all exponents of
%            the added homogenization variable in the retunred system are
%            zero.
%
%Homogenization of multivariate polynomial systems is often used in
%multivariate polynomial solving algorithms to better handle "solutions at
%infinite", which are generally not solutions that one cares about. For
%example, an application of homogenizing multivariate polynomials is in
%Section 8 of [1].
%
%EXAMPLE:
%Consider four equations in four unknowns.
% termMats=cell(4,1);
% termMats{1}=string2Terms('x1*x3-2*x2*x3-6*x1-x2-4*x4+12');
% termMats{2}=string2Terms('x1*x3-3*x2*x3-7*x1+12');
% termMats{3}=string2Terms('x1*x4^3-2*x1*x2-6*x4-x2+12');
% termMats{4}=string2Terms('x1^2+x2^2+x3^2+x4^2-6');
% [termMats1,alreadyHomogenized]=homogenizeMultiDimPolys(termMats)
% %The new system will correspond to
% %'x1*x3-2*x2*x3-6*x1*x5-x2*x5-4*x4*x5+12*x5^2'
% %'x2*x3-3*x3*x4-7*x1*x5+12*x5^2'
% %'x1*x4^3-2*x1*x2*x5^2-6*x4*x5^3-1*x2*x5^3+12*x5^4'
% %'x1^2+x2^2+x3^2+x4^2-6*x5^2'
% %The projective transformed solution is
% [termMats2,alreadyHomogenized]=homogenizeMultiDimPolys(termMats,true)
% %It is the same as the above solution except an extra equation consisting
% %of 
% % 'c1*x1+c2*x2+c3*x3+c4*x4+c5*x5'
% %has been added, where c1...c5 are random complex numbers (with both
% %components between 0 and 1).
%
%REFERENCES:
%[1] L. T. Watson, S. C. Billups, and A. P. Morgan, "Algorithm 652 HOMPACK:
%    A suit of codes for globally convergent homotopy algorithms," ACM
%    Transactions on Mathematical Software, vol. 13, no. 3, pp. 281-310,
%    Sep. 1987.
%
%March 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(makeProjTransform))
   makeProjTransform=false; 
end

isCell=true;
%If only a single matrix is passed.
if(~isa(termMats,'cell'))
    isCell=false;
    termMats={termMats};%Put it in a cell array.
end

numPoly=length(termMats);

%Determine the maximum number of variables.
n=0;
for curPoly=1:numPoly
   n=max([n,size(termMats{curPoly},1)-1]);
end

alreadyHomogenized=false;
for curPoly=1:numPoly
    termMat=termMats{curPoly}; 
    
    %The degrees of the terms.
    dTerms=sum(termMat(2:end,:),1);
    
    %The total degree of the polynomial.
    d=max(dTerms);
    
    numTerms=size(termMat,2);
    
    nCur=size(termMats,1)-1;
    
    %Augment to make have the same number of variables as the equation with
    %the most variables, plus one extra for the homogenization.
    termMat=[termMat;zeros(n-nCur,numTerms)];
        
    for i=1:numTerms
        %The degree of the term.
        deg=dTerms(i);
        
        %The homogenization makes all terms of degree d
        diff=d-deg;
        termMat(end,i)=diff;
        alreadyHomogenized=alreadyHomogenized||(d-deg~=0);
    end
   
    termMats{curPoly}=termMat;
end
alreadyHomogenized=~alreadyHomogenized;

%If only a single equation was passed as a matrix, not a cell array, then
%return a matrix.
if(isCell==false&&makeProjTransform==false)
    termMats=termMats{1};    
    return;
end

if(makeProjTransform)
    termMat=zeros(n+2,n+1);
    
    termMat(1,:)=rand(1,n+1)+1j*rand(1,n+1);
    termMat(2:end,:)=eye(n+1,n+1);
    
    termMats{end+1}=termMat;
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
