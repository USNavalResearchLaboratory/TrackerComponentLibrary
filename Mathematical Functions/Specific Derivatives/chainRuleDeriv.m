function [derivVal,derivStructure]=chainRuleDeriv(desDeriv,derivF,derivX,n1,derivStructure)
%%CHAINRULEDERIV For a function
%                f(x1(y1,...,yn2),x2(y1,...yn2),...,xn1(y1,...yn2))
%                Evaluate derivatives with respect to y1...yn2. For
%                example, one might want to evaluate d^3f/(dy1^2dy2). Here,
%                we say n1=numX and n2=numY.
%
%INPUTS: desDeriv An n2X1 or 1Xn2 vector containing the integer
%                 derivative numbers that one desires of the function f
%                 with respect to each of the y values. The order of the
%                 derivative is sum(desDeriv) and must be >=0.
%          derivF This can be a multidimensional or a function handle. As a
%                 function handle, this returns values of the partial
%                 derivatives of varying order of the function f with
%                 respect to x1...xn1. The function is called as
%                 derivF(xList), where xList is a 1Xn1 vector such that
%                 xList(i) lists the number of partial derivatives of f
%                 with respect to the ith x variable. For example, for
%                 n1=3, xList=[1;0;2] means that the third-order
%                 derivative d^3f/(dx1dx3^2) is desired. When a matrix is
%                 passed, then it is a matrix such that
%                 derivF(xList(1)+1,xList(2)+1,...,xList(n1)+1) returns the
%                 appropriate derivative as the function would. This can
%                 return multidimensional outputs if one wishes to compute
%                 multiple derivatives in paralllel. Its output
%                 dimensionality must match that of derivX or be scalar.
%          derivX This can be a multidimensional or a function handle. As a
%                 function handle, this returns values of the partial
%                 derivatives of the x values with respect to y variables.
%                 derivX(yListAug) returns the derivative of the x variable
%                 numbered from zero in yListAug(n2+1). The order of the
%                 derivative of each with respect to the ith y value is
%                 given by yListAug(i) and yListAug(n2=1) specified which
%                 of the n1 x values is selected. When a matrix is passed,
%                 then it is a matrix such that
%                 derivX(yListAug(1)+1,yListAug(2)+1,...,yListAug(n2)+1,yListAug(n2+1)+1)
%                 returns the appropriate derivative as the function would.
%                 This can return multidimensional outputs if one wishes to
%                 compute multiple derivatives in paralllel. Its output
%                 dimensionality must match that of derivF or be scalar.
%              n1 This is the integer number of x variables.
%  derivStructure A large part of the execution time of this function is
%                 determining how to form the necessary terms in a sum as
%                 products of values in derivX and derivF as a function of
%                 desDeriv and n1. If one is going to take derivatives of
%                 multiple functions with the same n1 and desDeriv values,
%                 then after a single function call, this structure can be
%                 saved and passed back during the next function call so as
%                 to reduce the execution time.
%
%OUTPUTS: derivVal The value of the desired derivative. If derivF and
%                  derivX returned matrices, then this is a matrix.
%   derivStructure A structure holding precomputed values that depend only
%                  on desDeriv and n1. This structure can be apssed on
%                  subsequent calls with different derivF and derivX but
%                  the same desDeriv and n1 in order to reduce the
%                  computation time.
%
%Each term in the sum for the desired derivative can be expressed as a
%collection of monomial terms. The monomial can be generated recursively.
%To see this, consider n2=2 and n1=3. In such an instance, if we wanted to
%find the derivative with respect to y1 (desDeriv=[1;0]), then from chain
%rule, one gets have
%df/dy1=(df/dx1)(dx1/dy1)+(df/dx2)(dx2/dy1)+(df/dx3)(dx3/dy1)
%If one wants higher-order derivatives, one then differentiates each
%product of terms and then differentiate the result and so on.
%Each time f or a derivative of f is differentiated, it spawns n1 terms.
%Each time x is differentiated, it spawns 1 term. However, with an
%increasing number of derivatives, one will end up with products of x
%terms. For example, one could have a monomial term of the form
%(d^2f/(dx1dx2))(dx1/dy1)(dx2/dy1)
%In such an instance, due to the chain rule, differentiating
%(dx1/dy1)(dx2/dy1) would spawn two terms.
%
%This function keeps track of lists of tuples representing derivatives with
%respect to f and the x terms with which those derivatives are multipled.
%In order to allocate enough memory up front, the maxNumTotalDerivTerms
%function is called to determine how many monomial terms of each order are
%needed. That is, how many are the product of f and one x value, f and two
%x values, etc.
%
%EXAMPLE:
%Consider an example that is related to differentiating the scalar complex
%response of an ideal isotropic antenna element with respect to azimuth and
%elevation.
% %We define the inner functions
% %x1=cos(y1)*sin(y2);
% %x2=sin(y1)*sin(y2);
% %x3=cos(y2);
% %We will designate a vector of the inner functions at a given set of y as
% %xVals. The outer function chosen is
% %f=exp(-1j*2*pi*(z.'*xVals);
% %Chosen parameters:
% z=[4;-10;7];
% y1=pi;
% y2=-1/10;
% %We would like to evaluate some third-order derivatives. For the
% %derivatives of x, we choose to build tables to pass to the
% %chainRuleDeriv function. the tables have to be include all partial
% %derivatives that are needed. We will build 4X4 tables. Such tables could
% %be used to compute all derivatives with sum(desDeriv)<=3 and some values
% %with sum(desDeriv)>3. Note that there is no need to fill entries of the
% %table higher than the maximum desired derivative level.
% n1=3;
% maxEls=3;
% %derivX has n2+1 columns.
% derivX=zeros(maxEls+1,maxEls+1,n1);
% %The last column selects the x variable.
% %Derivatives of cos(y1)*sin(y2), including the zeroth-order value.
% derivX(:,:,1)=[cos(y1)*sin(y2), cos(y1)*cos(y2),-cos(y1)*sin(y2),-cos(y1)*cos(y2);
%               -sin(y1)*sin(y2),-cos(y2)*sin(y1), sin(y1)*sin(y2), cos(y2)*sin(y1);
%               -cos(y1)*sin(y2),-cos(y1)*cos(y2), cos(y1)*sin(y2), cos(y1)*cos(y2);
%                sin(y1)*sin(y2), cos(y2)*sin(y1),-sin(y1)*sin(y2),-cos(y2)*sin(y1)];
% %Derivatives of sin(y1)*sin(y2), including the zeroth-order value.
% derivX(:,:,2)=[sin(y1)*sin(y2), cos(y2)*sin(y1),-sin(y1)*sin(y2),-cos(y2)*sin(y1);
%                cos(y1)*sin(y2), cos(y1)*cos(y2),-cos(y1)*sin(y2),-cos(y1)*cos(y2);
%               -sin(y1)*sin(y2),-cos(y2)*sin(y1), sin(y1)*sin(y2), cos(y2)*sin(y1);
%               -cos(y1)*sin(y2),-cos(y1)*cos(y2), cos(y1)*sin(y2), cos(y1)*cos(y2)];
% %Derivatives of cos(y2);
% derivX(:,:,3)=[cos(y2),-sin(y2), -cos(y2), sin(y2);
%           0,       0,        0,        0;
%           0,       0,        0,        0;
%           0,       0,        0,        0];
% %The values of xVals are the zeroth-order terms.
% xVals=[derivX(1,1,1);derivX(1,1,2);derivX(1,1,3)];
% %For derivF, we choose to pass a function handle instead of building a
% %large matrix.
% derivF=@(xDerivVal)(-1j*2*pi)^sum(xDerivVal)*prod(z(:).^xDerivVal(:))*exp(-1j*2*pi*z'*xVals);
% 
% desDeriv=[0,1];
% derivVal=chainRuleDeriv(desDeriv,derivF,derivX,n1)
% %The solution is about 15.5190 -13.5718i.
% desDeriv=[2,1];
% derivVal=chainRuleDeriv(desDeriv,derivF,derivX,n1)
% %The solution is about -1.1117e+03 - 9.8612e-01i.
% desDeriv=[4,2];
% derivVal=chainRuleDeriv(desDeriv,derivF,derivX,n1)
% %The solution is about -2.3163e+06 + 1.4360e+06i.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5)
    derivStructure=[];
end

if(~isa(derivF,'function_handle'))
    nDimsF=size(derivF);
else
    nDimsF=[];
end

if(~isa(derivX,'function_handle'))
    nDimsX=size(derivX);
else
    nDimsX=[];
end
k=sum(desDeriv);

desDerivOrig=desDeriv;

%Deal with the special case of no derivatives first.
if(k==0)
    fArg=zeros(1,n1);
    if(~isempty(nDimsF))
        idx=nDim2Index(nDimsF,fArg+1);
        derivVal=derivF(idx);
    else
        derivVal=derivF(fArg);
    end
    return
end

%Determine the maximum number of pairs, triples, etc. or derivatives
%needed.
numTermsList=maxNumTotalDerivTerms(k,n1);
n2=length(desDeriv);

%If the derivative structure was not provided, them compute it.
if(isempty(derivStructure))
    %Allocate space for the tuples of varying orders. We allocate twice as much
    %space as needed so that we can store tuples for the next derivative
    %without overwriting the ones for the current derivatives. Then we switch
    %next and current and do it again for one derivative higher until all
    %derivatives have been completed.

    %We separate the tuples for f from those for x, since the f tuples have n1
    %entires for derivatives (for each x) and the x derivatives have n2+1
    %entries, 1 to identify which x and then n2 to identify which derivatives
    %of y are taken.
    curFTuples=cell(k,1);
    nextFTuples=cell(k,1);
    curXTuples=cell(k,1);
    nextXTuples=cell(k,1);

    for i=1:k
        curMaxTerms=numTermsList(i);
        curFTuples{i}=zeros(n1,curMaxTerms);
        nextFTuples{i}=zeros(n1,curMaxTerms);
        %For the ith order, there are i x terms.
        curXTuples{i}=zeros(n2+1,i,curMaxTerms);
        nextXTuples{i}=zeros(n2+1,i,curMaxTerms); 
    end

    %Though we are preallocating the amount of space for all terms, only a
    %subset of the space is used most of the time. Thus, we need to have an
    %array indicating how many terms of each order are present. For each
    %derivative, there can be no terms higher than the derivative number.
    curNumMonomials=zeros(k,1);
    nextNumMonomials=zeros(k,1);

    %Insert the first-order terms by hand.
    %First, find the index with respect to which the derivative will be taken.
    derivIdx=find(desDeriv,1);
    desDeriv(derivIdx)=desDeriv(derivIdx)-1;

    %There are n1 first order monomials
    curNumMonomials(1)=n1;
    %Fill in the f and x terms.
    for curXVar=1:n1
        curFTuples{1}(curXVar,curXVar)=1;
        curXTuples{1}(derivIdx,1,curXVar)=1;
        curXTuples{1}(n2+1,1,curXVar)=curXVar-1;
    end

    %Now, compute all of the higher-order terms.
    for i=2:k%i is the derivative number.
        %We now go through all of the orders that are present and differentiate
        %the terms. Each derivative of f spans n1 terms. Each product of n x
        %terms spawns n terms.

        %First, find the index with respect to which the next derivative will
        %be taken.
        derivIdx=find(desDeriv,1);
        desDeriv(derivIdx)=desDeriv(derivIdx)-1;

        %Zero the count of the number of terms that have been placed into the
        %next batch.
        nextNumMonomials(1:i)=0;
        for order=1:(i-1)
            numMonomials=curNumMonomials(order);

            for curMonomial=1:numMonomials
                %First, differentiate the f term, resulting in n1 new terms
                for curXVar=1:n1
                    nextNumMonomials(order+1)=nextNumMonomials(order+1)+1;
                    curIdx=nextNumMonomials(order+1);

                    nextFTuples{order+1}(:,curIdx)=curFTuples{order}(:,curMonomial);
                    nextFTuples{order+1}(curXVar,curIdx)=nextFTuples{order+1}(curXVar,curIdx)+1;
                    %The f terms are multiplied by all of the existing x
                    %derivatives.
                    nextXTuples{order+1}(:,1:order,curIdx)=curXTuples{order}(:,1:order,curMonomial);
                    %And the f term is also multiplied by the new derivative of
                    %x with respect to the currently selected derivative index.
                    nextXTuples{order+1}(:,order+1,curIdx)=zeros(n2+1,1);
                    nextXTuples{order+1}(n2+1,order+1,curIdx)=curXVar-1;
                    nextXTuples{order+1}(derivIdx,order+1,curIdx)=1;
                end

                %Next, differentiate each of the x terms resulting in the same
                %number of new terms as there are x terms. The number of x
                %terms equals the order and the order of the result does not
                %change.
                for curNewX=1:order
                    nextNumMonomials(order)=nextNumMonomials(order)+1;
                    curIdx=nextNumMonomials(order);

                    %The f term remains unchanged.
                    nextFTuples{order}(:,curIdx)=curFTuples{order}(:,curMonomial);
                    %Increment the order of the x term with respect to y
                    nextXTuples{order}(:,1:order,curIdx)=curXTuples{order}(:,1:order,curMonomial);
                    nextXTuples{order}(derivIdx,curNewX,curIdx)=nextXTuples{order}(derivIdx,curNewX,curIdx)+1;
                end
            end
        end

        %Swap the order of the buffers to indicate that the next set is now
        %the current set. Matlab does not allow one to simply swap
        %pointers, so the commented out lines are just to remind one what
        %has to keep in mind when implementing this in C/Fortran/etc.
        %temp=curFTuples;
        curFTuples=nextFTuples;
        %nextFTuples=temp;

        %temp=curXTuples;
        curXTuples=nextXTuples;
        %nextXTuples=temp;

        %temp=curNumMonomials;
        curNumMonomials=nextNumMonomials;
        %nextNumMonomials=temp;
    end
else
    if(any(derivStructure.desDeriv(:)~=desDeriv(:)))
       error('The derivStructure is inconsistent with desDeriv.')
    end
    if(derivStructure.n1~=n1)
       error('The derivStructure is inconsistent with n1.')
    end

    curNumMonomials=derivStructure.curNumMonomials;
    curFTuples=derivStructure.curFTuples;
    curXTuples=derivStructure.curXTuples;
end
    
%When we get here, desDeriv should be all zeros and all of the indices for
%the monomial values are stored in curFTuples and curXTuples. Thus, the
%next step is to evaluate and sum the products of the terms.

derivVal=0;
for order=1:k
    numMonomials=curNumMonomials(order);
    
    for curMonomial=1:numMonomials
        monomialProdVal=1;
        %First, the f derivative term.
        fArg=curFTuples{order}(:,curMonomial);
        
        if(~isempty(nDimsF))
            idx=nDim2Index(nDimsF,fArg+1);
            termVal=derivF(idx);
        else
            termVal=derivF(fArg);
        end
        
        monomialProdVal=monomialProdVal.*termVal;
        
        %Now for the order number of derivatives of x
        for curX=1:order
            xArg=curXTuples{order}(:,curX,curMonomial);
            if(~isempty(nDimsX))
                idx=nDim2Index(nDimsX,xArg+1);

                termVal=derivX(idx);
            else
                termVal=derivX(xArg);
            end

            monomialProdVal=monomialProdVal.*termVal;
        end
        derivVal=derivVal+monomialProdVal;
    end
end

if(nargout>1)
    derivStructure.desDeriv=desDerivOrig;
    derivStructure.n1=n1;
    derivStructure.curNumMonomials=curNumMonomials;
    derivStructure.curFTuples=curFTuples;
    derivStructure.curXTuples=curXTuples;
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
