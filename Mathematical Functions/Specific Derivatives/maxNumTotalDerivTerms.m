function numTermsList=maxNumTotalDerivTerms(k,n1)
%%MAXTOTALNUMDERIVTERMS For a function
%                f(x1(y1,...,yn2),x2(y1,...yn2),...,xn1(y1,...yn2))
%                consider that someone wishes to evaluate a kth-order
%                partial derivative with respect to some combination of the
%                y terms. The total derivative can be evaluated by
%                performing a series of k sequential derivatives resulting
%                in a sum of products of partial derivatives of f and x.
%                This function returns the number of monomial partial
%                derivative terms of varying orders if one does not combine
%                like terms.
%                
%INPUTS:  k The order (k>=0) of the derivatives.
%        n1 The number of variables x that the function f depends on and
%           that themselves are functions of the y terms with which
%           derivatives are taken.
%
%OUTPUTS: numTermsList This is a kX1 vector such that numTermsList(i) is
%             the number of terms in the sum that are the products of i+1
%             derivatives.
%
%Consider n2=2 and n1=3. In such an instance, if we wanted to find
%the derivative with respect to y1, then from chain rule, we have
%df/dy1=(df/dx1)(dx1/dy1)+(df/dx2)(dx2/dy1)+(df/dx3)(dx3/dy1)
%In this instance, there are three terms. that consist of the product of
%two derivative terms. Thus,maxNumTotalDerivTerms(1,3)=3. Consider taking a
%second derivative. It could be with respect to y1 again or with respect to
%y2. the derivative must be applied to each of the n1 terms in df/dy1. For
%each term one applies the product rule. For the first term, the product
%rule tells us that
%(d/dy2){(df/dx1)(dx1/dy1)}=(d^2f/(dx1dy2))(dx1/dy1)+(df/dx1)(d^2x1/(dy1dy2))
%The term d^2f/(dx1dy2) will expand to n1 terms using the law of total
%probability. Each time f or a derivative of f is differentiated, it spawns
%n1 terms. Each time x is differentiated, it spawns 1 term. However, with
%an increasing number of derivatives, one will end up with products
%including one f term and an increasing number of x terms. For example, one
%could have a monomial term of the form
%(d^2f/(dx1dx2))(dx1/dy1)(dx2/dy1)
%This is the product of three derivatives. In such an instance, due to the
%chain rule, differentiating this would result in n1 terms that are the
%product of four derivatives, due to the differentiation of the f term plus
%2 terms that come from using the product rule and differentiating the x
%terms.
%
%In general, one can form a graph to help count how many individual
%monomial terms consisting of a derivative term of f times a number of
%derivatives of x:
%f-(n1)->pairs-(n1)->triplets-(n1)->quadruplets-(n1)->quintuplets->etc.
%        ^   |       ^     |        ^     |           ^    |
%        |   |       |     |        |     |           |    |
%        -(1)-        -(2)-          -(3)-             -(4)-
%
%Basically, start at f with the number 1 and follow all paths going k steps
%through the graph. In any path, multiply the values of all of the numbers
%in parantheses encountered together. For example, For k=4, one gets n^4 as
%a number of quintuplet terms to add by going straight from f to the right.
%However, one can also loop back. So, one path to quadruplets for k=4 is
%f-(n1)->pairs-(n1)>triplets-(n1)>quadruplets-(3)>quadruplets
%which results in 3*n1^3 terms to add. However, another path for k=4 is
%f-(n1)->pairs-(1)->pairs->(n1)->triplets-(n1)->quadruplets
%which results in n1^3 terms being added.
%
%This function traces all paths of depth k in such graphs and adds up how
%many values were found for pairs, triplets, etc., segregating the values
%in numtermsList. Thus, numTermsList(1) lists pairs, numTermsLists(2) lists
%triplets, etc. In each monomial, there is just one  derivative with
%respect to f. The rest of the terms are derivatives with respect to x
%values. Note that this counts something like (d^2f/(dx1^2))(dx1/dy1)^2 as
%three terms. Also, this is essentially an upper bound on the number of
%terms, because this method of working out the chain and product rules does
%not combine like terms. that result from going down different paths.
%
%EXAMPLE 1:
% numTermsList=maxNumTotalDerivTerms(3,2)
%One gets numTermsList=[2;12;8] for a total of 22 terms. If one were
%differentiating with respect to three different variables, then there
%would be no repeated terms. When differentiating with respect to repeated
%variables, then some terms will be repeated and can be combined, which
%means that the values in numTermsList are an upper bound.
%
%EXAMPLE 2:
% numTermsList=maxNumTotalDerivTerms(3,3)
%One will get numTermsList=[3;27;27] for a total of 57 terms.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numTermsList=zeros(k,1);

levelNumber=1;
goingDown=true;
orderNumber=zeros(k,1);
termCount=zeros(k,1);
takenRightBranch=false(k,1);

%Initialize for the first level.
orderNumber(1)=1;%We have pairs of terms.
termCount(1)=n1;

while(levelNumber>0)
    if(goingDown)
        if(levelNumber==k)
            %If we hit a leaf node, add in the accumulated number of terms.
            idx=orderNumber(levelNumber);
            numTermsList(idx)=numTermsList(idx)+termCount(levelNumber);

            goingDown=false;
            levelNumber=levelNumber-1;
            continue;
        else
            %If we are going down a level, then we go down another level,
            %taking the left (n1) branch first.
            takenRightBranch(levelNumber)=false;
            termCount(levelNumber+1)=n1*termCount(levelNumber);
            orderNumber(levelNumber+1)=orderNumber(levelNumber)+1;
            levelNumber=levelNumber+1;
        end
    else%We are going up. This means, that we must take the right branch to
        %go back down, unless we have taken both branches already.
        if(takenRightBranch(levelNumber))
            %If the right branch has been taken, then just go up another
            %level.
            levelNumber=levelNumber-1;
        else%Now, we will take the right branch.
            takenRightBranch(levelNumber)=true;
            goingDown=true;
            
            termCount(levelNumber+1)=orderNumber(levelNumber)*termCount(levelNumber);
            
            orderNumber(levelNumber+1)=orderNumber(levelNumber);
            levelNumber=levelNumber+1;
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
