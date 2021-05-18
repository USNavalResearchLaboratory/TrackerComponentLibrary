function areRPrime=areRelativelyPrime(a)
%%ARERELATIVELYPRIME Determine whether a set of numbers is relatively
%               prime (are coprime). A collections of numbers are
%               relatively prime if the greatest common divisor between all
%               pairs of the numbers is 1. Relatively prime numbers are
%               also known as coprime numbers.
%
%INPUTS: a An nXnumSets set of numSets sets of n-integers where one wishes
%          to determine whether the n numbers in each set are relatively
%          prime.
%
%OUTPUTS: areRPrime A numSetsX1 boolean array indicating whether the
%                   numbers in each set are coprime.
%
%Coprime numbers arise in a number of areas, such as in the Chinese
%remainder theorem. This function jst goes through all pairs of numbers and
%uses the gcd function to see whether the greatest common divisor between
%all pairs is 1.
%
%EXAMPLE:
% a=zeros(3,2);
% %This set of numbers is coprime even though 4 is not a prime number.
% a(:,1)=[7;9;4];
% %This set of numbers is not coprime
% a(:,2)=[3;2;9];
% areRPrime=areRelativelyPrime(a)
%The functions returns [1;0]
%
%September 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(a,1);
numSets=size(a,2);

areRPrime=true(numSets,1);

for curSet=1:numSets
    for n1=1:(N-1)
        for n2=(n1+1):N
            if(gcd(a(n1,curSet),a(n2,curSet))~=1)
                areRPrime(curSet)=false;
                break;
            end
        end
        
        if(areRPrime(curSet)~=true)
            break;
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
