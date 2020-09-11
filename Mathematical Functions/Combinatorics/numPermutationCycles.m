function [NCycle,signVal,permut]=numPermutationCycles(permut,option)
%%NUMPERMUTATIONCYCLES Determine the number of cycles in a given permutation
%                   and the sign of the permutation. The sign of the
%                   permutation is (-1)^(n-NCycle), where n is the
%                   number of items being permutation (==length(permut)).
%                   The sign of a permutation is also equal to
%                   (-1)^nInversions. An inversion is a pair of values
%                   whose order is reverse. Also, if desired, return either
%                   the original permutation, the inverse permutation, or
%                   the tagged permutation. An inverse permutation permInv
%                   is such that permut(permInv)=[1;2;3;...;n]. Tagged
%                   permutations can be used for computing inverse
%                   permutations as well as in algorithms for permuting
%                   matrix elements in place.
%
%INPUTS: permut A permutation of n elements. That is, an NX1 or 1XN vector
%              of numbers from 1 to n.
%       option An optional input specifying what the return value permut
%              should be. Possible values are:
%              0) (The default if omitted) permut on putput=permut on the
%                 input.
%              1) permut is the tagged permutation.
%             -1) permut is the inverse permutation.
%
%OUTPUTS: NCycle The number of cycles in the permutation permut.
%        signVal The sign of the permutation permut. This is also known as
%                the signature or parity of the permutation. A value of +1
%                indicates an even parity and -1 indicates an odd parity.
%         permut The value given by the input option.
%
%The algorithm is based on CYCLES from Chapter 16 of [1].
%
%REFERENCES:
%[1] A. Nijenhuis and H. S. Wilf, Combinatorial Algorithms for Computers
%    and Calculators, 2nd ed. New York: Academic press, 1978.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    option=0;
end

n=length(permut);

is=1;
NCycle=n;

for i=1:n
    i1=permut(i);
    while(i1>i)%Line 6
        NCycle=NCycle-1;
        i2=permut(i1);
        permut(i1)=-i2;
        i1=i2;
    end
    
    if(option~=0)
        is=-sign(permut(i));
    end
    
    permut(i)=is*abs(permut(i));%Line 5
end

signVal=1-2*mod(n-NCycle,2);
if(option>=0)
    return;
end

%We are only here if option=-1;
for i=1:n
    i1=-permut(i);
    if(i1<0)
        continue;
    end
    i0=i;
    while(1)
        i2=permut(i1);
        permut(i1)=i0;
        if(i2<0)
            break;
        end
        
        i0=i1;
        i1=i2;
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
