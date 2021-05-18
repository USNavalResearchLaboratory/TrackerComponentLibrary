function digVec=arbBaseInc(digVec,base)
%%ARBBASEINC Increment a vector whose elements represent digits in an
%            arbitrary base (radix) counting system. 
%
%INPUTS: digVec A vector of length n whose elements range in value from 0
%               to base-1 and which should be incremented by one.
%          base The base (radix) of the number represented by digVec, where
%               each element in the vector is a digit of radix base.
%OUTPUTS: digVec The vector incremented by one. If the incremented vector
%               requires more than n places to represent, then an empty
%               matrix will be returned. If any of the elements of digVec
%               are base or larger, then the output will have no meaning.
%
%This is just a simple function for computing the next n-tuple in a certain
%base. It is essentially a recursion for the method described in Chapter
%7.2.1.1 of [1].
%
%This function can be used in a loop to go through all possible groups of
%a certain number of values allowing for repeats when starting from a 
%vector of zeros.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 2:
%    Generating all Tuples and Permutations, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    n=length(digVec);

    curChoice=1;
    digVec(curChoice)=digVec(curChoice)+1;
    while(digVec(curChoice)>=base)
        %If we have reached the maximum value for this element, set this
        %element and those below it back to one. Then, increment the next
        %element.
        digVec(curChoice:-1:1)=0;
        curChoice=curChoice+1;

        if(curChoice>n)
            digVec=[];
            return
        else
            digVec(curChoice)=digVec(curChoice)+1;
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
