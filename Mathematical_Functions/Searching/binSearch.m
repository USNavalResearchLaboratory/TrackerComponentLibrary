function [val,idx]=binSearch(vec,key,choice,idxRange)
%%BINSEARCH Perform a binary search for the value key in the vector vec,
%           which has been sorted such that the elements are in ascending
%           order. The optional choice parameter sets what is returned if
%           key is not in vec.
%
%INPUTS: vec A vector with elements sorted in increasing order.
%            Alternatively, this can be a function handle to a function
%            that takes integers as inputs and is strictly increasing. 
%        key The value that one wishes to find in the vector vec.
%     choice An optional parameter that determines what is returned if key
%            is not found. If this parameter is omitted, then the default
%            is zero (the closest value).
%            0 means return the closest value.
%            1 means return the next lower value if there is one, otherwise
%              return the lowest value in vec.
%            2 means return the next higher value if there is one,
%              otherwise return the highest value in vec.
%   idxRange This is the 2X1 or 1X2 range of integers that are considered
%            as searchable indices in vec (or possible inputs to the
%            function vec). If vec is a function, then this must be
%            specified. Otherwise, the default if omitted or an empty
%            matrix is passed is [1;length(vec)].
%
%OUTPUTS: val Either key, if found, or a value nearby as determined by the
%             parameter closest.
%         idx The index of val in the vector vec.
%
%This is just a basic binary search. The search space is cut in half each
%time. In some cases, such as are elaborated in [1] the Fibonacci search
%can be faster. However, for most problems, there is little difference
%between a binary and Fibonacci search and in [2], it is shown that the
%binary search has a degree of optimality over the more complicated
%Fibonacci search in many instances.
%
%REFERENCES:
%[1] S. Nishihara and H. Nishino, "Binary search revisited: Another
%    advantage of Fibonacci search," IEEE Transactions on Computers, vol.
%    C-36, no. 9, pp. 1132-1135, Sep. 1987.
%[2] K. J. Overholt, "Optimal binary search methods," BIT Numerical
%    Mathematics, vol. 13, no. 1, pp. 84-91, 1973.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The default choice is to return the closest value.
if(nargin<3)
    choice=0;
end

if(nargin<4||isempty(idxRange))
    idxRange=[1;length(vec)];
end

UB=idxRange(2);
LB=idxRange(1);

UBVal=vec(UB);
if(UBVal<=key)
    val=UBVal;
    idx=UB;
    return;
end

LBVal=vec(LB);
if(LBVal>=key)
    val=vec(LB);
    idx=LB;
    return;
end

while(UB~=LB)
    mid=fix((UB+LB)/2);
    
    %If the search has reached the deepest level. 
    if(mid==UB||mid==LB)
        %First, check to see whether either one equals the key.
        if(vec(LB)==key)
            val=vec(LB);
            idx=LB;
            return;
        elseif(vec(UB)==key)
            val=vec(UB);
            idx=UB;
            return;
        end
        
        %If the key is not in vec, then the return value depends on the
        %choice parameter.
        switch(choice)
           case 1
               val=vec(LB);
               idx=LB;%Return the next lowest value.
               return;
           case 2
               val=vec(UB);
               idx=UB;%Return the next highest value.
               return;
           otherwise%Return the closest value.
               diff1=abs(key-vec(LB));
               diff2=abs(key-vec(UB));
               if(diff1<diff2)
                   idx=LB;
               else
                   idx=UB;
               end
               val=vec(idx);
               return;
        end
    end
 
    if(vec(mid)>key)
        UB=mid;
        continue;
    elseif(vec(mid)<key)
        LB=mid;
        continue;
    else
        idx=mid;
        val=vec(idx);
        return;
    end
end

idx=UB;
val=vec(UB);
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
