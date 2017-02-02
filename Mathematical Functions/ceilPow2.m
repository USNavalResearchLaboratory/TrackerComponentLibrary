function p2Vals=ceilPow2(nums)
%%CEILPOW2   Find the smallest powers of 2 that are greater than or equal
%            to the non-negative, real values in vals.
%
%INPUTS:  nums  An array or matrix of non-negative, real values.
%
%OUTPUTS:p2Vals Values corresponding to the power of 2 greater than or
%               equal to the values in nums having the same datatype as
%               nums (float, int32, etc.).
%
%This function checks for the data type of the passed value. If the
%value passed returns true from the function isfloat then an algorithm
%that manipulates the mantissa and exponent of a floating point number is
%used. If the data type is not a float, then it is assumed to be an integer
%and binary shift operations are used to find the next highest power of
%two.
%
%The implementation here could be easily translated in C, since the
%floating point operations correspond to functions in math.h and the
%bitshift function for integer values just implements the << and >>
%operators in C.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

assert(all(nums>=0))

if(isfloat(nums))
    %log2 with two return values is the frexp function in C.
    [F,E]=log2(nums);
    
    %pow2 with two arguments is the ldexp function in C.
    p2Vals=pow2(0.5*(F>0),E+(F>0.5));
else
    p2Vals=zeros(size(nums),'like',nums);
    numVals=size(nums(:));
    
    for curVal=1:numVals
        val=ones(1,1,'like',nums);
        temp=nums(curVal);
        
        if(temp==0)
            p2Vals(curVal)=0;
            continue; 
        end
        
        while(temp~=1)
            temp=bitshift(temp,-1); 
            val=bitshift(val,1); 
        end

        if(val==nums(curVal))
            p2Vals(curVal)=val;
        else
            p2Vals(curVal)=bitshift(val,1);
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
