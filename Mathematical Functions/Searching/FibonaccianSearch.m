function [val,idx]=FibonaccianSearch(vec,key,choice)
%%FIBONACCIANSEARCH Perform a Fibonaccian search for the value key in the
%           vector vec, which has been sorted such that the elements in it
%           are in ascending order.  The optional choice parameter sets
%           what is returned if key is not in vec.
%
%INPUTS: vec A 1XN or NX1 vector with elements sorted in increasing order.
%        key The value that one wishes to find in the vector vec.
%     choice An optional parameter that determines what is returned if key
%            is not found. If this parameter is omitted, then the default
%            is zero (the closest value).
%            0 means return the closest value.
%            1 means return the next lower value if there is one, otherwise
%              return the lowest value in vec.
%            2 means return the next higher value if there is one,
%              otherwise return the highest value in vec.
%
%OUTPUTS: val Either key, if found, or a value nearby as determined by the
%             parameter closest.
%         idx The index of val in the vector vec.
%
%This implements Algorithm F from Chapter 6.2.1 of [1]. The basic algorithm
%requires that the length of vec be a Fibonacci number -1. However, this
%has been modified as per the comments related to problem 14 in Chaopter
%6.2.1 of [1] to make the algorithm work for general-length vectors vec. It
%is necessary to include special cases for length 1 and 2 vectors.
%
%REFERENCES:
%[1] D. Knuth, The Art of Computer Programming: Sorting and Searching, 2nd
%    ed. Reading, MA: Addison-Wesley, 1998, vol. 3.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(choice))
   choice=0; 
end

N=length(vec);

%Special cases
if(N==1)
    val=vec;
    idx=1;
    return;
elseif(N==2)
    if(vec(1)==key)
        val=vec(1);
        idx=1;
        return;
    elseif(vec(2)==key)
        val=vec(2);
        idx=2;
        return;
    end
    
    switch(choice)
        case 1%Return the next lowest value.
            if(vec(2)<key)
                val=vec(2);
                idx=2;
            else
                val=vec(1);
                idx=1;
            end
            return;
        case 2%Return the next highest value.
            if(vec(1)>key)
                val=vec(1);
                idx=1;
            else
                val=vec(2);
                idx=2;
            end
            return;
        otherwise%Return the nearest value.
            if(abs(key-vec(2))<abs(key-vec(1)))
                val=vec(2);
                idx=2;
            else
                val=vec(1);
                idx=1;
            end
            return;
    end
end

k=FibonacciNumInv(N+1,2)-1;
%Offset from N being a Fibonacci number -1.
M=FibonacciNum(k+1)-1-N;

%Initialize, step F1.
i=FibonacciNum(k);
p=FibonacciNum(k-1);
q=FibonacciNum(k-2);

curStep=0;
while(1)
    curStep=curStep+1;
    
    %Compare, Step F2.
    if(key>vec(i))
        %Step F4
        
        %This is the correction suggested in the solution to Problem 14
        %that deals with vec having a length that is not a Fibonacci
        %number-1.
        if(curStep==1)
           i=i-M; 
        end

        if(p==1)
            %Algorithm ends without finding the value. i holds a value that
            %is too small.
            switch(choice)
                case 1%Return the next lowest value.
                    idx=i;
                    break;
                case 2%Return the next highest value.
                    idx=min(i+1,N);
                    break;
                otherwise%Return the nearest value.
                    if(i==N)
                        idx=N;
                        break;
                    end
                    iErr=abs(vec(i)-key);
                    ipErr=abs(vec(i+1)-key);
                    if(iErr<ipErr)
                        idx=i;
                    else
                        idx=i+1;
                    end
                    break;
            end
        end
        i=i+q;
        p=p-q;
        q=q-p;
    elseif(key<vec(i))
        %Step F3
        if(q==0)
            %Algorithm ends without finding the value. vec(i) holds a value
            %that is too large.
            switch(choice)
                case 1%Return the next lowest value.
                    idx=max(i-1,1);
                    break;
                case 2%Return the next highest value.
                    idx=i;
                    break;
                otherwise%Return the nearest value.
                    if(i==1)
                        idx=i;
                        break;
                    end
                    iErr=abs(vec(i)-key);
                    imErr=abs(vec(i-1)-key);
                    if(iErr<imErr)
                        idx=i;
                    else
                        idx=i-1;
                    end
                    break;
            end
        end
        i=i-q;
        
        diff=p-q;
        p=q;
        q=diff;
    else%(key==vec(i))
        idx=i;
        break;
    end
end

val=vec(idx);

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
