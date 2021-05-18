function tuple=getNextUnorderedTuple(tuple,maxVal,firstMostSig)
%%GETNEXTUNORDEREDTUPLE Unordered tuples are tuples where the ordering of
%           the elements in the tuple does not matter. Thus, the tuples
%           1,1,2,3 and 3,1,2,1 represent the same value. This function
%           takes a length-N tuple (starting with [0;0;...0]). and returns
%           the next tuple in the sequence. The sequence is arranged such
%           that the elements of the tuples are either all increasing or
%           all decreasing. If the final tuple in the sequence is passed,
%           then an empty matrix is returned.
%
%INPUTS: tuple A 1NN or NX1 vector holding the current tuple. The elements
%              in the tuple must be between 0 and maxVal and the components
%              must be ordered in the same manner as the output of this
%              function (see comments below).
%       maxVal The maximum value that each digit of the tuple can take
%              (>=0). This corresponds to the base of each digit -1.
%  firstIsMostSig This is a boolean variable indicating whether tuple(1) is
%              the most significant digit (or whether tuple(N) is the most
%              significant). This affects the ordering of the sequence that
%              produces the unordered tuples. The default if omitted or an
%              empty matrix is passed is false.
%
%OUTPUTS: tuple The next tuple in the sequence. This has the same
%               dimensions as tuple on the input unless tuple on the input
%               is the final tuple in the sequence, in which case this is
%               an empty matrix.
%
%If firstMostSig=true, then the tuples are generated such that
%tuple(1)>=tuple(2)...>=tuple(N)
%If firstMostSig=false, then the tuples are generated such that
%tuple(1)<=tuple(2)...<=tuple(N)
%The tuple provided on the input to this function must have the same
%ordering as that is specified by firstMostSig. The first tuple in both
%sequences is the all zero tuple.
%
%November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(firstMostSig))
    firstMostSig=false;
end

tupLength=length(tuple);

if(firstMostSig)
    curIdx=tupLength;
    while(1)
        tuple(curIdx)=tuple(curIdx)+1;
        
        if(tuple(curIdx)>maxVal&&curIdx==1)
            %If we have passed the final tuple.
            tuple=[];
            return;
        elseif(curIdx==1)
            %We have found a valid tuple.
            return; 
        end
        
        if(tuple(curIdx)>tuple(curIdx-1))
            tuple(curIdx)=0;
            curIdx=curIdx-1;
        else
            %We have found a valid tuple.
            return; 
        end
    end
else
    curIdx=1;
    while(1)
        tuple(curIdx)=tuple(curIdx)+1;
        
        if(tuple(curIdx)>maxVal&&curIdx==tupLength)
            %If we have passed the final tuple.
            tuple=[];
            return;
        elseif(curIdx==tupLength)
            %We have found a valid tuple.
            return; 
        end
        
        if(tuple(curIdx)>tuple(curIdx+1))
            tuple(curIdx)=0;
            curIdx=curIdx+1;
        else
            %We have found a valid tuple.
            return; 
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
