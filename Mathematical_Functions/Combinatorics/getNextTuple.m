function tuple=getNextTuple(param1,maxVals,firstIsMostSig)
%%GETNEXTTUPLE This function produces the next tuple in a counting series.
%           This can be used as a method of implementing nested loops. Each
%           index of the tuple counts from 0 to the value in maxVals and
%           resets. This is essentially a function for counting using
%           different bases for each digit.
%
%INPUTS: param1, maxVals If the first tuple in the series is desired, then
%                        param1 is the total number of digits (N) in the
%                        tuple. Otherwise, param1 is the previous tuple and
%                        maxVals is a vector whose elements correspond to
%                        the maximum value that each digit of the tuple can
%                        take (>=0). The elements of maxVals correspond to
%                        the base of each digit -1.
%         firstIsMostSig This is a boolean variable indicating whether
%                        tuple(1) is the most significant digit (or whether
%                        tuple(N) is the most significant). This affects
%                        the ordering of the tuples. The default if this
%                        parameter is omitted or an empty matrix is passed
%                        is true.
%
%OUTPUTS: tuple The first NX1 tuple in the series (all zeros) if only
%               param1 is provided as the number of digits in the tuple.
%               Otherwise, the next tuple in the series is returned. If one
%               is past the final tuple in the sequence, then an empty
%               matrix is returned.
%
%This function just implements the rules of integer counting with different
%bases for different digits.
%
%EXAMPLE:
% maxVals=[2;3;1];
% maxNumTuples=prod(maxVals+1);
% tuple(:,1)=getNextTuple(length(maxVals));
% for k=2:maxNumTuples
%     tuple(:,k)=getNextTuple(tuple(:,k-1),maxVals);
% end
%The matrix tuple now contains all of the tuples in the sequence. The first
%digit in each tuple is the most significant.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%If the first tuple is requested.
if(nargin==1)
   tuple=zeros(param1,1);
   return;
end

if(nargin<3||isempty(firstIsMostSig))
    firstIsMostSig=true; 
end

numDim=length(param1);
tuple=param1;%The previous tuple.

if(firstIsMostSig==false)
    curLevel=1;
    isAscending=true;
    while(curLevel<=numDim)
        if(curLevel<1)
            %If we have gotten a complete new tuple.
            return;
        elseif(isAscending)
            %Try incrementing the order at this level.
            tuple(curLevel)=tuple(curLevel)+1;
            if(tuple(curLevel)<=maxVals(curLevel))
                isAscending=false;
                curLevel=curLevel-1;
                continue;
            else
                %If the value is invalid, then just keep ascending.
                curLevel=curLevel+1;
                continue;
            end
        else%We are descending in the sums here.
            tuple(curLevel)=0;
            curLevel=curLevel-1;
        end
    end
else
    curLevel=numDim;
    isAscending=true;
    while(curLevel>=1)
        if(curLevel>numDim)
            %If we have gotten a complete new tuple.
            return;
        elseif(isAscending)
            %Try incrementing the order at this level.
            tuple(curLevel)=tuple(curLevel)+1;
            if(tuple(curLevel)<=maxVals(curLevel))
                isAscending=false;
                curLevel=curLevel+1;
                continue;
            else
                %If the value is invalid, then just keep ascending.
                curLevel=curLevel-1;
                continue;
            end
        else%We are descending in the sums here.
            tuple(curLevel)=0;
            curLevel=curLevel+1;
        end
    end
end

%If we get here, then we have gone past the final tuple.
tuple=[];
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
