function [arr,data]=getNextArrangement(data,m)
%%GETNEXTARRANGEMENT Return the next arrangement of n items put into m
%                    spaces. The arrangements are generated in
%                    lexicographic order. There is a total of
%                    binomial(n,m)*factorial(m) arrangements.
%
%INPUTS: data If the first arrangement is desired, then this is n, the
%             number of items to use and the second input must be given.
%             Otherwise, this is the data structure returned from the
%             previous call to this function and the second input must be
%             omitted or have an empty matrix passed.
%           m This input is only provided when obtaining the first
%             arrangement in the sequence. This is the number of slots to
%             use. m<=n or else no arrangements exist.
%
%OUTPUTS: arr The current mX1 arrangement vector. This the last arrangement
%             in the sequence has been passed, then this is an empty
%             matrix. Values can be from 1 to n.
%        data The data that should be passed to a subsequent call of this
%             function (without the m input) to get the next arrangement in
%             the sequence.
%
%Arragements are generating by going through all possible combinations of
%which of the n elements are in the m slots and then permuting the ordering
%of the elements into the slots. The getNextCombo function is used to
%generate the combinations and the getNextMultisetPermutation function is
%used to generate the permutations. As both produce lexicographic ordered
%sequences, the sequence produced by this function should be in
%lexicographic order.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin>1)
    %If we are getting the first arrangement.
    n=data;

    if(m>n||m==0)
        arr=[];
        return;
    end

    curCombo=0:(m-1);
    arr=getNextMultisetPermutation(curCombo+1,true);
    
    data=[];
    data.curCombo=curCombo;
    data.arr=arr;
    data.n=n;
    return;
end

%For subsequent arrangements, get the next permutation. If the final one
%has been reached, then get the next combination.

arr=getNextMultisetPermutation(data.arr);
if(~isempty(arr))
    data.arr=arr;
    return;
end

data.curCombo=getNextCombo(data.curCombo,data.n);

%If we have passed the final arrangement in the sequence.
if(isempty(data.curCombo))
    arr=[];
    return;
end

arr=getNextMultisetPermutation(data.curCombo+1,true);
data.arr=arr;
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
