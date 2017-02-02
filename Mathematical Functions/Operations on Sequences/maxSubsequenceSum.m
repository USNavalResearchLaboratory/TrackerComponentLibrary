function [maxSum,startIdx,endIdx]=maxSubsequenceSum(a)
%%MAXSUBSEQUENCESUM  Given a vector of numbers that might be positive or
%                    find the maximum value that can be obtained by summing
%                    a subsequence of the set a.
%
%INPUTS: a  A 1XN or NX1 vector of real numbers.
%
%OUTPUTS: maxSum The maximum value of sums possible by summing consecutive
%                values in a. If a only contains negative values, then this
%                will be zero as the minimum value is obtained by summing
%                nothing in a.
%startIdx,endIdx The beginning and ending indices of the maximum
%                subsequence of a. If a is all negative, these will both be
%                zero.
%
%The solution to the maximum subsequence sum is described in Chapter 2.4
%of [1]. Multiple approaches to the problem are described in [2]. The
%algorithm has O(N) complexity as noted in its presentation in [3].
%
%The algorithm has been slightly modified to handle arrays that contain all
%negative or zero elements.
%
%REFERENCES:
%[1] M.A.Weiss, Data Structures and Algorithm Analysis in C++, 2nd ed.
%    Reading, MA: Addison-Wesley, 1999.
%[2] J. Bentley, "Algorithm design techniques," Communications of the ACM,
%    vol. 27, no. 9, pp. 865-871, Sep. 1984.
%[3] D. Gries, "A note on the standard strategy for developing loop
%    invariants and loops," Department of Computer Science, Cornell
%    University, Ithaca, NY, Tech. Rep. TR 82-531, Oct. 1982.
%
%July 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numElA=length(a);

maxSum=0;
thisSum=0;
thisStartIdx=1;

startIdx=0;
endIdx=0;

maxVal=-Inf;
maxIdx=0;

for j=1:numElA
    %Finding the maximum value is needed so that something can be returned if
    %an array of only zero or negative numbers is passed.
    if(a(j)>maxVal)
        maxVal=a(j);
        maxIdx=j;
    end
    
    thisSum=thisSum+a(j);
    
    if(thisSum>maxSum)
        maxSum=thisSum;
        startIdx=thisStartIdx;
        endIdx=j;
    elseif(thisSum<0)
        thisSum=0;
        thisStartIdx=j+1;
    end
end

if(maxVal<=0)
    %If all of the entries are negative or zero, then return the index of
    %the maximum element.
    maxSum=maxVal;
    startIdx=maxIdx;
    endIdx=maxIdx;
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
