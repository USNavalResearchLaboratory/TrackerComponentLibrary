function [thePartition,theData]=getNextPartition(n,theData)
%%GETNEXTPARTITION Get the next partition of the integer n. A partition is
%                  a set of other integers that sum to n. This goes
%                  through all possible partitions in reverse
%                  lexicographic order. If this function is called only
%                  with n, then the first partition is returned. The
%                  values r, m, and d are needed to get subsequent
%                  partitions. If the last partition is passed, then an
%                  empty matrix is returned. The total number of
%                  partitions for a given n can be found using the
%                  numberOfPartitions function. Unlike the function
%                  getNextMPartition, this goes through partitions of n
%                  into all numbers of parts, whereas getNextMPartition
%                  only goes through partitions of m parts.
%
%INPUTS: n A positive integer that one wishes to partition. If this is the
%          only parameter passed (or theData is an empty matrix), then the 
%          first partition in reverse lexicographic order is returned.
%  theData After the first  partition has been returned, one must pass
%          theData retruend by the function back here to get subsequent
%          values. This is a structure with members:
%          r A parameter needed to get subsequent partitions that is
%            returned by this function upon generating a partition. The
%            first d entries are all of the digits in the partition,
%            without repeats. The other elements are needed to compute
%            future partitions.
%          m A vector specifying the number of times that the digits in r
%            are repeated. T
%          d The number of elements in r and m that contribute to the
%            current partition.
%
%OUTPUTS: thePartition The next partition, which is also encoded in the
%                      returned values of r, m, and d in data. The elements
%                      are in DECREASING order. This is a series of
%                      integers that sum to n. If n is the only parameter
%                      passed, then this is the first partition. If the
%                      parameters for the final partition were passed, then
%                      this is an empty matrix.
%              theData The updated values of r, m, and d in a structure
%                      that should be passed to get the next partition.
%                      Note that the entire partition is encoded in r, m,
%                      and d.
%
%The algorithm is based on NEXPAR in Chapter 9 of [1].
%
%The function is called as
%[thePartition,theData]=getNextPartition(n);
%to get the first parition and as
%[thePartition,theData]=getNextPartition(n,theData);
%to get subsequent partitions.
%
%REFERENCES:
%[1] A. Nijenhuis and H. S. Wilf, Combinatorial Algorithms for Computers
%    and Calculators, 2nd ed. New York: Academic press, 1978.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%If this is the first partition.
if(nargin==1||isempty(theData))
    %Step A
    theData.r(1)=n;
    theData.m(1)=1;
    theData.d=1;
    thePartition=n;
    return;
end

r=theData.r;
m=theData.m;
d=theData.d;

%If the last partition was passed.
if(m(d)==n)
    thePartition=[];
    return;
end

%Step B
if(r(d)==1)
    sigma=m(d)+1;
    d=d-1;
else
    sigma=1;
end

%Step C
f=r(d)-1;
if(m(d)~=1)
   m(d)=m(d)-1;
   d=d+1;
end

%Step D
r(d)=f;
m(d)=fix(sigma/f)+1;

%Step E
s=mod(sigma,f);
if(s~=0)
   d=d+1;
   r(d)=s;
   m(d)=1;
end

%Build the partition including the repeats
partSize=sum(m(1:d));
thePartition=zeros(partSize,1);

partIdx=1;
for idx=1:d
    numRep=m(idx);
    thePartition(partIdx:(partIdx+numRep-1))=r(idx);
    partIdx=partIdx+numRep;
end

theData.r=r;
theData.m=m;
theData.d=d;

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
