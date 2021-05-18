function [tuples,minCostVal]=assign3DBF(C,maximize)
%%ASSIGN3DBF Solve the operations research axial 3D assignment problem
%         using an inefficient brute-force algorithm. The branch-and-bound
%         solution in assign3DBB is usually faster, except in certain
%         problems (See the example below). Such problems are NP-hard and
%         thus cannot be solved in polynomial time, so this function can
%         only solve problems of a small size. The optimization problem
%         being solved is minimize (or maximize)
%         \sum_{i=1}^{n1}\sum_{j=1}^{n2}\sum_{k=1}^{n3}C_{i,j,k}*\rho_{i,j,k}
%         subject to
%         \sum_{i=1}^{n1}\sum_{j=1}^{n2}\rho_{i,j,k}<=1 for all k
%         \sum_{i=1}^{n1}\sum_{k=1}^{n3}\rho_{i,j,k}<=1 for all j
%         \sum_{j=1}^{n2}\sum_{k=1}^{n3}\rho_{i,j,k} =1 for all i
%         \rho_{i,j,k} = 0 or 1
%         assuming that n1<=n2<=n3, and C is an n1Xn2Xn3 cost matrix.
%         This is equivalent to the optimization problem
%         min (or max) sum_{i} C(i,phi_2(i),phi_3(i))
%         where phi_2 and phi_3 are length n1 arrangements of n2 and n3
%         items over which the minimization is performed.
%
%INPUTS: C An n1Xn2Xn3 cost hypermatrix with n1<=n2<=n3. C cannot contain
%          any NaNs and the largest finite element minus the smallest
%          element is a finite quantity (does not overflow) when performing
%          minimization and where the smallest finite element minus the
%          largest element is finite when performing maximization.
%          Forbidden assignments can be given costs of +Inf for
%          minimization and -Inf for maximization. During minimization,
%          there should be no elements with -Inf cost and no elements with
%          +Inf cost during maximization.
% maximize If true, the minimization problem is transformed into a
%          maximization problem. The default if this parameter is omitted
%          or an empty matrix is passed is false.
%
%OUTPUTS: tuples A 3Xn1 matrix where tuples(:,i) is the ith
%                assigned tuple in C as described above.
%     minCostVal The cost value of the optimal assignment found. This is
%                the sum of the values of the assigned elements in C.
%
%The brute-force enumeration algorithm goes through all arrangements of
%items in the second index to those in the first index. Conditioned on
%those assignments, the problem is reduced to a 2D assignment problem. The
%function assign2D is used to solve that problem in O(n^3) time. Thus, if a
%matrix of size nXnXn is passed, the algorithm takes O((n!)(n^3)) time as
%it has to run the 2D assignment algorithm n! times, one for each
%arrangement of the items in the second dimension with respect to the first
%(permutation in this instance).
%
%EXAMPLE:
%The algorithm is usually slower then branch-and bound, but an example of
%an exception is if n1 and n2 are small and n3 is large. For example:
% n1=5;
% n2=5;
% n3=100;
% C=randn(n1,n2,n3);
% maximize=false;
% tic
% [tuples,minCostVal]=assign3DBF(C,maximize)
% toc
% tic
% [tuples,minCostVal]=assign3DBB(C,maximize)
% toc
%The brute-force method will usually be faster than the branch-and-bound
%solution in this isntance due to the use of the 2D assignemnt algorithm.
%
%June 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(maximize))
   maximize=false; 
end

if(isempty(C))%The empty matrix special case.
    tuples=[];
    minCostVal=0;
    return
end

if(maximize)
    C=-C;
end

nVals=size(C);

%Deal with the special case of C being scalar.
if(isscalar(C))
    tuples=[1;1;1];
    minCostVal=C(1);
    return;
end

if(length(nVals)<3||~(nVals(1)<=nVals(2)&&nVals(2)<=nVals(3)))
    error('It is required that size(C,1)<=size(C,2)<=size(C,3)')
end

n1=nVals(1);
n2=nVals(2);
n3=nVals(3);

%Allocate space
CInner=zeros(n1,n3);
minCostVal=Inf;
col4Row=[];

[arrange1,data]=getNextArrangement(n2,n1);
while(~isempty(arrange1))
    for i=1:n1
        CInner(i,:)=reshape(C(i,arrange1(i),:),[1,n3]);
    end

    %col4row is n1X1.
    [col4Row, ~, costVal]=assign2D(CInner,false);
    
    %If it is feasible.
    if(~isempty(col4Row)&&costVal<minCostVal)
        minCostVal=costVal;
        phi2=arrange1;
        phi3=col4Row;
    end

    [arrange1,data]=getNextArrangement(data);
end

%If the problem was not feasible.
if(isempty(col4Row))
    tuples=[];
    return; 
end

%Adjust the gains for the case where the initial cost matrix is transformed
%so that maximization can be performed.
if(maximize)
    minCostVal=-minCostVal;
end

%Take the phi values and make the set of tuples that have been accepted.
tuples=[1:n1;phi2.';phi3.'];

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
