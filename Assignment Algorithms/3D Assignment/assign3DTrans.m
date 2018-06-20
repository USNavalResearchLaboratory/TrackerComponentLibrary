function [tuples,minCostVal,costGap,relGap]=assign3DTrans(C,n2Limits,n3Limits,maximize,maxIter,epsVal)
%%ASSIGN3DTRANS Approximate the solution to the axial 3D version of the
%         transportation problem using a Lagrangian relaxation technique.
%         The problem is NP-hard. The optimization problem being solved is:
%         minimize (or maximize)
%         \sum_{i=1}^{n1}\sum_{j=1}^{n2}\sum_{k=1}^{n3}C_{i,j,k}*\rho_{i,j,k}
%         subject to
%         \sum_{i=1}^{n1}\sum_{j=1}^{n2}\rho_{i,j,k}<=n3Limits(k) for all k
%         \sum_{i=1}^{n1}\sum_{k=1}^{n3}\rho_{i,j,k}<=n2Limits(j) for all j
%         \sum_{j=1}^{n2}\sum_{k=1}^{n3}\rho_{i,j,k} =1 for all i
%         \rho_{i,j,k} = 0 or 1
%         assuming that C is an n1Xn2Xn3 cost matrix with n1<=n2<=n3. The
%         axial 3D transportation problem differs from the data fusion
%         axial 3D assignment problen in that in the data fusion axial 3D
%         assignment problem n2Limits and n3Limits are all 1s.
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
% n2Limits An n2X1 vector that provides the limits on the sums over the
%          first and third indices for each item in the second index. All
%          elements of n2Limits must be integers >=1 and <=n2.
% n3Limits An n3X1 vector that provides the limits on the sums over the
%          first and second indices for each item in the third index. All
%         elements of n3Limits must be integers >=1 and <=n3.
% maximize If true, the minimization problem is transformed into a
%          maximization problem. The default if this parameter is omitted
%          or an empty matrix is passed is false.
%  maxIter This parameter is only used if algorithm=0 or 1. It is the
%          maximum number of iterations to perform. The default value if
%          this parameter is omitted is 50.
%   epsVal This parameter is the threshold to use for determining
%          convergence of the algorithm based on the relative duality gap.
%          The relative duality gap is (f-q)/abs(f), where f is the
%          smallest value of the primal function found and q is the largest
%          value of the dual function found. If omitted, the default value
%          is 0.05.
%
%OUTPUTS: tuples A 3Xn1 matrix where tuples(:,i) is the ith assigned tuple
%                in C as described above.
%     minCostVal The cost value of the optimal assignment found. This is
%                the sum of the values of the assigned elements in C.
%        costGap An upper bound for the difference between minCostVal and
%                the true global minimum (the duality gap). Note that
%                finite precision errors can push this to be negative when
%                at or very close to the optimum value.
%         relGap The relative duality gap that is used as the convergence
%                criterion with epsVal. Whereas costGap is the value of the
%                dual cost, relGap is (minCostVal-dualCost)/dualCost .
%                Note that finite precision errors can push this to be
%                negative when at or very close to the optimum value.
%
%This function implements the algorithm of [1]. The algorithm is a
%dual-primal Lagrangian relaxation algorithm.
%
%EXAMPLE:
%Here, we just solve a random problem:
% n1=4;
% n2=5;
% n3=6;
% C=randn(n1,n2,n3);
% n2Limits=2*ones(n2,1);
% n3Limits=ones(n3,1);
% maximize=false;
% [tuples,minCostVal,costGap,relGap]=assign3DTrans(C,n2Limits,n3Limits,maximize)
%
%REFERENCES:
%[1] A. M. Frieze and J. Yadegar, "An algorithm for solving 3-dimensional
%    assignment problems with application to scheduling a teaching
%    practice," The Journal of the Operational Research Society, vol. 32,
%    no. 11, pp. 989-995, Nov. 1981.
%
%June 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Equation numbers refer to those in [1].

if(nargin<4||isempty(maximize))
   maximize=false; 
end

if(nargin<5||isempty(maxIter))
    maxIter=50;
end

if(nargin<6||isempty(epsVal))
    epsVal=0;%0.05;
end

if(isempty(C))%The empty matrix special case.
    tuples=[];
    minCostVal=0;
    costGap=0;
    relGap=0;
    return
end

%This algorithm tries to maximize cost, so we will just flip the
%sign of the costs to use it to minimize cost.
if(maximize==false)
    C=-C;
end

nVals=size(C);

%Deal with the special case of C being scalar.
if(isscalar(C))
    %If the problem is feasible.
    if(isscalar(n2Limits)&&n2Limits==1&&isscalar(n3Limits)&&n3Limits==1)
        tuples=[1;1;1];
        minCostVal=C(1);
        costGap=0;
        relGap=0;
    else
        tuples=[];
        minCostVal=[];
        costGap=[];
        relGap=[];
    end
    return;
end

if(length(nVals)<3||~(nVals(1)<=nVals(2)&&nVals(2)<=nVals(3)))
    error('It is required that size(C,1)<=size(C,2)<=size(C,3)')
end

if(length(n2Limits)~=nVals(2))
    error('n2Limits has the wrong dimensonality.')
end

if(length(n3Limits)~=nVals(3))
    error('n3Limits has the wrong dimensonality.')
end

if(any(n2Limits<=0)||any(n2Limits>nVals(2))||any(n3Limits<=0)||any(n2Limits>nVals(3)))
    error('Elements of n2Limits and n3Limits must be greater than 0 and less than or equal to size(C,2) and size(C,3).')
end

n1=nVals(1);
n2=nVals(2);
n3=nVals(3);

u=zeros(n3,1);%Allocate and initialize the dual variables.

%Upper and lower bounds for the cost function encountered thus far.
lb=-Inf;
ub=Inf;
%Allocate space
dim2ForDim1Cur=zeros(n1,1);
dim3ForDim1Cur=zeros(n1,1);
for curIter=1:maxIter
    %Compute the w values in Equation 3.
    CDiff=bsxfun(@minus,C,reshape(u,[1,1,n3]));
    %The maxIdx is used in the computation of the subgradient.
    [w,maxIdx]=max(CDiff,[],3);

    %Perform the maximization in Equation 4 over xi. As opposed to
    %directly returning xi, the columns for rows (nonzero elements)
    %are returned.
    n1Limits=ones(n1,1);
    [X, maxVal]=solveTransportationProblem(w,n1Limits,n2Limits,true);
    if(isempty(X))%If the problem is not feasible.
        tuples=[];
        minCostVal=[];
        costGap=[];
        relGap=[];
        return;
    end
    
    for curRow=1:n1
       dim2ForDim1Cur(curRow)=find(X(curRow,:));
    end

    phi=maxVal+sum(u);%The dual cost.
    ub=min(ub,phi);%Adjust the upper bound of the dual cost.

    %Next, we have to get a primal solution.            
    %The unnumbered equation before Equation 18.
    a=zeros(n1,n3);
    for i=1:n1
        a(i,:)=reshape(C(i,dim2ForDim1Cur(i),:),[1,n3]);
    end

    %Perform the miximization in 8 to get eta so the primal
    %variables can be found.
    %If custom limits over the third dimension are provided.
    if(~isempty(n3Limits))
        n1Limits=ones(n1,1);
        [X, primalVal]=solveTransportationProblem(a,n1Limits,n3Limits,true);
        if(isempty(X))%If the problem is not feasible.
            tuples=[];
            minCostVal=[];
            costGap=[];
            relGap=[];
            return;
        end

        for curRow=1:n1
           dim3ForDim1Cur(curRow)=find(X(curRow,:));
        end
    else
        [dim3ForDim1Cur,~,primalVal]=assign2D(a,true);
    end

    if(primalVal>lb)
        lb=primalVal;

        %Save the best feasible solution found thus far so that it
        %can be returned.
        phi2=dim2ForDim1Cur;
        phi3=dim3ForDim1Cur;
        minCostVal=lb;
    end

    %Check for convergence as in Algorithm 0.
    costGap=ub-lb;
    relGap=costGap/abs(ub);
    %The isfinite condition is added because if costGap==0 and
    %ub==0, then the relative ratio will be a NaN and an
    %if-statement would have an error.
    if(costGap<=eps(lb)||(~isfinite(relGap)&&(relGap<epsVal)))
       break; 
    end

    %For the subgradient, we have to compute xhat from xi according
    %to the unnumbered equations in the middle of Equation 4. This
    %is done in the same manner as computing rho in Eq. 3.14b of
    %[1].
    x=zeros(n1,n2,n3);
    for i1=1:n1
        x(i1,dim2ForDim1Cur(i1),maxIdx(i1,dim2ForDim1Cur(i1)))=1;
    end

    %The direction of search from Eq. 10.
    v=1-reshape(sum(sum(x,1),2),[n3,1]);

    if(all(v==0))
        break;
    end

    %Compute the stepsize from Eq. 11.
    sigma=(phi-lb)/norm(v)^2;

    %Take the step.
    u=u+sigma*v;

    %The bounding of the step to enforce the inequality constraints
    %is described in the steps before Equation 3.
    u=max(0,u);
end    

%Adjust the gains for the case where the initial cost matrix is transformed
%so that minimization can be performed.
if(maximize==false)
    minCostVal=-minCostVal;
    costGap=-costGap;
end

%Take the phi values and make the set of tuples to return.
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
