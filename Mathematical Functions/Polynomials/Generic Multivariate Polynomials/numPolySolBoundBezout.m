function maxBound=numPolySolBoundBezout(termMats,method)
%%NUMPOLYSOLBOUNDBEZOUT Compute an upper bound on the number of solutions
%            to a system of n polynomial equations in n unknowns using one
%            of two types of Bézout bounds. Note that the bounds only
%            applies to systems that have a finite number of solutions
%            (non-degenerate systems).
%
%INPUTS: termMats  A cell array of matrices of terms of the polynomials.
%           Each term matrix is an (n+1)XnumTerms matrix such that
%           termMat(:,i)=[c,a1,a2,...,an] is a monomial term where c is the
%           value of of the monomial coefficient and the monomial is
%           x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1). For n variables,
%           there must be n term matrices. If in an equation only the first
%           nEq<n variables are used, then the term matrices still should
%           have space for the other coefficients, which are all set to
%           zero degree (as in the example below).
%       method An optional parameter indicating the method to use to
%           determine the upper bound. Possible values are:
%           0 (The default if omitted or an empty matrix is passed).
%             Compute the bound from Bézout's theorem. That is, compute the
%             total degree of the polynomials. This upper bound is equal to
%             the number of solutions of a homogenized system, but
%             generally greatly overestimates the number of solutions of
%             more common sparse systems. This bound is easy to compute.
%           1 Compute the multihomogenerous Bézout bound using brute force
%             to find the optiomal partition. This is significantly more
%             computationally intensive than 0. However, it is also
%             significantly better.
%           2 Compute an approximate multihomogeneous Bézout bound using a
%             method similar to that described in Section 3.4.2 of [1].
%
%OUTPUTS: maxBound The upper bound on the maximum number of solutions of
%                  the polynomial system.
%
%The bounds are described in Chapter 3 of [1]. Method 2 is similar to the
%heuristic algorithm in Section 3.4.3 of [1]. However, the description in
%[1] is inconsistent. Rather than going through all popssible partitions
%for finding the one that produces the minimum permanent of the degree
%matrix, as in the true multihomogenous Bézout bound, this function starts
%with everything in one partition and goes through the variables
%sequentially. A variable is either kept in the initial partition, moved
%to its own partition or put in another partition, depending on which has
%the least cost. This is a greedy algorithm.
%
%Though there may exist more efficient methods of computing the
%multihomogenous Bézout bound, it is shown in [2] that the multihomogeneous
%Bézout bound cannot be computed in polynomial time unlesss P=NP.
%
%EXAMPLE 1:
%The example here is that used in [1]. We have a system of eight equations.
%(Implicitly, this means all of these equations=0).
% termMats=cell(8,1);
% theEq='x1^2+x2^2-1';
% termMats{1}=string2Terms(theEq,8);
% theEq='x3^2+x4^2-1';
% termMats{2}=string2Terms(theEq,8);
% theEq='x5^2+x6^2-1';
% termMats{3}=string2Terms(theEq,8);
% theEq='x7^2+x8^2-1';
% termMats{4}=string2Terms(theEq,8);
% theEq='0.004731*x1*x3-0.3578*x2*x3-0.1238*x1-0.001637*x2-0.9338*x4+x7-0.3571';
% termMats{5}=string2Terms(theEq,8);
% theEq='0.2238*x1*x3+0.7623*x2*x3+0.2638*x1-0.07745*x2-0.6734*x4-0.6022';
% termMats{6}=string2Terms(theEq,8);
% theEq='x6*x8+0.3578*x1+0.004731*x2';
% termMats{7}=string2Terms(theEq,8);
% theEq='-0.7623*x1+0.2238*x2+0.3461';
% termMats{8}=string2Terms(theEq,8);
% maxBound=numPolySolBoundBezout(termMats,0)
% %Using the bound from Bézout's theorem, one gets 128, which is much
% %higher than the actual number of solutions.
% maxBound=numPolySolBoundBezout(termMats,1)
% %Though much slower, when using the multihomogenerous Bézout bound, one
% %gets a solution of 16, which coincides with the actual number of
% %solutions to the system.
% maxBound=numPolySolBoundBezout(termMats,2)
% %The approximate multihomogeneous Bézout bound gives 64, which is better
% %than the bound from Bézout's theorem, but is worse than the true
% multihomogeneous Bézout bound. On the other hand, the approximation is
% significantly faster than the exact multihomogeneous Bézout bound.
%
%EXAMPLE 2:
%This is the equation in example 2.3 of [1] with a=1.
% termMats=cell(3,1);
% theEq='0.5*x2^2+2*x2*x3+0.5*x3^2+x2^3*x3^2-1';
% termMats{1}=string2Terms(theEq,8);
% theEq='0.5*x3^2+2*x3*x1+0.5*x1^2+x3^3*x1^2-1';
% termMats{2}=string2Terms(theEq,8);
% theEq='0.5*x1^2+2*x1*x2+0.5*x2^2+x1^3*x2^2-1';
% termMats{3}=string2Terms(theEq,8);
% maxBound0=numPolySolBoundBezout(termMats,0)
% maxBound1=numPolySolBoundBezout(termMats,1)
% maxBound2=numPolySolBoundBezout(termMats,2)
%In this example, maxBound0 is 125, maxBound1 is 35, and maxBound2 is also
%35. In this instance, the approximate multihomogeneous Bézout bound equals
%the exact one.
%
%REFERENCES:
%[1] J. Verschelde, "Homotopy continuation methods for solving polynomial
%    systems," Ph.D. dissertation, Katholieke Universiteit Leuven, Leuven,
%    Belgium, May 1996.
%[2] G. Malajovich and K. Meer, "Computing multi-homogeneous B´ezout
%    numbers is hard," in Proceedings of the 22nd Annual Symposium on
%    Theoretical Aspects of Computer Science, Stuttgart, Germany, 24-26
%    Feb. 2005, pp. 244-255.
%
%March 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The number of variables/ number of equations.
n=length(termMats);

switch(method)
    case 0
        maxBound=1;
        for i=1:n
            curMat=termMats{i};
            d=max(sum(curMat(2:end,:),1));
            maxBound=maxBound*d;
        end
    case 1
        %The number of variables.= and equations
        maxBound=Inf;

        %Go through all bounds of the set into partitions.
        [q,recurVals]=getNextSetPartition(n);
        while(~isempty(q))
            %The number of items in the set.
            nc=recurVals.nc;
            p=recurVals.p;

            %Allocate space for the degree matrix.
            A=buildDegreeMatrix(termMats,q,nc);
            curBound=permGen(A,p(1:nc));
            maxBound=min(maxBound,curBound);

            [q,recurVals]=getNextSetPartition(q,recurVals);
        end
    case 2
        assignedPartition=ones(n,1);
        numInPartition=zeros(n,1);
        numInPartition(1)=n;
        numPartitions=1;
        
        A=buildDegreeMatrix(termMats,assignedPartition,numPartitions);
        maxBound=permGen(A,numInPartition(1:numPartitions));
        
        prevPartition.assignedPartition=assignedPartition;
        prevPartition.numInPartition=numInPartition;
        prevPartition.numPartitions=numPartitions;
        
        curMinPartition=prevPartition;
        
        for i=1:n
            %For each variable sequentially, we see if the bound is lower
            %if we keep the current partitioning scheme, or if we either
            %make the ith variable its own partition or move it to any of
            %the other partitions. Since the variables are being considered
            %sequentially, when starting, the ith variable is assigned to
            %the lowest numbered partition (1).
            
            %First, consider the case of the ith variable going into its
            %own partition. The test deals with the case for i=n if the ith
            %variable is already in its own partition.
            if(prevPartition.numInPartition(1)-1>0)
                numPartitions=prevPartition.numPartitions+1;
                assignedPartition=prevPartition.assignedPartition;
                assignedPartition(i)=numPartitions;
                numInPartition=prevPartition.numInPartition;
                numInPartition(1)=numInPartition(1)-1;
                numInPartition(numPartitions)=1;
                
                A=buildDegreeMatrix(termMats,assignedPartition,numPartitions);
                curBound=permGen(A,numInPartition(1:numPartitions));
                
                if(curBound<maxBound)
                    maxBound=curBound;
                    curMinPartition.assignedPartition=assignedPartition;
                    curMinPartition.numInPartition=numInPartition;
                    curMinPartition.numPartitions=numPartitions;
                end
            end
            
            %Next, we consider removing the ith variable from the first
            %partition and adding it to any of the other partitions. In
            %this instance, if i=n and there is only one item in the first
            %partition, we have to take into account that the first
            %partition will be going away. The prevPartition variable will
            %be modified so that this can be easily done in the loop below.
            if(prevPartition.numInPartition(1)-1==0)
                prevPartition.numPartitions=prevPartition.numPartitions-1;
                prevPartition.numInPartition=prevPartition.numInPartition(2:end);
                prevPartition.assignedPartition=prevPartition.assignedPartition-1;
                startPartition=1;
            else
                prevPartition.numInPartition(1)=prevPartition.numInPartition(1)-1;
                startPartition=2;
            end
            
            for curPartition=startPartition:prevPartition.numPartitions
                assignedPartition=prevPartition.assignedPartition;
                assignedPartition(i)=curPartition;
                numInPartition=prevPartition.numInPartition;
                numInPartition(curPartition)=numInPartition(curPartition)+1;
                
                A=buildDegreeMatrix(termMats,assignedPartition,prevPartition.numPartitions);
                curBound=permGen(A,numInPartition(1:prevPartition.numPartitions));
                
                if(curBound<maxBound)
                    maxBound=curBound;
                    curMinPartition.assignedPartition=assignedPartition;
                    curMinPartition.numInPartition=numInPartition;
                    curMinPartition.numPartitions=prevPartition.numPartitions;
                end
            end
            
            prevPartition=curMinPartition;
        end
    otherwise
        error('inValid method specified')
end

end

function A=buildDegreeMatrix(termMats,q,nc)

n=length(termMats);
for curEq=1:n
    terms=termMats{curEq};
    for curSet=1:nc            
        %Selected terms.
        idx=find(q==curSet)+1;
        sel=any(terms(idx,:)~=0,1);
        curDeg=max(sum(terms(idx,sel),1));

        if(~isempty(curDeg))
            A(curEq,curSet)=curDeg;
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
