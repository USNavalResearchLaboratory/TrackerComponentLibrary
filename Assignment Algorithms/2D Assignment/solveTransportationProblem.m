function [X,gain,exitCode]=solveTransportationProblem(C,s,d,maximize,maxIter)
%%SOLVETRANSPORTATIONPROBLEM Solve the transportation problem using a
%       strong polynomial-time algorithm. The transportation problem is a
%       generalization of the 2D assignment problem.
%
%INPUTS: C A numRowXnumCol cost matrix that does not contain any NaNs and
%          where the largest finite element minus the smallest element is a
%          finite quantity (does not overflow) when performing minimization
%          and where the smallest finite element minus the largest
%          element is finite when performing maximization.  Forbidden
%          assignments can be given costs of +Inf for minimization and
%          -Inf for maximization.
%        s A numRowX1 vector of positive supply for each row.
%        d A numColX1 vector of positive demands for each column.
% maximize If true, the minimization problem is transformed into a
%          maximization problem. The default if this parameter is omitted
%          or an empty matrix is passed is false.
%  maxIter An optional parameter giving the maximum number of iterations to
%          be performed in the minCostFlow subroutine. If omitted, the
%          default value in the minCostFlow function is used.
%
%OUTPUTS: X The flow matrix. This is explained with the further definition
%           of the problem below. If the minCostFlow subroutine failed due
%           to problem infeasbility, this will be an empty matrix.
%      gain The total cost of the minimization/ maximization. If the
%           minCostFlow subroutine failed due to problem infeasbility, this
%           will be an empty matrix.
%  exitCode An integer noting the success or failure of the minCostFlow
%           subroutine, and thus of the entire algorithm. See the comments
%           to the minCostflow function for more information.
%
%Depending on whether total supply (sum(s)) is greater than equal to or
%less than total demand (sum(d)), the exact problem being solved varies a
%bit. The transportation problem generally determines a matrix X with all
%X(i,j)>=0 so as to minimize the transportation cost sum(sum(C.*X)) subject
%to the constraints:
%1) Supply Constraints:          0<sum(X(i,:))<=s(i) for all i
%2) Minimum Demand Constraints:  0<sum(X(:,j))>=d(j) for all j
%For any feasible solution, sum(s)>=sum(d); That is, supply is at least as
%high as demand. This means that for any feasible instance, the second
%constraint becomes
%2) Fixed Demand Constraints:  0<sum(X(:,j))=d(j) for all j
%However, if sum(s)<sum(d), this function assumes a different set of
%constraints. Specifically:
%1) Supply Constraints:          0<sum(X(i,:))=s(i)  for all i
%2) Minimum Demand Constraints:  0<sum(X(:,j))<=d(j) for all j
%
%As noted in Chapters 7.4 and 7.5 of [1], the transportation problem is
%equivalent to the minimum cost flow problem. Thus, this function
%transforms the transportation problem into a minimum cost flow problem and
%uses the function minCostFlow to solve the problem.
%
%The transportation problem gets its name from the types of supply problems
%that often take its form. Consider the example from Chapter 5 of [2]:
%EXAMPLE 1: A balanced model
% C=[80, 215;
%    100,108;
%    102,68];
% s=[1000;1500;1200];
% d=[2300;1400];
% [X,gain,exitCode]=solveTransportationProblem(C,s,d,false)
%Here the minimum (cost) is 313200. In the above example, the supply and
%demand represented cars, the rows represented possible locations of
%factories, the columns cities where the cars would be sold and the costs
%transportation costs of moving a car from the factory to the city where it
%is to be sold. The optimal solution is
%X=[1000, 0;
%   1300, 200;
%   0, 1200];
%
%EXAMPLE 2:
%The first example was balanced. This example, also from [2], is a case
%where there is more demand than supply.
% C=[80, 215;
%    100,108;
%    102,68];
% s=[1000;1300;1200];
% d=[2300;1400];
% [X,gain,exitCode]=solveTransportationProblem(C,s,d,false)
%Here the optimal gain is 291600 and the assignment matrix is
%X =[1000, 0;
%    1300, 0;
%    0, 1200];
%
%EXAMPLE 3:
%On the other hand, if one wishes to run the company into the ground, one
%could maximize the cost instead. using the case for the first example,
%this is
% C=[80, 215;
%    100,108;
%    102,68];
% s=[1000;1500;1200];
% d=[2300;1400];
% [X,gain,exitCode]=solveTransportationProblem(C,s,d,true)
%where one gets a cost of 490600 and an assignment matrix of 
%X =[0,   1000;
%    1100,400;
%    1200,0];
%
%EXAMPLE 3: One can also consider the case of more supply than demand,
%whereby all of the demand is met and only some of the supply is used.
% C=[80, 215;
%    100,108;
%    102,68];
% s=[1000;1500;1200];
% d=[1900;1400];
% [X,gain,exitCode]=solveTransportationProblem(C,s,d,false)
%Here the cost is 273200 and the assignment is
%X =[1000, 0;
%    900,  200;
%    0,    1200];
%
%EXAMPLE 4: The algorithm will also work if there are state subsidies that
%lead to negative transportation costs. For example,
% C=[80, -215;
%    100,108;
%    102,68];
% s=[1000;1500;1200];
% d=[2300;1400];
% [X,gain,exitCode]=solveTransportationProblem(C,s,d,false)
%Now the cost is only 43800 and the solution is
% X =[0,   1000;
%     1500,0;
%      800,400];
%
%EXAMPLE 4:
%The transportation problem is a generalization of a 2D assignment problem.
%We can also consider a basic rectangular 2D assignment problem:
% C=[Inf,  2,   Inf,Inf,3;
%      7,  Inf, 23, Inf,Inf;
%     17,  24,  Inf,Inf,Inf;
%    Inf,  6,   13, 20, Inf];%2D cost matrix.
% s=[1;1;1;1];
% d=[1;1;1;1;1];
% [X,gain,exitCode]=solveTransportationProblem(C,s,d,false)
%Here, the minimum cost is 47 and the assignment matrix is 
%X =[0, 0, 0, 0, 1;
%    1, 0, 0, 0, 0;
%    0, 1, 0, 0, 0;
%    0, 0, 1, 0, 0];
%which is equivalent to the solution one would get passing the same
%assignment matrix to the assign2D algorithm.
%
%EXAMPLE 5:
%Note that the supply and demand do not have to be integers. A number of
%authors mention the theorem that as long as the supplies, demands and
%capacities are integer, then the solution will be integer. However,
%non-integer solutions are possible with non-integer supplies,capacities or
%demands. Consdier, the modified first example
% C=[80.1, 215.7;
%    100.2,108.5;
%    102.0,68.6];
% s=[1000.2;1500.3;1200];
% d=[2300;1400];
% [X,gain,exitCode]=solveTransportationProblem(C,s,d,false)
%In this instance, the gain is 314375.98 and the optimal association matrix
%is
% X =[1000.2, 0;
%     1299.8,200;
%     0   1200];
%
%REFERENCES:
%[1] C. H. Papadimitriou and K. Steiglitz, Combinatorial Optimization:
%    Algorithms and Complexity. Englewood Cliffs, NJ: Prentice-Hall Inc.,
%    1982.
%[2] H. A. Taha, Operations Research: An Introduction, 9th ed. Boston:
%    Prentice Hall, 2011.
%
%July 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(maxIter))
   maxIter=[];%Use the default in the minCostFlow function.
end

if(nargin<4||isempty(maximize))
   maximize=false; 
end

%The cost matrix must have all non-negative elements for the assignment
%algorithm to work. This forces all of the elements to be positive. The
%delta is added back in when computing the gain in the end.
if(maximize==true)
    CDelta=max(max(C));
    C=-C+CDelta;
else
    CDelta=min(min(C));
    C=C-CDelta;
end

totalSupplyUsed=min(sum(d),sum(s));

%Next, the problem will be reformulated as a minimum cost flow problem and
%solved using the minCostFlow function.

numRows=size(C,1);
numCols=size(C,2);
%Make each row a vertex and each column a vertex. an additional source and
%sink are added. The source is connected to all row vertices. The capacity
%of path from the source to the ith row vertex equals s(i). All column
%vertices are connected to the sink. The capacity path from each column 
%vertex to the sink equals d(i).
totalVertices=numRows+numCols+2;
sourceVertex=totalVertices-1;
sinkVertex=totalVertices;
AMat=zeros(totalVertices,totalVertices);%Cost matrix
CMat=zeros(totalVertices,totalVertices);%Capacity matrix
%The vertices connecting each row to each column get a capacity equal to
%the supply from that row and a cost equal to the given cost, assuming that
%an edge exists
for curRow=1:numRows
    for curCol=1:numCols
        if(C(curRow,curCol)~=Inf)
            AMat(curRow,numRows+curCol)=C(curRow,curCol);
            CMat(curRow,numRows+curCol)=min(s(curRow),d(curCol));
        end
    end
end

%Add in capacities for the edges associated with the source and sink nodes.
CMat(sourceVertex,1:numRows)=s(:);
CMat((numRows+1):(numRows+numCols),sinkVertex)=d(:);

%The demand vector is made assuming that the source can supply as much as
%possible to meet demand.
b=zeros(totalVertices,1);%The demand vector
b(sourceVertex)=totalSupplyUsed;
b(sinkVertex)=-totalSupplyUsed;

%Solve the optimization problem.
[F,gain,exitCode]=minCostFlow(AMat,CMat,b,maxIter);

%If a problem occurred so that minCostFlow failed.
if(isempty(F))
    X=[];
    gain=[];
    return;
end

%If F was obtained, then extract the X matrix.
colVertices=(numRows+1):(numRows+numCols);
X=F(1:numRows,colVertices);

%Adjust the gain for the initial offset of the cost matrix.
if(maximize==true)
    gain=-gain+CDelta*totalSupplyUsed;
else
    gain=gain+CDelta*totalSupplyUsed;
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
