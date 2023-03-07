function [x,M]=semidefQuadProgUnconst(B,a,b,A)
%%SEMIDEFQUADPROGUNCONST Solve an unconstrained quadratic programming
%        problem of the form: 
%          minimize_x  x'*B*x+2*a'*x
%        where B is a symmetric positive definite or a positive
%        semidefinite matrix. If B is positive semidefinite, multiple
%        solutions exist and one can choose a solution by providing b and
%        A inputs.
%
%INPUTS: B An nXn real, symmetric positive definite or positive
%          semidefinite matrix.
%        a A real nX1 vector.
%     b, A If B is semidefinite, then an infinite number of solutions
%          exist. These values choose a solution. If neither b nor A is
%          provided, then x is the solution minimizing norm(x). If just b
%          is provided, then x is the solution closest to b. If b and A are
%          provided, then x is the solution minimizing norm(A*x-b).
%
%OUTPUTS: x The nX1 solution.
%         M An nXn matrix. All other solutions are of the form
%           xOther=x+M*y, where y is a an nX1 real vector.
% 
%The solutions are described in Chapter 3 of [1].
%
%EXAMPLE:
%In this example, we solve an optimziation problem with a semidefinite B
%matrix. Since multiple solutions exist, we look at what one gets if 1) b
%and A are not provided, 2) only b is provided, and 3) if both b and A
%are provided. We show that the solution numbers that minimize the norm of
%x, is closest to b, and that minimize the distance between A*X and b are
%all as expected.
% %Create a positive semidefinite matrix. Here, we fix the eigenvalues but
% %choose random eigenvectors. Having at least one zero eigenvalue, with
% %all others positive, makes it semidefinite.
% V=randOrthoMat(4);
% D=diag([10;2;3;0]);
% B=V*D*inv(V);
% a=[1;2;-3;6];
% b=[-4; -7; 6; -3];
% A=[-23,  0,  6,   3;
%    -15,  5, -9,  -9;
%     21, 19, 15, -13;
%    -12, -8,  2,   3];
% 
% xMinNorm=semidefQuadProgUnconst(B,a);
% xNearB=semidefQuadProgUnconst(B,a,b);
% xAxNearB=semidefQuadProgUnconst(B,a,b,A);
% [~,smallestNormSol]=min([norm(xMinNorm),norm(xNearB),norm(xAxNearB)])
% [~,solClosest2b]=min([norm(xMinNorm-b),norm(xNearB-b),norm(xAxNearB-b)])
% [~,solAxNearestb]=min([norm(A*xMinNorm-b),norm(A*xNearB-b),norm(A*xAxNearB-b)])
%
%REFERENCES:
%[1] D. L. Nelson, "Quadratic programming techniques using matrix
%    pseudoinverses," Ph.D. dissertation, Texas Tech University, Lubbock,
%    TX, Dec. 1969.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4)
    A=[];
end

if(nargin<3)
    b=[];
end

numDim=length(a);
I=eye(numDim,numDim);

BPlus=pinv(B);
x=BPlus*a;
M=(I-BPlus*B);

if(rank(B)~=numDim)
    %Skip in the positive definite case.

    if(~isempty(b)&&isempty(A))
        %Choose the solution closes to b.
        %Case 3.3.1 in [1].
        x=x+M*b;
    elseif(~isempty(b)&&~isempty(A))
        %Choose the solution that minimizes the distance of the vector A*x
        %from b. Case 3.3.2 in [1]. If B is positive definite, if one were
        %to execute the following lines, one would probably get a bad
        %result, becausse M would be very small and y would be very big and
        %finite precision errors would dominate. Hence, it is best to test
        %the rank and skip this setp if B is full rank.
        y=pinv(A*M)*(b-A*BPlus*a);
        x=x+M*y;
    elseif(~isempty(A))
        error('If A is provided, then b must be provided.')
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
