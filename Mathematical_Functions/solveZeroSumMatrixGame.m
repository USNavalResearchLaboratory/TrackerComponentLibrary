function [optVal,x,y,exitCodes]=solveZeroSumMatrixGame(A,linProgOpts)
%%SOLVEZEROSUMMATRIXGAME Determine the optimal strategies of two players in
%               a zero-sum matrix game. This is a problem where player 1
%               needs to choose a y vector to minimize y'*A*x and player 2
%               such that sum(y)=1 and all(y>=0), and player 2 chooses a
%               x vector to maximize y'*A*x such that sum(x)=1 and
%               all(x>=0). The choice of a row of A can be through of as a
%               strategy choice for player 1 and the choice of a column of
%               A as a strategy choice for player 2.  If one player uses a
%               deterministic strategy, then the other player can usually
%               beat him easily. Thus, the optimal values of x and y
%               represent probabilities of choosing each strategy. This
%               type of game is often taught in introductory game theory
%               classes and can be solved using linear programming.
%
%INPUTS: A A numRowsXnumCols matrix of the payoffs. player 1 (row choice)
%          wants to minimize the payoff, player 2 (column choice) wants to
%          maximize the payoff. Payoffs can be negative.
% linProgOpts Linear programming is used to solve the problem. Thus, this
%          is an optional structure that passes parameters to the linear
%          programming algorithm. members of the structure can be:
%          algorithm This is either 0, 1, or 2. 0 means use
%                    linProgRevisedSimplex to solve the problem; 1 means
%                    use linProgInteriorPoint; 2 means use linProgSOCP to
%                    solve the problem. The default is 0.
%          maxIter, epsilon, alpha These three parameters are only used if
%                    algorithm=0 or 1. These are the same as the inputs in
%                    the aforementioned linear programming algorithms. See
%                    the comments to those functions for details. The
%                    defaults if omitted are the defaults of the respective
%                    linear programming algorithms.
%             params This parameter is only used if algorithm=2. It
%                    corresponds to the params input of the linProgSOCP
%                    function. See the comments to that function for more
%                    information. The defaults if omitted are the default
%                    of that function.
%
%OUTPUTS: optVal The value of y'*A*x obtained solving the zero-sum matrix
%           game.
%         x A numColsX1 vector solving max_x min_y y'*A*x such that
%           sum(y)=1, sum(x)=1, all(y>=0) and all(x>=0), which happens
%           to equal min_y max_x y'*A*x. Note that this only finds one
%           solution, even if multiple solutions exist.
%         y A numRowsX1 vector solving min_y max_x y'*A*x , which happens
%           to equal max_x min_y y'*A*x under the same constraints as
%           given for x above.
%  exitCode If algorithm=0 or 1, this is a 2X1 vector of the values of the
%           exitCode parameter from the chosen linear programming
%           algorithm. See the comments for the chosen algorithm for more
%           details. if algorithm=2, this is a 2X1 cell array holding the
%           info output of linProgSOCP. The algorithm is run twice on
%           different problems (once to solve for x, another time to
%           solve for y given x). Thus, there are two outputs.
% 
%The transformation of a zero-sum matrix game into a linear programming
%problem is given in Chapter 11 of [1].
%
%EXAMPLE 1:
%An often-used example of a zero-sum game (also given in Chapter 11 of
%[1]) is rock-paper-scissors. Each player can choose rock, paper, or
%scissors. If they choose the same thing, then it is a draw. For rock
%versus scissors, scissors lose. For rock versus paper, rock loses. For
%paper versus scissors, paper loses. Saying player 1 wins has a cost of -1
%and player 2 wins has a cost of 1, then the matrix of payoffs (strategies
%ordered rock-paper-scissors) is
% A=[0, 1, -1;
%   -1, 0,  1;
%    1, -1, 0];
% [optVal,x,y,exitCodes]=solveZeroSumMatrixGame(A)
%In this instance, the optimal strategy for both is x,y=[1/3;1/3;1/3],
%which means to just choose uniformly among rock, paper and scissors.
%optVal is also 0. 
%
%EXAMPLE 2:
%The French version of rock, paper, scissors is rock, paper, scissors well,
%where the rock and scissors lose to the well, but the paper wins. The
%matrix of payoffs in such an instance is
% A=[0, 1,-1, 1;
%   -1, 0, 1,-1;
%    1,-1, 0, 1;
%   -1, 1,-1, 0];
% [optVal,x,y,exitCodes]=solveZeroSumMatrixGame(A)
%Here, the optimal strategy for both players is [0;1/3;1/3;1/3].
%
%EXAMPLE: 3
%The previous two examples have games that are symmetric (The matrix is
%square and the off-diagonal elements are the same with flipped signs). An
%example of an asymmetric game is from Chapter 11 of [1]
% A=[0,    0,      1/6,   1/6;
%    0,   -1/6,    1/3,   1/6;
%    1/6,  1/6,   -1/6,  -1/6;
%    1/6,  0,      0,    -1/6;
%   -1/6,  1/3,    0,     1/2;
%   -1/6,  1/6,    1/6,   1/2;
%    0,    1/2,   -1/3,   1/6;
%    0,    1/3,   -1/6,   1/6];
% [optVal,x,y,exitCodes]=solveZeroSumMatrixGame(A)
%Here, one might get x=[1/3;1/3;1/3;0]; and y=[2/3;0;1/3;0;0;0;0;0]; These
%are not the same as the solution in [1], though optVal=1/18 is the same
%for both solutions.
%
%EXAMPLE 4:
%In this example, the convergence of a standard interior point algorithm
%is rather slow even though the problem is very simple. For small problems,
%the simplex algorithm (the default) is often the best despite its
%theoretical worst-case performance.
% A=[4, 1, -1;
%    3, 2, 5;
%    0, 1, 6];
% tic;[optVal,x,y,exitCodes]=solveZeroSumMatrixGame(A);toc
% linProgOpts.algorithm=1;
% tic;[optValI,xI,yI,exitCodesI]=solveZeroSumMatrixGame(A,linProgOpts);toc
%Here, one find the optimal value to be about 3000/1375.
%
%REFERENCES:
%[1] R. J. Vanderbei, Linear Programming: Foundations and Extensions, 4th
%    ed. New York: Springer, 2014.
%
%June 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    algorithm=0;
    maxIter=[];
    epsilon=[];
    alpha=[];
    params=[];
else
    if(isfield(linProgOpts,'algorithm'))
        algorithm=linProgOpts.algorithm;
    else
        algorithm=0;
    end
    
    if(isfield(linProgOpts,'maxIter'))
        maxIter=linProgOpts.maxIter;
    else
        maxIter=[];
    end
    
    if(isfield(linProgOpts,'epsilon'))
        epsilon=linProgOpts.epsilon;
    else
        epsilon=[];
    end
    
    if(isfield(linProgOpts,'alpha'))
        alpha=linProgOpts.alpha;
    else
        alpha=[];
    end

    if(isfield(linProgOpts,'params'))
        params=linProgOpts.params;
    else
        params=[];
    end
end

numRows=size(A,1);
numCols=size(A,2);

c=[zeros(numCols,1);1];
ALeq=[-A,ones(numRows,1)];
bLeq=zeros(numRows,1);

AEq=[ones(1,numCols),  0];
bEq=1;
maximize=true;

switch(algorithm)
    case 0
        exitCodes=zeros(2,1);
        [optVal,x,exitCodes(1)]=linProgRevisedSimplex(AEq,bEq,ALeq,bLeq,c,maximize,maxIter,epsilon);
    case 1
        exitCodes=zeros(2,1);
        [optVal,x,exitCodes(1)]=linProgInteriorPoint(AEq,bEq,ALeq,bLeq,c,maximize,maxIter,epsilon,alpha);
    case 2
        exitCodes=cell(2,1);
        params.alpha=alpha;
        params.eps=epsilon;
        params=[];
        [optVal,x,exitCodes{1}]=linProgSOCP(AEq,bEq,ALeq,bLeq,c,maximize,true,params);
    otherwise
        error('Unknown algorithm specified')
end
x=x(1:numCols);

%If desired, solve the reverse problem to get the optimal value of y.
if(nargout>2)
    c=[zeros(numRows,1);1];
    ALeq=[A',-ones(numCols,1)];
    bLeq=zeros(numCols,1);

    AEq=[ones(1,numRows),  0];
    bEq=1;
    maximize=false;
    switch(algorithm)
        case 0
            [optVal,y,exitCodes(2)]=linProgRevisedSimplex(AEq,bEq,ALeq,bLeq,c,maximize,maxIter,epsilon);
        case 1
            [optVal,y,exitCodes(2)]=linProgInteriorPoint(AEq,bEq,ALeq,bLeq,c,maximize,maxIter,epsilon,alpha);
        case 2
            [optVal,y,exitCodes{2}]=linProgSOCP(AEq,bEq,ALeq,bLeq,c,maximize,true,params);
        otherwise
            error('Unknown algorithm specified')
    end
    y=y(1:numRows);
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
