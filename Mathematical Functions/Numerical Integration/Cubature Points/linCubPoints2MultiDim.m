function [xi, w]=linCubPoints2MultiDim(d,k,algorithm,linRule)
%%LINCUBPOINTS2MULTIDIM Generate d-dimensional cubature points and weights
%                   from a function to perform 1D quadrature integration.
%                   This does not produce the  minimum number of points for
%                   the polynomial order of the result. This works when
%                   considering extending rules that are defined over
%                   infinite regions or over square/ cubic regions. This
%                   includes methods for integration over Gaussian PDFs
%                   as well as methods for integrating over a unit cube
%                   among others.
%
%INPUTS: d A positive integer specifying the dimensionality of the points
%          to be generated.
%        k A number specifying the maximum number of points from the 1D
%          integration rule to consider using in multiple dimensions. When
%          using Gauss-Hermite interpolation points from the
%          quadraturePoints1D(n) function (the default), the minimum
%          polynomial order of algorithm 0 (tensor product method) is
%          2*k-1. That of algorithm 1 (Smolyak's method) is 2*(k-d)+1.
% algorithm An optional value indicating which algorithm to use to
%          transform the 1D quadrature methods into multiple dimensions.
%          Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            tensor product rule described, for example, in Section III of
%            [1]. Though this method generates a very large number of
%            points, when used with a 1D quadrature method that generates
%            all positive weights, the weights returned by this function
%            will also be all positive.
%          1 Use Smolyak's algorithm, descibed in [2]. This tends to
%            require fewer points for the same polynomial accuracy as the
%            tensor product method, but it can result in negative weights
%            even if the weights from the 1D quadrature method are
%            positive.
%  linRule An optional function handle such that [xi,w]=linRule(n) provides
%          n cubature points and weights in 1D. The polynomial order of
%          the points presumably increases with the number of points. When
%          using algorithm 0, linRule only needs to work with k as the
%          argument. When using algorithm 1, linRule must work with values
%          from 1 to k-d+1. If omitted or an empty matrix is passed, then a
%          handle to the quadraturePoints1D(n) function, for integration
%          with a Gaussian PDF as the weighting function, is passed. An
%          alternative could be, for example, the GaussLegendrePoints1D
%          function for integrating over the region from -1 to 1.
%
%OUTPUTS: xi A d X numCubaturePoints matrix containing the cubature points
%            (Each "point" is a vector).
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%For more details on how to use these points when integrating over a
%Gaussian PDF, see the comments in the function fifthOrderCubPoints.m.
%
%REFERENCES:
%[1] K. Ito and K. Xiong,"Gaussian filters for non linear filtering
%    problems," IEEE Transactions on Automatic Control, vol. 45, no. 5,
%    p. 910, May 2000.
%[2] V. Kaarnioja, "Smolyak quadrature," Master's thesis, University of
%    Helsinki, Department of Mathematics and Statistics, Helsinki, Finland,
%    May 2013.
%    [Online]. Available: https://helda.helsinki.fi/bitstream/handle/10138/40159/thesis.pdf
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
   algorithm=0; 
end

if(nargin<4||isempty(linRule))
    linRule=@(i)quadraturePoints1D(i);
end

switch(algorithm)
    case 0
        [xi,w]=tensorProductQuad(d,k,linRule);%This function is below.
    case 1
        [xi,w]=SmolyakQuad(d,k,linRule);%This function is below.
    otherwise
       error('Invalid Algorithm Specified') 
end

end

function [xi,w]=tensorProductQuad(numDim,m,linRule)
%%TENSORPRODUCTQUAD Generate n-dimensional cubature points and weights of 
%                   polynomial order 2m-1. This does not produce the 
%                   minimum number of points for the selected polynomial
%                   order.
%
%INPUTS:    numDim A positive integer specifying the dimensionality of the
%                  points to be generated.
%           m      A number such that 2m-1 is the desired order of the
%                  cubature points.
%          linRule A function handle such that [xi,w]=linRule(m)
%                  provides m cubature points and weights in 1D.
%
%OUTPUTS:   xi      A numDim X numCubaturePoints matrix containing the
%                   cubature points. (Each "point" is a vector)
%           w       A numCubaturePoints X 1 vector of the weights
%                   associated with the cubature points.
%
%This function uses the quadrature rules mentioned in Section III of
%K. Ito and K. Xiong,"Gaussian filters for non linear filtering problems,"
%IEEE Transactions on Automatic Control, vol. 45, no. 5, p. 910, May 2000.
%to generate vectors and weights for the multidimensional Gauss-Hermite
%quadrature rule.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[xL,wL]=linRule(m);

numPoints=m^numDim;
xi=zeros(numDim,numPoints);
w=zeros(numPoints,1);

%We now have the quadrature points for a one-dimensional integration. We
%need to take all possible groups of n of these values, allowing repeats,
%to extend it to n dimensions. This means that there are m^n possible
%vectors.
curSel=1;
selVec=zeros(numDim,1);
while(curSel<=numPoints)
    xi(:,curSel)=xL(selVec+1);
    w(curSel)=prod(wL(selVec+1));
    
    %Move on to the next set of combinations.
    selVec=arbBaseInc(selVec,m);
    curSel=curSel+1;
end

end

function [xi,w]=SmolyakQuad(d,k,linRule)
%%SMOLYAKQUAD Generate n-dimensional cubature points and weights by
%             transforming a 1D quadrature function into multiple dimensions
%             using Smolyak's method. Though Smolyak's method results in
%             using fewer points than a simple tensor produce method, some
%             of the weights can be negative, which is undesirable when
%             performing covariance estimation.
%
%INPUTS:    d A positive integer >1 specifying the dimensionality of the
%             points to be generated.
%           k The integer order of the tensor product operator to use.
%             k>=d. Higher values of k mean higher orders of polynomial
%             estimation accuracy. When using Gauss-Hermite points as the
%             1D integration method, 2*(k-d)+1 is a lower bound for the
%             polynomial accuracy of the algorithm.
%     linRule An function handle such that [xi,w]=linRule(m) provides m
%             cubature points and weights in 1D for m=1 to k-d+1.
%
%OUTPUTS:   xi      A d X numCubaturePoints matrix containing the
%                   cubature points. (Each "point" is a vector)
%           w       A numCubaturePoints X 1 vector of the weights
%                   associated with the cubature points.
%
%This function implements Algorithm 1.12 in Chapter 1 of [1]. Unlike the
%algorithm given in Appendix C of [1], the total number of points needed
%is found before entering the loop so that everything can be preallocated.
%Also, no special checking is performed for the trivial case. 
%
%REFERENCES:
%[1] V. Kaarnioja, "Smolyak quadrature," Master's thesis, University of
%    Helsinki, Department of Mathematics and Statistics, Helsinki, Finland,
%    May 2013.
%    [Online]. Available: https://helda.helsinki.fi/bitstream/handle/10138/40159/thesis.pdf
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %Cell arrays with various sets of 1D cubature points.
    xiList=cell(k-d+1,1);
    wList=cell(k-d+1,1);
    num1DPoints=zeros(k-d+1,1);
    for i=1:(k-d+1)
        [xi,w]=linRule(i);
        xiList{i}=xi;
        wList{i}=w';
        num1DPoints(i)=length(w);
    end
    
    %Count the total number of cubature points that will be needed so that
    %sufficient memory can be allocated outside of the loop. This is based
    %on the rules on page 17 of [1].
    totalPoints=0;
    for L=max(d,k-d+1):k
        tuples=genAllTCompositions(L,d);
        totalPoints=totalPoints+sum(prod(num1DPoints(tuples),1));
    end
    
    %Initialize the nodes and weights to empty matrices; 
    xi=zeros(d,totalPoints);
    w=zeros(totalPoints,1);
    numAdded=0;
    for L=max(d,k-d+1):k
        %All tuples of dimension d whose 1-norm is equal to L.
        tuples=genAllTCompositions(L,d)';
        numTuples=size(tuples,1);

        for i=1:numTuples
           %Make all combinations of univariate quadrature nodes and
           %weights. 

            alpha=tuples(i,:);
            tempNodes=xiList{alpha(1)};
            tempWeights=wList{alpha(1)};

            for j=2:d
                temp1=xiList{alpha(j)};
                temp2=wList{alpha(j)};

                %Determine all vector combinations of the univariate
                %quadrature nodes and weights.
                tempNodes=combVec2D(tempNodes,temp1);
                tempWeights=combVec2D(tempWeights,temp2);
            end

            %Account for the Smolyak quadrature coefficient 
            tempWeights=(-1)^(k-L)*binomial(d-1,k-L)*prod(tempWeights);
            num2Add=size(tempWeights,2);
            newNumAdded=numAdded+num2Add;
            xi(:,(numAdded+1):newNumAdded)=tempNodes;
            w((numAdded+1):newNumAdded)=tempWeights;
            numAdded=newNumAdded;
        end
    end
    
    function combVec=combVec2D(vec1,vec2)
    %This subroutine is based on Appendix C.2 of [1]. It creates the vector
    %combinations on page 21.
        combVec=[kron(vec1,ones(1,length(vec2)));
                reshape(vec2'*ones(1,length(unique(vec1,'rows'))),1,[])];
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
