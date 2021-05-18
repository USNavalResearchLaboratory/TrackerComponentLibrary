function [intEst,totalErr,exitCode]=integrateUnifCub57Adaptive(f,lowerBounds,upperBounds,maxSearchReg,AbsTol,RelTol)
%%INTEGRATEUNIFCUB57ADAPTIVE Perform adaptive numerical integration over a
%                    function that takes a multidimensional input. The
%                    output of the function can be a scalar, a vector, or a
%                    matrix. This integration does not use a weighting
%                    function.
%
%INPUTS: f The handle to the function. f(x) takes a single numDimX1 vector,
%          where numDim>=2, and returns a scalar, vector or matrix output.
% lowerBounds,upperBounds The lower and upper bounds of integration as
%          numDimX1 vectors. It is assumed that all elements of lowerBounds
%          are <= those of upperBounds and that all bounds are finite.
% maxSearchReg The algorithm works by splitting space into increasingly
%          small regions. This optional parameter is the maximum number of
%          regions allowed. If omitted or an empty matrix is passed, the
%          default value of 500 is used. The number of function evaluations
%          to split each region is a constant that depends on the
%          dimensionality of the problem.
%   AbsTol The absolute error tolerance allowed,a positive scalar or a
%          vector/ matrix if the output is multidimensional and different
%          tolerances apply to different dimensions. If omitted or an empty
%          matrix is passed, the default value of 1e-12 is used.
%   RelTol The maximum relative error tolerance allowed, a positive scalar
%          or a vector/ matrix if the output is multidimensional and
%          different tolerances apply to different dimensions. If omitted
%          or an empty matrix is passed, the default value of 1e-9 is used.
%
%OUTPUTS: intEst The estimated value of the integral. The dimensionality of
%                this is the dimensionality of the output of f.
%       totalErr The estimated error. This has the same dimensionality as
%                intEst, representing the error in each dimension of the
%                output.
%       exitCode A value indicating the status of the algorithm on
%                termination. Possible values are
%                0 Termination occurred due to the absolute or relative
%                  error tolerances being fulfilled.
%                1 Termination occurred due to the search region having
%                  been split into the maximum number of subregions.
%
%The algorithm used is a slightly modified version of that of [1]. In [1],
%only functions with scalar outputs can be used. however, here, we allow
%multiple outputs. The algorithm of [1] is a modification of that of [2].
%Both utilize fifth and seventh-order cubature points, whereby the
%fifth-order points are a subset of the seventh-order ones. Integration
%over a region is performed and the dimension over which it is likely that
%the most benefit from splitting the region could provide is identified via
%numerical second derivatives. The error of the integral of the region
%comes from the difference between the fifth and seventh-order estimates.
%The algorithm thus, splits the region and keeps splitting the region with
%the highest error.
%
%In this implementation, a binary heap is used to store the regions instead
%of the semi-ordered tree used in [2]. Variants of this algorithm might
%consider the use of different sets of cubature points. However, a subset
%of the cubature points chosen here are used for the finite differences for
%
%EXAMPLE:
%Here, we consider a five-dimensional integral that is also used as an
%example in KorobovLatticeRules. It is worth noting how much better this
%function does than those rules. We consider the function sum(abs(x-1/2))
%in five dimensions with bounds as follows:
% f=@(x)sum(abs(x-1/2),1);
% lowerBounds=[-2;-1;0;4;3];
% upperBounds=[2;1;1;6;8];
% %The exact solution against which we will compare is easily found:
% exactSol=915;
% [intEst,totalErr,exitCode]=integrateUnifCub57Adaptive(f,lowerBounds,upperBounds)
%In this instance, intEst is equal to the exact solution.
%
%REFERENCES:
%[1] A. C. Genz and A. A. Malik, "Remarks on algorithm 006: An adaptive
%    algorithm for numerical integration over an n-dimensional rectangular
%    region," Journal of Computational and Applied Mathematics, vol. 6, no.
%    4, pp. 295-302, Dec. 1980.
%[2] P. van Dooren and L. de Ridder, "An adaptive algorithm for numerical
%    integration over an n-dimensional cube," Journal of Computational and
%    Applied Mathematics, vol. 2, no. 3, pp. 207-217, Sep. 1976.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=length(lowerBounds);

if(nargin<6||isempty(RelTol))
    RelTol=1e-9;
end

if(nargin<5||isempty(AbsTol))
    AbsTol=1e-12;
end

if(nargin<7||isempty(maxSearchReg))
    maxSearchReg=500;
end

%This binary heap stores the regions that are to be searched in order of
%the maximum error of the region.
searchRegions=BinaryHeap(maxSearchReg,true);

if(numDim<2)
   error('This function requires that the dimensionality of the integration region be >=2.') 
end

%The initial integration region is the whole space.
center=(lowerBounds+upperBounds)/2;
width=(upperBounds-lowerBounds)/2;

%If the integration volume is zero, then just retuen zero.
if(prod(width)==0)
   intEst=0;
   totalErr=0;
   exitCode=0;
   return;
end

%This will hold the integral estimate.
intEst=0;
%This will hold the sum of the errors of all of the subregions.
totalErr=0;

%The degree-7 BF rule in [1] is actually just a set of cubature points and
%weights. These lambda values go into making the cubature points. The
%degree-5 rule is just a subset of the degree 7 rule.
lambda2=sqrt(9/70);
lambda3=sqrt(9/10);
lambda4=sqrt(9/10);
lambda5=sqrt(9/19);

%cubature weights for the degree-seven cubature points.
%The 2^numDim term in the paper has been removed and is multiplied back
%with the region volume when using the cubature weights.
w1=(12824-9120*numDim+400*numDim^2)/19683;
w2=(980/6561);
w3=(1820-400*numDim)/19683;
w4=(200/19683);
w5=6859/19683/2^numDim;

%The first subset of cubature points (making up the fifth and seventh-order
%points). Just the origin.
xi1=zeros(numDim,1);

%Half of the second and third subsets of cubature points (making up the
%fifth and seventh-order points). The whole second set and + and - of these
%points, but we only store half of them as we need to keep track of each
%pair of function evaluations to perform finite differencing for
%determining which dimension should be split.
xi2Half=lambda2*eye(numDim);
xi3Half=lambda3*eye(numDim);
numXi23HalfPoints=numDim;

%The four subset of cubature points (making up the fifth and seventh-order
%points).
xi4=fullSymPerms([lambda4;lambda4;zeros(numDim-2,1)]);
numXi4Points=2*(numDim^2-numDim);

%The fifth subset of cubature points (making up the seventh-order points).
xi5=PMCombos(lambda5*ones(numDim,1));
numXi5Points=2^numDim;

%The degree-5 BF rule is similarly just a set of cubature points and
%weights given on the same page as the degree-7 rule. The points are a
%subset of the degree 7 points, but the weights are different.
wp1=(729-950*numDim+50*numDim^2)/729;
wp2=(245/486);
wp3=(265-100*numDim)/1458;
wp4=(25/729);

%This is false when the last time around the loop the top region in the
%searchRegions heap is halfway through dividing. That means that the
%results go into the thing already in the heap and the next time around, we
%investigate the other half of the split dimension.
shouldDivide=true;
while(1)
    %The volume of the current region being investigated.
    regionVol=2^numDim*prod(width);

    %The first cubature point --the one at the origin of the region. The
    %associated weight is w1.
    sum1=f(center+xi1.*width);
    
    %The next two sets of cubature points are fully symmetric sums of
    %f([lambda2;0;0;...;0]) and f([lambda3;0;0;...;0]), where the position
    %and sign of lambda2 and lambda4 change. For each pair of function
    %evaluations with the lambda in the same position but different signs,
    %the second derivative is numerically estimated. These numerical second
    %derivatives are then combined into a fourth-divided difference. T
    %The dimension with the largest fourth divided difference  will be used
    %for adaptively splitting the region.
    sum2=0;
    sum3=0;
    diffMax=0;
    maxIdx=1;
    for curPoint=1:numXi23HalfPoints
        %Cubature point/ function evaluation using lambda2 (associated with
        %w2).
        f1=f(center-xi2Half(:,curPoint).*width);
        f2=f(center+xi2Half(:,curPoint).*width);

        %Cubature point/ function evaluation using lambda3 (associated with
        %w3).
        f3=f(center-xi3Half(:,curPoint).*width);
        f4=f(center+xi3Half(:,curPoint).*width);
        
        sum2=sum2+f1+f2;
        sum3=sum3+f3+f4;
        deltaF1=f1+f2-2*sum1;
        deltaF2=f3+f4-2*sum1;

        %This is the DIF_i equation in [2]. This is essentially a constant
        %times an average of two finite difference estimates of the second
        %dervative along a particular dimension.
        diff=abs(deltaF1-(lambda2/lambda3)^2*deltaF2);
        %This lets us deal with the output of the function being
        %multidimensional.
        diff=max(diff(:));
        if(diff>diffMax)
            diffMax=diff;
            maxIdx=curPoint;
        end
    end

    %Next, we consider the set of cubature points that are fully symmetric
    %sums of the form f([lambda4;lambda4;0;...;0]).
    sum4=0;
    for curPoint=1:numXi4Points
        sum4=sum4+f(center+xi4(:,curPoint).*width);
    end
    
    %Next, go through all 2^numDim combinations that are +/- variations of 
    %f([lambda5;lambda5;lambda5;...;lambda5])
    sum5=0;
    for curPoint=1:numXi5Points
        sum5=sum5+f(center+xi5(:,curPoint).*width);
    end
    
    %Compute the values of the fifth-and seventh degree rules. The
    %seventh-order value is used as the "correct" value and the fifth-order
    %value is used for comparison to get an estimate of the error in this
    %region.
    fifthOrderEst=regionVol*(wp1*sum1+wp2*sum2+wp3*sum3+wp4*sum4);
    seventhOrderEst=regionVol*(w1*sum1+w2*sum2+w3*sum3+w4*sum4+w5*sum5);
    
    regionError=abs(seventhOrderEst-fifthOrderEst);
    
    %Add the value of this region to the integral estimate.
    intEst=intEst+seventhOrderEst;
    
    %Add the error to the estimate of the total error.
    totalErr=totalErr+regionError;
    
    if(shouldDivide==true)
        %Determine whether either of the error bounds has been met.
        if(all(totalErr(:)<AbsTol(:))||all(totalErr(:)./max(abs(intEst(:)))<RelTol))
            exitCode=0;
            return;
        end
        
        %If the heap of search regions is full, then terminate the
        %algorithm.
        if(searchRegions.heapSize()==maxSearchReg)
           exitCode=1;
           return;
        end

        %Otherwise, insert this region into the heap of search regions.
        regionParams.regionValue=seventhOrderEst;
        regionParams.center=center;
        regionParams.width=width;
        regionParams.splitIdx=maxIdx;
        regionParams.regionError=regionError;
        %When multidimensional, the key is the dimension with the highest
        %error.
        searchRegions.insert(max(regionError(:)),regionParams);
        
        %Now, we will divide the region in the heap with the highest error.
        topRegion=searchRegions.getTop();
        
        regionParams=topRegion.value;
        center=regionParams.center;
        width=regionParams.width;
        splitIdx=regionParams.splitIdx;
        
        width(splitIdx)=width(splitIdx)/2;
        center(splitIdx)=center(splitIdx)-width(splitIdx);
        
        %Since the region is dividing, we will remove its contribution to
        %the integral and the total error.
        intEst=intEst-regionParams.regionValue;
        totalErr=totalErr-regionParams.regionError;
        
        shouldDivide=false;
        continue;
    else%If here, then  we divided the current region last time in the
        %loop. Now, we have to investigate the other half.
        
        %Remove the now invalid old region.
        topRegion=searchRegions.deleteTop();
        splitIdx=topRegion.value.splitIdx;
        
        %Insert this half region into the heap of search regions in place
        %of the old one.
        regionParams.regionValue=seventhOrderEst;
        regionParams.center=center;
        regionParams.width=width;
        regionParams.splitIdx=maxIdx;
        regionParams.regionError=regionError;
        searchRegions.insert(max(regionError(:)),regionParams);
        
        %Adjust the center to the non-investigated half.
        center(splitIdx)=center(splitIdx)+2*width(splitIdx);
        
        shouldDivide=true;
        continue;
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
