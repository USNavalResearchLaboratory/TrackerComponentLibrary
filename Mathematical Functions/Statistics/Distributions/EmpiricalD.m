classdef EmpiricalD
%Functions to handle the empirical distribution for scalar samples. This
%is a discrete distribution whose values are determined by samples of data.
%Optionally, the samples can be weighted.
%Implemented methods are: mean, var, PMF, CDF, rand
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
function val=mean(empData,w)
%%MEAN  Obtain the mean of the empirical distribution for the specified
%       samples. Weights can be provided if the samples are weighted.
%
%INPUTS: empData An NX1 or 1XN vector of N scalar samples that make up the
%                empirical distribution. Values can be repeated.
%              w An optional NX1 or 1XN vector of weights associated with
%                the samples. These should be positive and sum to 1. If
%                omitted or an empty matrix is passed, the samples are
%                assumed to be uniformly weighted.
%
%OUTPUTS: val  The mean of the empirical distribution under consideration.
%
%The empirical distribution is discussed in many introductory statistics
%textbooks, such as in Chapter 4-2 of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    %Use uniform weighting if weighting is not specified.
    if(nargin<2||isempty(w))
       val=mean(empData(:));
       return;
    end

    val=sum(w(:).*empData(:));
end
    
function val=var(empData,w)
%%VAR   Obtain the variance of the empirical distribution for the specified
%       samples. Weights can be provided if the samples are weighted.
%
%INPUTS: empData An NX1 or 1XN vector of N scalar samples that make up the
%                empirical distribution. Values can be repeated.
%              w An optional NX1 or 1XN vector of weights associated with
%                the samples. These should be positive and sum to 1. If
%                omitted or an empty matrix is passed, the samples are
%                assumed to be uniformly weighted.
%
%OUTPUTS: val  The scalar variance of the empirical distribution under
%              consideration.
%
%The empirical distribution is discussed in many introductory statistics
%textbooks, such as in Chapter 4-2 of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    %Use uniform weighting if weighting is not specified.
    if(nargin<2||isempty(w))
        numEls=numel(empData);
        mu=mean(empData);
        diff=empData(:)-mu(:);
        val=sum(diff.*diff)/numEls;
        return;
    end

    mu=EmpiricalD.mean(empData,w);
    diff=empData-mu;
    val=sum((diff(:).*diff(:)).*w(:));
end

function val=PMF(x,empData,w,isSortedNoRepeats)
%%PMF Evaluate the empirical probability mass function (PMF) at given
%     points. Weights can be provided if the samples are weighted.
%
%INPUTS:    x The point(s) at which the empirical PMF is to be evaluated.
%             Values in x that are not in empData have a PMF value of 0.
%     empData An NX1 or 1XN vector of N scalar samples that make up the
%             empirical distribution. Values can be repeated.
%           w An optional NX1 or 1XN vector of weights associated with the
%             samples. These should be positive and sum to 1. If omitted or
%             an empty matrix is passed, the samples are assumed to be
%             uniformly weighted.
% isSortedNoRepeats If the points in empData are sorted in ascending order
%             and no points are repeated, then the algorithm is faster.
%             This indicates whether the values in empData are sorted in
%             ascending order without repeats. If this parameter is omitted
%             or an empty matrix is passed, a default value of false is
%             used.
%
%OUTPUTS:  val The value(s) of the empirical PMF evaluated at the given
%              points.
%
%The empirical distribution is discussed in many introductory statistics
%textbooks, such as in Chapter 4-2 of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<4||isempty(isSortedNoRepeats))
       isSortedNoRepeats=false; 
    end
    
    numEls=numel(empData);
    %Use the empirical distribution if weighting is not specified.
    if(nargin<3||isempty(w))
        w=ones(numEls,1)/numEls;
    end
    %If weighting is given, then we must figure out which of the discrete
    %samples x equals (if any).

    %If it must be sorted and repeats eliminated.
    if(isSortedNoRepeats==false)
        %First, sort the data.
        [empData,idx]=sort(empData,'ascend');
        w=w(idx);
        
        %Now, any repeats will be eliminated. However, we cannot just use
        %the unique command, because when multiple occurrences of a value
        %are combined, their weights must be added.
        empDataNew=zeros(numEls,1);
        wNew=zeros(numEls,1);
        
        empDataNew(1)=empData(1);
        wNew(1)=w(1);
        curUniqueEl=1;
        for curEl=2:numEls
            %If there is a repeat.
            if(empData(curEl)==empDataNew(curUniqueEl))
                wNew(curUniqueEl)=wNew(curUniqueEl)+w(curEl);
            else
                curUniqueEl=curUniqueEl+1;
                empDataNew(curUniqueEl)=empData(curEl);
                wNew(curUniqueEl)=w(curEl);
            end
        end
        numEls=curUniqueEl;
        empData=empDataNew(1:numEls);
        w=wNew(1:numEls);
    end

    numPoints=numel(x);
    
    val=zeros(size(x));
    for curEl=1:numPoints
        [foundVal,idx]=binSearch(empData, x(curEl));
        if(foundVal==x(curEl))
            val(curEl)=w(idx);
        end
    end
end

function [val,empData,w]=CDF(x,empData,w,isSortedNoRepeats)
%%CDF        Evaluate the cumulative distribution function of the
%            empirical distribution at desired points. Weights can be
%            provided if the samples are weighted.
%
%INPUTS:    x The point(s) at which the empirical CDF is to be evaluated.
%     empData An NX1 or 1XN vector of N scalar samples that make up the
%             empirical distribution. Values can be repeated.
%           w An optional NX1 or 1XN vector of weights associated with the
%             samples. These should be positive and sum to 1. If omitted or
%             an empty matrix is passed, the samples are assumed to be
%             uniformly weighted.
% isSortedNoRepeats If the points in empData are sorted in ascending order
%             and no points are repeated, then the algorithm is faster.
%             This indicates whether the values in empData are sorted in
%             ascending order without repeats. If this parameter is omitted
%             or an empty matrix is passed, a default value of false is
%             used.
%
%OUTPUTS: val The CDF of the empirical distribution evaluated at the
%             desired point(s).
%     empData The same as empData given on the input, except if
%             isSortedNoRepeats is false, empData has been sorted and any
%             repeats have been removed.
%           w The same as w on the input, except if isSortedNoRepeats is
%             false, then w has been sorted to reflect any change of
%             ordering in empData, and its values have been adjusted to
%             reflect the removal of any duplicate values.
%
%The empirical distribution is discussed in many introductory statistics
%textbooks, such as in Chapter 4-2 of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    numEls=numel(empData);
    if(nargin<3||isempty(w))
        w=ones(numEls,1)/numEls;
    end
    
    if(nargin<4||isempty(isSortedNoRepeats))
       isSortedNoRepeats=false; 
    end

    %If it must be sorted and repeats eliminated.
    if(isSortedNoRepeats==false)
        %First, sort the data.
        [empData,idx]=sort(empData,'ascend');
        w=w(idx);
        
        %Now, any repeats will be eliminated. However, we cannot just use
        %the unique command, because when multiple occurrences of a value
        %are combined, their weights must be added.
        empDataNew=zeros(numEls,1);
        wNew=zeros(numEls,1);
        
        empDataNew(1)=empData(1);
        wNew(1)=w(1);
        curUniqueEl=1;
        for curEl=2:numEls
            %If there is a repeat.
            if(empData(curEl)==empDataNew(curUniqueEl))
                wNew(curUniqueEl)=wNew(curUniqueEl)+w(curEl);
            else
                curUniqueEl=curUniqueEl+1;
                empDataNew(curUniqueEl)=empData(curEl);
                wNew(curUniqueEl)=w(curEl);
            end
        end
        numEls=curUniqueEl;
        empData=empDataNew(1:numEls);
        w=wNew(1:numEls);
    end
    minVal=empData(1);
    
    %Get the CDF values
    CDF=cumsum(w);
    CDF(end)=1;%Deal with any possible finite-precision issues.

    numPoints=numel(x);
    
    val=zeros(size(x));
    for curEl=1:numPoints
        if(x(curEl)<minVal)
            val(curEl)=0;
        else
            %Search for the point, or the next lowest value if the point is
            %not there.
            [~,idx]=binSearch(empData, x(curEl),1);
            val(curEl)=CDF(idx);
        end
    end
end

function [vals,empData,w]=invCDF(prob,empData,w,isSortedNoRepeats)
%%INVCDF Evaluate the inverse of the cumulative distribution function of
%        the empirical distribution. Weights can be provided if the samples
%        are weighted.
%
%INPUTS: prob The probability or probabilities (0<=prob<=1) at which the 
%             argument of the CDF is desired.
%     empData An NX1 or 1XN vector of N scalar samples that make up the
%             empirical distribution. Values can be repeated.
%           w An optional NX1 or 1XN vector of weights associated with the
%             samples. These should be positive and sum to 1. If omitted or
%             an empty matrix is passed, the samples are assumed to be
%             uniformly weighted.
% isSortedNoRepeats If the points in empData are sorted in ascending order
%             and no points are repeated, then the algorithm is faster.
%             This indicates whether the values in empData are sorted in
%             ascending order without repeats. If this parameter is omitted
%             or an empty matrix is passed, a default value of false is
%             used.
%
%OUTPUTS: val The argument(s) of the CDF that would give the probability or
%             probabilities in prob. When this function is given a
%             probability value that does not align with a value on the
%             discrete CDF, the function returns the next most probable
%             point.
%     empData The same as empData given on the input, except if
%             isSortedNoRepeats is false, empData has been sorted and any
%             repeats have been removed.
%           w The same as w on the input, except if isSortedNoRepeats is
%             false, then w has been sorted to reflect any change of
%             ordering in empData, and its values have been adjusted to
%             reflect the removal of any duplicate values.
%
%The empirical distribution is discussed in many introductory statistics
%textbooks, such as in Chapter 4-2 of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    numEls=numel(empData);
    if(nargin<3||isempty(w))
        w=ones(numEls,1)/numEls;
    end
    
    if(nargin<4||isempty(isSortedNoRepeats))
       isSortedNoRepeats=false; 
    end

    %If it is not sorted, then it must be sorted.
    if(isSortedNoRepeats==false)
        %First, sort the data.
        [empData,idx]=sort(empData,'ascend');
        w=w(idx);
        
        %Now, any repeats will be eliminated. However, we cannot just use
        %the unique command, because when multiple occurrences of a value
        %are combined, their weights must be added.
        empDataNew=zeros(numEls,1);
        wNew=zeros(numEls,1);
        
        empDataNew(1)=empData(1);
        wNew(1)=w(1);
        curUniqueEl=1;
        for curEl=2:numEls
            %If there is a repeat.
            if(empData(curEl)==empDataNew(curUniqueEl))
                wNew(curUniqueEl)=wNew(curUniqueEl)+w(curEl);
            else
                curUniqueEl=curUniqueEl+1;
                empDataNew(curUniqueEl)=empData(curEl);
                wNew(curUniqueEl)=w(curEl);
            end
        end
        numEls=curUniqueEl;
        empData=empDataNew(1:numEls);
        w=wNew(1:numEls);
    end
    
    CDF=cumsum(w);
    CDF(end)=1;%Deal with finite-precision issues.
    
    numPoints=numel(prob);
    
    vals=zeros(size(prob));
    for curEl=1:numPoints
        %Search for the point, or the next highest value if the point is
        %not there.
        [~,idx]=binSearch(CDF, prob(curEl),2);
        vals(curEl)=empData(idx);
    end
end

function [vals,empData,w]=rand(N,empData,w,isSortedNoRepeats)
%%RAND Generate empirically distributed random variables.
%
%INPUTS:    N If N is a scalar, then rand returns an NXN matrix of random
%             variables. If N=[M,N1] is a two-element row vector, then
%             rand returns an MXN1 matrix of  random variables.
%     empData An NX1 or 1XN vector of N scalar samples that make up the
%             empirical distribution. Values can be repeated.
%           w An optional NX1 or 1XN vector of weights associated with the
%             samples. These should sum to 1. If omitted or an empty matrix
%             is passed, the samples are assumed to be uniformly weighted.
% isSortedNoRepeats If the points in empData are sorted in ascending order
%             and no points are repeated, then the algorithm is faster.
%             This indicates whether the values in empData are sorted in
%             ascending order without repeats. If this parameter is omitted
%             or an empty matrix is passed, a default value of false is
%             used.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated empirical random variables.
%      empData The same as empData given on the input, except if
%              isSortedNoRepeats is false, empData has been sorted and any
%              repeats have been removed.
%            w The same as w on the input, except if isSortedNoRepeats is
%              false, then w has been sorted to reflect any change of
%              ordering in empData, and its values have been adjusted to
%              reflect the removal of any duplicate values.
%
%This is an implementation of the inverse transform algorithm of Chapter
%4.1 of [1].
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(isscalar(N))
        dims=[N, N];
    else
        dims=N;
    end

    numEls=numel(empData);
    if(nargin<3||isempty(w))
        w=ones(numEls,1)/numEls;
    end
    
    if(nargin<4||isempty(isSortedNoRepeats))
       isSortedNoRepeats=false; 
    end

    %If it is not sorted, then it must be sorted.
    if(isSortedNoRepeats==false)
        %First, sort the data.
        [empData,idx]=sort(empData,'ascend');
        w=w(idx);
        
        %Now, any repeats will be eliminated. However, we cannot just use
        %the unique command, because when multiple occurrences of a value
        %are combined, their weights must be added.
        empDataNew=zeros(numEls,1);
        wNew=zeros(numEls,1);
        
        empDataNew(1)=empData(1);
        wNew(1)=w(1);
        curUniqueEl=1;
        for curEl=2:numEls
            %If there is a repeat.
            if(empData(curEl)==empDataNew(curUniqueEl))
                wNew(curUniqueEl)=wNew(curUniqueEl)+w(curEl);
            else
                curUniqueEl=curUniqueEl+1;
                empDataNew(curUniqueEl)=empData(curEl);
                wNew(curUniqueEl)=w(curEl);
            end
        end
        numEls=curUniqueEl;
        empData=empDataNew(1:numEls);
        w=wNew(1:numEls);
    end
    
    CDF=[0;cumsum(w(1:(end-1)))];
    U=rand(dims);
    
    %Allocate space.
    vals=zeros(N);

    numPoints=numel(vals);
    for curPoint=1:numPoints
        [~,idx]=binSearch(CDF,U(curPoint),1);
        vals(curPoint)=empData(idx);
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
