function [idxPairs,ySolIdx,xValPairs,yValPairs]=findValueBrackets(xVals,yVals,numSolPerPoint,yDesired,maxNumSol,oneSided)
%%FINDVALUEBRACKETS Given a 1D function evaluated at multiple points,
%       sorted in order of the dependent variable, find all pairs of points
%       that bracket the value yDesired. This can be useful for bracketing
%       the roots of nonlinear functions that have been evaluated on a
%       grid. A local optimization algorithm can then be used on the
%       bracketed domain.
%
%INPUTS: xVals A length numPts vector of the dependent variable sorted in
%              increasing or decreasing order. If one does not wish to
%              return xValPairs, then this can be an empty matrix.
%        yVals If numSolPerPoint is omitted, this is a length numPts vector
%              of the independent values corresponding to the entries in
%              xVals. If numSolPerPoint is given, then there can be more
%              than one solution for each dependent variable. In such an
%              instance, this is the maxSolXnumPts matrix of the solutions
%              for each dependent variable. The order of the solutions at
%              each dependent variable matters, because bracketing is only
%              considered if yVals(i,k) and y(i,k+1) bracket the solution.
%              The number of solutions for each step need not be the same.
%              This multiple solution option can be useful for bracketing
%              solutions when ray tracing, whereby multiple solutions arise
%              from multiple bounces (and solutions from one bounce should
%              not be associated with those from a subsequent bounce). Note
%              that we assume that the regions from multiple solutions do
%              not overlap. The means that a maximum of only one bracket
%              per entry in xVals is taken.
% numSolPerPoint If only one solution is present for each entry in xVals,
%              then an empty matrix should be passed of this input.
%              Otherwise, this is a length numPts vector where
%              numSolPerPoint(i) indicates that
%              yVals(1:numSolPerPoint(i),i) holds i solutions for xVals(i).
%              If all entries in yVals have the same number of solutions,
%              then a scalar value can be passed for numSolPerPoint.
%     yDesired The scalar y value that should be bracketed.
%    maxNumSol The maximum number of solutions to find. The default if
%              omitted or an empty matrix is passed is 400.
%     oneSided If yDesired exactly equals y(i), then two brackets are
%              possible [i-1,i] and [i,i+1]. If this is true, then only a
%              single bracket will be returned, the [i-1,i] bracket. The
%              default if omitted or an empty matrix is passed is true.
%
%OUTPUTS: idxPairs A 2XnumSol matrix where idxPairs(1,i) and idxPairs(2,i)
%              are the indices of the values in xVals that bracket a value
%              of yDesired.
%      ySolIdx A numSolX1 vector selecting which of the solutions in yVals
%              are associated in the bracketing at each time. If
%              numSolPerPoint is an empty matrix, then this is just an
%              empty matrix, because there is no ambiguity (one could just
%              use a vector of all 1s).
%    xValPairs A 2XnumSol matrix of the values in xVals that correspond to
%              the pairs in idxPairs or an empty matrix is xVals was empty.
%    yValPairs A 2XnumSol matrix of the values in yVals that bracket the
%              solution. When multiple solutions are given per point, this
%              chooses the specific solution chosen.
%
%EXAMPLE:
%Here, we consider a sine wave and bracket where it equals 0.8. We display
%the sine wave and the bracketed points.
% numPts=1e3;
% x=linspace(0,5*pi,numPts);
% y=sin(x);
% yDesired=0.8;
% [~,~,xValPairs,yValPairs]=findValueBrackets(x,y,[],yDesired);
% 
% figure(1)
% clf
% hold on
% plot(x,y,'linewidth',2)
% scatter(xValPairs(:),yValPairs(:),300,'.')
%
%March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(oneSided))
    oneSided=true;
end

if(isempty(numSolPerPoint))
    numPts=length(yVals);
else
    numPts=size(yVals,2);
end

if(numPts<2)
    error('A minimum of two points must be given.');
end

if(nargin<5||isempty(maxNumSol))
    maxNumSol=200;
end

hasX=~isempty(xVals);

idxPairs=zeros(2,maxNumSol);
numBracketsFound=0;
if(isempty(numSolPerPoint))
    ySolIdx=[];
    for curPoint=2:numPts
       prevY=yVals(curPoint-1);
       curY=yVals(curPoint);

       if(oneSided&&curPoint<numPts)
           testVal=(curY<yDesired&&prevY>=yDesired)||(curY>yDesired&&prevY<=yDesired);
       else
           testVal=(curY<=yDesired&&prevY>=yDesired)||(curY>=yDesired&&prevY<=yDesired);
       end
       
       if(testVal)
            %The first bracketing case is
            %yVals(curPoint-1)>=yDesired>=yVals(curPoint).
            %The second bracketing case is
            %yVals(curPoint-1)<=yDesired<=yVals(curPoint).
            numBracketsFound=numBracketsFound+1;
            idxPairs(:,numBracketsFound)=[curPoint-1;curPoint];
           
            if(numBracketsFound==maxNumSol)
                break;
            end
        end
    end
else
    if(isscalar(numSolPerPoint))
        numSolPerPoint=numSolPerPoint*ones(numPts,1);
    end
    
    ySolIdx=ones(maxNumSol,1);

    for curPoint=2:numPts
        prevY=yVals(1:numSolPerPoint(curPoint-1),curPoint-1);
        curY=yVals(1:numSolPerPoint(curPoint),curPoint);

        for curSol=1:min(numSolPerPoint(curPoint-1),numSolPerPoint(curPoint))
            if(oneSided&&curPoint<numPts)
                testVal=(curY(curSol)<yDesired&&prevY(curSol)>=yDesired)||(curY(curSol)>yDesired&&prevY(curSol)<=yDesired);
            else
                testVal=(curY(curSol)<=yDesired&&prevY(curSol)>=yDesired)||(curY(curSol)>=yDesired&&prevY(curSol)<=yDesired);
            end

            if(testVal)
                %The first bracketing case is
                %yVals(curPoint-1)>=yDesired>=yVals(curPoint).
                %The second bracketing case is
                %yVals(curPoint-1)<=yDesired<=yVals(curPoint).
                numBracketsFound=numBracketsFound+1;
                idxPairs(:,numBracketsFound)=[curPoint-1;curPoint];
                ySolIdx(numBracketsFound)=curSol;
                
                %We only expect one solution per successful bracketing.
                break;
            end
        end
        if(numBracketsFound==maxNumSol)
            break;
        end
    end
    
    ySolIdx=ySolIdx(1:numBracketsFound);
end

%Shrink to fit.
idxPairs=idxPairs(:,1:numBracketsFound);

if(hasX&&nargout>2)
    xValPairs=zeros(2,numBracketsFound);
    for curBracket=1:numBracketsFound
        xValPairs(1,curBracket)=xVals(idxPairs(1,curBracket));
        xValPairs(2,curBracket)=xVals(idxPairs(2,curBracket));
    end
else
    xValPairs=[];
end

if(nargout>3)
    if(isempty(ySolIdx))
        yValPairs=zeros(2,numBracketsFound);
        for curBracket=1:numBracketsFound
            yValPairs(1,curBracket)=yVals(idxPairs(1,curBracket));
            yValPairs(2,curBracket)=yVals(idxPairs(2,curBracket));
        end
    else
        yValPairs=zeros(2,numBracketsFound);
        for curBracket=1:numBracketsFound
            yValPairs(1,curBracket)=yVals(ySolIdx(curBracket),idxPairs(1,curBracket));
            yValPairs(2,curBracket)=yVals(ySolIdx(curBracket),idxPairs(2,curBracket));
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
