function [peakVals,isMax,peakIdx]=findMaxMinPoints(data,allowEdges,get1DLedges)
%%FINDMAXMINPOINTS Given a grid of points in one or more dimensions, this
%    function finds all maxima and minima. A maximum/minimum is declared if
%    the point in question is higher/lower than all neighboring points
%    (including those diagonally offset), with slightly different behavior
%    in 1D when handling flat regions depending on get1DLedges. Points at
%    the edge of the matrix can be considered if desired, when considering
%    aliasing for neighbor determination.
%
%INPUTS: data An n1Xn2X...Xnk dimensional matrix.  Singleton dimensions
%             (ni=1) are ignored. Non-singleton dimensions must be > 2. to
%             consider points that are not on the edge. data must be a real
%             matrix.
%  allowEdges One has the option of finding detections at the edges (with
%             the "neighbors" being the points aliased at the other edge).
%             The default if omitted or an empty matrix is passed is false
%             (the edge points wil not be considered local maxima nor
%             minima).
% get1DLedges When data is 1D, if this input is true, then when there are
%             multiple equal values, the start of the region can be
%             declared a local maximum or a local minimum. Otherwise, the
%             criteria for local maxima and minima are just as specified
%             above. The default if omitted or an empty matrix is passed is
%             true.
%
%OUTPUTS: peakVals A numValsX1 vector of the maximum and minimum values
%                  found. If no values are found, this is an empty matrix. 
%            isMax A numValsX1 vector indicating whether each value in
%                  peakVals is a maximum or a minimum. 1 indicates maximum,
%                  0 indicated minimum.
%          peakIdx A numValsX1 vector of indices of the peaks such that
%                  data(peakIdx(i)) provides the ith peak value.
%
%With noisy or poorly discretized functions, one is going to want a more
%robust peak finding algorithm.
%
%EXAMPLE:
%This function is useful for finding maxima and minima of cost functions.
%Here, we find all maxima and minima of a surface.
% numPoints=100;
% pts=linspace(0,2*pi,numPoints);
% [X,Y]=meshgrid(pts,pts);
% Z=sin(2*X+3*Y)+cos(X-2*Y);
% %Display the surface whose peaks are desired.
% figure(1)
% clf
% hold on
% surface(X,Y,Z,'EdgeColor','None')
% [peakVals,isMax,peakIdx]=findMaxMinPoints(Z);
% scatter3(X(peakIdx),Y(peakIdx),Z(peakIdx),200,'.k')
%One sees that all 13 of the maxima and minima (-1 and 1)  that are found.
%Values on the edge are not counted for maxma nor minima.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(get1DLedges))
    get1DLedges=true;
end

if(nargin<2||isempty(allowEdges))
    allowEdges=false;
end

dims=size(data);

%Get rid of singleton dimensions.
dims=dims(dims~=1);
if(isempty(dims)||isscalar(dims)&&dims==1||any(dims==2))
   peakVals=[];
   peakIdx=[];
   isMax=[];
   return;
end

%A crude upper bound is just all the elements.
maxDetects=numel(data);

peakVals=zeros(maxDetects,1);
peakIdx=zeros(maxDetects,1);
isMax=zeros(maxDetects,1);

numDetects=0;
numDims=length(dims);
switch(numDims)
    case 1%Find peaks for linear values. 
        if(allowEdges==false)
            for idx1=2:(dims(1)-1)
                curVal=data(idx1);
    
                vals=[curVal-data(idx1-1);
                      curVal-data(idx1+1)];
                if(all(vals>0)||(get1DLedges&&vals(1)>0&&vals(2)==0))
                    numDetects=numDetects+1;
                    peakVals(numDetects)=curVal;
                    peakIdx(numDetects)=idx1;
                    isMax(numDetects)=1;
                elseif(all(vals<0)||(get1DLedges&&vals(1)<0&&vals(2)==0))
                    numDetects=numDetects+1;
                    peakVals(numDetects)=curVal;
                    peakIdx(numDetects)=idx1;
                    isMax(numDetects)=0;
                end
            end
        else
            for idx1=1:dims(1)
                curVal=data(idx1);
                
                idxPrev=mod(idx1-1-1,dims(1))+1;
                idxNext=mod(idx1+1-1,dims(1))+1;
    
                vals=[curVal-data(idxPrev);
                      curVal-data(idxNext)];
                if(all(vals>0)||(get1DLedges&&vals(1)>0&&vals(2)==0))
                    numDetects=numDetects+1;
                    peakVals(numDetects)=curVal;
                    peakIdx(numDetects)=idx1;
                    isMax(numDetects)=1;
                elseif(all(vals<0)||(get1DLedges&&vals(1)<0&&vals(2)==0))
                    numDetects=numDetects+1;
                    peakVals(numDetects)=curVal;
                    peakIdx(numDetects)=idx1;
                    isMax(numDetects)=0;
                end
            end
        end
    case 2
        data=reshape(data,dims(1),dims(2));

        if(allowEdges==false)
            for idx1=2:(dims(1)-1)
                for idx2=2:(dims(2)-1)
                    curVal=data(idx1,idx2);

                    vals=[curVal-data(idx1,idx2-1);
                          curVal-data(idx1-1,idx2-1);
                          curVal-data(idx1+1,idx2-1);
                          curVal-data(idx1,idx2+1);
                          curVal-data(idx1-1,idx2+1);
                          curVal-data(idx1+1,idx2+1);
                          curVal-data(idx1-1,idx2);
                          curVal-data(idx1+1,idx2)];

                    %If we have found a maximum or a minimum
                    if(all(vals>0))
                        numDetects=numDetects+1;
                        peakVals(numDetects)=curVal;
                        idx=sub2ind(dims,idx1,idx2);
                        peakIdx(numDetects)=idx;
                        isMax(numDetects)=1;
                    elseif(all(vals<0))
                        numDetects=numDetects+1;
                        peakVals(numDetects)=curVal;
                        idx=sub2ind(dims,idx1,idx2);
                        peakIdx(numDetects)=idx;
                        isMax(numDetects)=0;
                    end
                end
            end
        else
            for idx1=1:dims(1)
                idx1Prev=mod(idx1-1-1,dims(1))+1;
                idx1Next=mod(idx1+1-1,dims(1))+1;
                for idx2=1:dims(2)
                    curVal=data(idx1,idx2);

                    idx2Prev=mod(idx2-1-1,dims(2))+1;
                    idx2Next=mod(idx2+1-1,dims(2))+1;

                    vals=[curVal-data(idx1,idx2Prev);
                          curVal-data(idx1Prev,idx2Prev);
                          curVal-data(idx1Next,idx2Prev);
                          curVal-data(idx1,idx2Next);
                          curVal-data(idx1Prev,idx2Next);
                          curVal-data(idx1Next,idx2Next);
                          curVal-data(idx1Prev,idx2);
                          curVal-data(idx1Next,idx2)];

                    %If we have found a maximum or a minimum
                    if(all(vals>0))
                        numDetects=numDetects+1;
                        peakVals(numDetects)=curVal;
                        idx=sub2ind(dims,idx1,idx2);
                        peakIdx(numDetects)=idx;
                        isMax(numDetects)=1;
                    elseif(all(vals<0))
                        numDetects=numDetects+1;
                        peakVals(numDetects)=curVal;
                        idx=sub2ind(dims,idx1,idx2);
                        peakIdx(numDetects)=idx;
                        isMax(numDetects)=0;
                    end
                end
            end

        end
    otherwise%3 or more dimensions
        if(allowEdges==false)
            maxTupleVals=dims-3;
        else
            maxTupleVals=dims-1;
        end
        
        idxList=getNextTuple(numDims);
        while(~isempty(idxList))
            %Points that will be compared to their neighbors are never the
            %edge points if allowEdges is false. Thus, the minimum value
            %possible for an element is 2 and the maximum value is one less
            %than the number of things in that dimension.
            if(allowEdges==false)
                shiftIdxList=idxList+2;
            else
                shiftIdxList=idxList+1;
            end
            
            curIdx=nDim2Index(dims,shiftIdxList);
            curVal=data(curIdx);
            
            %Now, we look at all neighbors of the current point. This means
            %that all points go through all combinations of +1,-1, and 0,
            %except for the all zero case. We can get the values using
            %tuples.
            maxSignVals=2*ones(numDims,1);
            
            signList=getNextTuple(zeros(numDims,1),maxSignVals);
            
            isMaxCur=true;
            isMinCur=true;
            while(~isempty(signList)&&(isMaxCur==true||isMinCur==true))
                temp=signList;
                temp(temp==2)=-1;

                shiftIdxListCur=shiftIdxList+temp;

                if(allowEdges)
                    for curDim=1:length(dims)
                        shiftIdxListCur(curDim)=mod(shiftIdxListCur(curDim)-1,dims(curDim))+1;
                    end
                end

                shiftedIdx=nDim2Index(dims,shiftIdxListCur);
                compVal=data(shiftedIdx);
                
                if(compVal>=curVal)
                   isMaxCur=false; 
                end
                
                if(compVal<=curVal)
                   isMinCur=false; 
                end

                signList=getNextTuple(signList,maxSignVals);
            end
            
            if(~(isMaxCur&&isMinCur))
                if(isMaxCur)
                    numDetects=numDetects+1;
                    peakVals(numDetects)=curVal;
                    peakIdx(numDetects)=curIdx;
                    isMax(numDetects)=1;
                elseif(isMinCur)
                    numDetects=numDetects+1;
                    peakVals(numDetects)=curVal;
                    peakIdx(numDetects)=curIdx;
                    isMax(numDetects)=0;
                end
            end
            
            idxList=getNextTuple(idxList,maxTupleVals);
        end
end

%Size to fit.
peakVals=peakVals(1:numDetects);
peakIdx=peakIdx(1:numDetects);
isMax=isMax(1:numDetects);

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
