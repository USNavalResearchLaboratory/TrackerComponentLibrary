function [peakVals,isMax,peakIdx]=findMaxMinPoints(data)
%%FINDMAXMINPOINTS Given a grid of points in one or more dimensions, this
%                  function finds all maxima and minima. A maximum/minimum
%                  is declared if the point in question is higher/lower
%                  than all neighboring points (including those diagonally
%                  offset). Points at the edge of the matrix are not
%                  considered to be maxima or minima.
%
%INPUTS: data An n1Xn2X...Xnk dimensional matrix.  Singleton dimensions
%             (ni=1) are ignored. Non-singleton dimensions must be > 2. to
%             consider points that are not on the edge. data must be a real
%             matrix.
%
%OUTPUTS: peakVals A numValsX1 vector of the maximum and minimum values
%                  found. If no values are found, this is an empty matrix. 
%            isMax A numValsX1 vector indicating whether each value in
%                  peakVals is a maximum or a minimum. 1 indicates maximum,
%                  0 indicated minimum.
%          peakIdx A numValsX1 vector of indices of the peaks such that
%                  data(peakIdx(i)) provides the ith peak value.
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
% surface(X,Y,Z,'EdgeColor','None')
% [peakVals,isMax,peakIdx]=findMaxMinPoints(Z)
%One sees that all 13 of the maxma and minima (-1 and 1)  that are found.
%Values on the edge are not counted for maxma or minima.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

dims=size(data);

%Get rid of singleton dimensions.
dims=dims(dims~=1);
if(isempty(dims)||isscalar(dims)&&dims==1||any(dims==2))
   peakVals=[];
   peakIdx=[];
   isMax=[];
   return;
end

numDims=length(dims);
switch(numDims)
    case 1%Find peaks for linear values. 
        peakVals=[];
        peakIdx=[];
        isMax=[];
        for idx1=2:(dims(1)-1)
            curVal=data(idx1);
            
            vals=[curVal-data(idx1-1);
                  curVal-data(idx1+1)];
            if(all(vals>0))
                peakVals=[peakVals;curVal];
                peakIdx=[peakIdx;idx1];
                isMax=[isMax;1];
            elseif(all(vals<0))
                peakVals=[peakVals;curVal];
                peakIdx=[peakIdx;idx1];
                isMax=[isMax;0];
            end
        end
    case 2
        data=reshape(data,dims(1),dims(2));
        
        peakVals=[];
        peakIdx=[];
        isMax=[];
        
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
                    peakVals=[peakVals;curVal];
                    idx=sub2ind(dims,idx1,idx2);
                    peakIdx=[peakIdx;idx];
                    isMax=[isMax;1];
                elseif(all(vals<0))
                    peakVals=[peakVals;curVal];
                    idx=sub2ind(dims,idx1,idx2);
                    peakIdx=[peakIdx;idx];
                    isMax=[isMax;0];
                end
            end
        end
    otherwise%3 or more dimensions
        maxTupleVals=dims-3;
        
        peakVals=[];
        peakIdx=[];
        isMax=[];
        
        idxList=getNextTuple(numDims);
        while(~isempty(idxList))
            %Points that will be compared to their neighbors are never the
            %edge points. Thus, the minimum value possible for an element
            %is 2 and the maximum value is one less than the number of
            %things in that dimension.
            shiftIdxList=idxList+2;
            
            curIdx=nDim2Index(dims,shiftIdxList);
            curVal=data(curIdx);
            
            %Now, we look at all neighbors of the current point. This means
            %that all points go through all combinations of +1,-1, and 0,
            %except for the all zero case. We can get the values using
            %tuples.
            maxSignVals=2*ones(numDims,1);
            
            signList=getNextTuple(zeros(numDims,1),maxSignVals);
            
            isMax=true;
            isMin=true;
            while(~isempty(signList)&&(isMax==true||isMin==true))
                temp=signList;
                temp(temp==2)=-1;

                shiftedIdx=nDim2Index(dims,shiftIdxList+temp);
                compVal=data(shiftedIdx);
                
                if(compVal>=curVal)
                   isMax=false; 
                end
                
                if(compVal<=curVal)
                   isMin=false; 
                end

                signList=getNextTuple(signList,maxSignVals);
            end
            
            if(~(isMax&&isMin))
                if(isMax)
                    peakVals=[peakVals;curVal];
                    peakIdx=[peakIdx;curIdx];
                    isMax=[isMax;1];
                elseif(isMin)
                    peakVals=[peakVals;curVal];
                    peakIdx=[peakIdx;curIdx];
                    isMax=[isMax;0];
                end
            end
            
            idxList=getNextTuple(idxList,maxTupleVals);
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
