function hashVal=hashPoint2Grid(point,cellSizes,numCellsPerDim)
%%HASHPOINT2GRID Given a numDim-dimensional point, turn it into a scalar
%                value corresponding to a position in a finite grid of
%                points. This is useful for gating. One can make the grid
%                size slightly larger than the maximum allowable
%                uncertainty region for a target before it is declared
%                lost. Measurements are then hashed to discrete points and
%                the corners of the rectangular bounding box of each target
%                is hashed into the grid. The detections in each bin of the
%                grid into which the hashed bounding box falls can then be
%                used as detections that pass a first round of gating. For
%                a very large area, this can be much more efficient than
%                trying to gate all targets with all measurements.
%
%INPUTS: point A numDimXN set of N points that one wishes to hash onto a
%               grid. The grid starts from 0, so the points should be
%               appropriately shifted in all dimensions.
%     cellSizes A numDimX1 or 1XnumDim array of the widths of the cells in
%               each dimension of the grid.
% numCellsPerDim A numDimX1 or 1XnumDim array listing the number of cells
%               in each dimension.
%
%OUTPUTS: hashVal An NX1 array of the hash values for each of the points.
%                 This is an integer starting from 1.
%
%The idea of spatial hashing arises when performing collision detection in
%computer graphics, among many other applications. A simple explanation is
%in [1].
%
%Dimension of the point is gridded to a region of size cellSizes(i). This
%can be done using fix(point(i)/cellSizes(i)). However, a maximum number of
%cells in each dimension is allowed. Thus, we use
%mod(fix(point(i)/cellSizes(i)),numCellsPerDim(i)). Note that we do not
%have to specify whether the actual grid covers a specific range or
%negative values as the modulo operation just wraps values around in a
%particular dimension. Next, we want to handle multiple dimensions. For two
%dimensions, we can compute the index as above for the x-axis, and then use
%the same formula as for the x axis but times numCellsPerDim(1) for the y
%axis. This avoid collisions on the grid between axes. Similarly, for a
%z-axis, the multiplicative offset for higher dimensions is the product of
%numCellsPerDim for the lower dimensions. 1 is added to the final hash
%value so that it starts from 1 instead of 0.
%
%Note that the mod operation can produce an undesired behaviour if one is
%not careful. For example, suppose that point contains scalar values from
%minVal to maxVal. One might want to put the values in point into numBins
%bins between minVal and maxVal. One might initially try
% points=points-min(points);
% delta=max(points);
% cellSizes=delta/numBins;
% hashVal=hashPoint2Grid(points,cellSizes,numBins);
%However, this would alias the uppermost point into the same bin as the
%lowermost point, because delta*numBins is the non-inclusive upper bound on
%the binned values, it is not the center of a bin. One solution to this is
%to use
% hashVal=hashPoint2Grid(points,cellSizes+eps(cellSizes),numBins);
%so that the final bin ends just slightly after the final point. Another
%solution is to just add one one extra bin after the end fo the range.
%
%REFERENCES:
%[1] E. J. Hastings, J. Mesit, and R. K. Guha, "Optimization of large-
%    scale, real-time simulation by spatial hashing," in Proceedings of the
%    2005 Summer Computer Simulation Conference, vol. 37, no. 4,
%    Philadelphia, PA, 24-28 Jul. 2005.
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=length(cellSizes);

hashVal=mod(fix(point(1,:)/cellSizes(1)),numCellsPerDim(1));

prevCellProd=numCellsPerDim(1);
for curDim=2:numDim
    hashVal=hashVal+mod(fix(point(curDim,:)/cellSizes(curDim)),numCellsPerDim(curDim))*prevCellProd;

    prevCellProd=prevCellProd*numCellsPerDim(curDim);
end


%Make it so that the has values start from 1 and can be used as indices in
%Matlab.
hashVal=hashVal'+1;

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
