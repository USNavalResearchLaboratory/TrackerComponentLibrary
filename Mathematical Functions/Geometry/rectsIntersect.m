function val=rectsIntersect(rectMin1,rectMax1,rectMin2,rectMax2, containedIntersects)
%%RECTSINTERSECT  Determine whether two axis-aligned rectangles (or more
%                 general axis-aligned hyperrectangles if the number of
%                 dimensions is not 2) overlap at all. The handling of the
%                 case where one rectangle is completely contained in
%                 another without intersecting boundaries can be specified.
%
%INPUTS: rectMin1 A kX1 vector of the lower bounds of each of the
%                 dimensions of the first k-dimensional hyperrectangle.
%        rectMax1 A kX1 vector of the upper bounds of the first
%                 hyperectangle.
%        rectMin2 A kX1 vector of the lower bounds of the second
%                 hyperrectangle.
%        rectMax2 A kX1 vector of the upper bounds of the second
%                 hyperrectangle.
%containedIntersects An optional boolean parameter. If false, then two
%                 hyperrectangles where one is completely contained in
%                 another are not considered to intersect. Otherwise, they
%                 are considered to intersect.
%
%OUTPUTS: val A boolean value that is true if the hyperrectangles intersect
%             and false otherwise.
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<5)
        containedIntersects=false;
    end

    if(sum(rectMax1<rectMin2)>0||sum(rectMax2<rectMin1)>0)
        val=false;
    else
        val=true;
    end
    
    %If they do not intersect, check whether one is contained in another,
    %if containedIntersects is true.
    if(val==false&&containedIntersects)
        val=rectContainedInRect(rectMin1,rectMax1,rectMin2,rectMax2);
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
