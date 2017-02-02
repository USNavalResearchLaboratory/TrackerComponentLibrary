function val=rectContainedInRect(rectMin1,rectMax1,rectMin2,rectMax2)
%RECTCONTAINEDINRECT Determine whether an axis-aligned  rectangle (or a
%                    more general axis-aligned hyperrectangle if the number
%                    of dimensions is not two) is completely engulfed by
%                    another rectangle. This does not indicate which
%                    rectangle is contained in the other.
%
%INPUTS: rectMin1 A kX1 vector of the lower bounds of each of the
%                 dimensions of the first k-dimensional hyperrectangle.
%        rectMax1 A kX1 vector of the upper bounds of the first
%                 hyperectangle.
%        rectMin2 A kX1 vector of the lower bounds of the second
%                 hyperrectangle.
%        rectMax2 A kX1 vector of the upper bounds of the second
%                 hyperrectangle.
%
%OUTPUTS: val A boolean value that is true if one hyperrectangle is
%             completely contained in the other and is false otherwise.
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

k=length(rectMin1);

for curIdx=1:k
    if(rectMin1(curIdx)<rectMin2(curIdx))
        val=false;
        return
    end
    
    if(rectMax1(curIdx)>rectMax2(curIdx))
       val=false;
       return
    end
end

val=true;
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
