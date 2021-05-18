function wrapVals=wrapRange(vals,minBound,maxBound,mirrorWrap)
%%WRAPRANGE Bound values to a specific range, wrapping around if the value
%           is outside the range. The parameter mirrorWrap determines how
%           the wrapping is performed. For example, if mirrorWrap=false,
%           minBound=-pi and maxBound=pi, then a value that is some eps
%           above pi will be wrapped to a value some eps above -pi.
%           Similarly, a value some eps below -pi is wrapped to a value
%           some eps below pi. On the other hand, if mirrorWrap==true, then
%           leaving the region of [minBound,maxBound] goes the other way.
%           For example, if mirrorWrap=true, minBound=-pi/2, and
%           maxBound=pi/2, then a value that is some eps above pi/2 will be
%           mapped to a point some eps below pi/2 and a value some eps
%           below -pi/2 will be mapped to a point some eps above -pi/2.
%
%INPUTS: vals A vector or matrix of real values that should be wrapped
%             to the range minBound->maxBound.
%    minBound The lower scalar bound of the output parameters.
%    maxBound A value > minBound that is the upper bound to the allowable
%             range of output values.
%  mirrorWrap A value that determines whether a values outside the bounds
%             maps back into the valid region offset to the bound that it
%             passed or offset to the other boundary. If this parameter is
%             omitted, then mirrorWrap=false is assumed and the function
%             behaves similar to a shifted modulo function.
%
%OUTPUTS: wrapVals The set of vals wrapped such that
%                  minBound<=wrapVals<maxBound.
%
%If mirrowWrap=false, then the values are shifted to correspond to a region
%from zero to maxBound-minBound and the mod function is used to wrap vals
%to that region. The results are then shifted back by minBound to get
%values in the desired region.
%
%If mirrorWrap=true, then if the bounds are truly -pi/2->pi/2, the results
%are analogous to what one gets by taking the sine and then inverse sine
%of the angular value. On the other hand, if the bounds are still pi apart
%but shifted, then one just shifts the value back to zero evaluates
%asin(sin(*)) of the value and shifts everything to the original position.
%In the more general case, there will be both shifting and
%scaling, but a combination of sine and inverse since will still produce
%the desired result
%
%When considering mirrorWrap=false, the wrapping is not as simple as one
%might think. One would expect that values within the interval of minBound
%to maxBound would not change. However, one cannot just use a shifted
%modulo operation, because mod(vals-minBound,maxBound-minBound) can return
%zero if vals-minBound is numerically equal to maxBound-minBound. For
%example, consider maxBound=pi, minBound=-pi and vals=pi-eps(pi). Thus, we
%first do a check to see whether the value is in the primary range. If so,
%then the value is just returned directly.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<4)
       mirrorWrap=false; 
    end

    if(maxBound<=minBound)
        error('The maximum bound must be less than the minimum bound')
    end
    
    if(mirrorWrap==false)
        wrapVals=mod(vals-minBound,maxBound-minBound)+minBound;
        
        %This takes care of numerical precision issues where something very
        %near the upper edge of the primary interval will be wrapped into
        %the bottom of the interval, because vals-minBound is numerically
        %equal to maxBound-minBound.
        sel=vals<maxBound&vals>=minBound;
        wrapVals(sel)=vals(sel);
    else%Otherwise
        spreadVal=maxBound-minBound;
        
        %First, shift the values so that the region of wrapping is centered
        %symmetrically around zero.
        vals=vals-minBound-spreadVal/2;
        
        %Now, scale the values so that  the bounds reduce to (+/-)pi/2,
        %meaning that the spreadVal becomes pi.
        vals=(pi/spreadVal)*vals;
        
        %Now, take the tangent and two-quadrant inverse tangent to wrap
        %everything back where it should be. Note that Matlab properly
        %handles infinite values in the inverse tangent.
        vals=asin(sin(vals));
        
        %Scale everything back to the original size
        vals=(spreadVal/pi)*vals;
        
        %Shift the origin back to where it should be.
        wrapVals=vals+minBound+spreadVal/2;
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
