function val=aPoly(x,t,numDim)
%%APOLY     The drift function for a linear continuous-time motion model
%           with a given order in a specified number of Cartesian
%           dimensions. The order of the linear filter, that is the number
%           of moments of position, does not ned to be explicitely
%           specified.
%
%INPUTS:    x The xDimXN state vector of N targets in the order of
%             [position;velocity;acceleration;etc] for however many
%             derivatives of position there are.
%           t An unused time component so that aLinear can be used with
%             Runge-Kutta methods that expect the function to take two
%             parameters.
%      numDim The number of dimensions of the simulation problem. If
%             the numDim parameter is omitted, then numDim=3 (3D motion) is
%             assumed. The dimensionality of the state must be an integer
%             multiple of numDim.
%
%OUTPUTS: val The time-derivative of the N state vectors under the linear
%             motion model.
%
%The drift function corresponds to the state transition given in
%discrete-time by the function FPolyKal.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<3)
        numDim=3; 
    end
    
    numTar=size(x,2);

    val=[x((numDim+1):end,:);zeros(numDim,numTar)];
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
