function xTimePos=LorentzTrans(xTimePos,v,vectorType,c)
%%LORENTZTRANS Perform a Lorentz transformation of a time-space vector
%     between two inertial coordinate systems. The Lorentz transform
%     embodies the coordinate system transformation inherent to special
%     relativity theory and plays a role in transforming space-time
%     vectors. The matrix can be obtained for the case where an interval is
%     specified as x=[c*Delta t;Delta x']' (the real case), where c is the
%     speed of light in a vacuum, Delta t is a time interval, and x is a 3D
%     distance vector, or for the case where x=[Delta t;Delta x']' (the
%     real, asymmetric case), or for the case where an interval is
%     specified as x=[1i*c*Delta t;Delta x']' where 1i=sqrt(-1).
%
%INPUTS: xTimePos The 4XN [time offset;position offset] values relative to
%            a common event, as viewed by the reference inertial coordinate
%            system. The exact format depends on the vectorType input.
%          v The 3XN velocity vectors of the origin of the second inertial
%            coordinate system measured with respect to the first inertial
%            coordinate system. Note that norm(vVec)<c, where c is the
%            speed of light in a vacuum. If all fo the vectors are the
%            same, then a single 3X1 vector can be passed.
% vectorType A string specifying the type of Lorentz transform matrix to
%            obtain. This can be
%            'Real' (The default if this parameter is omitted or an empty
%                   matrix is passed). All of the entries in the Lorentz
%                   transform matrix are real and the interval being
%                   transformed is assumed to have the real form
%                   xTimePos=[c*Delta t;Delta x']'
%            'RealAsymmetric' All of the entries in the Lorentz transform
%                   matrix are real and the interval being transformed is
%                   assumed to have the real form
%                   xTimePos=[Delta t;Delta x']'
%            'Complex' Some of the entries in Lorentz matrix can be complex
%                   and the interval being transformed is assumed to have a
%                   complex time component. That is,
%                   xTimePos=[1i*c*Delta t;Delta x']' The advantage of this
%                   system is that L is now a 4D rotation matrix.
%          c The speed of light in a vacuum. If omitted or an empty matrix
%            is passed,t he default of c=Constants.speedOfLight, which has
%            united of meters per second, is used.
%
%OUTPUTS: xTimePos A 4XN set of [time;position offset] values relative to
%                  the moving (non-reference) inertial coordinate system.
%
%This function multiplied xTimePos by the outputs of LorentzTransMatrix.
%
%A derivation of the real (symmetric and asymmetric) forms of the Lorentz
%transformation can be found in Chapter 1.2.2 of [1]. The complex form of
%the Lorentz transformation and its relation to rotation matrices is
%provided in Appendix D of [2].
%
%Note that events that are simultaneous in one inertial coordinate system
%are not necessarily simultaneous in another.
%
%In accelerating systems, space-time intervals are generally obtained by
%integrating over a series of instantaneous Lorentz transformations.
%
%REFERENCES:
%[1] G. Ludyk, Einstein in Matrix Form: Exact Derivation of the Theory of
%    Special and General Relativity without Tensors. Heidelberg: Springer,
%    2013.
%[2] M. D. Shuster, "A survey of attitude representations," The Journal of
%    Astronautical Sciences, vol. 41, no. 4, pp. 439-517, Oct.-Dec. 1993.
%
%February 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<4||isempty(c))
        c=Constants.speedOfLight;
    end
    
    if(nargin<3||isempty(vectorType))
        vectorType='Real';
    end

    if(size(v,2)==1)
        xTimePos=LorentzTransMatrix(v,vectorType,c)*xTimePos;
    else
        N=size(v,2);
        for k=1:N
            xTimePos(:,k)=LorentzTransMatrix(v(:,k),vectorType,c)*xTimePos(:,k);
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
