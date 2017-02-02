function L=LorentzTransMatrix(vVec,vectorType)
%%LORENTZTRANSMATRIX  Get the matrix for the Lorentz transformation of
%                     time-space vectors between two inertial coordinate
%                     systems. The Lorentz transform embodies the
%                     coordinate system transformation inherent to special
%                     relativity theory and plays a role in transformaing
%                     space-time vectors. The matrix can be obtained for
%                     either the case where an interval is specified as
%                     z=[c*Delta t;Delta x']', where c is the speed of
%                     light, Delta t is a time interval, and x is a
%                     3D distance vector or for the case where an interval
%                     is specified as z=[1i*c*Delta t;Delta x']' where
%                     1i=sqrt(-1).
%
%INPUTS:vvec The 3X1 velocity vector of the origin of the second inertial
%            coordinate system measured with respect to the first inertial
%            coordinate system. Note that
%            norm(vVec)<=Constants.speedOfLight.
% vectorType A string specifying the type of Lorentz transform matrix to
%            obtain. This can be
%            'Real' The default if this parameter is omitted. All of the
%                   entries in the Lorentz transforms matrix are real and
%                   the interval being transformed is assumed to have the
%                   real form z=[c*Delta t;Delta x']'
%         'Complex' Some of the entries in Lorentz matrix can be complex
%                   and the interval being transformed is assumed to have a
%                   complex time component. That is,
%                   z=[1i*c*Delta t;Delta x']' The advantage of this system
%                   is that L is now a 4D rotation matrix.
%
%OUTPUTS: L If vectorType='Real', L is the symmetric Lorentz
%           transformation matrix such that if z=[c*Delta t;Delta x']' is a
%           4-dimensional vector of duration Delta t in time and Delta x in
%           3D space (c is the speed of light), in the first inertial
%           coordinate system, then L*z is the equivalent vector in the
%           second inertial coordinate system. If vectorType='Complex',
%           then L is the Hermitian matrix such that is
%           z=[1i*c*Delta t;Delta x']' in the first inertial coordinate
%           system, then L*z is the equivalent vector in the second
%           inertial coordinate system. Note that regardless of the
%           vectorType, inv(L)=LorentzTransMatrix(-v), because the
%           direction of the velocity of the first coordinate system
%           measured in the second coordinate system is -v. Note that
%           det(L)=1. If vectorType='Complex', the L*transpose(L)=eye(4).
%
%A derivation of the real form of the Lorentz transformation can be found
%in Chapter 1.2.2 of [1]. The complex form of the Lorentz transformation
%and its relation to rotation matrices is provided in Appendix D of [2].
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
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Check the size of v.
if(size(vVec,1)~=3||size(vVec,2)~=1)
   error('The vector v is not 3X1.') 
end

if(nargin<2)
    vectorType='Real';
end

c=Constants.speedOfLight;

%The min operation deals with precision problems if norm(vVec) is just a
%smidge over c that would normally cause complex numbers to pop up. It is
%still assumed that within numerical precision, norm(vVec)<=c.
v=min(norm(vVec),c);

%Equation 1.20 in Ludyk
gamma=1/sqrt(1-v^2/c^2);
L11=gamma;

switch(vectorType)
    case 'Real'
        L21=-(gamma/c)*vVec;

        %Deal with the possible limit when v->c.
        L21(isnan(L21))=0;
        
        temp=(gamma-1)*(vVec*vVec')/v^2;

        %If v is very small, then uVec will contain garbage and temp should be
        %replaced with zeros (the asymptotic value of v->0.
        if(any(~isfinite(temp)))
            temp=zeros(3,3);
        end

        L22=eye(3)+temp;

        L=[L11, L21';
           L21, L22];
    case 'Complex'
        beta=vVec/c;
        betaHat=vVec/norm(vVec);
        
        L22=eye(3)+(gamma-1)*(betaHat*betaHat');
        
        %Equation D10 in Shuster, with the rows and columns reversed.
        L=[L11, -1i*gamma*beta';
           1i*gamma*beta, L22];
    otherwise
        error('Invalid vectorType given');
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
