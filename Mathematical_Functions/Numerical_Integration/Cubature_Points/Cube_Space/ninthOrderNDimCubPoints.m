function [xi,w]=ninthOrderNDimCubPoints(numDim,algorithm)
%%NINTHORDERNDIMCUBPOINTS Generate ninth-order cubature points for
%               integration over a numDim-dimensional cube with bounds in
%               coordinates of (-1,1).
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated. numDim>1.
%     algorithm An optional parameter specifying the algorithm to be used
%               to generate the points. Possible values are:
%               0 (The default if omitted or an empty matrix is passed)
%                 Formula Cn 9-1 in [1], pg. 236,
%                 4*(numDim^4-5*numDim^3+14*numDim^2-7*numDim+3)/3 points,
%                 variant 1, numDim>=4.
%               1 (The default if omitted or an empty matrix is passed)
%                 Formula Cn 9-1 in [1], pg. 236,
%                 4*(numDim^4-5*numDim^3+14*numDim^2-7*numDim+3)/3 points,
%                 variant 1, numDim>=4.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    algorithm=0; 
end

switch(algorithm)
    case 0%Cn 9-1 in [1], pg. 236,
          %4*(numDim^4-5*numDim^3+14*numDim^2-7*numDim+3)/3 points, variant
          %1.
        if(numDim<4)
           error('This formula is for numDim>=4') 
        end
        
        V=2^numDim;
        
        %The monomial integral values:
        I2=2^numDim/(2+1);
        I4=2^numDim/(4+1);
        I6=2^numDim/(6+1);
        I8=2^numDim/(8+1);
        I22=2^numDim/((2+1)*(2+1));
        I42=2^numDim/((4+1)*(2+1));
        I44=2^numDim/((4+1)*(4+1));
        I62=2^numDim/((6+1)*(2+1));
        I222=2^numDim/((2+1)^3);
        I422=2^numDim/((4+1)*(2+1)*(2+1));
        I2222=2^numDim/((2+1)^4);
    
        v=sqrt((I2*I8-I4*I6+sqrt(I2^2*I8^2+4*I4^3*I8+4*I2*I6^3-6*I2*I4*I6*I8-3*I4^2*I6^2))/(2*(I2*I6-I4^2)));
        u=sqrt((I2*I8-I4*I6-sqrt(I2^2*I8^2+4*I4^3*I8+4*I2*I6^3-6*I2*I4*I6*I8-3*I4^2*I6^2))/(2*(I2*I6-I4^2)));
        
        F=(I62-I44)/(4*u^2*v^2*(u^2-v^2)^2);
        H=(I422-(numDim-3)*I2222)/(8*v^8);
        I=(I422-v^2*I222)/(16*(numDim-3)*u^6*(u^2-v^2));
        J=(I2222-16*u^8*I)/(16*v^8);
        E=(u^2*I22-I42)/(4*v^4*(u^2-v^2))-F*u^2/(v^2)-2*(numDim-2)*(H+(numDim-3)*J);
        D=(I42-v^2*I22)/(4*u^4*(u^2-v^2))-F*v^2/(u^2)-2*(numDim-2)*(numDim-3)*I;
        C=(u^2*I2-I4)/(2*v^2*(u^2-v^2))-2*(numDim-1)*(E+F+(numDim-2)*(H+(2/3)*(numDim-3)*J));
        B=(I4-v^2*I2)/(2*u^2*(u^2-v^2))-2*(numDim-1)*(D+F+(2/3)*(numDim-2)*(numDim-3)*I);
        A=V-2*numDim*(B+C+(numDim-1)*(D+E+2*F+(1/3)*(numDim-2)*(2*H+(numDim-3)*(I+J))));
        
        %The formula for the number of points as given in the text is
        %incorrect.
        numPoints=4*(numDim^4-5*numDim^3+14*numDim^2-7*numDim+3)/3;

        xi=zeros(numDim,numPoints);
        w=zeros(numPoints,1);
        
        %xi(:,1) is all zeros.
        w(1)=A;
        curStartIdx=2;
        
        xiCur=fullSymPerms([u;zeros(numDim-1,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStartIdx:(curStartIdx+num2Add-1))=xiCur;
        w(curStartIdx:(curStartIdx+num2Add-1))=B;
        curStartIdx=curStartIdx+num2Add;
        
        xiCur=fullSymPerms([v;zeros(numDim-1,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStartIdx:(curStartIdx+num2Add-1))=xiCur;
        w(curStartIdx:(curStartIdx+num2Add-1))=C;
        curStartIdx=curStartIdx+num2Add;
 
        xiCur=fullSymPerms([u;u;zeros(numDim-2,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStartIdx:(curStartIdx+num2Add-1))=xiCur;
        w(curStartIdx:(curStartIdx+num2Add-1))=D;
        curStartIdx=curStartIdx+num2Add;
        
        xiCur=fullSymPerms([v;v;zeros(numDim-2,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStartIdx:(curStartIdx+num2Add-1))=xiCur;
        w(curStartIdx:(curStartIdx+num2Add-1))=E;
        curStartIdx=curStartIdx+num2Add;
        
        xiCur=fullSymPerms([u;v;zeros(numDim-2,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStartIdx:(curStartIdx+num2Add-1))=xiCur;
        w(curStartIdx:(curStartIdx+num2Add-1))=F;
        curStartIdx=curStartIdx+num2Add;
        
        xiCur=fullSymPerms([v;v;v;zeros(numDim-3,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStartIdx:(curStartIdx+num2Add-1))=xiCur;
        w(curStartIdx:(curStartIdx+num2Add-1))=H;
        curStartIdx=curStartIdx+num2Add;
        
        xiCur=fullSymPerms([u;u;u;u;zeros(numDim-4,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStartIdx:(curStartIdx+num2Add-1))=xiCur;
        w(curStartIdx:(curStartIdx+num2Add-1))=I;
        curStartIdx=curStartIdx+num2Add;
        
        xiCur=fullSymPerms([v;v;v;v;zeros(numDim-4,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStartIdx:(curStartIdx+num2Add-1))=xiCur;
        w(curStartIdx:(curStartIdx+num2Add-1))=J;
    case 1%Cn 9-1 in [1], pg. 236,
          %4*(numDim^4-5*numDim^3+14*numDim^2-7*numDim+3)/3 points, variant
          %2
                  if(numDim<4)
           error('This formula is for numDim>=4') 
        end
        
        V=2^numDim;
        
        %The monomial integral values:
        I2=2^numDim/(2+1);
        I4=2^numDim/(4+1);
        I6=2^numDim/(6+1);
        I8=2^numDim/(8+1);
        I22=2^numDim/((2+1)*(2+1));
        I42=2^numDim/((4+1)*(2+1));
        I44=2^numDim/((4+1)*(4+1));
        I62=2^numDim/((6+1)*(2+1));
        I222=2^numDim/((2+1)^3);
        I422=2^numDim/((4+1)*(2+1)*(2+1));
        I2222=2^numDim/((2+1)^4);
    
        u=sqrt((I2*I8-I4*I6+sqrt(I2^2*I8^2+4*I4^3*I8+4*I2*I6^3-6*I2*I4*I6*I8-3*I4^2*I6^2))/(2*(I2*I6-I4^2)));
        v=sqrt((I2*I8-I4*I6-sqrt(I2^2*I8^2+4*I4^3*I8+4*I2*I6^3-6*I2*I4*I6*I8-3*I4^2*I6^2))/(2*(I2*I6-I4^2)));
        
        F=(I62-I44)/(4*u^2*v^2*(u^2-v^2)^2);
        H=(I422-(numDim-3)*I2222)/(8*v^8);
        I=(I422-v^2*I222)/(16*(numDim-3)*u^6*(u^2-v^2));
        J=(I2222-16*u^8*I)/(16*v^8);
        E=(u^2*I22-I42)/(4*v^4*(u^2-v^2))-F*u^2/(v^2)-2*(numDim-2)*(H+(numDim-3)*J);
        D=(I42-v^2*I22)/(4*u^4*(u^2-v^2))-F*v^2/(u^2)-2*(numDim-2)*(numDim-3)*I;
        C=(u^2*I2-I4)/(2*v^2*(u^2-v^2))-2*(numDim-1)*(E+F+(numDim-2)*(H+(2/3)*(numDim-3)*J));
        B=(I4-v^2*I2)/(2*u^2*(u^2-v^2))-2*(numDim-1)*(D+F+(2/3)*(numDim-2)*(numDim-3)*I);
        A=V-2*numDim*(B+C+(numDim-1)*(D+E+2*F+(1/3)*(numDim-2)*(2*H+(numDim-3)*(I+J))));
        
        %The formula for the number of points as given in the text is
        %incorrect.
        numPoints=4*(numDim^4-5*numDim^3+14*numDim^2-7*numDim+3)/3;

        xi=zeros(numDim,numPoints);
        w=zeros(numPoints,1);
        
        %xi(:,1) is all zeros.
        w(1)=A;
        curStartIdx=2;
        
        xiCur=fullSymPerms([u;zeros(numDim-1,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStartIdx:(curStartIdx+num2Add-1))=xiCur;
        w(curStartIdx:(curStartIdx+num2Add-1))=B;
        curStartIdx=curStartIdx+num2Add;
        
        xiCur=fullSymPerms([v;zeros(numDim-1,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStartIdx:(curStartIdx+num2Add-1))=xiCur;
        w(curStartIdx:(curStartIdx+num2Add-1))=C;
        curStartIdx=curStartIdx+num2Add;
 
        xiCur=fullSymPerms([u;u;zeros(numDim-2,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStartIdx:(curStartIdx+num2Add-1))=xiCur;
        w(curStartIdx:(curStartIdx+num2Add-1))=D;
        curStartIdx=curStartIdx+num2Add;
        
        xiCur=fullSymPerms([v;v;zeros(numDim-2,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStartIdx:(curStartIdx+num2Add-1))=xiCur;
        w(curStartIdx:(curStartIdx+num2Add-1))=E;
        curStartIdx=curStartIdx+num2Add;
        
        xiCur=fullSymPerms([u;v;zeros(numDim-2,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStartIdx:(curStartIdx+num2Add-1))=xiCur;
        w(curStartIdx:(curStartIdx+num2Add-1))=F;
        curStartIdx=curStartIdx+num2Add;
        
        xiCur=fullSymPerms([v;v;v;zeros(numDim-3,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStartIdx:(curStartIdx+num2Add-1))=xiCur;
        w(curStartIdx:(curStartIdx+num2Add-1))=H;
        curStartIdx=curStartIdx+num2Add;
        
        xiCur=fullSymPerms([u;u;u;u;zeros(numDim-4,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStartIdx:(curStartIdx+num2Add-1))=xiCur;
        w(curStartIdx:(curStartIdx+num2Add-1))=I;
        curStartIdx=curStartIdx+num2Add;
        
        xiCur=fullSymPerms([v;v;v;v;zeros(numDim-4,1)]);
        num2Add=size(xiCur,2);
        xi(:,curStartIdx:(curStartIdx+num2Add-1))=xiCur;
        w(curStartIdx:(curStartIdx+num2Add-1))=J;
    otherwise
        error('Unknown algorithm specified');
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
