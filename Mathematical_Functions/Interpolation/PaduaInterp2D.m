function [interpVals,C]=PaduaInterp2D(deg,f,xEval,C)
%%PADUAINTERP2D Perform 2D interpolation on the square (-1,1)X(-1,1) using
%               Padua interpolation. Given a function handle, f, the
%               function is evaluated on a select number of points and then
%               values at other points can be interpolated. Alternatively,
%               this function can compute a matrix C so that subsequent
%               calls for interpolation of the same function do not have to
%               evaluate f.
%
%INPUTS: deg The integer polynomial degree of interpolation (deg>0).
%          f A handle to a function whose outputs should be interpolated. f
%            must take a 2XN matrix of points with elements ranging from -1
%            to 1 and return an NX1 or 1XN vector of values. If C is
%            provided and is not an empty matrix, then an empty matrix can
%            be passed for f.
%      xEval A 2XnumPoints set of numPoints points where the function
%            should be interpolated. If one just wants to compute C, then
%            this parameter can be omitted or an empty matrix passed.
%          C An optional parameter. If interpolation is meant to be
%            performed with the same function f over more than one call of
%            this function, then this is the output matrix passed back in
%            to simplify computation. If C is not available, then this
%            parameter can be omitted or an empty matrix passed.
%
%OUTPUTS: interpVals A 2XnumPoints set of the values of the function
%            interpolated at the points in xEval. If not points were
%            provided (only C is desired), then this will be an empty
%            matrix.
%          C A matrix that can be passed back to this function to
%            accelerate its execution in subsequent calls. It is Equation
%            13 from [1] computed as in Section 3.2.
%
%This function implements the algorithm described in [1].
%
%EXAMPLE:
%Here, we interpolate two points of a nonlinear function.
% f=@(x)sin(pi*(x(1,:)+2*x(2,:)));
% xEval=[1/4,   0;
%   -8/10, 0];
% deg=23;
% interpVals=PaduaInterp2D(deg,f,xEval)
% trueVals=f(xEval)
%One can see that the interpolated values are very close to the true
%values.
%
%REFERENCES:
%[1] M. Caliari, S. de Marchi, and M. Vianello, "Algorithm 886: Padua2D -
%   Lagrange interpolation at Padua points on bivariate domains," ACM
%   Transactions on Mathematical Software, vol. 35, no. 3, Oct. 2008,
%   Article 21.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(C))
    numPoints=(deg+1)*(deg+2)/2;%This equals binomial(deg+2,2).

    xi=zeros(2,numPoints);
    w=zeros(numPoints,1);

    curPoint=1;
    for j=0:(deg+1)
        cosVal2=-cos(j*pi/(deg+1));

        for k=mod(j,2):2:deg
            %These are putting in the values from Equations 1 and 2.
            xi(1,curPoint)=-cos(k*pi/deg);
            xi(2,curPoint)=cosVal2;
            %Equation 7 for an interior point. The follow two conditional
            %statements deal with vertex and edge points.
            w(curPoint)=2/(deg*(deg+1));

            %j=0,k=0
            if(k==0||k==deg)
                w(curPoint)=w(curPoint)/2;
            end

            if(j==0||j==deg+1)
                w(curPoint)=w(curPoint)/2;
            end

            curPoint=curPoint+1;
        end
    end

    fVals=f(xi);

    C=findC(deg,fVals,w);
    %We will transform C into the C0 in Equation 13
    C(end,1)=C(end,1)/2;
    for curRow=2:(deg+1)
        C(curRow,(deg-curRow+3):(deg+1))=0;
    end
end

if(nargin>2&&~isempty(xEval))
    numxPoints=size(xEval,2);
    
    interpVals=zeros(numxPoints,1);
    
    for k=1:numxPoints
        T1=ChebyshevPoly(xEval(1,k),deg);
        T1(2:end)=T1(2:end)*sqrt(2);%Normalize
        T2=ChebyshevPoly(xEval(2,k),deg);
        T2(2:end)=T2(2:end)*sqrt(2);%Normalize
        %Equation 14
        interpVals(k)=T1'*C*T2;
    end
else
    interpVals=[];
end
end


function C=findC(deg,fVals,wVals)
%Section 3.2

%Compute G(f) as in Equation 16
G=zeros(deg+1,deg+2);
curPoint=1;
for s=0:(deg+1)
    for r=0:deg
        if(mod(s+r,2)==0) 
            G(r+1,s+1)=fVals(curPoint)*wVals(curPoint);
            curPoint=curPoint+1;
        end
    end
end

%Compute P2
P2=zeros(deg+1,deg+2);
for i=0:(deg+1)
    tau=-cos(i*pi/(deg+1));
    %The sqrt(2) term comes from the assumption that deg!=0 and from the
    %deginition og the normalized Chebyshev polynomial given between
    %Equations 6 and 7.
    P2(:,i+1)=sqrt(2)*ChebyshevPoly(tau,deg);
    P2(1,i+1)=P2(1,i+1)/sqrt(2);
end

%Compute P1
P1=zeros(deg+1,deg+1);
for i=0:deg
    tau=-cos(i*pi/deg);
    %The sqrt(2) term comes from the assumption that deg!=0 and from the
    %deginition og the normalized Chebyshev polynomial given between
    %Equations 6 and 7.
    P1(:,i+1)=sqrt(2)*ChebyshevPoly(tau,deg);
    P1(1,i+1)=P1(1,i+1)/sqrt(2);
end

%Equation 17
C=P1*G*P2';

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
