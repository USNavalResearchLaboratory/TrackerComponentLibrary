function [interpPolyA,interpPolyC]=RKInterpPolys(x0,t0,x1,t1,f,order,solutionChoice,k,interpMainOrder)
%%RKINTERPPOLYS Get polynomials to interpolate between two steps of a
%               Runge-Kutta method or of a Rosenbrock method. Depending on
%               the method used, either a low-computational complexity
%               explicit solution from a  direct formula is used, or a
%               Hermite interpolating polynomial is obtained by taking a
%               sufficient number of small steps between the two
%               Runge-Kutta solutions.
%
%INPUTS: (x0,t0), (x1,t1) The values of the function at time t0 and t1.
%               Presumably, (x0,t0) were inputs to the RungeKStep function
%               and with a step duration of deltaT=t1-t0, and x1 was the
%               main output. t0 and t1 are scalar. x0 and x1 can be xDimX1
%               column vectors.
%            f  f(xVal,curT) returns the derivative of xVal with
%               respect to time taken at time curT.
% (order,solutionChoice) The pair of order and solutionChoice identify
%               which algorithm was used to perform the Runge-Kutta step
%               and thus indicate what the best interpolation routine
%               should be. If a Rosenbrock method is used, set
%               solutionChoice=-1;
%            k  The values of the derivatives f evaluated at various
%               points as determined by the algorithm used for the
%               RungeKStep function. If a Rosenbrock method is used, this
%               parameter can be omitted or an empty matrix can be passed.
% interpMainOrder An optional parameter indicating whether the interpolation
%               should be the same as the main order of the Runge-Kutta
%               pair or whether it should be the same as the secondary
%               order. If omitted, the default value is true.
%
%OUTPUTS: interpPolysA,interpPolysC An xDimX(polyOrder+1) and
%               xDimX(polyOrder) matrices that hold the interpolating
%               coefficients (A) and  control points (C) in Newton's form
%               for interpolation in each dimension. The argument of the
%               polynomials, sigma, is a value from 0->1 indicating how far
%               between the endpoints the value is. That is, to evaluate
%               the function at a time t0<=t<=t1, then compute
%               z=(t-t0)/deltaT
%               and use
%               evalPolyNewton(interpPolysA(i,:),interpPolysC(i,:),z)
%               to find the value in the ith dimensions of x.
%               Though the coefficients can be converted to a standard
%               polynomial form using the convertPolynomialForm function, a
%               loss of precision might occur.
%
%The method of interpolation varies depending on the Runge-Kutta formula
%used. For orders 1->3, cubic Hermite interpolation is performed using the
%two points. For order four, if the default algorithm (solutionChoice=0) is
%used, a fourth-order interpolant given in [1] is used. For order 5, if the
%default formula is used, then a fourth-order interpolant given in [2] is
%used. For other orders/ solutionChoice values, a bunch of intermediate
%Runge-Kutta steps are taken and the HermiteInterpPoly function is used to
%find a matching Hermite interpolation polynomial whose order is order.
%
%When performing Hermite interpolation, it is checked whether the last
%value in k is f evaluated at x1,t1 (the isFSAL parameter), if so, then an
%extra evaluation of f is avoided. 
%
%REFERENCES:
%[1] M. K. Horn, "Fourth- and fifth-order, scaled Runge-Kutta algorithms
%    for treating dense output," SIAM Journal on Numerical Analysis, vol.
%    20, no. 3, pp. 558-568, 1983.
%[2] J. R. Dormand, Numerical Methods for Differential Equations. Boca
%    Raton: CRC Press, 1996.
%
%February 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<9)
    interpMainOrder=true;
end

deltaT=t1-t0;
xDim=size(x0,1);

%If a Rosenbrock method was used, then just find the Hermite interpolating
%polynomial between the two points.
if(solutionChoice==-1)
        interpPolyA=zeros(xDim,4);
        interpPolyC=zeros(xDim,3);
        dxdt0=f(x0,t0);
        dxdt1=f(x1,t1);
        for curDim=1:xDim
            %The output should be scaled over a range from 0->1, so the
            %derivatives must be scaled according to the change of
            %variables. The substitution is sigma=(t-t0)/deltaT, so
            %dSigma=(1/deltaT)*dt Thus, the dxdt terms have to be
            %multiplied by deltaT.
            [interpPolyA(curDim,:),interpPolyC(curDim,:)]=HermiteInterpPoly([0;1],[[x0(curDim);deltaT*dxdt0(curDim)],[x1(curDim);deltaT*dxdt1(curDim)]]);
        end
        return;
end

%Otherwise, assume that a Runge-Kutta method is used.
dxdt0=k(:,1);
switch(order)
    case 1
        [interpPolyA,interpPolyC]=RK2PtHermiteInterp(x0,deltaT,dxdt0,x1,t1,order,solutionChoice,k,f);
    case 2%Find the Hermite interpolating polynomial from the two points.
        [interpPolyA,interpPolyC]=RK2PtHermiteInterp(x0,deltaT,dxdt0,x1,t1,order,solutionChoice,k,f);
    case 3%Find the Hermite interpolating polynomial from the two points.
        [interpPolyA,interpPolyC]=RK2PtHermiteInterp(x0,deltaT,dxdt0,x1,t1,order,solutionChoice,k,f);
    case 4
        switch(solutionChoice)
            case 0%Fehlberg's RK4(5)6 formula
                if(interpMainOrder==true)
                    %The fourth-order interpolation routine is taken from
                    %[1].
                    B= [-311/360,       269/108,        -301/120,   1,  0;
                        0,              0,              0,          0,  0;
                        14848/4275,     -4096/513,      7168/1425,  0,  0;
                        -371293/75240,  199927/22572,   -28561/8360,0,  0;
                        42/25,          -3,             57/50,      0,  0;
                        -102/55,        40/11,          -96/55,     0,  0;
                        5/2,            -4,             3/2,        0,  0];

                    %Horn's method requires an extra function evaluation.
                    k(:,7)=f(x0+deltaT*((1/6)*k(:,1)+(1/6)*k(:,5)+(2/3)*k(:,6)),t1);

                    %Build the polynomials for each dimension in terms of
                    %an offset sigma from 0->1 that takes the time from t0
                    %to t1:
                    interpPolyA=(deltaT*k*B);
                    interpPolyA(:,end)=x0;

                    %The polynomial is given in a non-Hermite form, so the
                    %control points are all zero.
                    interpPolyC=zeros(xDim,size(interpPolyA,2)-1);
                    %Also, the ordering of the coefficients in interpPolyA
                    %must be reversed.
                    interpPolyA=fliplr(interpPolyA);
                else
                    %Use a computationally-expensive method that works for
                    %any order to get the fifth-order coefficients.
                    [interpPolyA,interpPolyC]=RKStepsForHermite(x0,t0,x1,t1,f,order,solutionChoice,interpMainOrder,k); 
                end
            otherwise
                %Use a computationally-expensive method that works for any
                %order.
                [interpPolyA,interpPolyC]=RKStepsForHermite(x0,t0,x1,t1,f,order,solutionChoice,interpMainOrder,k);
        end
    case 5%Use an explicit formula, if available
        switch(solutionChoice)
            case 0%RK5(4)FM
                if(interpMainOrder==false)
                    %The fourth-order interpolation routine is taken from
                    %Chapter 6.5 of [2].
                    B=zeros(7,5);
                    B(1,:)=[435,    -1184,  1098,   -384, 0]*(-1/384);
                    B(2,:)=0;
                    B(3,:)=[6,      -14,    9,      0,    0]*(500/1113);
                    B(4,:)=[9,      -16,    6,      0,    0]*(-125/192);
                    B(5,:)=[35,     -64,    26,     0,    0]*(729/6784);
                    B(6,:)=[15,     -28,    12,     0,    0]*(-11/84);
                    B(7,:)=[5,      -8,     3,      0,    0]*(1/2);
                else
                    %The fifth-order interpolation routine is taken from
                    %Chapter 6.5 of [2].

                    %There are two more rows in the A matrix.
                    a8=[5207/48000, 0, 92/795, -79/960, 53217/848000, -11/300, 4/125, 0, 0];
                    a9=[613/6144, 0, 125/318, -125/3072, 8019/108544, -11/192, 1/32, 0, 0];
                    c8=1/5;
                    c9=1/2;
                    %Expand the k vector by two columns corresponding to
                    %the values in the two extra rows of the A matrix.
                    k=[k,zeros(xDim,2)];
                    k(:,8)=f(x0+deltaT*(k*a8'),t0+c8*deltaT);
                    k(:,9)=f(x0+deltaT*(k*a9'),t0+c9*deltaT);

                    B=zeros(9,6);
                    B(1,:)=[696,    -2439,  3104,   -1710,  384,    0]*(1/384);
                    B(2,:)=0;
                    B(3,:)=[24,     -51,    32,     -6,     0,  0]*(-500/1113);
                    B(4,:)=[24,     -51,    32,     -6,     0,  0]*(-125/192);
                    B(5,:)=[24,     -51,    32,     -6,     0,  0]*(2187/6784);
                    B(6,:)=[24,     -51,    32,     -6,     0,  0]*(-11/84);
                    B(7,:)=[32,     -63,    38,     -7,     0,  0]*(1/8);
                    B(8,:)=[0,      1,      -2,     1,      0,  0]*(125/24);
                    B(9,:)=[3,      -7,     5,      -1,     0,  0]*(16/3);
                end
                %Build the polynomials for each dimension in terms of an
                %offset sigma from 0->1 that takes the time from t0 to t1:
                interpPolyA=(deltaT*k*B);
                interpPolyA(:,end)=x0;
                
                %The polynomial is given in a non-Hermite form, so the
                %control points are all zero.
                interpPolyC=zeros(xDim,size(interpPolyA,2)-1);
                %Also, the ordering of the coefficients in interpPolyA must
                %be reversed.
                interpPolyA=fliplr(interpPolyA);
            otherwise
                %Use a computationally-expensive method that works for any
                %order.
                [interpPolyA,interpPolyC]=RKStepsForHermite(x0,t0,x1,t1,f,order,solutionChoice,interpMainOrder,k);
        end
    otherwise
    %If no explicit interpolation formula is available, then just take
    %small Runge-Kutta steps between the endpoints and find the Hermite
    %interpolating polynomial for each dimension of the state. This is a
    %computationally-expensive solution.
    [interpPolyA,interpPolyC]=RKStepsForHermite(x0,t0,x1,t1,f,order,solutionChoice,interpMainOrder,k);
end

end

function [interpPolyA,interpPolyC]=RK2PtHermiteInterp(x0,deltaT,dxdt0,x1,t1,order,solutionChoice,k,f)
%Find the Hermite interpolating polynomial from the two points.
    xDim=length(x0);
    interpPolyA=zeros(xDim,4);
    interpPolyC=zeros(xDim,3);

    [~,isFSAL]=RungeKStep(order,solutionChoice);
    if(isFSAL)
        dxdt1=k(:,end);
    else
        dxdt1=f(x1,t1);
    end
    for curDim=1:xDim
        %The output should be scaled over a range from 0->1, so the
        %derivatives must be scaled according to the change of
        %variables. The substitution is sigma=(t-t0)/deltaT, so
        %dSigma=(1/deltaT)*dt Thus, the dxdt terms have to be
        %multiplied by deltaT.
        [interpPolyA(curDim,:),interpPolyC(curDim,:)]=HermiteInterpPoly([0;1],[[x0(curDim);deltaT*dxdt0(curDim)],[x1(curDim);deltaT*dxdt1(curDim)]]);
    end
end


function [interpPolyA,interpPolyC]=RKStepsForHermite(x0,t0,x1,t1,f,order,solutionChoice,interpMainOrder,k)
%%RKSTEPSFORHERMITE Given estimates x0 and x1 at times t0 and t1, create a
%                   set of interpolating polynomials for each of the
%                   dimensions of x across the timespan by performing a
%                   series of fixed-length Runge-Kutta steps to get enough
%                   estimates of x and its derivative at different times t
%                   so that the interpolating polynomial has order order.
%                   This is an inefficient method of obtaining an
%                   interpolating polynomial compared to techniques that
%                   have been crafted for specified RUnge-Kutta formulae.
%
%INPUTS: (x0,t0) The xDimX1 state and the 1X1 time at the initial time.
%        (x1,t1) The xDimX1 state and the 1X1 time at the final step.
%           f    f(x,t) returns the derivative of x with respect to time
%                taken at time t.
% order, solutionChoice Two values (see the function RungeKStep for
%                details) that specify the order of the polynomial
%                approximation and the algorithm used for the Runge-Kutta
%                step. These are also used to determine whether the last
%                element in k is the derivative used in the next step.
%interpMainOrder A parameter indicating whether the interpolation
%               should be the same as the main order of the Runge-Kutta
%               pair or whether it should be the same as the secondary
%               order.
%            k  The values of the derivatives f evaluated at various
%               points as determined by the algorithm used for the
%               RungeKStep function. This is used so as to avoid needless
%               extra evaluation of f.
%
%OUTPUTS: interpPolysA,interpPolysC An xDimX(polyOrder+1) and
%               xDimX(polyOrder) matrices that hold the interpolating
%               coefficients (A) and  control points (C) in Newton's form
%               for interpolation in each dimension in the range (t0,t1). 
%
%The function just performs an appropriate number of uniform Runge-Kutta
%steps between t0 and t1 so that the function HermiteInterpPoly can be used
%to get an interpolating polynomial.
%
%February 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

xDim=size(x0,1);
deltaT=t1-t0;

[orders,isFSAL]=RungeKStep(order,solutionChoice);
if(interpMainOrder==true)
    %To get a Hermite interpolating polynomials (one for each dimensions)
    %of the given order matching points and derivatives, one needs
    %ceil((1+order)/2) points.
    numPoints=ceil((1+order)/2);
else
    %Interpolate to match the subsidiary interpolating order.
    numPoints=ceil((1+orders(2))/2);
end

deltaTStep=deltaT/(numPoints-1);
tVals=linspace(t0,t1,numPoints);
x=zeros(xDim,numPoints);
dxdt=zeros(xDim,numPoints);

x(:,1)=x0;
x(:,end)=x1;
fCur=k(:,1);
dxdt(:,1)=fCur;

%Fill in the intermediate values using short Runge-Kutta steps.
curT=t0;
for curStep=2:(numPoints-1)
    [x(:,curStep),~,kCur]=RungeKStep(x(:,curStep-1),curT,f,deltaTStep,fCur,order,solutionChoice);
    curT=curT+deltaTStep;
    if(isFSAL)
        fCur=kCur(:,end);
    else
        fCur=f(x(:,curStep),curT);
    end
    dxdt(:,curStep)=fCur;
end

if(isFSAL)
    dxdt(:,end)=k(:,end);
else
    dxdt(:,end)=f(x1,t1);
end

%Allocate space for the interpolating polynomials.
interpPolyA=zeros(xDim,2*numPoints);
interpPolyC=zeros(xDim,2*numPoints-1);
for curDim=1:xDim
    %The scaling is from a change of variables so that the range of
    %interpolation parameters goes from 0->1 and not from t0 to t1.
    [interpPolyA(curDim,:),interpPolyC(curDim,:)]=HermiteInterpPoly((tVals-t0)/deltaT,[x(curDim,:);deltaT*dxdt(curDim,:)]);
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
