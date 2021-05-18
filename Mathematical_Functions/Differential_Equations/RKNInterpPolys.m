function [interpPolyA,interpPolyC]=RKNInterpPolys(x0,t0,x1,t1,df,probIsGeneral,order,solutionChoice,g,interpMainOrder)
%%RKNINTERPPOLYS Get polynomials to interpolate between two steps of a
%                general or special Runge-Kutta Nyström method. Depending
%                on the method used, either a low-computational complexity
%                explicit solution from a direct formula is used, or a
%                Hermite interpolation polynomial is obtained by taking a
%                sufficient number of small steps betwen the two
%                Runge-Kutta solutions.
%
%INPUTS: (x0,t0), (x1,t1) The values of the function at time t0 and t1.
%               Presumably, (x0,t0) were inputs to the RungeKNystroemGStep
%               or the RungeKNystroemSStep function and with a step
%               duration of deltaT=t1-t0, and x1 was the main output. t0
%               and t1 are scalar. x0 and x1 are xDimX1 column vectors,
%               where xDim is a multiple of two (values and their
%               derivatives).
%            df df(xVec,t) returns the second dervative of derivative of x
%               with respect to time taken at time t. The output is half
%               the dimensionality of x0.
% probIsGeneral Indicates whether the problem being solved is the general
%               problem, where df depends on the frist derivatives, of the
%               special problem, where df does not depend on the first
%               derivatives. The default if omitted or am empty matrix is
%               passed is true.
% order,solutionChoice  A pair of optional parameters that specify the
%               highest order of the embedded Runge-Kutta-Nyström pair
%               (special or general) to use as well as the specific
%               algorithm that was used to get x1 and t1. Details are given
%               in the comments to the RungeKNystroemGStep function, for
%               the general problem, and in the RungeKNystroemSStep
%               function for the special problem. If omitted or empty
%               matrices are passed, the default order of 5 is used and the
%               default solutionChoice of 0 is used.
%             g The values of the derivatives df evaluated at various
%               points as determined by the algorithm used for the
%               RungeKNystroemGStep or the RungeKNystroemSStep function.
% interpMainOrder An optional parameter indicating whether the
%               interpolation should be the same as the main order of the
%               Runge-Kutta pair or whether it should be the same as the
%               secondary order. If omitted, the default value is false.
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
%two points. For orders four, five and six, some of the special methods
%have explicit interpolation solutions. Otherwise, a bunch of intermediate
%Runge-Kutta steps are taken and the HermiteInterpPoly function is used to
%find a matching Hermite interpolation polynomial whose order is order.
%
%When performing Hermite interpolation, it is checked whether the last
%value in g is df evaluated at x1,t1 (the isFSAL parameter), if so, then an
%extra evaluation of df is avoided. 
%
%The explicit algorithm used for (order, solutionChoice) pairs of (4,0),
%(5,0), (6,0) and (6,1) of the special problem are given in Table 1.4 of
%[1] and in [2].
%
%REFERENCES:
%[1] J. R. Dormand, Numerical Methods for Differential Equations.
%     Raton: CRC Press, 1996.
%[2] J. R. Dormand and P. J. Prince, "Runge-Kutta-Nystrom triples,"
%    Computers and Mathematics with Applications, vol. 13, no. 12, pp.
%    937-949, 1987.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<10)
    interpMainOrder=false;
end

deltaT=t1-t0;
xDim=length(x0);
xPosDim=xDim/2;

d2xdt20=g(:,1);

switch(order)
    case 1
        [interpPolyA,interpPolyC]=RKN2PtHermiteInterp(x0,deltaT,d2xdt20,x1,t1,probIsGeneral,order,solutionChoice,g);
    case 2%Find the Hermite interpolating polynomial from the two points.
        [interpPolyA,interpPolyC]=RKN2PtHermiteInterp(x0,deltaT,d2xdt20,x1,t1,probIsGeneral,order,solutionChoice,g);
    case 3%Find the Hermite interpolating polynomial from the two points.
        [interpPolyA,interpPolyC]=RKN2PtHermiteInterp(x0,deltaT,d2xdt20,x1,t1,probIsGeneral,order,solutionChoice,g);
    case 4
        switch(solutionChoice)
            case 0
                if(probIsGeneral)
                    [interpPolyA,interpPolyC]=RKNStepsForHermite(x0,t0,x1,t1,df,probIsGeneral,order,solutionChoice,interpMainOrder,g);
                else
                    %RKN4(3)4FM, the interpolation routine is from [2] and
                    %is order 3.
                    if(interpMainOrder==false)
                        %Interpolant for the values
                        BBar=zeros(4,6);
                        BBar(1,:)=[4,   -13, 15,    -7, 0,  0]*(-1/14);
                        BBar(2,:)=[6,   -17, 14,    0,  0,  0]*(8/81);
                        BBar(3,:)=[12,  -25, 10,    0,  0,  0]*(-25/567);
                        BBar(4,:)=[12,  -19, 7,     0,  0,  0]*(1/54);

                        %Interpolant for the first derivative
                        B=zeros(4,5);
                        B(1,:)=[960,    -2484,	2135,   -658,   0]*(-1/658);
                        B(2,:)=[360,    -814,	501,    0,      0]*(32/3807);
                        B(3,:)=[288,    -482,   147,    0,      0]*(-250/26649);
                        B(4,:)=[2880,   -3692,	1047,   0,      0]*(1/2538);

                        %Build the polynomials for the values in each
                        %dimension in terms of an offset sigma from 0->1
                        %that takes the time from t0 to t1:
                        interpPolyA=(deltaT^2*g*BBar);
                        interpPolyA(:,end-1)=deltaT*x0((xPosDim+1):end);
                        interpPolyA(:,end)=x0(1:xPosDim);

                        %Next, build the polynomial for the derivatives
                        interpPolyADeriv=(deltaT*g*B);
                        interpPolyADeriv(:,end)=x0((xPosDim+1):end);

                        interpPolyA=[interpPolyA;[zeros(xPosDim,1),interpPolyADeriv]];

                        %The polynomial is given in a non-Hermite form, so
                        %the control points are all zero.
                        interpPolyC=zeros(xDim,size(interpPolyA,2)-1);
                        %Also, the ordering of the coefficients in
                        %interpPolyA must be reversed.
                        interpPolyA=fliplr(interpPolyA);
                    else
                        [interpPolyA,interpPolyC]=RKNStepsForHermite(x0,t0,x1,t1,df,probIsGeneral,order,solutionChoice,interpMainOrder,g);
                    end
                end
            otherwise
                %Use a computationally-expensive method that works for any
                %order.
                [interpPolyA,interpPolyC]=RKNStepsForHermite(x0,t0,x1,t1,df,probIsGeneral,order,solutionChoice,interpMainOrder,g);
        end
    case 5
        switch(solutionChoice)
            case 0
                if(probIsGeneral)
                    [interpPolyA,interpPolyC]=RKNStepsForHermite(x0,t0,x1,t1,df,probIsGeneral,order,solutionChoice,interpMainOrder,g);
                else
                    %RKN5(4)4, the interpolation routine is from Table 14.8
                      %of [1] and is order 4.
                    if(interpMainOrder==false)
                        %Interpolant for the values
                        BBar=zeros(5,6);
                        BBar(1,:)=[147, -464,   510,    -216,   0, 0]*(-1/432);
                        BBar(2,:)=[9,   -25,    20,     0,      0, 0]*(25/376);
                        BBar(3,:)=[3,   -5,     0,      0,      0, 0]*(-125/1512);
                        BBar(4,:)=[33,  -76,    42,     0,      0, 0]*(-2197/142128);
                        BBar(5,:)=[1,   -2,     1,      0,      0, 0]*(1/2);

                        %Interpolant for the first derivative
                        B=zeros(5,5);                 
                        B(1,:)=[0, 441,  -928,   510,   0]*(-1/432);
                        B(2,:)=[0, 27,   -50,    20,    0]*(25/376);
                        B(3,:)=[0, 9,   -10,     0,     0]*(-125/1512);
                        B(4,:)=[0, 99,  -152,    42,    0]*(2197/142128);
                        B(5,:)=[0, 3,    -4,     1,     0]*(1/2);

                        c5=1;
                        a5=[23/432, 25/94, 125/756, 2197/142128, 0];
                        g=[g,zeros(xPosDim,1)];
                        g(:,5)=df(x0(1:xPosDim)+deltaT*(g*a5'),t0+c5*deltaT);

                        %Build the polynomials for the values in each
                        %dimension in terms of an offset sigma from 0->1
                        %that takes the time from t0 to t1:
                        interpPolyA=(deltaT^2*g*BBar);
                        interpPolyA(:,end-1)=deltaT*x0((xPosDim+1):end);
                        interpPolyA(:,end)=x0(1:xPosDim);

                        %Next, build the polynomial for the derivatives
                        interpPolyADeriv=(deltaT*g*B);
                        interpPolyADeriv(:,end)=x0((xPosDim+1):end);

                        interpPolyA=[interpPolyA;[zeros(xPosDim,1),interpPolyADeriv]];

                        %The polynomial is given in a non-Hermite form, so
                        %the control points are all zero.
                        interpPolyC=zeros(xDim,size(interpPolyA,2)-1);
                        %Also, the ordering of the coefficients in
                        %interpPolyA must be reversed.
                        interpPolyA=fliplr(interpPolyA);
                    else
                        [interpPolyA,interpPolyC]=RKNStepsForHermite(x0,t0,x1,t1,f,probIsGeneral,order,solutionChoice,interpMainOrder,g);
                    end
                end
            otherwise
                %Use a computationally-expensive method that works for any
                %order.
                [interpPolyA,interpPolyC]=RKNStepsForHermite(x0,t0,x1,t1,f,probIsGeneral,order,solutionChoice,interpMainOrder,g);
        end
    case 6
        switch(solutionChoice)
            case 0
                if(probIsGeneral)
                    [interpPolyA,interpPolyC]=RKNStepsForHermite(x0,t0,x1,t1,df,probIsGeneral,order,solutionChoice,interpMainOrder,g);
                else
                    %RKN6(4)6FM, the interpolant is from [2] and has order
                    %5.
                    if(interpMainOrder==false)
                        c7=1/5;
                        a7=[-13176257/721586250, 175419/7815500, 133933/8379000, -44981/1029000, 639035/15146439, 456/336875,0];

                        %Interpolant for the values
                        BBar=zeros(7,7);
                        BBar(1,:)=[6785,    -20973, 22360,  -9092,  1071,   0,  0]*(1/2142);
                        BBar(2,:)=[2050,    -6174,  6215,   -2100,  0,      0,  0]*(-5/1044);
                        BBar(3,:)=[1010,    -3162,  3415,   -1340,  0,      0,  0]*(-5/1368);
                        BBar(4,:)=[554,     -1530,  1387,   -444,   0,      0,  0]*(-5/504);
                        BBar(5,:)=[195,     -558,   530,    -175,   0,      0,  0]*(3125/112404);
                        BBar(6,:)=[5,       -9,     5,      -1,     0,      0,  0]*(1/24);
                        BBar(7,:)=[1,       -3,     3,      -1,     0,      0,  0]*(245/24);

                        %Interpolant for the first derivative
                        B=zeros(7,6);
                        B(1,:)=[40710,  -104865,    89440,  -27276, 2142,   0]*(1/2142);
                        B(2,:)=[1230,   -3087,      2486,   -630,   0,      0]*(-25/522);
                        B(3,:)=[606,    -1581,      1366,   -402,   0,      0]*(-25/684);
                        B(4,:)=[1662,   -3825,      2774,   -666,   0,      0]*(-5/252);
                        B(5,:)=[234,    -558,       424,    -105,   0,      0]*(15625/112404);
                        B(6,:)=[30,     -45,        20,     -3,     0,      0]*(1/24);
                        B(7,:)=[2,      -5,         4,      -1,     0,      0]*(245/8);

                        %Expand the g vector by 1 columns corresponding to
                        %the values in the extra row of the A matrix.
                        g=[g,zeros(xPosDim,1)];
                        g(:,7)=df(x0(1:xPosDim)+deltaT*(g*a7'),t0+c7*deltaT);

                        %Build the polynomials for the values in each
                        %dimension in terms of an offset sigma from 0->1
                        %that takes the time from t0 to t1:
                        interpPolyA=(deltaT^2*g*BBar);
                        interpPolyA(:,end-1)=deltaT*x0((xPosDim+1):end);
                        interpPolyA(:,end)=x0(1:xPosDim);

                        %Next, build the polynomial for the derivatives
                        interpPolyADeriv=(deltaT*g*B);
                        interpPolyADeriv(:,end)=x0((xPosDim+1):end);

                        interpPolyA=[interpPolyA;[zeros(xPosDim,1),interpPolyADeriv]];

                        %The polynomial is given in a non-Hermite form, so
                        %the control points are all zero.
                        interpPolyC=zeros(xDim,size(interpPolyA,2)-1);
                        %Also, the ordering of the coefficients in
                        %interpPolyA must be reversed.
                        interpPolyA=fliplr(interpPolyA);
                    else
                        [interpPolyA,interpPolyC]=RKNStepsForHermite(x0,t0,x1,t1,f,probIsGeneral,order,solutionChoice,interpMainOrder,g);
                    end
                end
            case 1
                if(probIsGeneral)
                    [interpPolyA,interpPolyC]=RKNStepsForHermite(x0,t0,x1,t1,df,probIsGeneral,order,solutionChoice,interpMainOrder,g);
                else
                    if(interpMainOrder==false)
                        %RKN6(4)6FD, the interpolant is from [2] and has
                        %order 5.
                        R=sqrt(8501);

                        BBar=zeros(6,7);
                        BBar(1,:)=[900,                     -3819,                      6386,                       -5244,                      2106,   0,  0]*(1/4212);
                        BBar(2,:)=0;
                        BBar(3,:)=[1800*(5860823-152228*R), -6*(4929647204-156109769*R),(22190560391-1109665151*R), 18*(81356461+25954829*R),   0,      0,  0]/22529243880;
                        BBar(4,:)=[1800*(5860823+152228*R), -6*(4929647204+156109769*R),(22190560391+1109665151*R), 18*(81356461-25954829*R),   0,      0,  0]/22529243880;
                        BBar(5,:)=[225,                     -651,                       620,                        -195,                       0,      0,  0]*(-200/17901);
                        BBar(6,:)=[300,                     -823,                       757,                        -234,                       0,      0,  0]*(1/220);

                        B=zeros(6,6);
                        B(1,:)=[5400,                       -19095,                         25544,                          -15732,                     4212,   0]*(1/4212);
                        BBar(2,:)=0;
                        B(3,:)=[5400*(5860823-152228*R),    -15*(4929647204-156109769*R),   2*(22190560391-1109665151*R),   27*(81356461+25954829*R),   0,      0]*(1/11264621940);
                        B(4,:)=[5400*(5860823+152228*R),    -15*(4929647204+156109769*R),   2*(22190560391+1109665151*R),   27*(81356461-25954829*R),   0,      0]*(1/11264621940);
                        B(5,:)=[270,                        -651,                           496,                            -117,                       0,      0]*(-1000/17901);
                        B(6,:)=[1800,                       -4115,                          3028,                           -702,                       0,      0]*(1/220);

                        %Build the polynomials for the values in each
                        %dimension in terms of an offset sigma from 0->1
                        %that takes the time from t0 to t1:
                        interpPolyA=(deltaT^2*g*BBar);
                        interpPolyA(:,end-1)=deltaT*x0((xPosDim+1):end);
                        interpPolyA(:,end)=x0(1:xPosDim);

                        %Next, build the polynomial for the derivatives
                        interpPolyADeriv=(deltaT*g*B);
                        interpPolyADeriv(:,end)=x0((xPosDim+1):end);

                        interpPolyA=[interpPolyA;[zeros(xPosDim,1),interpPolyADeriv]];

                        %The polynomial is given in a non-Hermite form, so
                        %the control points are all zero.
                        interpPolyC=zeros(xDim,size(interpPolyA,2)-1);
                        %Also, the ordering of the coefficients in
                        %interpPolyA must be reversed.
                        interpPolyA=fliplr(interpPolyA);
                    else
                        [interpPolyA,interpPolyC]=RKNStepsForHermite(x0,t0,x1,t1,f,probIsGeneral,order,solutionChoice,interpMainOrder,g);
                    end
                end
            otherwise
                %Use a computationally-expensive method that works for any
                %order.
                [interpPolyA,interpPolyC]=RKNStepsForHermite(x0,t0,x1,t1,df,probIsGeneral,order,solutionChoice,interpMainOrder,g);
        end
    otherwise
    %If no explicit interpolation formula is available, then just take
    %small Runge-Kutta steps between the endpoints and find the Hermite
    %interpolating polynomial for each dimension of the state. This is a
    %computationally-expensive solution.
    [interpPolyA,interpPolyC]=RKNStepsForHermite(x0,t0,x1,t1,df,probIsGeneral,order,solutionChoice,interpMainOrder,g);
end
end

function [interpPolyA,interpPolyC]=RKN2PtHermiteInterp(x0,deltaT,d2xdt20,x1,t1,probIsGeneral,order,solutionChoice,g)
%Find the Hermite interpolating polynomial from the two points.
    xDim=length(x0);
    xPosDim=xDim/2;
    interpPolyA=zeros(xDim,4);
    interpPolyC=zeros(xDim,3);

    if(probIsGeneral)
        [~,isFSAL]=RungeKNystroemGStep(order,solutionChoice);
    else
        [~,isFSAL]=RungeKNystroemSStep(order,solutionChoice);
    end

    if(isFSAL)
        d2xdt21=g(:,end);
    else
        if(probIsGeneral)
            d2xdt21=f(x1,t1);
        else
            d2xdt21=f(x1(1:xPosDim),t1);
        end
    end

    for curDim=1:xPosDim
        %The output should be scaled over a range from 0->1, so the
        %derivatives must be scaled according to the change of
        %variables. The substitution is sigma=(t-t0)/deltaT, so
        %dSigma=(1/deltaT)*dt Thus, the dxdt terms have to be
        %multiplied by deltaT.

        %The derivatives for the first xPosDim terms are just the
        %second set of terms.
        [interpPolyA(curDim,:),interpPolyC(curDim,:)]=HermiteInterpPoly([0;1],[[x0(curDim);deltaT*x0(xPosDim+curDim)],[x1(curDim);deltaT*x1(xPosDim+curDim)]]);
    end


    for curDim=(xPosDim+1):xDim
        %The derivatives of the second half of the terms are the second
        %derivative terms.
        [interpPolyA(curDim,:),interpPolyC(curDim,:)]=HermiteInterpPoly([0;1],[[x0(curDim);deltaT*d2xdt20(curDim-xPosDim)],[x1(curDim);deltaT*d2xdt21(curDim-xPosDim)]]);
    end
end


function [interpPolyA,interpPolyC]=RKNStepsForHermite(x0,t0,x1,t1,df,probIsGeneral,order,solutionChoice,interpMainOrder,g)
%%RKNSTEPSFORHERMITE Given estimates x0 and x1 at times t0 and t1, create a
%                   set of interpolating polynomials for each of the
%                   dimensions of x across the timespan by performing a
%                   series of fixed-length Runge-Kutta Nyström steps to get
%                   enough
%                   estimates of x and its derivative at different times t
%                   so that the interpolating polynomial has order order.
%                   This is an inefficient method of obtaining an
%                   interpolating polynomial compared to techniques that
%                   have been crafted for specified RUnge-Kutta formulae.
%
%INPUTS: (x0,t0) The xDimX1 state and the 1X1 time at the initial time.
%        (x1,t1) The xDimX1 state and the 1X1 time at the final step.
%        df      The function df returns the second derivatuve d2xdt2. If 
%                the general problem is being solved, then
%                probIsGeneral=true and df is called as df(xVec,t) where
%                the full vector of N values and their N first derivatives
%                is passed. If probIsGeneral=false, then df is called as
%                df(x,t), where x is just the values without their
%                derivatives.
%  probIsGeneral Indicates whether the problem being solved is the general
%                problem, where df depends on the frist derivatives, of the
%                special problem, where df does not depend on the first
%                derivatives.
% order,solutionChoice  A pair of optional parameters that specify the
%                highest order of the embedded Runge-Kutta pair to use as
%                well as the specific algorithm to use. Details are given
%                in the comments to the RungeKNystroemGStep function, for
%                the general problem, and in the RungeKNystroemSStep
%                function for the special problem.
%interpMainOrder A parameter indicating whether the interpolation
%               should be the same as the main order of the Runge-Kutta
%               pair. If false, then the order of the subsidiary formula is
%               used.
%             g The values of the second derivatives df evaluated at
%               various points as determined by the selected
%               algorithm.
%
%OUTPUTS: interpPolysA,interpPolysC An xDimX(polyOrder+1) and
%               xDimX(polyOrder) matrices that hold the interpolating
%               coefficients (A) and  control points (C) in Newton's form
%               for interpolation in each dimension in the range (t0,t1). 
%
%The function just performs an appropriate number of uniform
%Runge-Kutta-Nyström steps using the appropriate method (general or
%special) between t0 and t1 so that the function HermiteInterpPoly can be
%used to get an interpolating polynomial.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

xDim=size(x0,1);
xPosDim=xDim/2;

deltaT=t1-t0;

if(probIsGeneral)
    [orders,isFSAL]=RungeKNystroemGStep(order,solutionChoice);
else
    [orders,isFSAL]=RungeKNystroemSStep(order,solutionChoice);
end

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
dxVecdt=zeros(xDim,numPoints);

x(:,1)=x0;
x(:,end)=x1;
dfCur=g(:,1);
%The first half of the derivative of the state vector at a particular point
%is the second half of the state--from the definition of the Nyström
%problem.
dxVecdt(:,1)=[x0((xPosDim+1):end); dfCur];

%Fill in the intermediate values using short Runge-Kutta Nyström steps.
curT=t0;
for curStep=2:(numPoints-1)
    if(probIsGeneral)
        [x(:,curStep),~,gCur]=RungeKNystroemGStep(x(:,curStep-1),curT,df,deltaTStep,dfCur,order,solutionChoice);
    else
        [x(:,curStep),~,gCur]=RungeKNystroemSStep(x(:,curStep-1),curT,df,deltaTStep,dfCur,order,solutionChoice);
    end
    
    curT=curT+deltaTStep;
    if(isFSAL)
        dfCur=gCur(:,end);
    else
        if(probIsGeneral)
            dfCur=df(x(:,curStep),curT);
        else
            dfCur=df(x(1:xPosDim,curStep),curT);
        end
    end
    dxVecdt(:,curStep)=[x((xPosDim+1):end,curStep);dfCur];
end
if(isFSAL)
    dxVecdt(:,end)=[x((xPosDim+1):end,curStep);g(:,end)];
else
    if(probIsGeneral)
        dxVecdt(:,end)=[x((xPosDim+1):end,curStep);df(x1,t1)];
    else
        dxVecdt(:,end)=[x((xPosDim+1):end,curStep);df(x1(1:xPosDim),t1)];
    end
end

%Allocate space for the interpolating polynomials.
interpPolyA=zeros(xDim,2*numPoints);
interpPolyC=zeros(xDim,2*numPoints-1);
for curDim=1:xDim
    %The scaling is from a change of variables so that the range of
    %interpolation parameters goes from 0->1 and not from t0 to t1.
    [interpPolyA(curDim,:),interpPolyC(curDim,:)]=HermiteInterpPoly((tVals-t0)/deltaT,[x(curDim,:);deltaT*dxVecdt(curDim,:)]);
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
