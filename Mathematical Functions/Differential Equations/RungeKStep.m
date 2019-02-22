function [xPredMain,xPredSubsid,k,orders,isFSAL]=RungeKStep(xVal,curT,f,deltaT,fCur,order,solutionChoice)
%%RUNGEKSTEP Perform a single step of an explicit Runge-Kutta method whose
%            main integrating order is given by order. All implemented
%            Runge-Kutta routines are embedded routines, meaning that a
%            secondary result of a different order, either one order higher
%            or one lower, can also be returned along with the order of the
%            second output. This second output can be used in an algorithm
%            that has an adaptive step size (deltaT is automatically
%            varied), such as in the function RKAdaptiveOverRange.
%            Runge-Kutta methods are derivative-free techniques
%            for solving ordinary differential equations. That is,
%            integrating dx/dt=f(x,t) given initial conditions (x0,t0).
%            This function propagates xVal from time curT to  (curT+deltaT)
%            in one step. 32 different methods can be selected. To
%            automatically propagate a state over multiple steps with the
%            same step size, the functions RungeKSteps or RungeKAtTimes can
%            be used. To propagate a state using adaptive stepsize
%            correction, the function RKAdaptiveOverRange can be used.
%            Those functions call this function as a subroutine.
%
%INPUTS: The function can either be called with two inputs as
%        [orders,isFSAL]=RungeKStep(order,solutionChoice)
%        to just get the orders of a particular embedded Runge-Kutta method
%        and to figure out whether k(:,end) can be passed in place of fCur
%        (see below for details), or it can be run as
%        [xPredMain,xPredSubsid,k,orders,isFSAL]=RungeKStep(xVal,curT,f,deltaT,fCur,order,solutionChoice)
%        to actually do a full step and provide the orders. All of the
%        inputs in the full step are:
%              xVal The value of the (scalar or vector) state over which
%                   integration is being performed.
%              curT The time at which xVal is taken.
%                 f f(x,t) returns the derivative of x with respect to time
%                   taken at time t.
%            deltaT The size of the single (time) step over which the 
%                   Runge-Kutta integration is performed.
%              fCur The value f(xVal,curT). This is requested so that
%                   methods that are FSAL can pass k(:,end) on subsequent
%                   steps instead of having to perform a redundant
%                   evaluation of f. If omitted or an empty matrix is
%                   passed, the function f(xVal,curT) is evaluated to get
%                   fCur.
%             order The main integration order of the Runge-Kutta method.
%                   If this parameter is omitted, then the default order of
%                   5 is used. Order can range from 1 to 8.
%    solutionChoice Different Runge-Kutta formulae exist for each order. If
%                   this parameter is provided, then the selected formula
%                   for the given order is used. Otherwise the default
%                   (solutionChoice=0) formula is used. The algorithms
%                   chosen by pairs of (order, solutionChoice) are:
%                   (1,0) Fehlberg's RK1(2)3F method, Table XIV in [3].
%                   (1,1) The Euler-Cauchy RK1(2)2F method as given by
%                         Fehlberg in Table XV  of [3].
%                   (2,0) Bettis' RK2(3)4 method from Table II in [12].
%                   (2,1) Fehlberg's RK2(3)4 method from Table XI in [3].
%                   (2,2) Fehlberg's RK2(3)3 method from Table XII in [3].
%                   (3,0) Bogacki and Shampine's RK3(2)4F formula of [8].
%                   (3,1) Fehlberg's RK3(4)5F method (Formula 1), Table VII
%                         of [3].
%                   (3,2) Fehlberg's RK3(4)5F method (Formula 2), Table
%                         VIII of [3].
%                   (3,3) The Cash-Karp RK3(2)4 algorithm of [10].
%                   (4,0) Fehlberg's RK4(5)6 formula, Table III in [3].
%                   (4,1) The Cash-Karp RK4(3)6 algorithm of [10].
%                   (4,2) Fehlberg's RK4(5)6 formula, Table II of [3].
%                   (4,3) Sarafyan's RK4(5)6 method as given by Fehlberg in 
%                         Table IV of [3].
%                   (5,0) RK5(4)7FM from Dormand and Prince [5].
%                   (5,1) RK5(4)7FEq3 from Higham and Hall in [9].
%                   (5,2) RK5(4)7FEq1 from Higham and Hall in [9].
%                   (5,3) RK5(4)7FEq2 from Higham and Hall in [9].
%                   (5,4) The Cash-Karp RK5(4)6 algorithm of [10].
%                   (5,5) RK5(4)7FC from Dormand and Prince [7].
%                   (5,6) RK5(4)7FS from Dormand and Prince [5].
%                   (5,7) RK5(4)6M from Dormand and Prince [5].
%                   (5,8) Verner's RK5(6)8 procedure of Table 5 in [4].
%                   (5,9) Fehlberg's RK5(6)8 formula of Table II of [2].
%                   (6,0) RK6(5)8S from Dormand and Prince [7].
%                   (6,1) RK6(5)8C from Dormand and Prince [7].
%                   (6,2) RK6(5)8M from Prince and Dormand [7].
%                   (6,3) Verner's RK6(7)10 solution in [4].
%                   (6,4) Fehlberg's's RK6(7)10 formula of [2].
%                   (7,0) Verner's RK7(8)13 algorithm from [4].
%                   (7,1) Fehlberg's RK7(8)13 formula of [2].
%                   (8,0) RK8(7)13M from Prince and Dormand [7].
%                   (8,1) Verner's RK8(9)16 algorithm in [4].
%
%OUTPUTS: If the function is run with two inputs as
%         [orders,isFSAL]=RungeKStep(order,solutionChoice)
%         Then there are just two outputs, orders and isFSAL, described
%         below.
%         Otherwise, all of the outputs are
%         xPredMain The state propagated forward an interval of deltaT at
%                   the order of precision given by the input order. This
%                   is the main integration order.
%         xPredSubsid The subsidiary estimate. This is a different order
%                   than the main integration order can can be used in
%                   algorithms that adaptively adjust the stepsize.
%                 k The values of the derivatives f evaluated at various
%                   points as determined by the selected algorithm. These
%                   can sometimes be reused in interpolation routines to
%                   reduce the number of computations required. In all of
%                   the methods, k(:,1) is equal to f(xVal,curT). In FSAL
%                   methods, k(:,end) is f(xPredMain,curT+deltaT)
%            orders A 2X1 vector whereby order(1)=order, the input
%                   parameter, and orders(2) is the order of xPredSubsid.
%                   The value of orders(2) depends on the algorithm chosen
%                   by the combination of the order and solutionChoice
%                   inputs.
%            isFSAL Indicates whether the first evaluation of the f of the
%                   next step is equal to k(:,end). If so, this can be
%                   passed to the function on the next step rather than
%                   having to make an additional evaluation of f.
%
%The formulae are given in this file in the form of the components of their
%Butcher table. The Butcher table, which was introduced in [1], is a way of
%expressing Runge-Kutta formulae. A very basic introduction to Runge-Kutta
%integration is given in Chapter 5.5 of [11].
%
%The formulae were referred to above using a notation such as RK4(5)6. In
%general, the normenclature is
%RKq(p)s[F]X
%where
%q is the order of the main integrating formula
%p is the order of the subsidiary formula
%s is the number of stages (derivative function evaluations)
%The presence of an F indicates that the function is first same as last
%(FSAL). That is, f evaluated at the next step is k(:,end) and thus
%k(:,end) can be passed for fCur for the next step. The absence of an F
%indicates that the function is not FSAL.
%The other values X depend on whether the author of the paper deriving the
%formula abides by general conventions/ added. Common last letters of M, S,
%C, and G, mean that the formula was derived..
%M minimizing the error norm
%S with an enlarged stability region
%C as a compromise between 'M' and 'S',
%The Eq1 to Eq3 suffixes with Higham and Hall's paper just identify the
%formula.
%
%REFERENCES:
%[1] J. Butcher, "On Runge-Kutta processes of high order," Journal of the
%    Australian Mathematical Society, vol. 4, no. 2, pp. 179-194, May 1964.
%[2] E. Fehlberg, "Classical fifth-, sixth-, seventh-, and eighth-order
%    Runge-Kutta formulas with stepsize control," George C. Marshall Space
%    Flight Center, National Aeronautics and Space Administration,
%    Marshall, AL, Tech. Rep. NASA TR R-287, Oct. 1968.
%[3] E. Fehlberg, "Low-order classical Runge-Kutta formulas with stepsize
%    control and their application to some heat transfer problems," George
%    C. Marshall Space Flight Center, National Aeronautics and Space
%    Administration, Marshall, AL, Tech. Rep. NASA TR R-315, Jul. 1969.
%[4] J. H. Verner, "Explicit Runge-Kutta methods with estimates of the
%    local truncation error," SIAM Journal on Numerical Analysis, vol. 15,
%    no. 4, pp. 772-790, 1978.
%[5] J. R. Dormand and P. J. Prince, "A family of embedded Runge-Kutta
%    formulae," Journal of Computational and Applied Mathematics, vol. 6,
%    no. 1, pp. 19-26, Mar. 1980.
%[6] P. J. Prince and J. R. Dormand, "High order embedded Runge-Kutta
%    formulae," Journal of Computational and Applied Mathematics, vol. 7,
%    no. 1, pp. 67-75, Mar. 1981.
%[7] J. R. Dormand and P. J. Prince, "A reconsideration of some embedded
%    Runge-Kutta formulae," Journal of Computational and Applied
%    Mathematics, vol. 15, no. 2, pp. 203-211, Jun. 1986.
%[8] P. Bogacki and L. F. Shampine, "A 3(2) pair of Runge-Kutta formulas,"
%    Applied Mathematics Letters, vol. 2, no. 4, pp. 321-325, 1989.
%[9] D. J. Higham and G. Hall, "Embedded Runge-Kutta formulae with stable
%    equilibrium states," Journal of Computational and Applied Mathematics,
%    vol. 29, no. 1, pp. 25-33, Jan. 1990.
%[10] J. R. Cash and A. H. Karp, "A variable order Runge-Kutta method for
%    initial value problems with rapidly varying right-hand sides," ACM
%    Transactions on Mathematical Software, vol. 16, no. 3, pp. 201-222,
%    Sep. 1990.
%[11] R. L. Burden and J. D. Faires, Numerical Analysis, 9th ed. Boston,
%    MA: Brooks/ Cole, 2011.
%[12] M. K. Horn, "Scaled Runge-Kutta algorithms for treating the problem
%    of dense output," National Aeronautics and Space Administration,
%    Houston, TX, Tech. Rep. 58239, 1982.
%
%February 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin==2)
    %If the function was just called as
    %[orders,isFSAL]=RungeKStep(order,solutionChoice);
    %So that one can get the orders.
    order=xVal;
    solutionChoice=curT;
else
    if(nargin<7)
        solutionChoice=0;
    end

    if(nargin<6)
        order=5;
    end
    
    if(nargin<5||isempty(fCur))
        fCur=f(xVal,curT);
    end
end

%The isFSAL variable indicates whether the last value of k is equal to f
%evaluated at the updated value.
switch(order)
    case 1
        switch(solutionChoice)
            case 0
                %Fehlberg's RK1(2)3F method, Table XIV in [3].
                A=[0,       0,          0;
                   1/2,     0,          0;
                   1/256,   255/256,    0];
                bOrderMain=[1/256; 255/256; 0];
                bOrderSubsid=[1/512; 255/256; 1/512];
                c=[0; 1/2; 1];
                orders=[1;2];
                isFSAL=true;
            case 1
                %The Euler-Cauchy RK1(2)2F method as given by Fehlberg in
                %Table XV  of [3]. Table XVI of [3] report shows this
                %method to be less efficient and less accurate than the
                %algorithm of Table XIV.
                A=[0, 0;
                   1, 0];
                c=[0; 1];
                bOrderMain=[1;0];
                bOrderSubsid=[1/2;1/2];
                orders=[1;2];
                isFSAL=true;
            otherwise
                error('Unknown Solution Choice');
        end
    case 2
        switch(solutionChoice)
            case 0%Bettis' RK2(3)4 method from Table II in [12]
                A=[0,   0,          0,          0;
                   1/2, 0,          0,          0;
                   1/4, 1/4,        0,          0;
                   0,   -143/144,   287/144,    0];
               c=[0; 1/2;   1/2;    1];
               bOrderMain=[0;   0;  1;  0];
               bOrderSubsid=[1/6;   0;  2/3;    1/6];
               orders=[2;3];
               isFSAL=false;
            case 1
                %Fehlberg's RK2(3)4 method from Table XI in [3].
                A=[0,       0,      0,      0;
                   1/4,     0,      0,      0;
                   -189/800,729/800,0,      0;
                   214/891, 1/33,   650/891,0];
                c=[0; 1/4; 27/40; 1];
                bOrderMain=[214/891; 1/33; 650/891; 0];
                bOrderSubsid=[533/2106; 0; 800/1053; -1/78];
                orders=[2;3];
                isFSAL=false;
            case 2
                %Fehlberg's RK2(3)3 method from Table XII in [3]. Table XVI
                %of Fehlberg's report shows this method to be less
                %efficient and less accurate than the algorithm of Table
                %XI.
                A=[0,   0,  0;
                   1,   0,  0;
                   1/4, 1/4,0];
                c=[0; 1; 1/2];
                bOrderMain=[1/2;1/2;0];
                bOrderSubsid=[1/6;1/6;2/3];
                orders=[2;3];
                isFSAL=false;
            otherwise
                error('Unknown Solution Choice');
        end
    case 3
        switch(solutionChoice)
            case 0
                %Bogacki and Shampine's RK3(2)4F formula of [8].
                A=[0,   0,      0,  0;
                   1/2, 0,      0,  0;
                   0,   3/4     0,  0;
                   2/9, 1/3,    4/9,0];
                c=[0; 1/2; 3/4; 1];
                bOrderMain=[2/9;1/3;4/9; 0];
                bOrderSubsid=[7/24; 1/4; 1/3; 1/8];
                orders=[3;2];
                isFSAL=true;
            case 1
                %Fehlberg's RK3(4)5F method (Formula 1), Table VII of [3].
                A=[0,       0,          0,          0,      0;
                   1/4,     0,          0,          0,      0;
                   4/81,    32/81,      0,          0,      0;
                   57/98,   -432/343,   1053/686,   0,      0;
                   1/6,     0,          27/52,      49/156, 0];
                c=[0; 1/4; 4/9; 6/7; 1];
                bOrderMain=[1/6; 0; 27/52; 49/156; 0];
                bOrderSubsid=[43/288; 0; 243/416; 343/1872; 1/12];
                orders=[3;4];
                isFSAL=true;
            case 2
                %Fehlberg's RK3(4)5F method (Formula 2), Table VIII  of [3].
                A=[0,       0,              0,          0,          0;
                   2/7,     0,              0,          0,          0;
                   77/900,  343/900,        0,          0,          0;
                   805/1444,-77175/54872,   97125/54872,0,          0;
                   79/490,  0,              2175/3626,  2166/9065,  0];
                c=[0; 2/7; 7/15; 35/38; 1];
                bOrderMain=[79/490; 0; 2175/3626; 2166/9065; 0];
                bOrderSubsid=[229/1470; 0; 1125/1813; 13718/ 81585; 1/18];
                orders=[3;4];
                isFSAL=true;
            case 3
                %The Cash-Karp RK3(2)4 algorithm of [10], which was
                %specifically designed for order 3, not the truncated order
                %5 version.
                A=[0,   0,      0,  0;
                   1/5, 0,      0,  0;
                   3/40,9/40,   0,  0;
                   3/10,-9/10,  6/5,0];
                c=[0; 1/5; 3/10; 3/5];
                bOrderMain=[19/54; 0; -10/27; 55/54];
                bOrderSubsid=[-3/2; 5/2; 0; 0];
                orders=[3;2];
                isFSAL=false;
            otherwise
                error('Unknown Solution Choice');
        end
    case 4
        switch(solutionChoice)
             case 0
                %Fehlberg's RK4(5)6 formula (Formula 2), Table III in [3].
                A=[0,           0,          0,          0,          0,      0;
                   1/4,         0,          0,          0,          0,      0;
                   3/32,        9/32,       0,          0,          0,      0;
                   1932/2197,   -7200/2197, 7296/2197,  0,          0,      0;
                   439/216,     -8,         3680/513,   -845/4104,  0,      0;
                   -8/27,       2,          -3544/2565, 1859/4104,  -11/40, 0];
                c=[0; 1/4; 3/8; 12/13; 1; 1/2];
                bOrderMain=[25/216; 0; 1408/2565;  2197/4104;   -1/5;  0];
                bOrderSubsid=[16/135; 0; 6656/12825; 28561/56430; -9/50; 2/55];
                orders=[4;5];
                isFSAL=false;
            case 1
                %The Cash-Karp RK4(3)6 algorithm of [10] (which is a 
                %truncated version of their RK5(4) algorithm.
                A=[0,           0,          0,          0,              0,          0;
                   1/5,         0,          0,          0,              0,          0;
                   3/40,        9/40,       0,          0,              0,          0;
                   3/10,        -9/10,      6/5,        0,              0,          0;
                   -11/54,      5/2,        -70/27,     35/27,          0,          0;
                   1631/55296,  175/512,    575/13824,  44275/110592,   253/4096,   0];
                c=[0;    1/5;    3/10;   3/5;    1;  7/8];
                bOrderMain=[2825/27648; 0; 18575/48384; 13525/55296; 277/14336; 1/4];
                bOrderSubsid=[19/54; 0;  -10/27; 55/54;  0;  0];
                orders=[4;3];
                isFSAL=false;
            case 2
                %Fehlberg's RK4(5)6 formula (Formula 1) of Table II of [3].
                %Table XVI of Fehlberg's report shows this method to
                %be less efficient with comparable accuracy to the
                %algorithm of Table III.
                A=[0,       0,          0,          0,      0,      0;
                   2/9,     0,          0,          0,      0,      0;
                   1/12,    1/4,        0,          0,      0,      0;
                   69/128,  -243/128,   135/64,     0,      0,      0;
                   -17/12,  27/4,       -27/5,      16/15,  0,      0;
                   65/432,  -5/16,      13/16,      4/27,   5/144,  0];
                c=[0; 2/9; 1/3; 3/4; 1; 5/6];
                bOrderMain=[1/9; 0; 9/20; 16/45; 1/12; 0];
                bOrderSubsid=[47/450; 0; 12/25; 32/225; 1/30; 6/25];
                orders=[4;5];
                isFSAL=false;
            case 3
                %Sarafyan's RK4(5)6 method as given by Fehlberg in [3],
                %Table IV. Table XVI of Fehlberg's report shows this method
                %to be significantly less efficient than either of 
                %Fehlberg's methods, though it might have a higher
                %accuracy.
                A=[0,       0,      0,      0,      0,          0;
                   1/2,     0,      0,      0,      0,          0;
                   1/4,     1/4,    0,      0,      0,          0;
                   0,       -1,     2,      0,      0,          0;
                   7/27,   10/27,   0,      1/27,   0,          0;
                   28/625, -1/5,    546/625,54/625, -378/625,   0];
                c=[0; 1/2; 1/2; 1; 2/3; 1/5];
                bOrderMain=[1/6; 0; 2/3; 1/6; 0; 0];
                bOrderSubsid=[1/24; 0; 0; 5/48; 27/56; 125/336];
                orders=[4;5];
                isFSAL=false;
            otherwise
                error('Unknown Solution Choice');
        end
    case 5
        switch(solutionChoice)
            case 0
                %RK5(4)7FM from Dormand and Prince [5].
                A = [0,         0,          0,          0,          0,              0,    0;
                     1/5,       0,          0,          0,          0,              0,    0;
                     3/40,      9/40,       0,          0,          0,              0,    0;
                     44/45,     -56/15,     32/9,       0,          0,              0,    0;
                     19372/6561,-25360/2187,64448/6561, -212/729,   0,              0,    0; 
                     9017/3168, -355/33,    46732/5247, 49/176,     -5103/18656,    0,    0;
                     35/384,    0,          500/1113,   125/192,    -2187/6784,     11/84,0];
                c = [0; 1/5; 3/10; 4/5; 8/9; 1; 1];
                bOrderMain = [35/384; 0; 500/1113; 125/192; -2187/6784; 11/84; 0];
                bOrderSubsid = [5179/57600; 0; 7571/16695; 393/640; -92097/339200; 187/2100; 1/40];
                orders=[5;4];
                isFSAL=true;
            case 1
                %RK5(4)7FEq3 from Higham and Hall in [9].
                A = [0,             0,              0,              0,                  0,          0,      0;
                     11/45,         0,              0,              0,                  0,          0,      0;
                     11/120,        11/40,          0,              0,                  0,          0,      0;
                     106865/87808,  -408375/87808,  193875/43904,   0,                  0,          0,      0;
                     79503/121000,  -1053/440,      147753/56870,   27048/710875,       0,          0,      0;
                     89303/78045,   -2025/473,      994650/244541,  -2547216/28122215,  475/2967,   0,      0;
                     1247/10890,    0,              57375/108053,   -1229312/1962015,   125/207,    43/114, 0];

                c=[0;11/45;11/30;55/56;9/10;1;1];
                bOrderMain=[1247/10890;0;57375/108053;-1229312/1962015;125/207;43/114;0];
                bOrderSubsid=[21487/185130;0;963225/1836901;-39864832/33354255;2575/3519;4472/4845;-1/10];
                orders=[5;4];
                isFSAL=true;
            case 2
                %RK5(4)7FEq1 from Higham and Hall in [9].
                A=[0,       0,      0,      0,      0,      0,      0;
                   2/9,     0,      0,      0,      0,      0,      0;
                   1/12,    1/4,    0,      0,      0,      0,      0;
                   1/8,     0,      3/8,    0,      0,      0,      0;
                   91/500,  -27/100,78/125, 8/125,  0,      0,      0;
                   -11/20,  27/20,  12/5,   -36/5,  5,      0,      0;
                   1/12,    0,      27/32,  -4/3,   125/96, 5/48,   0];
                c=[0; 2/9; 1/3; 1/2; 3/5; 1; 1];
                bOrderMain=[1/12; 0; 27/32; -4/3; 125/96; 5/48; 0];
                bOrderSubsid=[2/15; 0; 27/80; -2/15; 25/48; 1/24; 1/10];
                orders=[5;4];
                isFSAL=true;
            case 3
                %RK5(4)7FEq2 from Higham and Hall in [9].
                A=[0,           0,          0,              0,              0,              0,      0;
                   2/13,        0,          0,              0,              0,              0,      0;
                   3/52,        9/52,       0,              0,              0,              0,      0;
                   12955/26244, -15925/8748,12350/6561,     0,              0,              0,      0;
                   -10383/52480,13923/10496,-176553/199424, 505197/997120,  0,              0,      0;
                   1403/7236,   -429/268,   733330/309339,  -7884/8911,     104960/113967,  0,      0;
                   181/2700,    0,          656903/1846800, 19683/106400,   34112/110565,   67/800, 0];
                c=[0; 2/13; 3/13; 5/9; 3/4; 1; 1];
                bOrderMain=[181/2700; 0; 656903/1846800; 19683/106400; 34112/110565; 67/800; 0];
                bOrderSubsid=[11377/154575; 0; 35378291/105729300; 343359/1522850; 535952/1947645; 134/17175; 1/12];
                orders=[5;4];
                isFSAL=true;
            case 4
                %The Cash-Karp RK5(4)6 algorithm of [10].
                A=[0,           0,          0,          0,              0,          0;
                   1/5,         0,          0,          0,              0,          0;
                   3/40,        9/40,       0,          0,              0,          0;
                   3/10,        -9/10,      6/5,        0,              0,          0;
                   -11/54,      5/2,        -70/27,     35/27,          0,          0;
                   1631/55296,  175/512,    575/13824,  44275/110592,   253/4096,   0];
                c=[0;    1/5;    3/10;   3/5;    1;  7/8];
                bOrderMain=[2825/27648; 0; 18575/48384; 13525/55296; 277/14336; 1/4];  
                bOrderSubsid=[37/348; 0; 250/621; 125/594; 0; 512/1771];
                orders=[5;4];
                isFSAL=false;
            case 5
                %RK5(4)7FC from Dorman and Prince [7].
                A = [0,         0,          0,          0,              0,          0,      0;
                    1/5,        0,          0,          0,              0,          0,      0;
                    3/40,       9/40,       0,          0,              0,          0,      0;
                    264/2197,   -90/2197,   840/2197,   0,              0,          0,      0;
                    932/3645,   -14/27,     3256/5103,  7436/25515,     0,          0,      0;
                    -367/513,   30/19,      9940/5643,  -29575/8208,    6615/3344,  0,      0;
                    35/432,     0,          8500/14553, -28561/84672,   405/704,  19/196,   0];
                c=[0; 1/5; 3/10; 6/13; 2/3; 1; 1];
                bOrderMain=[35/432; 0; 8500/14553; -28561/84672; 405/704; 19/196; 0];
                bOrderSubsid=[11/108; 0; 6250;14553; 2197/21168; 81/176; 171/1960; 1/40];
                orders=[5;4];
                isFSAL=true;
            case 6
                %RK5(4)7FS from Dorman and Prince [5]. This is supposed to 
                %have a large stability region, even if it is not as 
                %efficient as other methods.
                A=[0,       0,      0,      0,          0,      0,      0;
                   2/9,     0,      0,      0,          0,  	0,      0;
                   1/12,    1/4,    0,      0,          0,      0,      0;
                   55/325,  -25/108,50/81,  0,          0,      0,      0;
                   83/330,  -13/22, 61/66,  9/100,      0,      0,      0;
                   -19/28,  9/4,    1/7,    -27/7,      22/7,   0,      0;
                   19/200,  0,      3/5,    -243/400,   33/40,  7/80,   0];

                c=[0; 2/9; 1/3; 5/9; 2/3; 1; 1];
                bOrderMain=[19/200; 0; 3/5; -243/400; 33/40; 7/80; 0];
                bOrderSubsid=[431/5000; 0; 333/500; -7857/10000; 957/1000; 193/2000; -1/50];
                orders=[5;4];
                isFSAL=true;
            case 7
                %RK5(4)6M from Dormand and Prince [5].
                A=[0,       0,      0,          0,      0,      0;
                   1/5,     0,      0,          0,      0,      0;
                   3/40,    9/40,   0,          0,      0,      0;
                   3/10,    -9/10,  6/5,        0,      0,      0;
                   226/729, -25/27, 880/729,    55/729, 0,      0;
                   -181/270,5/2,    -266/297,   -91/27, 189/55, 0];
                c=[0; 1/5; 3/10; 3/5; 2/3; 1];
                bOrderMain=[19/216; 0; 1000/2079; -125/216; 81/88; 5/56];
                bOrderSubsid=[31/540; 0; 190/297; -145/108; 351/220; 1/20];
                orders=[5;4];
                isFSAL=false;
            case 8
                %Verner's RK5(6)8 procedure of Table 5 in [4].
                A=[0,           0,          0,          0,          0,          0,  0,          0;
                   1/18,        0,          0,          0,          0,          0,  0,          0;
                   -1/12,       1/4,        0,          0,          0,          0,  0,          0;
                   -2/81,       4/27,       8/81,       0,          0,          0,  0,          0;
                   40/33,       -4/11,      -56/11,     54/11,      0,          0,  0,          0;
                   -369/73,     72/73,      5380/219,   -12285/584, 2695/1752,  0,  0,          0;
                   -8716/891,   656/297,    39520/891,  -416/11,    52/27,      0,  0,          0;
                   3015/256,    -9/4,       -4219/78,   5985/128,   -539/384,   0,  693/3328,   0];
                c=[0; 1/18; 1/6; 2/9; 2/3; 1; 8/9; 1];
                bOrderMain=[3/80; 0; 4/25; 243/1120; 77/160; 73/700; 0; 0];
                bOrderSubsid=[57/640; 0; -16/65; 1377/2240; 121/320; 0; 891/8320; 2/35];
                orders=[5;6];
                isFSAL=false;
            case 9
                %Fehlberg's RK5(6)8 formula of Table II of [2].
                A=[0,       0,      0,      0,      0,      0,  0,  0;
                   1/6,     0,      0,      0,      0,      0,  0,  0;
                   4/75,    16/75,  0,      0,      0,      0,  0,  0;
                   5/6,     -8/3,   5/2,    0,      0,      0,  0,  0;
                   -8/5,    144/25, -4,     16/25,  0,      0,  0,  0;
                   361/320, -18/5,  407/128,-11/80, 55/128, 0,  0,  0;
                   -11/640, 0,      11/256, -11/160,11/256, 0,  0,  0;
                   93/640,  -18/5,  803/256,-11/160,99/256, 0,  1,  0];
                c=[0; 1/6; 4/15; 2/3; 4/5; 1; 0; 1];
                bOrderMain=[31/384; 0; 1125/2816; 9/32; 125/768; 5/66; 0; 0];
                bOrderSubsid=[7/1408; 0; 1125/2816; 9/32; 125/768; 0; 5/66; 5/66];
                orders=[5;6];
                isFSAL=false;
            otherwise
                error('Unknown Solution Choice');
        end
    case 6
        switch(solutionChoice)
            case 0
                %RK6(5)8S from Dormand and Prince [7].
                A=[0,       0,          0,              0,              0,          0,      0,  0;
                   1/4,     0,          0,              0,              0,          0,      0,  0;
                   3/25,    9/50,       0,              0,              0,          0,      0,  0;
                   102/343, -1368/343,  1560/343,       0,              0,          0,      0,  0;
                   -3/100,  36/25,      -12/13,         147/1300,       0,          0,      0,  0;
                   37/225,  -48/25,     872/351,        49/1053,        2/81,       0,      0,  0;
                   11/648,  14/3,       -10193/2106,    -30331/50544,   1025/1944,  59/48,  0,  0;
                   796/1701,-352/63,    134093/22113,   -78281/75816,   -9425/20412,781/504,0,  0];
                c=[0; 1/4; 3/10; 6/7; 3/5; 4/5; 1; 1];
                bOrderMain=[29/324; 0; 3400/7371; -16807/25272; -125/1944; 25/24; 1/84; 1/8];
                bOrderSubsid=[2041/21600; 0; 748/1755; -2401/46800; 11/108; 59/160; 3/50; 0];
                orders=[6;5];
                isFSAL=false;
            case 1
                %RK6(5)8C from Dormand and Prince [7].
                A=[0,               0,          0,              0,                  0,                  0,                  0,  0;
                   1/10,            0,          0,              0,                  0,                  0,                  0,  0;
                   1/36,            5/36,       0,              0,                  0,                  0,                  0,  0;
                   10/243,          20/243,     8/81,           0,                  0,                  0,                  0,  0;
                   4047/5500,       -18/55,     -4212/1375,     17901/5500,         0,                  0,                  0,  0;
                   -5587/4125,      24/55,      9576/1375,      -140049/23375,      38/51,              0,                  0,  0;
                   12961/2376,      -35/33,     -160845/5434,   1067565/38896,      -103375/47736,      32875/35568,        0,  0;
                   702799/199584,   -1865/2772, -2891375/152152,19332955/1089088,   -5356375/4009824,   2207875/2987712,    0,  0];

                c=[0; 1/10; 1/6; 2/9; 3/5; 4/5; 1; 1];
                bOrderMain=[1/12; 0; -216/1235; 6561/12376; 1375/5304; 1375/5928; -5/168; 1/10];
                bOrderSubsid=[163/1440; 0; -2628/6175; 13851/17680; 1525/7956; 6575/23712; 3/50; 0];
                orders=[6;5];
                isFSAL=false;
            case 2
                %RK6(5)8M from Prince and Dormand [7].
                A=[0,               0,          0,                  0,                  0,              0,              0,  0;
                   1/10,            0,          0,                  0,                  0,              0,              0,  0;
                   -2/81,           20/81,      0,                  0,                  0,              0,              0,  0;
                   615/1372,        -270/343,   1053/1372,          0,                  0,              0,              0,  0;
                   3243/5500,       -54/55,     50949/71500,        4998/17875,         0,              0,              0,  0;
                   -26492/37125,    72/55,      2808/23375,         -24206/37125,       338/459,        0,              0,  0;
                   5561/2376,       -35/11,     -24117/31603,       899983/200772,      -5225/1836,     3925/4056,      0,  0;
                   465467/266112,   -2945/1232, -5610201/14158144,  10513573/3212352,   -424325/205632, 376225/454272,  0,  0];
                c=[0; 1/10; 2/9; 3/7; 3/5; 4/5; 1; 1];
                bOrderMain=[61/864; 0; 98415/321776; 16807/146016; 1375/7344; 1375/5408; -37/1120; 1/10];
                bOrderSubsid=[821/10800; 0; 19683/71825; 175273/912600; 395/3672; 785/2704; 3/50; 0];
                orders=[6;5];
                isFSAL=false;
            case 3
                %Verner's RK6(7)10 solution in [4].
                A=[0,               0,  0,              0,              0,              0,                  0,          0,  0,          0;
                   1/12,            0,  0,              0,              0,              0,                  0,          0,  0,          0;
                   0,               1/6,0,              0,              0,              0,                  0,          0,  0,          0;
                   1/16,            0,  3/16,           0,              0,              0,                  0,          0,  0,          0;
                   21/16,           0,  -81/16,         9/2,            0,              0,                  0,          0,  0,          0;
                   1344688/250563,  0,  -1709184/83521, 1365632/83521,  -78208/250563,  0,                  0,          0,  0,          0;
                   -559/384,        0,  6,              -204/47,        14/39,          -4913/78208,        0,          0,  0,          0;
                   -625/224,        0,  12,             -456/47,        48/91,          14739/136864,       6/7,        0,  0,          0;
                   -12253/99144,    0,  16/27,          16/459,         29072/161109,   -2023/75816,        112/12393,  0,  0,          0;
                   30517/2512,      0,  -7296/157,      268728/7379,    2472/2041,      -3522621/10743824,  132/157,    0,  -12393/4396,0];
                c=[0; 1/12; 1/6; 1/4; 3/4; 16/17; 1/2; 1; 2/3; 1];
                bOrderMain=[7/90; 0; 0; 16/45; 16/45; 0; 2/15; 7/90; 0; 0];
                bOrderSubsid=[2881/40320; 0; 0; 1216/2961; -2624/4095; 24137569/57482880; -4/21; 0; 4131/3920; -157/1260];
                orders=[6;7];
                isFSAL=false;
            case 4
                %Fehlberg's's RK6(7)10 formula of [2].
                A=[0,           0,      0,          0,              0,          0,          0,          0,  0,  0;
                   2/33,        0,      0,          0,              0,          0,          0,          0,  0,  0;
                   0,           4/33,   0,          0,              0,          0,          0,          0,  0,  0;
                   1/22,        0,      3/22,       0,              0,          0,          0,          0,  0,  0;
                   43/64,       0,      -165/64,    77/32,          0,          0,          0,          0,  0,  0;
                   -2383/486,   0,      1067/54,    -26312/1701,    2176/1701,  0,          0,          0,  0,  0;
                   10077/4802,  0,      -5643/686,  116259/16807,   -6420/16807,1053/2401,  0,          0,  0,  0;
                   -733/176,    0,      141/8,      -335763/23296,  216/77,     -4617/2816, 7203/9152,  0,  0,  0;
                   15/352,      0,      0,          -5445/46592,    18/77,      -1215/5632, 1029/18304, 0,  0,  0;
                   -1833/352,   0,      141/8,      -51237/3584,    18/7,       -729/512,   1029/1408,  0,  1,  0]; 
                c=[0; 2/33; 4/33; 2/11; 1/2; 2/3; 6/7; 1; 0; 1];
                bOrderMain=[77/1440;   0;  0;  1771561/6289920; 32/105; 243/2560; 16807/74880; 11/270; 0; 0];
                bOrderSubsid=[11/864;    0;  0;  1771561/6289920; 32/105; 243/2560; 16807/74880; 0; 11/270; 11/270];
                orders=[6;7];
                isFSAL=false;
            otherwise
                error('Unknown Solution Choice');
        end
    case 7
        switch(solutionChoice)
            case 0
                %Verner's RK7(8)13 algorithm from [4].
                A=[0,               0,      0,          0,              0,                  0,              0,              0,          0,          0,          0,  0,      0;
                   1/4,             0,      0,          0,              0,                  0,              0,              0,          0,          0,          0,  0,      0;
                   5/72,            1/72,   0,          0,              0,                  0,              0,              0,          0,          0,          0,  0,      0;
                   1/32,            0,      3/32,       0,              0,                  0,              0,              0,          0,          0,          0,  0,      0;
                   106/125,         0,      -408/125,   352/125,        0,                  0,              0,              0,          0,          0,          0,  0,      0;
                   1/48,            0,      0,          8/33,           125/528,            0,              0,              0,          0,          0,          0,  0,      0;
                   -1263/2401,      0,      0,          39936/26411,    -64125/26411,       5520/2401,      0,              0,          0,          0,          0,  0,      0;
                   37/392,          0,      0,          0,              1625/9408,          -2/15,          61/6720,        0,          0,          0,          0,  0,      0;
                   17176/25515,     0,      0,          -47104/25515,   1325/504,           -41792/25515,   20237/145800,   4312/6075,  0,          0,          0,  0,      0;
                   -23834/180075,   0,      0,          -77824/1980825, -636635/633864,     254048/300125,  -183/7000,      8/11,       -324/3773,  0,          0,  0,      0;
                   12733/7600,      0,      0,          -20032/5225,    456485/80256,       -42599/7125,    339227/912000,  -1029/4180, 1701/1408,  5145/2432,  0,  0,      0;
                   -27061/204120,   0,      0,          40448/280665,   -1353775/1197504,   17662/25515,    -71687/1166400, 98/225,     1/16,       3773/11664, 0,  0,      0;
                   11203/8680,      0,      0,          -38144/11935,  2354425/458304,      -84046/16275,   673309/1636800, 4704/8525,  9477/10912,  -1029/992,  0,  729/341, 0];
                c=[0; 1/4; 1/12; 1/8; 2/5; 1/2; 6/7; 1/7; 2/3; 2/7; 1; 1/3; 1];
                bOrderMain=[13/288;  0; 0; 0; 0; 32/125; 31213/144000; 2401/12375;  1701/14080; 2401/19200; 19/950; 0;         0];
                bOrderSubsid=[31/720; 0; 0; 0; 0; 16/75;  16807/79200;  16807/79200; 243/1760;   0;          0;      243/1760; 31/720];
                orders=[7;8];
                isFSAL=false;
            case 1
                %Fehlberg's RK7(8)13 formula of [2].
                A=[0,           0,      0,      0,          0,          0,      0,          0,      0,      0,      0,  0,  0;
                   2/27,        0,      0,      0,          0,          0,      0,          0,      0,      0,      0,  0,  0;
                   1/36,        1/12,   0,      0,          0,          0,      0,          0,      0,      0,      0,  0,  0;
                   1/24,        0,      1/8,    0,          0,          0,      0,          0,      0,      0,      0,  0,  0;
                   5/12,        0,      -25/16, 25/16,      0,          0,      0,          0,      0,      0,      0,  0,  0;
                   1/20,        0,      0,      1/4,        1/5,        0,      0,          0,      0,      0,      0,  0,  0;
                   -25/108,     0,      0,      125/108,    -65/27,     125/54, 0,          0,      0,      0,      0,  0,  0;
                   31/300,      0,      0,      0,          61/225,     -2/9,   13/900,     0,      0,      0,      0,  0,  0;
                   2,           0,      0,      -53/6,      704/45,     -107/9, 67/90,      3,      0,      0,      0,  0,  0;
                   -91/108,     0,      0,      23/108,     -976/135    311/54, -19/60,     17/6,   -1/12,  0,      0,  0,  0;
                   2383/4100,   0,      0,      -341/164,   4496/1025,  -301/82, 2133/4100, 45/82,  45/164, 18/41,  0,  0,  0;
                   3/205,       0,      0,      0,          0,          -6/41,   -3/205,    -3/41,  3/41,   6/41,   0,  0,  0;
                   1777/4100,   0,      0,      -341/164,   4496/1025,  -289/82, 2193/4100  51/82,  33/164, 12/41,  0,  1,  0];

                c=[0; 2/27; 1/9; 1/6; 5/12; 1/2; 5/6; 1/6; 2/3; 1/3; 1; 0; 1];
                bOrderMain=[41/840; 0; 0; 0; 0; 34/105; 9/35; 9/35; 9/280; 9/280; 41/840; 0; 0];
                bOrderSubsid=[0; 0; 0; 0; 0; 34/105; 9/35; 9/35; 9/280; 9/280; 0; 41/840; 41/840];
                orders=[7;8];
                isFSAL=false;
            otherwise
                error('Unknown Solution Choice');
        end
    case 8
        switch(solutionChoice)
            case 0
                %RK8(7)13M from Prince and Dormand [7].
                A=[0,                       0,      0,      0,                          0,                      0,                      0,                      0,                      0,                      0,                      0,                      0,  0;
                   1/18,                    0,      0,      0,                          0,                      0,                      0,                      0,                      0,                      0,                      0,                      0,  0;
                   1/48,                    1/16,   0,      0,                          0,                      0,                      0,                      0,                      0,                      0,                      0,                      0,  0;
                   1/32,                    0,      3/32,   0,                          0,                      0,                      0,                      0,                      0,                      0,                      0,                      0,  0;
                   5/16,                    0,      -75/64, 75/64,                      0,                      0,                      0,                      0,                      0,                      0,                      0,                      0,  0;     
                   3/80,                    0,      0,      3/16,                       3/20,                   0,                      0,                      0,                      0,                      0,                      0,                      0,  0;
                   29443841/614563906,      0,      0,      77736538/692538347,         -28693883/1125000000,   23124283/1800000000,    0,                      0,                      0,                      0,                      0,                      0,  0;
                   16016141/946692911,      0,      0,      61564180/158732637,         22789713/633445777,     545815736/2771057229,   -180193667/1043307555,  0,                      0,                      0,                      0,                      0,  0;
                   39632708/573591083,      0,      0,      -433636366/683701615,       -421739975/2616292301,  100302831/723423059,    790204164/839813087,    800635310/3783071287,   0,                      0,                      0,                      0,  0;
                   246121993/1340847787,    0,      0,      -37695042795/15268766246,   -309121744/1061227803,  -12992083/490766935,    6005943493/2108947869,  393006217/1396673457,   123872331/1001029789,   0,                      0,                      0,  0;
                   -1028468189/846180014,   0,      0,      8478235783/508512852,       1311729495/1432422823,  -10304129995/1701304382,-48777925059/3047939560,15336726248/1032824649, -45442868181/3398467696,3065993473/597172653,   0,                      0,  0;
                   185892177/718116043,     0,      0,      -3185094517/667107341,      -477755414/1098053517,  -703635378/230739211,   5731566787/1027545527,  5232866602/850066563,   -4093664535/808688257,  3962137247/1805957418,  65686358/487910083,     0,  0;
                   403863854/491063109,     0,      0,      -5068492393/434740067,      -411421997/543043805,   652783627/914296604,    11173962825/925320556,  -13158990841/6184727034,3936647629/1978049680,  -160528059/685178525,   248638103/1413531060,   0,  0];
                c=[0; 1/18; 1/12; 1/8; 5/16; 3/8; 59/400; 93/200; 5490023248/9719169821; 13/20; 1201146811/1299019798; 1; 1];
                bOrderMain=[14005451/335480064; 0; 0; 0; 0; -59238493/1068277825; 181606767/758867731; 561292985/797845732; -1041891430/1371343529; 760417239/1151165299; 118820643/751138087; -528747749/2220607170; 1/4];
                bOrderSubsid=[13451932/455176623; 0; 0; 0; 0; -808719846/976000145; 1757004468/5645159321; 656045339/265891186; -3867574721/1518517206; 465885868/322736535; 53011238/667516719;  2/45;0];
                orders=[8;7];
                isFSAL=false;
            case 1
                %Verner's RK8(9)16 algorithm.
                A=[0,                           0,      0,              0,                      0,                  0,                      0,                          0,                      0,                          0,                      0,              0,          0,          0,  0,          0;
                   1/12,                        0,      0,              0,                      0,                  0,                      0,                          0,                      0,                          0,                      0,              0,          0,          0,  0,          0;
                   1/27,                        2/27,   0,              0,                      0,                  0,                      0,                          0,                      0,                          0,                      0,              0,          0,          0,  0,          0;
                   1/24,                        0,      1/8,            0,                      0,                  0,                      0,                          0,                      0,                          0,                      0,              0,          0,          0,  0,          0;
                   vf(4,94)/375,                0,      vf(-94,-84)/125,vf(328,208)/375,        0,                  0,                      0,                          0,                      0,                          0,                      0,              0,          0,          0,  0,          0;
                   vf(9,-1)/150,                0,      0,              vf(312,32)/1425,        vf(69,29)/570,      0,                      0,                          0,                      0,                          0,                      0,              0,          0,          0,  0,          0;
                   vf(927,-347)/1250,           0,      0,              vf(-16248,7328)/9375,   vf(-489,179)/3750,  vf(14268,-5798)/9375,   0,                          0,                      0,                          0,                      0,              0,          0,          0,  0,          0;
                   2/27,                        0,      0,              0,                      0,                  vf(16,-1)/54,           vf(16,1)/54,                0,                      0,                          0,                      0,              0,          0,          0,  0,          0;
                   19/256,                      0,      0,              0,                      0,                  vf(118,-23)/512,        vf(118,23)/512,             -9/256,                 0,                          0,                      0,              0,          0,          0,  0,          0;
                   11/144,                      0,      0,              0,                      0,                  vf(266,-1)/864,         vf(266,1)/864,              -1/16,                  -8/27,                      0,                      0,              0,          0,          0,  0,          0;
                   vf(5034,-271)/61440,         0,      0,              0,                      0,                  0,                      vf(7859,-1626)/10240,       vf(-2232,813)/20480,    vf(-594,271)/960,           vf(657,-813)/5120,      0,              0,          0,          0,  0,          0;
                   vf(5996,-3794)/405,          0,      0,              0,                      0,                  vf(-4342,-338)/9,       vf(154922,-40458)/135,      vf(-4176,3794)/45,      vf(-340864,242816)/405,     vf(26304,-15176)/45,    -26624/81,      0,          0,          0,  0,          0;
                   vf(3793,2168)/103680,        0,      0,              0,                      0,                  vf(4042,2263)/13824,    vf(-231278,40717)/69120,    vf(7947,-2168)/11520,   vf(1048,-542)/405,          vf(-1383,542)/720,      2624/1053,      3/1664,     0,          0,  0,          0;
                   -137/1296,                   0,      0,              0,                      0,                  vf(5642,-337)/864,      vf(5642,337)/864,           -299/48,                184/81,                     -44/9,                  -5120/1053,     -11/468,    16/9,       0,  0,          0;
                   vf(33617,-2168)/518400,      0,      0,              0,                      0,                  vf(-3846,31)/13824,     vf(155338,-52807)/345600,   vf(-12537,2168)/57600,  vf(92,542)/2025,            vf(-1797,-542)/3600,    320/567,        -1/1920,    4/105,      0,  0,          0;
                   vf(-36487,-30352)/279600,    0,      0,              0,                      0,                  vf(-29666,-4499)/7456,  vf(2779182,-615973)/186400, vf(-94329,91056)/93200, vf(-232192,121408)/17475,   vf(101226,-22764)/5825, -169984/9087,   -87/30290,  492/1165,   0,  1260/233,   0];

                c=[0; 1/12; 1/9; 1/6; vf(2,2)/15; vf(6,1)/15; vf(6,-1)/15; 2/3; 1/2; 1/3; 1/4; 4/3; 5/6; 1; 1/6; 1];
                bOrderMain=[103/1680; 0; 0; 0; 0; 0; 0; -27/140;  76/105; -201/280; 1024/1365;  3/7280;   12/35;  9/280; 0;    0];
                bOrderSubsid=[23/525; 0; 0; 0; 0; 0; 0; 171/1400; 86/525; 93/280;   -2048/6825; -3/18200; 39/175; 0;     9/25; 233/4200];
                orders=[8;9];
                isFSAL=false;
            otherwise
                error('Unknown Solution Choice');
        end
    otherwise
        error('An invalid order was entered')
end

if(nargin==2)
    %If the function was just called with two input parameters of the form 
    %[orders,isFSAL]=RungeKStep(order,solutionChoice);
    %Then the first output should just be orders, the second isFSAL, and
    %there are no other outputs.
    xPredMain=orders;
    xPredSubsid=isFSAL;
else
    %Use the values from the selected Butcher table.
    k=funcValsFromButcherTableau(xVal,curT,f,deltaT,fCur,A,c);
    xPredMain=xVal+deltaT*sum(bsxfun(@times,bOrderMain,k'),1)';
    if(nargout>1)
        xPredSubsid=xVal+deltaT*sum(bsxfun(@times,bOrderSubsid,k'),1)';
    end
end

end

function val=vf(a,b)
%%VF This function implements 1-b*sqrt(6) and is used to keep the notation
%    simple when implementing Verner's order 8(9) algorithm from
%    J. H. Verner, "Explicit Runge-Kutta methods with estimates of the
%    local truncation error," SIAM Journal on Numerical Analysis, vol. 15,
%    no. 4, pp. 772-790, 1978.
%
%February 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=a+b*sqrt(6);
end

function k=funcValsFromButcherTableau(x,t,f,deltaT,fCur,A,c)
%%FUNCVALSFROMBUTCHERTABLEAU A Butcher tableau is a way of storing
%                         coefficients used in Runge-Kutta algorithms. This
%                         function computes the values of the function f at
%                         the points in terms of the state x and the time t
%                         as implied by the time-increment deltaT and the
%                         components A and c of the Butcher tableau. This
%                         function only works for explicit Runge-Kutta
%                         methods, so A is a lower-triangular square matrix
%                         of the coefficients that go into determining the
%                         offsets in the state x, where the first row and
%                         the last column of A are all zeros. c is the set
%                         of coefficients that are multiplied by the deltaT
%                         term; f is the differential equation function.
%                         fCur is f(x,t). The diagonal of A must be zero.
%
%The idea of arranging coefficients for Runge-Kutta methods into a table
%comes from
%J. Butcher, "On Runge-Kutta processes of high order," Journal of the
%Australian Mathematical Society, vol. 4, no. 2, pp. 179-194, May 1964.
%
%February 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    xDim=size(x,1);
    
    numVals=size(A,1);
    k=zeros(xDim,numVals);
    k(:,1)=fCur;
    for curVal=2:numVals
        fVal=f(x+deltaT*(k*A(curVal,:)'),t+c(curVal)*deltaT);
        k(:,curVal)=fVal;
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
