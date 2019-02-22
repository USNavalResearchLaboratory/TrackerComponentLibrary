function [xPredMain,xPredSubsid,g,orders,isFSAL,subsidType]=RungeKNystroemSStep(xVec,t,df,deltaT,dfCur,order,solutionChoice)
%%RUNGEKNYSTROMSSTEP Perform one step of a special (non-general) explicit
%                   Runge-Kutta Nyström method. Such methods are used to
%                   efficiently implement second order differential
%                   equations of the form d^2xdt^2=df(x,t), where the
%                   function df does not depend on dxdt. Runge-Kutta- 
%                   Nyström methods can be made to a higher order than
%                   general Runge-Kutta methods while having the same
%                   number of stages. All of the methods implemented here
%                   are embedded methods, meaning that a secondary result
%                   of a different order can also be returned along with
%                   the order of the second output. This second output
%                   can be used in an algorithm that has an adaptive step
%                   size (deltaT is automatically varied). Note that the
%                   function RungeKNystroemGStep handles the generalized
%                   problem where d2xdt2 depends on dxdt (the first
%                   derivative term).
%
%INPUTS: The function can either be called with two inputs as
%        [orders,isFSAL,subsidType]=RungeKNystroemSStep(order,solutionChoice)
%        to just get the orders of a particular general embedded Runge-
%        Kutta-Nyström method, to figure out whether g(:end) can be
%        passed in place of dfCur during a subsequent step, and to figure
%        out what type of information is contained in xPredSubsid, or it
%        can be run as
%        [xPredMain,xPredSubsid,g,orders,isFSAL,subsidType]=RungeKNystroemSStep(xVec,t,df,deltaT,dfCur,order,solutionChoice)
%        to actually do a full step. All of the inputs in the full step are:
%           xVec The value of the vector state over which integration is
%                being performed. The dimensionality must be a multiple of
%                two with the second half of the elements being the
%                derivatives of the first half of the elements (e.g.
%                velocity is the second half and position is the first
%                half).
%              t The time at which xVal is taken.
%             df df(x,t) returns the second dervative of derivative of x
%                with respect to time taken at time t. The output is half
%                the dimensionality of xVec. Note that the first input is
%                just the first half of xVec, not all of xVec as df is not
%                supposed to depend on the derivative components.
%         deltaT The size of the single (time) step over which the 
%                Runge-Kutta-Nyström integration is performed.
%          dfCur The value df(xVal,t). This is requested so that
%                methods that are FSAL can pass g(:,end) on subsequent
%                steps instead of having to perform a redundant
%                evaluation of df. If omitted or an empty matrix is
%                passed, the function df(xVal,curT) is evaluated to get
%                dfCur.
%          order The main integration order of the Runge-Kutta method.
%                If this parameter is omitted, then the default order of
%                5 is used. Order can range from 4 to 12.
% solutionChoice Different Runge-Kutta formulae exist for some orders. If
%                this parameter is provided, then the selected formula
%                for the given order is used. Otherwise the default
%                (solutionChoice=0) formula is used. The algorithms
%                chosen by pairs of (order, solutionChoice) are:
%                 (4,0) RKN4(3)4FM from Table 3 of [8].
%                 (4,1) RKN4(5)5F Table 15 from [1]
%                 (5,0) RKN5(4)4 from Table 14.6 in Chapter 14.4 of [14].
%                 (5,1) RKN5(4)5F from [3]
%                 (5,2) RKN5(6)7F from Table 12 in [2]
%                 (5,3) RKN5(6)7F from Table 12 in [1]
%                 (5,4) RKN5(4) from [13]
%                 (6,0) RKN6(4)6FM from Table 5 of [8]
%                 (6,1) RKN6(4)6FD from Table 2 of [7]
%                 (6,2) RKN6(7)9F from Table 9 of [2]
%                 (6,3) RKN6(5)8F from Table 9 in [1]
%                 (7,0) RKN7(6)9FT from Table III of [4]
%                 (7,1) RKN7(8)11F from Table 6 of [2]
%                 (7,2) RKN7(8)10F from Table 6 in [1]
%                 (8,0) RKN8(6)9FM from Table 1 of [9] with corrections
%                       from [11]
%                 (8,1) RKN8(9)14F from Table 2 in [2]
%                 (8,2) RKN8(9)12F from Table 3 in [1]
%                 (9,0) RKN9(8)12F from Table 2 in [6]
%                 (10,0) RKN10(11)18F from [5].
%                 (11,0) RKN11(12)21F from [0] as RKN 11(12). The actual
%                        coefficients are not given in the paper and were
%                        taken from [12]. (This is RKN 11(12) - 20, not not
%                        RKN 11(12) - 20 opt)
%                 (11,1) RKN11(12)21F, which is mentioned in [6] as
%                        RKN 11(12) - 20 opt, but whose coefficients are
%                        not directly provided. The coefficients were taken
%                        from [12].
%                 (11,2) RKN11(10)18F from [6]. These are the RKN11(10)-17
%                        coefficients (not the RKN11(10-17 opt
%                        coefficients) in Table 1. The actual coefficients
%                        are not given in the paper and were taken from
%                        [12].
%                 (11,3) RKN11(10)18F from [6]. These are the RKN11(10)-17
%                        opt coefficients in Table 1. The actual
%                        coefficients are not given in the paper and were
%                        taken from [12].
%                 (12,0) RKN12(10)17 from [9]. However, [9] just says that
%                        the authors could be contacted for the
%                        coefficients. The coefficients were thus taken (in
%                        rational form) from the comments in the function
%                        in [10], as the comments said that those were the
%                        coefficients provided by the authors of [9].
%
%OUTPUTS: If the function is run with two inputs as
%         [orders,isFSAL,subsidType]=RungeKNystroemSStep(order,solutionChoice)
%         Then there are just three outputs, orders, isFSAL and subsidType, 
%         described below.
%         Otherwise, all of the outputs are
%         xPredMain The entire state vector state propagated forward an
%                   interval of deltaT at the order of precision given by
%                   the input order. This is the main integration order.
%         xPredSubsid The subsidiary estimate. This is a different order
%                   than the main integration order can can be used in
%                   algorithms that adaptively adjust the stepsize.
%                   Depending on the combination of (order,
%                   solutionChoice), this is either a full state vector or
%                   just the first or second (first derivatives) half or
%                   the state vector. The output depends on the
%                   (order, solutionChoice) pair and the output subsidType
%                   indicates the type of subsidary output.
%                 g The values of the second derivatives df evaluated at
%                   various points as determined by the selected
%                   algorithm. These can sometimes be reused in
%                   interpolation routines to reduce the number of
%                   computations required. In all of the methods, g(:,1) is
%                   equal to df(xVal,t). In FSAL methods, g(:,end) is
%                   df(xPredMain,t+deltaT)
%            orders A 2X1 vector whereby order(1)=order, the input
%                   parameter, and orders(2) is the order of xPredSubsid.
%                   The value of orders(2) depends on the algorithm chosen
%                   by the combination of the order and solutionChoice
%                   inputs.
%            isFSAL Indicates whether the first evaluation of the f of the
%                   next step is equal to g(:,end). If so, this can be
%                   passed to the function on the next step rather than
%                   having to make an additional evaluation of f.
%        subsidType Indicates the nature of the output xPredSubsid.
%                   Possible values are:
%                   0 xPredSubsid is a complete state vector, having the
%                     same dimensionality as xVec.
%                   1 xPredSubsid does not contain any first derivatives
%                     and thus is half the dimensionality of xVec.
%                   2 xPredSubsid only contains first derivatives and is
%                     thus half the dimensionality of xVec.
%
%The general formula used for performing the step is discussed in all of
%the references.
%
%The formulae were referred to above using a notation such as RKN4(5)6. In
%general, the normenclature is
%RKNq(p)s[F]X
%where
%q is the order of the main integrating formula
%p is the order of the subsidiary formula
%s is the number of stages (derivative function evaluations)
%The presence of an F indicates that the function is first same as last
%(FSAL). That is, f evaluated at the next step is g(:,end) and thus
%g(:,end) can be passed for fCur for the next step. The absence of an F
%indicates that the function is not FSAL.
%The other values X depend on whether the author of the paper deriving the
%formula abides by general conventions/ added. Common last letters of M, S,
%C, and G, mean that the formula was derived..
%M minimizing the error norm
%S with an enlarged stability region
%C as a compromise between 'M' and 'S',
%Other last letters depend on the authors.
%Many authors do not include an F and just reduce the number of stages
%counted by 1. In the text above (except when differentiation between
%algorithms using the names given by the authors), the algorithms were
%specified with s being the total number of stages, and F indicating
%whether it is FSAL.
%
%REFERENCES:
%[0] S. Filippi and J. Gräf, "Ein Runge-Kutta-Nyström Formelpaar der
%    Ordnung 11(12) für Differentialgleichungen der Form y'' = f(x,y),"
%    Computing, vol. 34, no. 3, pp. 271-282, 1985.
%[1] E. Fehlberg, "Classical eighth- and lower order Runge-Kutta-Nystrom
%    formulas with stepsize control for special second-order differential
%    equations," National Aeronautics and Space Administration, Marshall Space
%    Flight Center, AL, Tech. Rep. NASA TR R-381, Mar. 1972.
%[2] E. Fehlberg, "Classical eighth- and lower-order Runge-Kuttta-Nyström
%    formulas with a new stepsize control procedure for special second-order
%    differential equations," National Aeronautics and Space Administration,
%    Marshall Space Flight Center, AL, Tech. Rep. NASA TR R-410, Jun. 1973.
%[3] D. G.Bettis, "A Runge-Kutta Nyström Algorithm," Celestial Mechanics, 
%    vol. 8, no. 2, pp. 229-233, Sep. 1973.
%[4] J. R. Dormand and P. J. Prince, "New Runge-Kutta algorithms for
%    numerical simulation in dynamical astronomy," Celestial Mechanics,
%    vol. 18, no. 3, pp. 223-232, Oct. 1978.
%[5] E. Fehlberg, S. Filippi, and J. Gräf "Ein Runge-Kutta-Nyström-
%    Formelpaar der Ordnung 10(11) für Differentialgleichungen der Form
%    y''= f(x,y)," Zeitschrift für Angewandte Mathematik und Mechanik, 
%    vol. 66, no. 7, pp. 265-270, 1986.
%[6] S. Filippi and J. Gräf "New Runge-Kutta-Nyström formula-pairs of
%     order 8(7), 9(8), 10(9) and 11(10) for differential equations of the
%     form y'' = f (x, y)," Journal of Computational and Applied
%     Mathematics, vol. 14, no. 3, pp. 361-370, Mar. 1986.
%[7] J. R. Dormand and P. J. Prince, "Runge-Kutta-Nystrom triples,"
%    Computers
%    and Mathematics with Applications, vol. 13, no. 12, pp. 937-949, 1987.
%[8] J. R. Dormand, M. E. A. El-Mikkawy, and P. J. Prince, "Families of
%    Runge-Kutta-Nystrom formulae," IMA Journal of Numerical Analysis, vol.
%    7, no. 2, pp. 235-250, Apr. 1987.
%[9] J. R. Dormand, M. E. A. El-Mikkawy, and P. J. Prince, "High-order
%    embedded Runge-Kutta-Nystrom formulae," IMA Journal of Numerical
%    Analysis, vol. 7, no. 2, pp. 423-430, Apr. 1987.
%[10] R. W. Brankin, I. Gladwell, J. R. Dormand, P. J. Prince, and W. L.
%     Seward, "Algorithm 670: A Runge-Kutta-Nyström code," ACM Transactions
%     on Mathematical Software, vol. 15, no. 1, pp. 31-40, Mar. 1989.
%[11] J. R. Dormand, M. E. A. El-Mikkawy, and P. J. Prince,
%     "Corrigendum: High-order embedded Runge-Kutta-Nystrom formulae," IMA
%     Journal of Numerical Analysis, vol. 11, no. 2, p. 297, Apr. 1991.
%[12] J. Gräf. (2006, 15 Jan.) Bäume, RKN-type methods by Filippi and Gräf
%     (www.josef-graef.de)" [Online]. Available:
%     http://www.josef-graef.de/baeume/rkni.html
%[13] M. Mohamad, N. Senu, M. Suleiman,and F. Ismail,"An embedded 5(4)
%     explicit Runge-Kutta-Nyström method with dissipation of high order,"
%     in Proceedings of the 20th National Symposium on Mathematical
%     Sciences, vol. 1522, Putrajaya, Malaysia, 18-20 Dec. 2012.
%[14] J. R. Dormand, Numerical Methods for Differential Equations.
%     Raton: CRC Press, 1996.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin==2)
    %If the function was just called as
    %[orders,isFSAL]=RungeKNystroemGStep(order,solutionChoice);
    %So that one can get the orders.
    order=xVec;
    solutionChoice=t;
else
    if(nargin<7)
        solutionChoice=0;
    end

    if(nargin<6)
        order=5;
    end

    xDim=length(xVec)/2;
    x=xVec(1:xDim,1);%The position
    dx=xVec((xDim+1):end);%The first derivative 
    
    if(nargin<5||isempty(dfCur))
        dfCur=df(x,t);
    end
end

switch(order)
    case 4
        switch(solutionChoice)
            case 0
                %RKN4(3)4FM from Table 3 of [8].
                A=[0,       0,          0,      0;
                   1/32,    0,          0,      0;
                   7/1000,  119/500,    0,      0;
                   1/14,    8/27,       25/189, 0];
                c=   [ 0;      1/4;    7/10;    1];
                bHat= [1/14;   8/27;   25/189;  0];
                bPHat=[1/14;   32/81;  250/567; 5/54];
                b=    [-7/150; 67/150; 3/20;    -1/20];
                bP=    [13/21; -20/27; 275/189; -1/3];
                isFSAL=true;
                orders=[4;3];
            case 1%RKN4(5)5F Table 15 from [1]
                A=[0,       0,      0,      0,      0;
                   1/18,    0,      0,      0,      0;
                   0,       2/9,    0,      0,      0;
                   1/3,     0,      1/6,    0,      0;
                   13/120,  3/10,   3/40,   1/60,   0];
                c=[0; 1/3; 2/3; 1; 1];
                bHat=[13/120; 3/10; 3/40; 1/60; 0];
                bPHat=[1/8; 3/8; 3/8; 1/8; 0];
                
                b=[bHat(1:3);0;bHat(4)];
                bP=[];
                isFSAL=true;
                orders=[4;5];
        otherwise
                error('Unknown Solution Choice');
        end
    case 5
        switch(solutionChoice)
            case 0
                %RKN5(4)4 from Table 14.6 in Chapter 14.4 of [14].
                A=[0,           0,          0,          0;
                   1/50,        0,          0,          0;
                   -3/250,      24/125,     0,          0;
                   3708/28561,  3525/28561, 4935/28561, 0];
               c=[0; 1/5; 3/5; 12/13];
               bHat=[23/432; 25/94; 125/756; 2197/142128];
               bPHat=[23/432; 125/376; 625/1512; 28561/142128];
               b=[3/1352; 1475/4056; 125/1352; 169/4056];
               bP=[];
               isFSAL=false;
               orders=[5;4];
            case 1
                %RKN5(4)5F from [3]
                A=[0,       0,      0,      0,      0,      0;
                   1/128,   0,      0,      0,      0,      0;
                   1/96,    1/48,   0,      0,      0,      0;
                   1/24,    0,      1/12,   0,      0,      0;
                   9/128,   0,      9/64,   9/128,  0,      0;
                   7/90,    0,      4/15,   1/15,   4/45,   0];
                c=[0;1/8;1/4;1/2;3/4;1];

                bHat=[7/90; 0; 4/15; 1/15; 4/45; 0];
                bPHat=[7/90; 0; 16/45; 2/15; 16/45; 7/90];
                b=[1/6; 0; 0; 1/3; 0; 0];
                bP=[0; 0; 2/3; -1/3; 2/3; 0];
                isFSAL=true;
                orders=[5;4];
            case 2%RKN5(6)7F from Table 12 in [2]
                s5=sqrt(5);
                
                gamma51=1/300;
                A=zeros(7,7);
                A(2,1)= 1/2;
                A(3,1)= (35+11*s5)/300;
                A(3,2)= (5+2*s5)/150;
                A(4,1)= (25-3*s5)/600;
                A(4,2)=-s5/300;
                A(4,3)= (13-5*s5)/120;
                A(5,1)= (25+3*s5)/600;
                A(5,2)= s5/300;
                A(5,3)= 0;
                A(5,4)= (13+5*s5)/120;
                A(6,1)= 1/12-gamma51;
                A(6,2)= gamma51;
                A(6,3)= 0;
                A(6,4)= (5+s5)/24+s5*gamma51;
                A(6,5)= (5-s5)/24-s5*gamma51;
                A(7,1)= 1/12;
                A(7,2)= 0;
                A(7,3)= 0;
                A(7,4)= (5+s5)/24;
                A(7,5)= (5-s5)/24;
                A(7,6)= 0;
              
                c=zeros(7,1);
                c(1)=0;
                c(2)=1;
                c(3)=(5+s5)/10;
                c(4)=(5-s5)/10;
                c(5)=(5+s5)/10;
                c(6)=1;
                c(7)=1;
                
                bHat=A(end,:)';
                bPHat=zeros(7,1);
                bPHat(1)=1/12;
                bPHat(2)=0;
                bPHat(3)=0;
                bPHat(4)=5/12;
                bPHat(5)=5/12;
                bPHat(6)=1/12;
                bPHat(7)=0;
                
                b=[];
                bP=[bPHat(1:5);0;bPHat(6)];
                isFSAL=true;
                orders=[5;6];
                
            case 3
                %RKN5(6)7F from Table 12 in [1]
                A=[0,           0,      0,      0,      0,          0,      0;
                   1/288,       0,      0,      0,      0,          0,      0;
                   1/216,       1/108,  0,      0,      0,          0,      0;
                   0,           0,      1/8,    0,      0,          0,      0;
                   16/125,      0,      4/125,  4/25,   0,          0,      0;
                   -247/1152,   0,      12/19,  7/432,  4375/65664, 0,      0;
                   11/240,      0,      108/475,8/45,   125/2736,   1/300,  0];
                
                c=[0; 1/12; 1/6; 1/2; 4/5; 1; 1];
                bHat=[11/240; 0; 108/475; 8/45; 125/2736; 1/300;0];
                bPHat=[1/24; 0; 27/95; 1/3; 125/456; 1/15;0];
                
                b=[bHat(1:5);0;bHat(6)];
                bP=[];
                isFSAL=true;
                orders=[5;6];
            case 4
                %RKN5(4) from [13]. Their use of b and bHat is opposite
                %that here.
                A=zeros(5,5);
                A(2,1)=0.02;
                A(3,1)=0.0231572886463073748489489742954;
                A(3,2)=0.0568427113536926251510510257046;
                A(4,1)=0.038935257798435815732212478526;
                A(4,2)=0.08724018833332508230340009768;
                A(4,3)=0.053824553868239101964387423802;
                A(5,1)=0.060247084262281349359489566347;
                A(5,2)=0.13367601834127956487747098916;
                A(5,3)=0.07871700038420723371804669813;
                A(5,4)=0.047359897012231852044992746376;

                c=[0; 0.2; 0.4; 0.6; 0.8];

                bHat=zeros(5,1);
                bHat(1)=0.131944444444444444444444444444;
                bHat(2)=-0.0555555555555555555555555555555;
                bHat(3)=0.5;
                bHat(4)=-0.194444444444444444444444444444;
                bHat(5)=0.118055555555555555555555555556;
        
                bPHat=zeros(5,1);
                bPHat(1)=0.131944444444444444444444444444;
                bPHat(2)=-0.0694444444444444444444444444444;
                bPHat(3)=0.833333333333333333333333333333;
                bPHat(4)=-0.486111111111111111111111111111;
                bPHat(5)=0.590277777777777777777777777778;

                b=zeros(5,1);
                b(1)=0.129813652373927633573147287493;
                b(2)=-0.021269231003914335250309119279;
                b(3)=0.409925778768177204312043632881;
                b(4)=-0.108631807353655270499082391230;
                b(5)=0.0901616072154647678642005901345;

                bP=zeros(5,1);
                bP(1)=0.131944444444444444444444444444;
                bP(2)=-0.0694444444444444444444444444444;
                bP(3)=0.833333333333333333333333333333;
                bP(4)=-0.486111111111111111111111111111;
                bP(5)=0.590277777777777777777777777778;
                isFSAL=false;
                orders=[5;4];
            otherwise
                error('Unknown Solution Choice');
        end
    case 6
        switch(solutionChoice)
            case 0
                %RKN6(4)6FM from Table 5 of [8] (Has a dense formula)
                A=[0,               0,              0,              0,              0,              0;
                   1/200,           0,              0,              0,              0,              0;
                   -1/2200,         1/22,           0,              0,              0,              0;
                   637/6600,        -7/110,         7/33,           0,              0,              0;
                   225437/1968750,  -30073/281250,  65569/281250,   -9367/984375,   0,              0;
                   151/2142,        5/116,          385/1368,       55/168,         -6250/28101,    0];
                c=    [0;           1/10;       3/10;          7/10;           17/25;         1];
                bHat= [151/2142;    5/116;      385/1368;      55/168;         -6250/28101;   0];
                bPHat=[151/2142;    25/522;     275/684;       275/252;        -78125/112404; 1/12];
                b=    [1349/157500; 7873/50000; 192199/900000; 521683/2100000; -16/125;       0];
                bP=   [1349/157500; 7873/45000; 27457/90000;   521683/630000;  -2/5;          1/12];
                isFSAL=true;
                orders=[6;4];
            case 1
                %RKN6(4)6FD from Table 2 of [7]
                R = sqrt(8581);
                A=zeros(6,6);
                A(2,1) = (26131-209*R)/810000;
                A(3,1) = (26131-209*R)/607500;
                A(3,2) = (26131-209*R)/303750;
                A(4,1) = (980403512254+7781688431*R)/11694469921875;
                A(4,2) = -(1262884486208+15385481287*R)/11694469921875;
                A(4,3) = (7166233891441+78694563299*R)/46777879687500;
                A(5,1) = -9*(329260+3181*R)/27040000;
                A(5,2) = 27*(35129+3331*R)/13520000;
                A(5,3) = -27*(554358343+31040327*R)/464060480000;
                A(5,4) = 153*(8555257-67973*R)/2745920000;
                A(6,1) = 329/4212;
                A(6,2) = 0;
                A(6,3) = (84119543+366727*R)/409622616;
                A(6,4) = (84119543-366727*R)/409622616;
                A(6,5) = 200/17901;
                c=[0; (209-R)/900; (209-R)/450; (209+R)/450; 9/10; 0];

                bHat=[329/4212; 0; (84119543+366727*R)/409622616; (84119543-366727*R)/409622616;  200/17901; 0];
                bPHat=[329/4212; 0; (389225579+96856*R)/1024056540; (389225579-96856*R)/1024056540; 2000/17901; 1/20];
                b=[(2701+23*R)/4563; -(9829+131*R)/9126; 5*(1798+17*R)/9126; 0; 0; 0];
                bP=[115/2106; 0; (84119543+366727*R)/256014135; (84119543-366727*R)/256014135; 6950/17901; -1/10];
                isFSAL=true;
                orders=[6;4];
            case 2%RKN6(7)9F from Table 9 of [2]
                s21=sqrt(21);
                
                gamma73=1/100;
                A=zeros(9,9);
                A(2,1)= 1/50;
                A(3,1)= 2/75;
                A(3,2)= 4/75;
                A(4,1)= 277/31104;
                A(4,2)= 95/15552;
                A(4,3)=-35/31104;
                A(5,1)= 5/192;
                A(5,2)= 0;
                A(5,3)= 25/1344;
                A(5,4)= 9/112;
                A(6,1)= (56-5*s21)/4116;
                A(6,2)= 0;
                A(6,3)=-25*(77-16*s21)/28812;
                A(6,4)= 9*(133-27*s21)/9604;
                A(6,5)= (441-95*s21)/4116;
                A(7,1)= (781+103*s21)/11760;
                A(7,2)= 0;
                A(7,3)=-25*(1369+599*s21)/279888;
                A(7,4)=-9*(2389+513*s21)/6860;
                A(7,5)= (315+127*s21)/2940;
                A(7,6)= 69*(225+49*s21)/4760;
                A(8,1)= 1/20-1/54*gamma73;
                A(8,2)= 0;
                A(8,3)= 0;
                A(8,4)= gamma73;
                A(8,5)= 8/45+(1/81)*gamma73;
                A(8,6)= 7*(7+s21)/360-7*(23+5*s21)/324*gamma73;
                A(8,7)= 7*(7-s21)/360-7*(23-5*s21)/324*gamma73;
                A(9,1)= 1/20;
                A(9,2)= 0;
                A(9,3)= 0;
                A(9,4)= 0;
                A(9,5)= 8/45;
                A(9,6)= 7*(7+s21)/360;
                A(9,7)= 7*(7-s21)/360;
                A(9,8)= 0;
              
                c=zeros(9,1);
                c(1)=0;
                c(2)=1/5;
                c(3)=2/5;
                c(4)=1/6;
                c(5)=1/2;
                c(6)=(7-s21)/14;
                c(7)=(7+s21)/14;
                c(8)=1;
                c(9)=1;
                
                bHat=A(end,:)';
                bPHat=zeros(7,1);
                bPHat(1)=1/20;
                bPHat(2)=0;
                bPHat(3)=0;
                bPHat(4)=0;
                bPHat(5)=16/45;
                bPHat(6)=49/180;
                bPHat(7)=49/180;
                bPHat(8)=1/20;
                bPHat(9)=0;
                
                b=[];
                bP=[bPHat(1:7);0;bPHat(8)];
                isFSAL=true;
                orders=[5;6];
            case 3
                %RKN6(5)8F from Table 9 in [1]
                A=[0,           0,      0,          0,          0,          0,          0,          0;
                    1/200,      0,      0,          0,          0,          0,          0,          0;
                    1/150,      1/75,   0,          0,          0,          0,          0,          0;
                    2/75,       0,      4/75,       0,          0,          0,          0,          0;
                    9/200,      0,      9/100,      9/200,      0,          0,          0,          0;
                    199/3600,   -19/150,47/120,     -119/1200,  89/900,     0,          0,          0;
                    -179/1824,  17/38,  0,          -37/152,    219/456,    -157/1824,  0,          0;
                    61/1008,    0,      475/2016,   25/504,     125/1008,   25/1008,    11/2016,    0];
                c=[0; 1/10; 1/5; 2/5; 3/5; 4/5; 1; 1];
                bHat=[61/1008; 0; 475/2016; 25/504; 125/1008; 25/1008; 11/2016; 0];
                bPHat=[19/288; 0; 25/96; 25/144; 25/144; 25/96; 19/288; 0];
                b=[bHat(1:6);0;bHat(7)];
                bP=[];
                isFSAL=true;
                orders=[6;5];
        otherwise
                error('Unknown Solution Choice');
        end
    case 7
        switch(solutionChoice)
            case 0
                %RKN7(6)9FT from Table III of [4]
                s21=sqrt(21);
                lambda=1/20;
                A=[0,                   0,                      0,                      0,                              0,                      0,                      0,              0,  0;
                   1/200,               0,                      0,                      0,                              0,                      0,                      0,              0,  0;
                   1/150,               1/75,                   0,                      0,                              0,                      0,                      0,              0,  0;
                   171/8192,            45/4096,                315/8192,               0,                              0,                      0,                      0,              0,  0;
                   5/288,               25/528,                 25/672,                 16/693,                         0,                      0,                      0,              0,  0;
                   (1003-205*s21)/12348,-25*(751-173*s21)/90552,25*(624-137*s21)/43218, -128*(361-79*s21)/237699,       (3411-745*s21)/24696,   0,                      0,              0,  0;
                   (793+187*s21)/12348, -25*(331+113*s21)/90552,25*(1044+247*s21)/43218,-128*(14885+3779*s21)/9745659,  (3327+797*s21)/24696,   -(581+127*s21)/1722,    0,              0,  0;
                   -(157-3*s21)/378,    25*(143-10*s21)/2772,   -25*(876+55*s21)/3969,  1280*(913+18*s21)/596673,       -(1353+26*s21)/2268,    7*(1777+377*s21)/4428,  7*(5-s21)/36,   0,  0;
                   1/20,                0,                      0,                      0,                              8/45,                   7*(7+s21)/360,          7*(7-s21)/360,  0,  0];


                c=    [0; 1/10; 1/5; 3/8; 1/2; (7-s21)/14; (7+s21)/14; 1; 1];
                bHat= [1/20; 0; 0; 0; 8/45;  7*(7+s21)/360; 7*(7-s21)/360;  0;      0];
                bPHat=[1/20; 0; 0; 0; 16/45; 49/180;        49/180;         1/20;   0];
                b=    [1/20; 0; 0; 0; 8/45;  7*(7+s21)/360; 7*(7-s21)/360; -lambda; lambda];
                bP=[];
                isFSAL=true;
                orders=[7;6];
            case 1%RKN7(8)11F from Table 6 of [2]
                s21=sqrt(21);

                gamma93=1/10;
                
                A=zeros(11,11);
                A(2,1)=25/882;
                
                A(3,1)=50/1323;
                A(3,2)=100/1323;
                
                A(4,1)=73/9800;
                A(4,2)=17/4900;
                A(4,3)=-1/1400;
                
                A(5,1)=25/1176;
                A(5,2)=0;
                A(5,3)=225/2744;
                A(5,4)=625/4116;
                
                A(6,1)=(661+73*s21)/42000;
                A(6,2)=0;
                A(6,3)=9*(637+141*s21)/98000;
                A(6,4)=(1127+226*s21)/11760;
                A(6,5)=(357+76*s21)/42000;
                
                A(7,1)=(531-73*s21)/19200;
                A(7,2)=0;
                A(7,3)=0;
                A(7,4)=-7*(105-73*s21)/15360;
                A(7,5)=7*(399+73*s21)/230400;
                A(7,6)=73*(21-5*s21)/11520;
                
                A(8,1)=(789-89*s21)/29400;
                A(8,2)=0;
                A(8,3)=0;
                A(8,4)=-(105-17*s21)/8400;
                A(8,5)=-(77+103*s21)/25200;
                A(8,6)=13*(147-29*s21)/17640;
                A(8,7)=2*(325-51*s21)/11025;
                A(9,1)=(579+235*s21)/29400;
                A(9,2)=0;
                A(9,3)=0;
                A(9,4)=(315-191*s21)/4200;
                A(9,5)=-(161-15*s21)/12600;
                A(9,6)=-3*(21-5*s21)/490;
                A(9,7)=2*(395+51*s21)/11025;
                A(9,8)=(43+9*s21)/280;
                A(10,1)=1/20-(4/49)*gamma93;
                A(10,2)=0;
                A(10,3)=0;
                A(10,4)=gamma93;
                A(10,5)=-1/9*gamma93;
                A(10,6)=0;
                A(10,7)=8/45+(64/441)*gamma93;
                A(10,8)=7*(7+s21)/360-10*(21+5*s21)/441*gamma93;
                A(10,9)=7*(7-s21)/360-10*(21-5*s21)/441*gamma93;
                A(11,1)=1/20;
                A(11,2)=0;
                A(11,3)=0;
                A(11,4)=0;
                A(11,5)=0;
                A(11,6)=0;
                A(11,7)=8/45;
                A(11,8)=7*(7+s21)/360;
                A(11,9)=7*(7-s21)/360;
                A(11,10)=0;
                
                c=[0; 5/21; 10/21; 1/7; 5/7; (7+s21)/14; 1/2; (7-s21)/14; (7+s21)/14; 1; 1];
                bHat=A(end,:)';
                bPHat=[1/20; 0; 0; 0; 0; 0; 16/45; 49/180; 49/180; 1/20; 0];
                 
                b=[];
                bP=[bPHat(1:9);0;bPHat(10)];
                isFSAL=true;
                orders=[7;8];
            case 2
                %RKN7(8)10F from Table 6 in [1]
                A=[0,                           0,              0,                      0,                  0,                          0,                      0,                      0,                      0,      0;
                   361/281250,                  0,              0,                      0,                  0,                          0,                      0,                      0,                      0,      0;
                   10437/7600,                  -343/304,       0,                      0,                  0,                          0,                      0,                      0,                      0,      0;
                   547/319200,                  1125/342304,    -1/4729200,             0,                  0,                          0,                      0,                      0,                      0,      0;
                   74/9975,                     -1125/791578,   -1/157640,              311/22200,          0,                          0,                      0,                      0,                      0,      0;
                   1028/29925,                  -6375/1583156,  -55/319221,             -13/1665,           467/8100,                   0,                      0,                      0,                      0,      0;
                   148349/19254600,             6375/1583156    0,                      1299964/14060925,   4783/253350,                173101/3040200,         0,                      0,                      0,      0;
                   116719112/18953746875,       1125/791578,    1680359/2992696875,     51962281/585871875, 104130509/855056250,        1995658/47503125,       15029/253125,           0,                      0,      0;
                   604055892451/4935014784000,  0,              -206360699/115664409000,0,                  32963694031/528751584000,   9676095011/39166784000, 1641775937/176250528000,2851784579/47000140800, 0,      0;
                   67/2016,                     0,              0,                      440/3969,           25/252,                     425/3024,               5/72,                   625/14112,              11/4536,0];
                c=[0;19/375; -7/10; 1/10; 1/5; 2/5; 3/5; 4/5; 1; 1];
                bHat=[67/2016; 0; 0; 440/3969; 25/252; 425/3024; 5/72; 625/14112; 11/4536;0];
                bPHat=[23/2016; 0; 0; 880/3969; -25/2016; 1075/3024; 65/1008; 4225/14112; 1087/18144;0];
                b=[bHat(1:8);0;bHat(9)];
                bP=[];
                isFSAL=true;
                orders=[7;8];
            otherwise
                error('Unknown Solution Choice');
        end
    case 8
        switch(solutionChoice)
            case 0
                %RKN8(6)9FM from Table 1 of [9] with corrections from [11].
                A=[0,                       0,                  0,                      0,                  0,          0,          0,          0,  0;
                   1/800,                   0,                  0,                      0,                  0,          0,          0,          0,  0;
                   1/600,                   1/300,              0,                      0,                  0,          0,          0,          0,  0;
                   9/200,                   -9/100,             9/100,                  0,                  0,          0,          0,          0,  0;
                   -66701/197352,           28325/32892,        -2665/5482,             2170/24669,         0,          0,          0,          0,  0;
                   227015747/304251000,     -54897451/30425100, 12942349/10141700,      -9499/304251,       539/9250,   0,          0,          0,  0;
                   -1131891597/901789000,   41964921/12882700,  -6663147/3220675,       270954/644135,      -108/5875,  114/1645,   0,          0,  0;
                   13836959/3667458,        -17731450/1833729,  1063919505/156478208,   -33213845/39119552, 13335/28544,-705/14272, 1645/57088, 0,  0;
                   223/7938,                0,                  1175/8064,              925/6048,           41/448,     925/14112,  1175/72576, 0,  0];
                c=[0; 1/20; 1/10; 3/10; 1/2; 7/10; 9/10; 1; 1];
                bHat= [223/7938;          0; 1175/8064;        925/6048;          41/448;           925/14112;         1175/72576;        0;                 0];
                bPHat=[223/7938;          0; 5875/36288;       4625/21168;        41/224;           4625/21168;        5875/36288;        223/7938;          0];
                b=    [7987313/109941300; 0; 1610737/44674560; 10023263/33505920; -497221/12409600; 10023263/78180480; 1610737/402071040; 0;                 0];
                bP=   [7987313/109941300; 0; 1610737/40207104; 10023263/23454144; -497221/6204800;  10023263/23454144; 1610737/40207104; -4251941/54970650; 3/20];
                isFSAL=true;
                orders=[8;6];
            case 1%RKN8(9)14F from Table 2 in [2]
                s15=sqrt(15);
                
                gamma126=5/2;
                A=zeros(14,14);
                A(2,1)=1/18;
                
                A(3,1)=2/27;
                A(3,2)=4/27;
                
                A(4,1)=7/128;
                A(4,2)=5/64;
                A(4,3)=-1/128;
                
                A(5,1)=89/3240;
                A(5,2)=31/540;
                A(5,3)=11/1080;
                A(5,4)=-16/405;
                
                A(6,1)=11/120;
                A(6,2)=0;
                A(6,3)=9/40;
                A(6,4)=-4/15;
                A(6,5)=9/20;
                
                A(7,1)=33259/7085880;
                A(7,2)=0;
                A(7,3)=343/157464;
                A(7,4)=-4708/885735;
                A(7,5)=1879/393660;
                A(7,6)=-139/885735;
                
                A(8,1)=29/1920;
                A(8,2)=0;
                A(8,3)=0;
                A(8,4)=0;
                A(8,5)=99/2560;
                A(8,6)=1/30720;
                A(8,7)=729/10240;
                
                A(9,1)=13/1215;
                A(9,2)=0;
                A(9,3)=0;
                A(9,4)=0;
                A(9,5)=1/144;
                A(9,6)=1/77760;
                A(9,7)=87/2240;
                A(9,8)=-8/8505;
                
                A(10,1)=22/1215;
                A(10,2)=0;
                A(10,3)=0;
                A(10,4)=0;
                A(10,5)=0;
                A(10,6)=1/4860;
                A(10,7)=3/28;
                A(10,8)=256/8505;
                A(10,9)=1/15;
                
                A(11,1)=(7561-1454*s15)/420000;
                A(11,2)=0;
                A(11,3)=0;
                A(11,4)=0;
                A(11,5)=-9*(1373+45*s15)/800000;
                A(11,6)=(379-145*s15)/1334000;
                A(11,7)=729*(6997-1791*s15)/78400000;
                A(11,8)=(999-473*s15)/183750;
                A(11,9)=27*(19407-3865*s15)/5600000;
                A(11,10)=297*(78-19*s15)/700000;
                
                A(12,1)=(12647-2413*s15)/840000;
                A(12,2)=0;
                A(12,3)=0;
                A(12,4)=0;
                A(12,5)=-9*(1373-45*s15)/800000;
                A(12,6)=-(29-61*s15)/336000;
                A(12,7)=729*(14743+3789*s15)/19600000;
                A(12,8)=(999+143*s15)/183750;
                A(12,9)=27*(20157+4315*s15)/5600000;
                A(12,10)=27*(1641+463*s15)/1400000;
                A(12,11)=-(1/56)*(27+7*s15);
                
                A(13,1)=9/280-(35/6561)*gamma126;
                A(13,2)=0;
                A(13,3)=0;
                A(13,4)=0;
                A(13,5)=0;
                A(13,6)=0;
                A(13,7)=gamma126;
                A(13,8)=16/315-(160/19683)*gamma126;
                A(13,9)=243/1540+(35/2673)*gamma126;
                A(13,10)=243/3080+(7/2673)*gamma126;
                A(13,11)=25*(5+s15)/1386-3500*(31+8*s15)/216513*gamma126;
                A(13,12)=25*(5-s15)/1386-3500*(31-8*s15)/216513*gamma126;
                
                A(14,1)=9/280;
                A(14,2)=0;
                A(14,3)=0;
                A(14,4)=0;
                A(14,5)=0;
                A(14,6)=0;
                A(14,7)=0;
                A(14,8)=16/315;
                A(14,9)=243/1540;
                A(14,10)=243/3080;
                A(14,11)=25*(5+s15)/1386;
                A(14,12)=25*(5-s15)/1386;
                A(14,13)=0;
                
                c=[0; 1/3; 2/3; 1/2; 1/3; 1; 1/9; 1/2; 1/3; 2/3; (5-s15)/10; (5+s15)/10; 1; 1];

                bHat=A(end,:)';
                bPHat=[9/280; 0; 0; 0; 0; 0; 0; 32/315; 729/3080; 729/3080; 125/693; 125/693; 9/280; 0];
                b=[];
                bP=[bPHat(1:12);0;bPHat(13)];
                isFSAL=true;
                orders=[8;9];
            case 2
                %RKN8(9)12F from Table 3 in [1]
                A=[0,                       0,          0,                  0,              0,                      0,                      0,                      0,                      0,              0,              0,      0;
                   49/12800,                0,          0,                  0,              0,                      0,                      0,                      0,                      0,              0,              0,      0;
                   49/96900,                49/4800,    0,                  0,              0,                      0,                      0,                      0,                      0,              0,              0,      0;
                   16825/381024,            -625/11907, 18125/190512,       0,              0,                      0,                      0,                      0,                      0,              0,              0,      0;
                   23/840,                  0,          50/609,             9/580,          0,                      0,                      0,                      0,                      0,              0,              0,      0;
                   533/68040,               0,          5050/641277,        -19/5220,       23/12636,               0,                      0,                      0,                      0,              0,              0,      0;
                   -4469/85050,             0,          -2384000/641277,    3869/19575,     -1451/15795,            502/135,                0,                      0,                      0,              0,              0,      0;
                   694/10125,               0,          0,                  -5504/10125,    424/2025,               -104/2025,              364/675,                0,                      0,              0,              0,      0;
                   30203/691200,            0,          0,                  0,              9797/172800,            79391/518400,           20609/345600,           70609/2073600,          0,              0,              0,      0;
                   1040381917/14863564800,  0,          548042275/109444608,242737/5345280, 569927617/6900940800,   -2559686731/530841600,  -127250389/353894400,   -53056229/2123366400,   23/5120,        0,              0,      0;
                   -33213637/179088000,     0,          604400/324597,      63826/445875,   0,                      -6399863/2558400,       110723/511680,          559511/35817600,        372449/7675200, 756604/839475,  0,      0;
                   121/4200,                0,          0,                  0,              43/525,                 33/350,                 17/140,                 3/56,                   31/1050,        512/5775,       1/550,  0];
                
                c=[0; 7/80; 7/40; 5/12; 1/2; 1/6; 1/3; 2/3; 5/6; 1/12; 1; 1];
                bHat=[121/4200; 0; 0; 0; 43/525; 33/350; 17/140; 3/56; 31/1050; 512/5775; 1/550; 0];
                bPHat=[41/840; 0; 0; 0; 34/105; 9/35; 9/280; 9/280; 9/35; 0; 41/840; 0];
                b=[bHat(1:10);0;bHat(11)];
                bP=[];
                isFSAL=true;
                orders=[8;9];
            otherwise
                error('Unknown Solution Choice');
        end
    case 9
        switch(solutionChoice)
            case 0
                %RKN9(8)12F from Table 2 in [6].
                A=zeros(12,12);
                A(2,1)=  0.1998028426208654875176403333e-4;
                A(3,1)=  0.2664037901611539833568537778e-4;
                A(3,2)=  0.5328075803223079667137075556e-4;
                A(4,1)=  0.2865365237083083248442544836e2;
                A(4,2)= -0.5904090661728807966777976253e2;
                A(4,3)=  0.3047405980201280273890986972e2;
                A(5,1)=  0.1065903941831243410656110659e-3;
                A(5,2)=  0;
                A(5,3)=  0.2171588709672456801997590475e-3;
                A(5,4)=  0.1307385389198734629886539400e-8;
                A(6,1)=  0.3030980653858376916419326961e-1;
                A(6,2)=  0;
                A(6,3)= -0.7126134845184682180630126302e-1;
                A(6,4)=  0.1568175392509382439415139758e-4;
                A(6,5)=  0.4928662616150566473771384202e-1;
                A(7,1)= -0.1477808541165541247980136312;
                A(7,2)=  0;
                A(7,3)=  0.3682723455773522839746690025;
                A(7,4)=  0.4848245662307982600517077975e-3;
                A(7,5)= -0.2155027916119219160929380291;
                A(7,6)=  0.3865397359925407390123095003;
                A(8,1)= -0.1089954028613558927819286465e-1;
                A(8,2)=  0;
                A(8,3)=  0;
                A(8,4)=  0.7635032467467202877246657217e-2;
                A(8,5)=  0.4547118895489745578797952498e-1;
                A(8,6)=  0.4826544540148205388198081045e-1;
                A(8,7)=  0.3452787346228887673098587201e-1;
                A(9,1)=  0.5456430942591812422112950282;
                A(9,2)=  0;
                A(9,3)= -0.8991472409701229621282353299;
                A(9,4)= -0.1676560070320166624537773864e-1;
                A(9,5)=  0.3402552136626613562636251327;
                A(9,6)=  0.1681846222333179210929658770;
                A(9,7)=  0.4235220058882386425502896119e-1;
                A(9,8)=  0.6652778464370135979569806931e-1;
                A(10,1)=-0.1113406602941552830476082277e1;
                A(10,2)= 0;
                A(10,3)= 0.7546972259090522109972731093;
                A(10,4)= 0.1118600256491523208841810745e-1;
                A(10,5)= 0.7937353844859527762086562061;
                A(10,6)=-0.4656855788249036634796925437;
                A(10,7)= 0.3836619386801684522200183160;
                A(10,8)=-0.2737580302105589161215079526e-1;
                A(10,9)= 0.6449737368017565619020116575e-1;
                A(11,1)= 0.3114812407323702239935577607e1;
                A(11,2)= 0;
                A(11,3)= 0;
                A(11,4)= 0;
                A(11,5)=-0.4620858770988909664725067967e1;
                A(11,6)= 0.2485087828585910565589259449e1;
                A(11,7)=-0.9881138857312002169169839148;
                A(11,8)= 0.5265010973962881009401492752;
                A(11,9)=-0.4664170309918856325664075129e-1;
                A(11,10)=0.2921302651339753843370630214e-1;
                A(12,1)=-0.3518240655555007629763408588e-2;
                A(12,2)= 0;
                A(12,3)= 0;
                A(12,4)= 0;
                A(12,5)= 0.6879703880602338697982077081e-1;
                A(12,6)= 0.1177879603832509115294976608;
                A(12,7)= 0.1377362395616285338725155827;
                A(12,8)= 0.1008501654461988983215135519;
                A(12,9)= 0.6135631323226062775059274661e-1;
                A(12,10)=0.1699052322619264917582309574e-1;
                A(12,11)=0;

                c(1)=  0;
                c(2)=  0.6321437219823755989424902073e-2;
                c(3)=  0.1264287443964751197884980415e-1;
                c(4)=  0.4166666666666666666666666666;
                c(5)=  0.2544604380000000000000000000e-1;
                c(6)=  0.1292344072000000000000000000;
                c(7)=  0.2970774243000000000000000000;
                c(8)=  0.5;
                c(9)=  0.7029225757000000000000000000;
                c(10)=  0.8958905519456625480019182180;
                c(11)= 0.1e-1;
                c(12)= 0.1e-1;
        
                bHat=A(end,:)';
                
                bPHat=zeros(12,1);
                bPHat(1)= -0.3518240655555007629763408588e-2;
                bPHat(2)=  0;
                bPHat(3)=  0;
                bPHat(4)=  0;
                bPHat(5)=  0.7059336055058270664872682481e-1;
                bPHat(6)=  0.1352694242367759659250257652;
                bPHat(7)=  0.1959479526240353954141973686;
                bPHat(8)=  0.2017003308923977966430271039;
                bPHat(9)=  0.2065330725713466061937737980;
                bPHat(10)= 0.1631986677839733565676051585;
                bPHat(11)= 0.3027543199644318023740738961e-1;
                bPHat(12)= 0;
                lambda=1e-5;
                b=[bHat(1:10);-lambda;lambda];
                bP=[];
                isFSAL=true;
                orders=[9;8];
            otherwise
                error('Unknown Solution Choice');
        end
    case 10
        switch(solutionChoice)
            case 0
                %RKN10(11)18F from [5].
                A=zeros(18,18);
                A(2,1)=1/1152;

                A(3,1)=1/864;
                A(3,2)=1/432;

                A(4,1)=13/1250;
                A(4,2)=-8/625;
                A(4,3)=14/625;

                A(5,1)=5/768;
                A(5,2)=0;
                A(5,3)=9/448;
                A(5,4)=25/5376;

                A(6,1)=-3039616/4084101;
                A(6,2)=0;
                A(6,3)=2539520/1058841;
                A(6,4)=-122828800/28588707;
                A(6,5)=11976704/4084101;

                A(7,1)=647/77760;
                A(7,2)=0;
                A(7,3)=1696/53865;
                A(7,4)=3125/401436;
                A(7,5)=416/52245;
                A(7,6)=2401/1249421760;

                A(8,1)=241283/75937500;
                A(8,2)=0;
                A(8,3)=6036992/841640625;
                A(8,4)=-10282/2508975;
                A(8,5)=2571776/816328125;
                A(8,6)=4895639/3660415312500;
                A(8,7)=-9778/18984375;

                A(9,1)=1183/61440;
                A(9,2)=0;
                A(9,3)=0;
                A(9,4)=0;
                A(9,5)=-76/4515;
                A(9,6)=2401/29061120;
                A(9,7)=83/1920;
                A(9,8)=3125/39424;

                A(10,1)=38725/995328;
                A(10,2)=0;
                A(10,3)=0;
                A(10,4)=0;
                A(10,5)=8500/73143;
                A(10,6)=66327625/5178691584;
                A(10,7)=-875/31104;
                A(10,8)=390625/3902976;
                A(10,9)=12625/117612;

                A(11,1)=1220213/530841600;
                A(11,2)=0;
                A(11,3)=0;
                A(11,4)=0;
                A(11,5)=0;
                A(11,6)=-87153899/346851901440;
                A(11,7)=-88751/89579520;
                A(11,8)=4980625/2890432512;
                A(11,9)=56803/100362240;
                A(11,10)=4721/37324800;

                A(12,1)=7363/1075200;
                A(12,2)=0;
                A(12,3)=0;
                A(12,4)=0;
                A(12,5)=0;
                A(12,6)=0;
                A(12,7)=13/10752;
                A(12,8)=23125/3311616;
                A(12,9)=-1213/5913600;
                A(12,10)=71/11289600;
                A(12,11)=827/50400;

                A(13,1)=28785017/18748422000;
                A(13,2)=0;
                A(13,3)=0;
                A(13,4)=0;
                A(13,5)=32358656/916116075;
                A(13,6)=57712484053/25083781628400;
                A(13,7)=23367299/8436789900;
                A(13,8)=-53812922875/1333912728456;
                A(13,9)=-131961188/28356988275;
                A(13,10)=-201142364/172251127125;
                A(13,11)=5438368/78118425;
                A(13,12)=-9801248/997548615;

                A(14,1)=-14302315753/247479170400;
                A(14,2)=0;
                A(14,3)=0;
                A(14,4)=0;
                A(14,5)=-35512576/549669645;
                A(14,6)=-264924469621/30100537954080;
                A(14,7)=0;
                A(14,8)=-1600322764375/2000869092684;
                A(14,9)=7825742534/85070964825;
                A(14,10)=391181138/103350676275;
                A(14,11)=1348797728/2109197475;
                A(14,12)=31510496/41757849;
                A(14,13)=-60857638/180788355;

                A(15,1)=751967018007/17275912192000;
                A(15,2)=0;
                A(15,3)=0;
                A(15,4)=0;
                A(15,5)=56879/1256675;
                A(15,6)=458244129931/187916176588800;
                A(15,7)=-156110412347/971770060800;
                A(15,8)=90915884904125/307286650159104;
                A(15,9)=76173973209/3266227148800;
                A(15,10)=-74357445229/59520916224000;
                A(15,11)=-753288691/6642959400;
                A(15,12)=-105043369/516674620;
                A(15,13)=1484948603/4417136640;
                A(15,14)=200345391/15116423168;

                A(16,1)=1711560397019/171738272563200;
                A(16,2)=0;
                A(16,3)=0;
                A(16,4)=0;
                A(16,5)=-9704761/549669645;
                A(16,6)=-1231127701027/5604172884541440;
                A(16,7)=35652964503/268341050880;
                A(16,8)=-2445542604468125/15273543270408192;
                A(16,9)=401638932149833/3181988181335040;
                A(16,10)=2657236327733/394461344793600;
                A(16,11)=75523039283/352197629280;
                A(16,12)=10005789083/34241436180;
                A(16,13)=-2064099266483/9660277831680;
                A(16,14)=5723201/1110896640;
                A(16,15)=353/14580;

                A(17,1)=2463899/150843000;
                A(17,2)=0;
                A(17,3)=0;
                A(17,4)=0;
                A(17,5)=0;
                A(17,6)=0;
                A(17,7)=0;
                A(17,8)=0;
                A(17,9)=96676/1714125;
                A(17,10)=-306176/11998875;
                A(17,11)=1769048/11998875;
                A(17,12)=167336/3999625;
                A(17,13)=333323/2285500;
                A(17,14)=1108557/31997000;
                A(17,15)=5304/81625;
                A(17,16)=648/35915;

                A(18,1)=35873/1524600;
                A(18,2)=0;
                A(18,3)=0;
                A(18,4)=0;
                A(18,5)=0;
                A(18,6)=0;
                A(18,7)=0;
                A(18,8)=0;
                A(18,9)=634/5775;
                A(18,10)=6/1925;
                A(18,11)=2616/21175;
                A(18,12)=248/1925;
                A(18,13)=153/3850;
                A(18,14)=351/15400;
                A(18,15)=664/17325;
                A(18,16)=216/21175;
                A(18,17)=1/7260;

                c=[0; 1/24; 1/12; 1/5; 1/4; 16/21; 1/3; 2/15; 1/2; 5/6; 1/12; 1/4; 1/3; 2/3; 3/4; 11/12; 1; 1];
                bHat=A(end,:)';
                bPHat=[653/27720; 0; 0; 0; 0; 0; 0; 0; 584/2625; 0; 1296/9625; 272/1575; 81/1400; 81/1400; 272/1575; 1296/9625; 653/27720; 0];
                b=[bHat(1:16);0;bHat(17)];
                bP=[];
                isFSAL=true;
                orders=[10;11];
            otherwise
                error('Unknown Solution Choice');
        end
    case 11
        switch(solutionChoice)
            case 0%RKN11(12)21F from [0] as RKN 11(12). The actual
                %coefficients are not given in the paper and were taken
                %from [12]. (This is RKN 11(12) - 20, not not RKN 11(12) -
                %20 opt)

                A=zeros(21,21);
                A(2,1)=0.5425347222222222222222222225e-2;
                A(3,1)=0.7233796296296296296296296301e-2;
                A(3,2)=0.1446759259259259259259259258e-1;
                A(4,1)=0.1033340000000000000000000000;
                A(4,2)=-0.1703680000000000000000000001;
                A(4,3)=0.2182840000000000000000000001;
                A(5,1)=0.3257575757575757575757575757e-1;
                A(5,2)=0;
                A(5,3)=0.8780487804878048780487804878e-1;
                A(5,4)=0.4619364375461936437546193648e-2;
                A(6,1)=0.4748507367563208000018089973e-1;
                A(6,2)=0;
                A(6,3)=0.1717457745568296060188724074;
                A(6,4)=0.5999520621410557300028385132e-1;
                A(6,5)=0.6253892679998000604987985799e-2;
                A(7,1)=0.6897742477766473591191027092e-2;
                A(7,2)=0;
                A(7,3)=0.6890192907847598728171600460e-2;
                A(7,4)=0.1163520086778733422170956371e-1;
                A(7,5)=-0.1344647159854360748127587603e-1;
                A(7,6)=-0.7266646548577990597963152304e-3;
                A(8,1)=0.1349737761438706440043676626e-2;
                A(8,2)=0;
                A(8,3)=0;
                A(8,4)=0.6075294762105980384386224708e-3;
                A(8,5)=-0.6793734869158527204348147999e-3;
                A(8,6)=-0.4099894989642249543680965575e-4;
                A(8,7)=0.4639442565473848303431810095e-3;
                A(9,1)=0.3523183546767263585786907763e-2;
                A(9,2)=0;
                A(9,3)=0;
                A(9,4)=0;
                A(9,5)=-0.2584392204058368115482471755e-5;
                A(9,6)=0.2758330832880140022000070680e-6;
                A(9,7)=0.4071837479676998987495697937e-2;
                A(9,8)=0.1240728753267650778083067676e-1;
                A(10,1)=-0.1340488944178098066953190204e-1;
                A(10,2)=0;
                A(10,3)=0;
                A(10,4)=0;
                A(10,5)=0.4748167290240606363971037112e-1;
                A(10,6)=0.6978075811258520497345866824e-3;
                A(10,7)=-0.2104272720250956294017709271;
                A(10,8)=0.1479810168188958303177052879;
                A(10,9)=0.2638347286869889791697278318;
                A(11,1)=0.4338886679167641664845535142e-2;
                A(11,2)=0;
                A(11,3)=0;
                A(11,4)=0;
                A(11,5)=0;
                A(11,6)=-0.1244276150937331764200609062e-5;
                A(11,7)=0.7318827258766195401581730488e-2;
                A(11,8)=0.1702903672564807515932491560e-1;
                A(11,9)=0.2562063489400243489046551726e-2;
                A(11,10)=0.2430123168781616965467654675e-5;
                A(12,1)=0.1213181442859181818082415173e-2;
                A(12,2)=0;
                A(12,3)=0;
                A(12,4)=0;
                A(12,5)=0;
                A(12,6)=-0.8052083719028318452003960397e-6;
                A(12,7)=-0.1270548819575905357299570659e-2;
                A(12,8)=0.2106288690010669033987958908e-2;
                A(12,9)=0.1062863271197933834077206495e-2;
                A(12,10)=0.1608786873238876934044295570e-5;
                A(12,11)=-0.3000881629932153739368538180e-3;
                A(13,1)=0.4791457562121135376714969290e-2;
                A(13,2)=0;
                A(13,3)=0;
                A(13,4)=0;
                A(13,5)=0;
                A(13,6)=0;
                A(13,7)=0.1997283495783110806313165151e-1;
                A(13,8)=0.3263466633670545485582605094e-1;
                A(13,9)=-0.1997857874597377833235256280e-2;
                A(13,10)=0.5702678234802111502863498282e-7;
                A(13,11)=0.5809882374787801958204563035e-2;
                A(13,12)=-0.1316104038363047044175700713e-1;
                A(14,1)=0.2680759937173601150132530305e-2;
                A(14,2)=0;
                A(14,2)=0;
                A(14,4)=0;
                A(14,5)=0;
                A(14,6)=0;
                A(14,7)=0;
                A(14,8)=0.9931732719629747807246879088e-2;
                A(14,9)=0.1935996027139212373939522857e-2;
                A(14,10)=-0.1715138004351496028992556713e-6;
                A(14,11)=-0.1146222506237529255095036749e-2;
                A(14,12)=-0.8065227820301908673948257037e-3;
                A(14,13)=0.2044281181255939407738294573e-3;
                A(15,1)=0.7450337675216960918357469139e-2;
                A(15,2)=0;
                A(15,3)=0;
                A(15,4)=0;
                A(15,5)=0;
                A(15,6)=0;
                A(15,7)=0;
                A(15,8)=0;
                A(15,9)=0.2826065109541195073618795087;
                A(15,10)=0.2441581049489756027618007528e-4;
                A(15,11)=-0.1686732143404799746723882008;
                A(15,12)=0.6570612744839619885192244320e-1;
                A(15,13)=0.6623464602149089266097689541e-1;
                A(15,14)=-0.1475488235692384826810242958;
                A(16,1)=0.3509222722230964288525276387e-1;
                A(16,2)=0;
                A(16,3)=0;
                A(16,4)=0;
                A(16,5)=0;
                A(16,6)=0.1838947430224798509482261328e-4;
                A(16,7)=-0.8175298816697407045626030708e-1;
                A(16,8)=-0.2421643610089377664864385391;
                A(16,9)=-0.5664019757562306013172434991;
                A(16,10)=0.1771874628203468525014600127e-3;
                A(16,11)=0.5014423068259850996056210242;
                A(16,12)=0.3064012127817816729679168466;
                A(16,13)=-0.1276823788051560205864585339;
                A(16,14)=0.3247356420074775307952922239;
                A(16,15)=0.3618473796262191775472173802e-1;
                A(17,1)=0.1424997987538662936856380379e-2;
                A(17,2)=0;
                A(17,3)=0;
                A(17,4)=0;
                A(17,5)=0;
                A(17,6)=-0.1702826215079596625334655900e-4;
                A(17,7)=0.7380344520119479801877385713e-1;
                A(17,8)=0.2242802589974297123749555342;
                A(17,9)=0;
                A(17,10)=0.6886064175777859968642845433e-2;
                A(17,11)=0.1949958766431660483377549679e-1;
                A(17,12)=-0.1884013399697137820365353780;
                A(17,13)=0.4419137084831260467574116494e-1;
                A(17,14)=0.3818089550488958391368765979e-1;
                A(17,15)=0.5548020401907779252418285064e-1;
                A(17,16)=0.1347154383332695875617293519e-1;
                A(18,1)=0.3109485239098494273513625450e-1;
                A(18,2)=0;
                A(18,3)=0;
                A(18,4)=0;
                A(18,5)=0;
                A(18,6)=-0.4593156369223773119051692264e-4;
                A(18,7)=0.1216757909683287036676281540;
                A(18,8)=0.4415347466907301164835429459;
                A(18,9)=-0.8175066254337973834117229592;
                A(18,10)=0.1105937580434339623330850743;
                A(18,11)=-0.6993517935619971070857051264;
                A(18,12)=-0.7109440534947157331118419857;
                A(18,13)=0.9046040009912160496690560557;
                A(18,14)=0.1167087382162618864114174391e1;
                A(18,15)=-0.2208250617018881371815421672;
                A(18,16)=0.2197668861551819694073021920;
                A(18,17)=0.3475356444046391477824783619e-1;
                A(19,1)=0.1980871288612999895697641722e-1;
                A(19,2)=0;
                A(19,3)=0;
                A(19,4)=0;
                A(19,5)=0;
                A(19,6)=0.8215177736719382423031156570e-5;
                A(19,7)=-0.3591138465638208710224814228e-1;
                A(19,8)=-0.1102176949862629516954052482;
                A(19,9)=0.1174009534262675892260977919;
                A(19,10)=-0.9315711242536572430298479335e-2;
                A(19,11)=0.2059214179779386484105583443;
                A(19,12)=0.2751729989969760041839703778;
                A(19,13)=-0.7252234352249861849631508218e-1;
                A(19,14)=-0.1434854210237120617066495151;
                A(19,15)=0.9817452069416899555388063316e-1;
                A(19,16)=0.4788142501027539102933573843e-1;
                A(19,17)=0.2528961144767536430942173711e-1;
                A(19,18)=0.4994699814223580378253569976e-2;
                A(20,1)=0.5433685416605002611128713384e-1;
                A(20,2)=0;
                A(20,3)=0;
                A(20,4)=0;
                A(20,5)=0;
                A(20,6)=0;
                A(20,7)=0;
                A(20,8)=0;
                A(20,9)=0;
                A(20,10)=0;
                A(20,11)=-0.7028196177287693599155993281;
                A(20,12)=-0.7574471441516858729426845068e-1;
                A(20,13)=0.5970091564181668639770945761;
                A(20,14)=0.5077998923905797158831410566;
                A(20,15)=0.7803513721630428611707284145e-2;
                A(20,16)=0.4910953868318853268158721640e-1;
                A(20,17)=0.6805137650188406090883857916e-1;
                A(20,18)=-0.2240631594717180596069817908e-1;
                A(20,19)=0.1686031620961012499691011166e-1;
                A(21,1)=0.2162496580462148648465375811e-1;
                A(21,2)=0;
                A(21,3)=0;
                A(21,4)=0;
                A(21,5)=0;
                A(21,6)=0;
                A(21,7)=0;
                A(21,8)=0;
                A(21,9)=0;
                A(21,10)=0;
                A(21,11)=0.1146455196485693981118106667;
                A(21,12)=0.1079474883241069775546326492;
                A(21,13)=0.2689620095914987721580753858e-1;
                A(21,14)=0.2788652071946772011310243241e-1;
                A(21,15)=0.9963462823471806922021940733e-1;
                A(21,16)=0.4982293074498374663642074395e-1;
                A(21,17)=0.4009037653486573601359440580e-1;
                A(21,18)=0.1796700744724135163951682987e-2;
                A(21,19)=0.9561341361896604598689526565e-2;
                A(21,20)=0.9332692289624888711718839168e-4;
                
                c(1)=0;
                c(2)=0.1041666666666666666666666667;
                c(3)=0.2083333333333333333333333333;
                c(4)=0.5500000000000000000000000000;
                c(5)=0.5000000000000000000000000000;
                c(6)=0.7556188816150179549777576320;
                c(7)=0.1500000000000000000000000000;
                c(8)=0.5832390688876070299554460664e-1;
                c(9)=0.2000000000000000000000000000;
                c(10)=0.6872598700965161331199471409;
                c(11)=0.2500000000000000000000000000;
                c(12)=0.7500000000000000000000000000e-1;
                c(13)=0.3100000000000000000000000000;
                c(14)=0.1600000000000000000000000000;
                c(15)=0.4600000000000000000000000000;
                c(16)=0.6100000000000000000000000000;
                c(17)=0.7600000000000000000000000000;
                c(18)=0.8500000000000000000000000000;
                c(19)=0.9200000000000000000000012761;
                c(20)=0.1000000000000000000000000000e1;
                c(21)=0.1000000000000000000000000000e1;
                
                bHat=A(end,:)';
                
                bPHat=zeros(21,1);
                bPHat(1)=0.2147820901127829614676264646e-1;
                bPHat(2)=0;
                bPHat(3)=0;
                bPHat(4)=0;
                bPHat(5)=0;
                bPHat(6)=0;
                bPHat(7)=0;
                bPHat(8)=0;
                bPHat(9)=0;
                bPHat(10)=0;
                bPHat(11)=0.1639110625076501029623128826e+00;
                bPHat(12)=0.1177981233399404688143633333e00;
                bPHat(13)=0.2804873935493439037262171756e-01;
                bPHat(14)=0.2906195615716010883438936279e-01;
                bPHat(15)=0.1902857565710350238475328245e00;
                bPHat(16)=0.1221149275129523396892898086e00;
                bPHat(17)=0.1755497879515386479192666401e00;
                bPHat(18)=0.3293567770604707380921199778e-03;
                bPHat(19)=0.1287384223535162591549020064e00;
                bPHat(20)=0.2268365846293389152046665775e-01;
                bPHat(21)=0;
                
                b=[bHat(1:19);0;bHat(20)];
                bP=[];
                isFSAL=true;
                orders=[11;12];
            case 1
                %RKN11(12)21F, which is mentioned in [6] as RKN 11(12)
                %- 20 opt, but whose coefficients are not directly
                %provided. The coefficients were taken from [12].
                A=zeros(21,21);
                A(2,1)=0.5425347222222222222222222225e-2;
                A(3,1)=0.7233796296296296296296296301e-2;
                A(3,2)=0.1446759259259259259259259258e-1;
                A(4,1)=0.1033340000000000000000000000;
                A(4,2)=-0.1703680000000000000000000001;
                A(4,3)=0.2182840000000000000000000001;
                A(5,1)=0.3257575757575757575757575757e-1;
                A(5,2)=0;
                A(5,3)=0.8780487804878048780487804878e-1;
                A(5,4)=0.4619364375461936437546193648e-2;
                A(6,1)=0.4748507367563208000018089973e-1;
                A(6,2)=0;
                A(6,3)=0.1717457745568296060188724074;
                A(6,4)=0.5999520621410557300028385132e-1;
                A(6,5)=0.6253892679998000604987985799e-2;
                A(7,1)=0.6897742477766473591191027092e-2;
                A(7,2)=0;
                A(7,3)=0.6890192907847598728171600460e-2;
                A(7,4)=0.1163520086778733422170956371e-1;
                A(7,5)=-0.1344647159854360748127587603e-1;
                A(7,6)=-0.7266646548577990597963152304e-3;
                A(8,1)=0.1349737761438706440043676626e-2;
                A(8,2)=0;
                A(8,3)=0;
                A(8,4)=0.6075294762105980384386224708e-3;
                A(8,5)=-0.6793734869158527204348147999e-3;
                A(8,6)=-0.4099894989642249543680965575e-4;
                A(8,7)=0.4639442565473848303431810095e-3;
                A(9,1)=0.3523183546767263585786907763e-2;
                A(9,2)=0;
                A(9,3)=0;
                A(9,4)=0;
                A(9,5)=-0.2584392204058368115482471755e-5;
                A(9,6)=0.2758330832880140022000070680e-6;
                A(9,7)=0.4071837479676998987495697937e-2;
                A(9,8)=0.1240728753267650778083067676e-1;
                A(10,1)=-0.1340488944178098066953190204e-1;
                A(10,2)=0;
                A(10,3)=0;
                A(10,4)=0;
                A(10,5)=0.4748167290240606363971037112e-1;
                A(10,6)=0.6978075811258520497345866824e-3;
                A(10,7)=-0.2104272720250956294017709271;
                A(10,8)=0.1479810168188958303177052879;
                A(10,9)=0.2638347286869889791697278318;
                A(11,1)=0.4338886679167641664845535142e-2;
                A(11,2)=0;
                A(11,3)=0;
                A(11,4)=0;
                A(11,5)=0;
                A(11,6)=-0.1244276150937331764200609062e-5;
                A(11,7)=0.7318827258766195401581730488e-2;
                A(11,8)=0.1702903672564807515932491560e-1;
                A(11,9)=0.2562063489400243489046551726e-2;
                A(11,10)=0.2430123168781616965467654675e-5;
                A(12,1)=0.1315981034561655557920249494e-2;
                A(12,2)=0;
                A(12,3)=0;
                A(12,4)=0;
                A(12,5)=0;
                A(12,6)=-0.8782527297984828178544800505e-6;
                A(12,7)=-0.1379773004506861094363558123e-2;
                A(12,8)=0.2445615867912388681251959553e-2;
                A(12,9)=0.1156657207478131509996721934e-2;
                A(12,10)=0.1754639209365541673290930900e-5;
                A(12,11)=-0.3268452848936317136608093302e-3;
                A(13,1)=0.4757648391153826839907752946e-2;
                A(13,2)=0;
                A(13,3)=0;
                A(13,4)=0;
                A(13,5)=0;
                A(13,6)=0;
                A(13,7)=0.1549791354107748594550575393e-1;
                A(13,8)=0.2578358836615907329049723352e-1;
                A(13,9)=0.7723584402647932338238513404e-3;
                A(13,10)=-0.5453046730326065812618845430e-7;
                A(13,11)=0.3225703225768342511932988988e-2;
                A(13,12)=-0.6708837121456218561009455083e-2;
                A(14,1)=0.2062453723473214735559719306e-2;
                A(14,2)=0;
                A(14,3)=0;
                A(14,4)=0;
                A(14,5)=0;
                A(14,6)=0;
                A(14,7)=0;
                A(14,8)=0.8391108840414372379428007718e-2;
                A(14,9)=0.1082586893526879718172321785e-2;
                A(14,10)=-0.9552204475152784971565217179e-7;
                A(14,11)=-0.8100311594783424692710533987e-3;
                A(14,12)=-0.2552456978531004853224335309e-2;
                A(14,13)=0.2156065659208820171850554315e-3;
                A(15,1)=0.4612793314617264068992923362e-2;
                A(15,2)=0;
                A(15,3)=0;
                A(15,4)=0;
                A(15,5)=0;
                A(15,6)=0;
                A(15,7)=0;
                A(15,8)=0;
                A(15,9)=0.3671972962570833488917211322;
                A(15,10)=0.4624739571820681985355479553e-4;
                A(15,11)=-0.3442546230850344075615885905;
                A(15,12)=0.1304094665996243603039865640;
                A(15,13)=0.1564155198841677601881860482;
                A(15,14)=-0.2022811803466452827111516326;
                A(16,1)=0.4764722488136250883238518887e-1;
                A(16,2)=0;
                A(16,3)=0;
                A(16,4)=0;
                A(16,5)=0;
                A(16,6)=0.2469440106146505625377840436e-4;
                A(16,7)=-0.1042862143549249441082095027;
                A(16,8)=-0.2981571590755317696349688954;
                A(16,9)=-0.4329399367490009012615423678;
                A(16,10)=-0.1878322210652231000279915795e-3;
                A(16,11)=0.5314840674401284108321048956;
                A(16,12)=0.3397410273320287108334763297;
                A(16,13)=-0.1598998723085592532962583681;
                A(16,14)=0.2238129751121937855575210677;
                A(16,15)=0.2259623062043221028926586297e-1;
                A(17,1)=0.1902451058984853527346206363e-1;
                A(17,2)=0;
                A(17,3)=0;
                A(17,4)=0;
                A(17,5)=0;
                A(17,6)=-0.4882846979701356315661132055e-7;
                A(17,7)=0.9958369503541465773976755947e-2;
                A(17,8)=0.2277751445458477254054119562e-1;
                A(17,9)=0;
                A(17,10)=0.2996011280298835998018263725e-2;
                A(17,11)=-0.3658645910207106377190192858e-1;
                A(17,12)=0.1322941269080730874394828249e-1;
                A(17,13)=0.1130307620345566171348268026;
                A(17,14)=0.7159961376794738694817346241e-1;
                A(17,15)=0.3153763895783865533949080181e-1;
                A(17,16)=0.1683350472924228303302745352e-1;
                A(18,1)=0.1029170267943476427286491449e-1;
                A(18,2)=0;
                A(18,3)=0;
                A(18,4)=0;
                A(18,5)=0;
                A(18,6)=-0.3578990654680073931152431606e-4;
                A(18,7)=0.1201493680572643049749038619;
                A(18,8)=0.3737402617729644461656779293;
                A(18,9)=-0.1395209500373147796992155667e1;
                A(18,10)=-0.4824305334489776287397577572e-2;
                A(18,11)=0.1378982579792852870857948043e1;
                A(18,12)=-0.7213193572545734105413799248;
                A(18,13)=-0.4435071523221466577698068378;
                A(18,14)=0.8886457915442507835748030669;
                A(18,15)=0.5526124056477397137931770684e-1;
                A(18,16)=0.2686835085031771552727496965e-1;
                A(18,17)=0.7286509636076835577261033251e-2;
                A(19,1)=0.1526550484230059201602058577e-1;
                A(19,2)=0;
                A(19,3)=0;
                A(19,4)=0;
                A(19,5)=0;
                A(19,6)=0.3324351624050014581690792469e-5;
                A(19,7)=-0.1633568188026327026286314180e-1;
                A(19,8)=-0.4739891788157655468203173752e-1;
                A(19,9)=0.5726143210435682891937291916e-1;
                A(19,10)=-0.2118780632671606033608191069e-2;
                A(19,11)=0.1527606615230892001844132357;
                A(19,12)=0.2473451593183186070325554611;
                A(19,13)=0.2433033586985477852062515545e-2;
                A(19,14)=-0.1308274895314965385646731789;
                A(19,15)=0.4539155318997185313360896864e-1;
                A(19,16)=0.7975518417920696991874667819e-1;
                A(19,17)=-0.3631391480459251475363844572e-1;
                A(19,18)=0.5129790733334246421751407628e-1;
                A(20,1)=0.3754662379446297418280769111e-1;
                A(20,2)=0;
                A(20,3)=0;
                A(20,4)=0;
                A(20,5)=0;
                A(20,6)=0;
                A(20,7)=0;
                A(20,8)=0;
                A(20,9)=0;
                A(20,10)=0;
                A(20,11)=0.4929894742182617608855792266;
                A(20,12)=0.4526697482923057917557427810e-1;
                A(20,13)=-0.4593369440456942747900428666;
                A(20,14)=0.4759722032869840752638237447e-1;
                A(20,15)=0.4654889251258994259593464956;
                A(20,16)=-0.3321755294689402941095099898;
                A(20,17)=0.3952785978121922231177107745;
                A(20,18)=-0.2123660239786777322823065681;
                A(20,19)=0.1971068138456693033445858405e-1;
                A(21,1)=0.2285890998625555379695120949e-1;
                A(21,2)=0;
                A(21,3)=0;
                A(21,4)=0;
                A(21,5)=0;
                A(21,6)=0;
                A(21,7)=0;
                A(21,8)=0;
                A(21,9)=0;
                A(21,10)=0;
                A(21,11)=0.1203126822802527505218717499;
                A(21,12)=0.1155833456919026151325826139;
                A(21,13)=0.3502768234280783458331654724e-1;
                A(21,14)=0.9350324233933777658031475452e-2;
                A(21,15)=0.1001712618244591513505119362;
                A(21,16)=0.3318719494747818080400948715e-1;
                A(21,17)=0.3540771699207406021731391914e-1;
                A(21,18)=0.1686603296267794771166345990e-1;
                A(21,19)=0.1114152181526187933663041330e-1;
                A(21,20)=0.9332692289624888711718839168e-4;
                
                c(1)=0;
                c(2)=0.1041666666666666666666666667;
                c(3)=0.2083333333333333333333333333;
                c(4)=0.5500000000000000000000000000;
                c(5)=0.5000000000000000000000000000;
                c(6)=0.7556188816150179549777576320;
                c(7)=0.1500000000000000000000000000;
                c(8)=0.5832390688876070299554460664e-1;
                c(9)=0.2000000000000000000000000000;
                c(10)=0.6872598700965161331199471409;
                c(11)=0.2500000000000000000000000000;
                c(12)=0.8015624999999999999999999974e-1;
                c(13)=0.2943749999999999999999999981;
                c(14)=0.1295312499999999999999999991;
                c(15)=0.4735937499999999999999999989;
                c(16)=0.5828124999999999999999999964;
                c(17)=0.7271874999999999999999999954;
                c(18)=0.7698437499999999999999999929;
                c(19)=0.9148977819391580085721647998;
                c(20)=0.1000000000000000000000000000e1;
                c(21)=0.1000000000000000000000000000e1;
                
                bHat=A(end,:)';
                
                bPHat=zeros(21,1);
                bPHat(1)=0.2242929404397600530131470701e-1;
                bPHat(2)=0;
                bPHat(3)=0;
                bPHat(4)=0;
                bPHat(5)=0;
                bPHat(6)=0;
                bPHat(7)=0;
                bPHat(8)=0;
                bPHat(9)=0;
                bPHat(10)=0;
                bPHat(11)=0.1900303404131709886214564320;
                bPHat(12)=0.1308259238734620016022237099;
                bPHat(13)=0.1903103706904367636603518091e-1;
                bPHat(14)=0.4536441087947228435282166600e-3;
                bPHat(15)=0.2082960450051360233092459511;
                bPHat(16)=0.5819210373332354680001930773e-1;
                bPHat(17)=0.1664102268922875908112086647;
                bPHat(18)=0.4287072842584152852324451740e-1;
                bPHat(19)=0.1374085143361514078379615741;
                bPHat(20)=0.2405214209881250798376173855e-1;
                bPHat(21)=0;
                
                b=[bHat(1:19);0;bHat(20)];
                bP=[];
                isFSAL=true;
                orders=[11;12];
            case 2
                %RKN11(10)18F from [6]. These are the RKN11(10)-17
                %coefficients (not the RKN11(10-17 opt coefficients) in
                %Table 1. The actual coefficients are not given in the
                %paper and were taken from [12].

                A=zeros(18,18); 

                A(2,1)=0.8680555555555555555555555555e-3;
                A(3,1)=0.1157407407407407407407407407e-2;
                A(3,2)=0.2314814814814814814814814815e-2;
                A(4,1)=0.1040000000000000000000000000e-1;
                A(4,2)=-0.1280000000000000000000000000e-1;
                A(4,3)=0.2240000000000000000000000000e-1;
                A(5,1)=0.6510416666666666666666666668e-2;
                A(5,2)=0;
                A(5,3)=0.2008928571428571428571428571e-1;
                A(5,4)=0.4650297619047619047619047619e-2;
                A(6,1)=-0.7442558350050598650718971132;
                A(6,2)=0;
                A(6,3)=0.2398395982021852194994338203e1;
                A(6,4)=-0.4296409767675047353488214866e1;
                A(6,5)=0.2932519053764830987284594640e1;
                A(7,1)=0.8320473251028806584362139914e-2;
                A(7,2)=0;
                A(7,3)=0.3148612271419288963148612272e-1;
                A(7,4)=0.7784553453103358941400372646e-2;
                A(7,5)=0.7962484448272561967652406931e-2;
                A(7,6)=0.1921688957938430654513332561e-5;
                A(8,1)=0.3177389300411522633744855977e-2;
                A(8,2)=0;
                A(8,3)=0.7172885695720783440081685729e-2;
                A(8,4)=-0.4098087864566207315736505945e-2;
                A(8,5)=0.3150419446837017896449420999e-2;
                A(8,6)=0.1337454518694045049020649864e-5;
                A(8,7)=-0.5150551440329218106995884776e-3;
                A(9,1)=0.1925455729166666666666666669e-1;
                A(9,2)=0;
                A(9,3)=0;
                A(9,4)=0;
                A(9,5)=-0.1683277962347729789590254748e-1;
                A(9,6)=0.8261897683227625088090204146e-4;
                A(9,7)=0.4322916666666666666666666683e-1;
                A(9,8)=0.7926643668831168831168831192e-1;
                A(10,1)=0.3890677244084362139917695514e-1;
                A(10,2)=0;
                A(10,3)=0;
                A(10,4)=0;
                A(10,5)=0.1162107105259560039921797094;
                A(10,6)=0.1280779593149063653526890507e-1;
                A(10,7)=-0.2813143004115226337448560306e-1;
                A(10,8)=0.1000838847074642529187983710;
                A(10,9)=0.1073444886576199707512838846;
                A(11,1)=0.2298638614607445987654320989e-2;
                A(11,2)=0;
                A(11,3)=0;
                A(11,4)=0;
                A(11,5)=0;
                A(11,6)=-0.2512712158652423387447237099e-3;
                A(11,7)=-0.9907510109453589391860996819e-3;
                A(11,8)=0.1723141771801382228570780760e-2;
                A(11,9)=0.5659797947913478216508519521e-3;
                A(11,10)=0.1264842678326474622770919127e-3;
                A(12,1)=0.6848028273809523809523809520e-2;
                A(12,2)=0;
                A(12,3)=0;
                A(12,4)=0;
                A(12,5)=0;
                A(12,6)=0;
                A(12,7)=0.1209077380952380952380952379e-2;
                A(12,8)=0.6982995613017934446505875037e-2;
                A(12,9)=-0.2051204004329004329004328998e-3;
                A(12,10)=0.6288973922902494331065759621e-5;
                A(12,11)=0.1640873015873015873015873020e-1;
                A(13,1)=-0.7759862989766188419387152553e-2;
                A(13,2)=0;
                A(13,3)=0;
                A(13,4)=0;
                A(13,5)=0.1224804770337618196047355743e-1;
                A(13,6)=0.2077286426548597692800167757e-2;
                A(13,7)=-0.3084352111917326351227812139e-1;
                A(13,8)=-0.1322398920232751620391476668;
                A(13,9)=-0.5694717797969194984675336503e-3;
                A(13,10)=-0.1144721852146542255214994914e-2;
                A(13,11)=0.1366017573669854049407166350;
                A(13,12)=0.7718593382280344668606066469e-1;
                A(14,1)=0.3147136427366806555520298306e-1;
                A(14,2)=0;
                A(14,3)=0;
                A(14,4)=0;
                A(14,5)=0.1569719156364494947817982402;
                A(14,6)=-0.6654987695858527013495559480e-2;
                A(14,7)=0;
                A(14,8)=0.8269612470709131479069575139e-1;
                A(14,9)=0.5277046384601007307143946853e-1;
                A(14,10)=0.3564062192959960179721795389e-2;
                A(14,11)=-0.3782761972062383585113875547e-2;
                A(14,12)=-0.8098383997314998286826664921e-1;
                A(14,13)=-0.1383011879288579268975993216e-1;
                A(15,1)=-0.8551325588087360429217853897e-2;
                A(15,2)=0;
                A(15,3)=0;
                A(15,4)=0;
                A(15,5)=-0.8401264079085973583990782598e-1;
                A(15,6)=0.1186337961569032036179412516e-2;
                A(15,7)=-0.1135641886457673424136874024;
                A(15,8)=-0.2190092537441172654013264721;
                A(15,8)=0.4620369354790417202229644247e-1;
                A(15,10)=-0.1120372383953845501878892251e-2;
                A(15,11)=0.2618995077705879099608244632;
                A(15,12)=0.2841919156005766259624988398;
                A(15,13)=0.1007728344392311617920989094;
                A(15,14)=0.1325349183291664781212037926e-1;
                A(16,1)=0.3227323796009178105251112707e-1;
                A(16,2)=0;
                A(16,3)=0;
                A(16,4)=0;
                A(16,5)=-0.6552253393177109965556635085e-3;
                A(16,6)=-0.2053435081848825845464654784e-2;
                A(16,7)=0.8183285936474925294914041188e-1;
                A(16,8)=0.9687035581635444102020453817e-1;
                A(16,9)=0.1009741982188193601657085464;
                A(16,10)=0.7943878492642055987892133169e-2;
                A(16,11)=0.3988879155234310326668935987e-1;
                A(16,12)=0.6260072946369460370714950309e-1;
                A(16,13)=-0.2982691974303517060774392676e-1;
                A(16,14)=0.1160314062805673150276068512e-1;
                A(16,15)=0.2566041805220703528164152583e-1;
                A(17,1)=0.1598229814960049517404649275e-1;
                A(17,2)=0;
                A(17,3)=0;
                A(17,4)=0;
                A(17,5)=0;
                A(17,6)=0;
                A(17,7)=0;
                A(17,8)=0;
                A(17,9)=0.5685521470453431953596264635e-1;
                A(17,10)=-0.2342953776266563711445070295e-1;
                A(17,11)=0.1485105390278212809150227904;
                A(17,12)=0.3889673830278987492323711070e-1;
                A(17,13)=0.1485093897141690039134706881;
                A(17,14)=0.3093315321482061128468761163e-1;
                A(17,15)=0.6806447423483744354174394038e-1;
                A(17,16)=0.1567773041409260782627942270e-1;
                A(18,1)=0.2352719229768410096278945558e-1;
                A(18,2)=0;
                A(18,3)=0;
                A(18,4)=0;
                A(18,5)=0;
                A(18,6)=0;
                A(18,7)=0;
                A(18,8)=0;
                A(18,9)=0.1095238095238095238095226633;
                A(18,10)=0.5714285714285714285717006159e-2;
                A(18,11)=0.1235521235521235521235522436;
                A(18,12)=0.1287533440342429106474043912;
                A(18,13)=0.3989010989010989010989112569e-1;
                A(18,14)=0.2382352941176470588235556148e-1;
                A(18,15)=0.3621808143547273982056220803e-1;
                A(18,16)=0.8997524140506862358205344942e-2;
                A(18,11)=0;
                
                c(1)=0;
                c(2)=0.4166666666666666666666666667e-1;
                c(3)=0.8333333333333333333333333333e-1;
                c(4)=0.2000000000000000000000000000;
                c(5)=0.2500000000000000000000000000;
                c(6)=0.7619047619047619047619047647;
                c(7)=0.3333333333333333333333333333;
                c(8)=0.1333333333333333333333333337;
                c(9)=0.5000000000000000000000000000;
                c(10)=0.8333333333333333333333333333;
                c(11)=0.8333333333333333333333333333e-1;
                c(12)=0.2500000000000000000000000000;
                c(13)=0.3333333333333333333333333333;
                c(14)=0.6666666666666666666666666666;
                c(15)=0.7500000000000000000000000000;
                c(16)=0.9242424242424242424242460752;
                c(17)=0.1000000000000000000000000000e+1;
                c(18)=0.1000000000000000000000000000e+1;

                bHat=A(end,:)';

                bPHat=zeros(18,1);
                bPHat(1)=0.2352719229768410096278945558e-1;
                bPHat(2)=0;
                bPHat(3)=0;
                bPHat(4)=0;
                bPHat(5)=0;
                bPHat(6)=0;
                bPHat(7)=0;
                bPHat(8)=0;
                bPHat(9)=0.2190476190476190476190453266;
                bPHat(10)=0.3428571428571428571430203695e-1;
                bPHat(11)=0.1347841347841347841347842657;
                bPHat(12)=0.1716711253789905475298725216;
                bPHat(13)=0.5983516483516483516483668854e-1;
                bPHat(14)=0.7147058823529411764706668444e-1;
                bPHat(15)=0.1448723257418909592822488321;
                bPHat(16)=0.1187673186546905831283162769;
                bPHat(17)=0.2173881673881673881673791154e-1;
                bPHat(18)=0;

                lambda=1.4e-5;
                b=[bHat(1:16);-lambda;lambda];
                bP=[];
                isFSAL=true;
                orders=[11;10];
            case 3
                %RKN11(10)18F from [6]. These are the RKN11(10)-17 opt
                %coefficients in Table 1. The actual coefficients are
                %not given in the paper and were taken from [12].

                A=zeros(18,18);

                A(2,1)=0.8627885207850662210258997910e-3;
                A(3,1)=0.1150384694380088294701199722e-2;
                A(3,2)=0.2300769388760176589402399443e-2;
                A(4,1)=0.1048813009577021740354573069e-1;
                A(4,2)=-0.1307378621495554451132432640e-1;
                A(4,3)=0.2258565611918532710777859571e-1;
                A(5,1)=0.6540727920973800634409298201e-2;
                A(5,2)=0;
                A(5,3)=0.2032812414440663315636725431e-1;
                A(5,4)=0.4905792595569273240473447117e-2;
                A(6,1)=-0.9866296281968030800021038325;
                A(6,2)=0;
                A(6,3)=0.3094647192454822576423231408e1;
                A(6,4)=-0.5345954041360107668766936002e1;
                A(6,5)=0.3554526041401407323265167534e1;
                A(7,1)=0.8305872179116544972367559928e-2;
                A(7,2)=0;
                A(7,3)=0.3137039735660881318054124472e-1;
                A(7,4)=0.8393242819482349481524813340e-2;
                A(7,5)=0.7484610697242405536405592985e-2;
                A(7,6)=0.1432503105442384716344571743e-5;
                A(8,1)=0.3326256442868476167267206871e-2;
                A(8,2)=0;
                A(8,3)=0.7743460616381346943027737308e-2;
                A(8,4)=-0.3941298661638303825193706970e-2;
                A(8,5)=0.3040306178972026742079685521e-2;
                A(8,6)=0.1041535422216837899256070373e-5;
                A(8,7)=-0.5280958053014040723291822062e-3;
                A(9,1)=0.2005356082630297417517719541e-1;
                A(9,2)=0;
                A(9,3)=0;
                A(9,4)=0;
                A(9,5)=-0.2615084559917886386194315847e-1;
                A(9,6)=0.7037244725789190110573832314e-4;
                A(9,7)=0.4880358378547006117565645252e-1;
                A(9,8)=0.8503189863902488973500376967e-1;
                A(10,1)=0.4460073080576536587846203080e-1;
                A(10,2)=0;
                A(10,3)=0;
                A(10,4)=0;
                A(10,5)=0.1871119657446980285366273356;
                A(10,6)=0.1327438600652787132655931398e-1;
                A(10,7)=-0.9125776086410149316508222463e-1;
                A(10,8)=0.8389592974628709573178294907e-1;
                A(10,9)=0.1351312854472502692784561413;
                A(11,1)=0.2325938867166884468088530959e-2;
                A(11,2)=0;
                A(11,3)=0;
                A(11,4)=0;
                A(11,5)=0;
                A(11,6)=-0.2259845307860368767743334584e-3;
                A(11,7)=-0.1003643483120693922361997315e-2;
                A(11,8)=0.1729593142949256358697473742e-2;
                A(11,9)=0.5237186040042921580080116021e-3;
                A(11,10)=0.1225996220085200365645366915e-3;
                A(12,1)=0.6706133017061727106482743292e-2;
                A(12,2)=0;
                A(12,3)=0;
                A(12,4)=0;
                A(12,5)=0;
                A(12,6)=0;
                A(12,7)=0.1064819794414139802620133553e-2;
                A(12,8)=0.6034723610379959516770029935e-2;
                A(12,9)=-0.1700013721173865454683816653e-3;
                A(12,10)=0.4600910262781625002278745485e-5;
                A(12,11)=0.1698301524635864177584319607e-1;
                A(13,1)=-0.4045295089969598263407026818e-2;
                A(13,2)=0;
                A(13,3)=0;
                A(13,4)=0;
                A(13,5)=0.2724713861791203266764527045e-2;
                A(13,6)=0.2038629948138949746024409902e-2;
                A(13,7)=-0.1658035279760125107209717027e-1;
                A(13,8)=-0.1023528784753465524938428817;
                A(13,9)=-0.1954120844057010823475271332e-2;
                A(13,10)=-0.1163323402579112625631942833e-2;
                A(13,11)=0.1179292783494758994718313003;
                A(13,12)=0.7614367853166331480772294322e-1;
                A(14,1)=0.2928470788832902336077398654e-1;
                A(14,2)=0;
                A(14,3)=0;
                A(14,4)=0;
                A(14,5)=0.2557604301751744221950094250e-1;
                A(14,6)=-0.2770490798872626353718504215e-2;
                A(14,7)=0;
                A(14,8)=0.1214828801961999047764661964;
                A(14,9)=0.1014842148793759653727359216e-1;
                A(14,10)=0.1554260948206398650240168286e-2;
                A(14,11)=-0.2876086635474601440456614479e-1;
                A(14,12)=-0.2609805932570869996876447678e-1;
                A(14,13)=0.3118070841327781719668312337e-1;
                A(15,1)=0.1094868564814643322352823048e-2;
                A(15,2)=0;
                A(15,3)=0;
                A(15,4)=0;
                A(15,5)=-0.1513004798311404998761949729e-1;
                A(15,6)=0.8033888332100436127218323929e-3;
                A(15,7)=-0.2921656340391943713085827453e-1;
                A(15,8)=-0.1354774693869972349812600372;
                A(15,9)=-0.2956294173068199055515575719e-1;
                A(15,10)=-0.9865464273698914825160869666e-4;
                A(15,11)=0.1927563437672329668442225360;
                A(15,12)=0.1739470522627136604341763792;
                A(15,13)=0.6679125168274304414563517938e-1;
                A(15,14)=0.6067050107115917156903642081e-1;
                A(16,1)=0.3618969021059104097769243441e-1;
                A(16,2)=0;
                A(16,3)=0;
                A(16,4)=0;
                A(16,5)=-0.1031559446790956758888425790e-1;
                A(16,6)=-0.1591154266286366850436734751e-2;
                A(16,7)=0.3465993450035860634024176283e-1;
                A(16,8)=0.8393633290039513855306722695e-1;
                A(16,9)=0.2659698264204566928627519117;
                A(16,10)=0.7575916354995080153158336955e-3;
                A(16,11)=0.3586548416043737832991305347e-1;
                A(16,12)=0.1281974574530537903485720644;
                A(16,13)=-0.8743695425560904559912489751e-1;
                A(16,14)=-0.1052135460954427771075420252;
                A(16,15)=0.4837379673467808096434532545e-1;
                A(17,1)=-0.1751456120510270062862652183e-2;
                A(17,2)=0;
                A(17,3)=0;
                A(17,4)=0;
                A(17,5)=0;
                A(17,6)=0;
                A(17,7)=0;
                A(17,8)=0;
                A(17,9)=-0.8210098851189696378716630381;
                A(17,10)=-0.1140654283891927400308961378e-2;
                A(17,11)=0.2075616695240093175504560878;
                A(17,12)=-0.1249838280156769622567987062;
                A(17,13)=0.6128417654283759930709135422;
                A(17,14)=0.6137087696615115366457015345;
                A(17,15)=-0.2436848586879919336178716828e-3;
                A(17,16)=0.1501730378383994225818006495e-1;
                A(18,1)=0.2366173330336702092688069602e-1;
                A(18,2)=0;
                A(18,3)=0;
                A(18,4)=0;
                A(18,5)=0;
                A(18,6)=0;
                A(18,7)=0;
                A(18,8)=0;
                A(18,9)=0.4647282849123642792055088353e-1;
                A(18,10)=0.2147842164498174139028692258e-2;
                A(18,11)=0.1227540335672074472801687518;
                A(18,12)=0.1319559360176842028133017261;
                A(18,13)=0.6242703419929078584193419875e-1;
                A(18,14)=0.5589976511302376555244459167e-1;
                A(18,15)=0.4617524523272118896242712402e-1;
                A(18,16)=0.8505581910970986563263335804e-2;
                A(18,17)=0;

                c(1)=0;
                c(2)=0.4154006549790373696025353238e-1;
                c(3)=0.8308013099580747392050706476e-1;
                c(4)=0.2000000000000000000000000000;
                c(5)=0.2520898437499999999999999985;
                c(6)=0.7957255359724471616290815220;
                c(7)=0.3333333333333333333333333333;
                c(8)=0.1388644685058374375856553886;
                c(9)=0.5055859374999999999999999950;
                c(10)=0.8634309895833333333333333224;
                c(11)=0.8333333333333333333333333333e-1;
                c(12)=0.2474804687499999999999999997;
                c(13)=0.3814192708333333333333333298;
                c(14)=0.5685026041666666666666666572;
                c(15)=0.7570703124999999999999999947;
                c(16)=0.9267069277071608028412563026;
                c(17)=0.1000000000000000000000000000e1;
                c(18)=0.1000000000000000000000000000e1;

                bHat=A(end,:)';

                bPHat=zeros(18,1);
                bPHat(1)=0.2366173330336702092688069602e-1;
                bPHat(2)=0;
                bPHat(3)=0;
                bPHat(4)=0;
                bPHat(5)=0;
                bPHat(6)=0;
                bPHat(7)=0;
                bPHat(8)=0;
                bPHat(9)=0.9399576592997175908715355979e-1;
                bPHat(10)=0.1572715624097437897482037968e-1;
                bPHat(11)=0.1339134911642263061238204566;
                bPHat(12)=0.1753521743129962154232149388;
                bPHat(13)=0.1009197850107464658399929787;
                bPHat(14)=0.1295483255584122995391455550;
                bPHat(15)=0.1900765843315102809837810498;
                bPHat(16)=0.1160489203807327640439274289;
                bPHat(17)=0.2075606376706250905726295692e-1;
                bPHat(18)=0;
                
                lambda=1.4e-6;
                b=[bHat(1:16);-lambda;lambda];
                bP=[];
                isFSAL=true;
                orders=[11;10];
            otherwise
                error('Unknown Solution Choice');
        end
    case 12
        switch(solutionChoice)
            case 0
                %RKN12(10)17 from [9]. However, [9] just says that the
                %authors could be contacted for the coefficients. The
                %coefficients were thus taken (in rational form) from the
                %comments in the function in [10], as the comments said
                %that those were the coefficients provided by the authors
                %of [9].

                A=zeros(17,17);
                A(2,1)=1/5000;

                A(3,1)=1/3750;
                A(3,2)=1/1875;

                A(4,1)=7/2400;
                A(4,2)=-1/240;
                A(4,3)=1/160;

                A(5,1)=2/1215;
                A(5,2)=0;
                A(5,3)=4/729;
                A(5,4)=32/18225;

                A(6,1)=152/78125;
                A(6,2)=0;
                A(6,3)=1408/196875;
                A(6,4)=2048/703125;
                A(6,5)=432/546875;

                A(7,1)=29/51200;
                A(7,2)=0;
                A(7,3)=341/387072;
                A(7,4)=-151/345600;
                A(7,5)=243/716800;
                A(7,6)=-11/110592;

                A(8,1)=37/12000;
                A(8,2)=0;
                A(8,3)=0;
                A(8,4)=2/1125;
                A(8,5)=27/10000;
                A(8,6)=5/3168;
                A(8,7)=224/20625;

                A(9,1)=100467472123373/27511470744477696;
                A(9,2)=0;
                A(9,3)=101066550784375/25488568483854336;
                A(9,4)=49478218404275/15475202293768704;
                A(9,5)=21990175014231/2674726322379776;
                A(9,6)=-3576386017671875/2723635603703291904;
                A(9,7)=16163228153/1654104722787;
                A(9,8)=38747524076705/10316801529179136;

                A(10,1)=62178936641284701329/16772293867250014666848;
                A(10,2)=0;
                A(10,3)=46108564356250/9072835168325103;
                A(10,4)=1522561724950/1296119309760729;
                A(10,5)=-45978886013453735443/2174186242050927827184;
                A(10,6)=299403512366617849203125/4981371278573254356053856;
                A(10,7)=15571226634087127616/774466927638876610083;
                A(10,8)=-133736375367792139885/4717207650164066625051;
                A(10,9)=7461389216/501451974639;

                A(11,1)=501256914705531962342417557181/14270506505142656332600844507392;
                A(11,2)=0;
                A(11,3)=-1143766215625/132752960853408;
                A(11,4)=-6864570325/1185294293334;
                A(11,5)=194348369382310456605879163404183/99893545535998594328205911551744;
                A(11,6)=-94634958447010580589908066176109375/27549212808177898050085930321520256;
                A(11,7)=-17006472665356285286219618514/155584463413110817059022733377;
                A(11,8)=33530528814694461893884349656345/14270506505142656332600844507392;
                A(11,9)=-13439782155791134368/17777268379678341919;
                A(11,10)=1441341768767571/13159456712985856;

                A(12,1)=105854110734231079069010159870911189747853/5156624149476760916008179453333467046288864;
                A(12,2)=0;
                A(12,3)=-144579793509250000/19842290513127000261;
                A(12,4)=-101935644099967250/48188419817594143491;
                A(12,5)=1585474394319811696785932424388196965/1709257457318830856936350991091849456;
                A(12,6)=-843499776333774172853009613469456309715703125/510505790798199330684809765880013237582597536;
                A(12,7)=-15057703799298260121553794369056896088480/714327132646734138085088291809720015274157;
                A(12,8)=1749840442221344572962864758990584360232600/1450300542040339007627300471250037606768743;
                A(12,9)=-11255775246405733991656178432768/27206626483067760480757659602193;
                A(12,10)=669010348769579696/7368057640845834597;
                A(12,11)=4598083098752/858563707934367;

                A(13,1)=-1639758773684715326849438048667467886824967397/11447568726280607813664651120965112496134881280;
                A(13,2)=0;
                A(13,3)=3942453384375/314673684985856;
                A(13,4)=11737114158175/1719466921529856;
                A(13,5)=-23710715033675876683332701739887457/4940189888325748664958546898558976;
                A(13,6)=498150575499633273684774666731162498301909124515625/87415924307623977386706008889913792042985180430336;
                A(13,7)=64881557768202140428371179540010005713998551/85896810580242200654071863296887242202224768;
                A(13,8)=-2336309182318568698279006266321563486172654055/18316109962048972501863441793544179993815810048;
                A(13,9)=-493399374030747471036018890494175/251658285736841065236836942273664;
                A(13,10)=418285003077108927126515545155/455369916679568501838710898688;
                A(13,11)=-15171723902781457/63532954684873728;
                A(13,12)=1501203688494867/9434957026426880;

                A(14,1)=34188549803371802849576690267872548602326398788953/42496542183406636759747616530102745233754251202880;
                A(14,2)=0;
                A(14,3)=-18971246281693750/1138830954584356089;
                A(14,4)=-59230464334542700/2765732318276293359;
                A(14,5)=5147939981309774383134903239728881770043/305929030949718561059100251282184099064;
                A(14,6)=-362572021355026772337065830211467821556305840522907812/324512095420929759624784749347170583153994213035432256;
                A(14,7)=-60305503318319653518547439098565661266182518307816/17856872599361492097414471889911176856851308259643;
                A(14,8)=-1036461878759982363277481306266144563833492657780645/67994467493450618815596186448164392374006801924608;
                A(14,9)=128398681100219349205889126776607047000/7473801441221286756994805323613917077;
                A(14,10)=-49156374556350058671822606102117/9039888303968618912866414995904;
                A(14,11)=12253036339964386945/8828680926314891943;
                A(14,12)=-647188390508758231059/1092148506009694282240;
                A(14,13)=10915833599872/368729913707897;

                A(15,1)=-4939337286263213195547765488387521892799075623007291241961609516532/5408250052307451520718178852915698257207815452080611897685945761264;
                A(15,2)=0;
                A(15,3)=7588799849596321243074032368290625/3147217749590114939838670370597819616;
                A(15,4)=16870665568420512953501332587233725/955405388268427749593882076788623812;
                A(15,5)=-80864251591837801485030858227147601466956843757908779606/54447992506702009927986632715967769032585338753056786562;
                A(15,6)=4610328329649866588704236006423149172472141907645890762410296050212/2135428689710103309390449198881479603148467934048051598947383737508;
                A(15,7)=4159963831215576225909381034291748993887819834160487158570788681/1040533184037697645660563795162185415624171583014576682740416336;
                A(15,8)=738139214212435127943380193414870655354213707189052136566460666444958/259596002510757672994472584939953516345975141699869371088925396540699;
                A(15,9)=-333683433458405281346882867597135977469443722954786270692/132102862435303266640535426836147775872819092781208127980;
                A(15,10)=42661937996741208687503901295747546613008142604821349179/55162410119399855550108207148248549410926885937244965785;
                A(15,11)=-630755628691078947314733435975762542732598947/333503232300511886435069380727586592765317456;
                A(15,12)=1522350657470125698997653827133798314909646891/1520094067152619944607524353149267399623188480;
                A(15,13)=305575414262755427083262606101825880/65839748482572312891297405431209259829;
                A(15,14)=256624643108055110568255672032710477795/22874609758516552135947898572671559986304;

                A(16,1)=-571597862947184314270186718640978947715678864684269066846/207705506488030390761613596901272001190776700439774478634;
                A(16,2)=0;
                A(16,3)=66981514290625/1829501741761029;
                A(16,4)=43495576635800/4443075658562499;
                A(16,5)=-127865248353371207265315478623656127/10401415428935853634424440540325344;
                A(16,6)=131656514265807573955723157408023481433806699348396032656/92668695535091962564795912774190176478892159517481612467;
                A(16,7)=3881494143728609118531066904799685950051960514138645179820/2446349095978358868919950548516272963929118212742344026549;
                A(16,8)=16292266704968075585259245375842819400619822954470178684291/66288722243155885736983218667976563740242178853010092663614;
                A(16,9)=-43986024977384568043684084266385512680544563954/4922783599524658241955780540171948284522386185;
                A(16,10)=285912200202585226675651763671663063668290787/65371192072964016939690070594254881767827200;
                A(16,11)=-6776815256667778089672518929/3693654613173093729492918708;
                A(16,12)=398946554885847045598775476868169/344154261237450078839899047372800;
                A(16,13)=-76630698033396272/4432017119727044925;
                A(16,14)=28401702316003037/1469612686944417840;
                A(16,15)=66049942462586341419969330578128801/12691068622536592094919763114637498325;

                A(17,1)=83940754497395557520874219603241359529066454343054832302344735/64192596456995578553872477759926464976144474354415663868673233;
                A(17,2)=0;
                A(17,3)=892543892035485503125/51401651664490002607536;
                A(17,4)=-12732238157949399705325/686579204375687891972088;
                A(17,5)=5290376174838819557032232941734928484252549/357179779572898187570048915214361602000384;
                A(17,6)=2687322933801750693719999180471745666665021538793817303193221/2863980005760296740624015421425947092438943496681472214589916;
                A(17,7)=-197649786681880330585741729796159873563741413724149351549277865/378029217824623393200881653405474359138017953416246216408422692;
                A(17,8)=-100286075630483975704018828319990067604207336241794360144098685695/20486915674765670626893195919603679319429068544972409068469849579;
                A(17,9)=8739866119696575810411768434844068608106287881671139259/2282122412587168891929052689609009868137678763277087160;
                A(17,10)=-7922242431969626895355493632206885458496418610471389/748272134517487495468365669337985635214015258726400;
                A(17,11)=2777643183645212014464950387658055285/1141545470045611737197667093465955392;
                A(17,12)=-1372659703515496442825084239977218110461/1313121960368535725613950174847107891200;
                A(17,13)=6144417902699179309851023/85608793932459282773805825;
                A(17,14)=140294243355138853053241/64884622846351585391642880;
                A(17,15)=168671028523891369934964082754523881107337/24062875279623260368388427013982199424119600;
                A(17,16)=0;

                c=[0; 1/50; 1/25; 1/10; 2/15; 4/25; 1/20; 1/5; 1/4; 1/3; 1/2; 5/9; 3/4; 6/7; 8437/8926; 1; 1];
                bHat(1,1)=63818747/5262156900;
                bHat(2,1)=0;
                bHat(3,1)=0;
                bHat(4,1)=0;
                bHat(5,1)=0;
                bHat(6,1)=0;
                bHat(7,1)=22555300000000/261366897038247;
                bHat(8,1)=1696514453125/6717619827072;
                bHat(9,1)=-45359872/229764843;
                bHat(10,1)=19174962087/94371046000;
                bHat(11,1)=-19310468/929468925;
                bHat(12,1)=16089185487681/146694672924800;
                bHat(13,1)=1592709632/41841694125;
                bHat(14,1)=52675701958271/4527711056573100;
                bHat(15,1)=12540904472870916741199505796420811396/2692319557780977037279406889319526430375;
                bHat(16,1)=0;
                bHat(17,1)=0;

                bPHat(1,1)=63818747/5262156900;
                bPHat(2,1)=0;
                bPHat(3,1)=0;
                bPHat(4,1)=0;
                bPHat(5,1)=0;
                bPHat(6,1)=0;
                bPHat(7,1)=451106000000000/4965971043726693;
                bPHat(8,1)=8482572265625/26870479308288;
                bPHat(9,1)=-181439488/689294529;
                bPHat(10,1)=57524886261/188742092000;
                bPHat(11,1)=-38620936/929468925;
                bPHat(12,1)=144802669389129/586778691699200;
                bPHat(13,1)=6370838528/41841694125;
                bPHat(14,1)=368729913707897/4527711056573100;
                bPHat(15,1)=111940113324845802831946788738852162520696/1316544263754897771229629968877248424453375;
                bPHat(16,1)=-113178587/12362232960;
                bPHat(17,1)=1/40;

                b(1,1)=27121957/1594593000;
                b(2,1)=0;
                b(3,1)=0;
                b(4,1)=0;
                b(5,1)=0;
                b(6,1)=0;
                b(7,1)=4006163300000/55441463008113;
                b(8,1)=9466403125/25445529648;
                b(9,1)=-163199648/406149975;
                b(10,1)=23359833/69636250;
                b(11,1)=-18491714/140828625;
                b(12,1)=11052304606701/58344472186000;
                b(13,1)=1191129152/44377554375;
                b(14,1)=2033811086741/124730332137000;
                b(15,1)=3616943474975740389660406409450169802/951830146690244407118982233597812374375;
                b(16,1)=0;
                b(17,1)=0;

                bP(1,1)=27121957/1594593000;
                bP(2,1)=0;
                bP(3,1)=0;
                bP(4,1)=0;
                bP(5,1)=0;
                bP(6,1)=0;
                bP(7,1)=4217014000000/55441463008113;
                bP(8,1)=47332015625/101782118592;
                bP(9,1)=-652798592/1218449925;
                bP(10,1)=70079499/139272500;
                bP(11,1)=-36983428/140828625;
                bP(12,1)=99470741460309/233377888744000;
                bP(13,1)=4764516608/44377554375;
                bP(14,1)=14236677607187/124730332137000;
                bP(15,1)=198066487470143918516004831967805004004/2855490440070733221356946700793437123125;
                bP(16,1)=1/50;
                bP(17,1)=0;
                isFSAL=false;
                orders=[12;10];
            otherwise
                error('Unknown Solution Choice');
        end
    otherwise
        error('An invalid order was entered')
end


%Determine subsidiary type.
if(isempty(b))
    subsidType=2;
elseif(isempty(bP))
    subsidType=1;
else
    subsidType=0;
end

if(nargin==2)
    %If the function was just called with two input parameters of the form 
    %Then the first output should just be orders, the second isFSAL, the
    %third is subsidType and there are no other outputs.
    xPredMain=orders;
    xPredSubsid=isFSAL;
    g=subsidType;
else
    numStages=length(c);
    g=zeros(xDim,numStages);

    g(:,1)=dfCur;
    for curStage=2:numStages
        xCur=x+c(curStage)*deltaT*dx+deltaT^2*(g*A(curStage,:)');

        g(:,curStage)=df(xCur,t+c(curStage)*deltaT);
    end

    xNewMain =x+deltaT*dx+deltaT^2*sum(bsxfun(@times,bHat,g'),1)';
    dxNewMain=dx+deltaT*sum(bsxfun(@times,bPHat,g'),1)';
    xPredMain=[xNewMain;dxNewMain];
    
    if(~isempty(b))
        %If coefficients for s subsidiary position were provided.
        xNewSubsid=x+deltaT*dx+deltaT^2*sum(bsxfun(@times,b,g'),1)';
    else
        xNewSubsid=[];
    end

    if(~isempty(bP))
        %If coefficients for a subsidiary velocity were provided.
        dxNewSubsid=dx+deltaT*sum(bsxfun(@times,bP,g'),1)';
    else
        dxNewSubsid=[];
    end
    
    xPredSubsid=[xNewSubsid;dxNewSubsid];
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
