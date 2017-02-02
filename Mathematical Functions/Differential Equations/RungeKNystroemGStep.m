function [xPredMain,xPredSubsid,g,orders,isFSAL,subsidType]=RungeKNystroemGStep(xVec,t,df,deltaT,dfCur,order,solutionChoice)
%%RUNGEKNYSTROEMGSTEP Perform one step of a general explicit Runge-Kutta-
%                     Nyström method. Such methods are used to efficiently
%                     implement second order differential equations of the
%                     form d^2xdt^2=df(xVec,t), where xVec consists of the
%                     stacked values x and dxdt. Runge-Kutta-Nyström
%                     methods can be made to a higher order than general
%                     Runge-Kutta methods while having the same number of
%                     stages. All of the methods implemented here are
%                     embedded methods, meaning that a secondary result of
%                     a different order can also be returned along with
%                     the order of the second output. This second output
%                     can be used in an algorithm that has an adaptive step
%                     size (deltaT is automatically varied). Note that the
%                     function RungeKNystroemStep is made for the
%                     specialized problem where df does not depend on the
%                     first derivative terms.
%
%INPUTS: The function can either be called with two inputs as
%        [orders,isFSAL,subsidType]=RungeKNystroemGStep(order,solutionChoice)
%        to just get the orders of a particular general embedded Runge-
%        Kutta-Nyström method and to figure out whether g(:end) can be
%        passed in place of dfCur during a subsequent step, and to figure
%        out what type of information is contained in xPredSubsid, or it
%        can be run as
%        [xPredMain,xPredSubsid,g,orders,isFSAL]=RungeKNystroemGStep(xVec,t,df,deltaT,dfCur,order,solutionChoice)
%        to actually do a full step and provide the orders. All of the
%        inputs in the full step are:
%        xVec    The value of the vector state over which integration is
%                being performed. The dimensionality must be a multiple of
%                two with the second half of the elements being the
%                derivatives of the first half of the elements [x; dxdt]
%                (e.g. velocity is the second half and position is the
%                first half).
%        t       The time at which xVal is taken.
%        df      df(xVec,t) returns the second dervative of derivative of x
%                with respect to time taken at time t. The output is half
%                the dimensionality of xVec.
%        deltaT  The size of the single (time) step over which the 
%                Runge-Kutta-Nyström integration is performed.
%        dfCur   The value df(xVal,t). This is requested so that
%                methods that are FSAL can pass g(:,end) on subsequent
%                steps instead of having to perform a redundant
%                evaluation of df. If omitted or an empty matrix is
%                passed, the function df(xVal,curT) is evaluated to get
%                dfCur.
%         order  The main integration order of the Runge-Kutta method.
%                If this parameter is omitted, then the default order of
%                5 is used. Order can range from 5 to 7.
%solutionChoice  Different Runge-Kutta formulae exist for some orders. If
%                this parameter is provided, then the selected formula
%                for the given order is used. Otherwise the default
%                (solutionChoice=0) formula is used. The algorithms
%                chosen by pairs of (order, solutionChoice) are:
%                (5,0) RKNG5(4)7FM of Table 14.4 of [2].
%                (5,1) RKNG5(6)9F of Table 7 of [1].
%                (6,0) RKNG6(7)11F of Table 5 of [1].
%                (7,0) RKNG7(8)14F of table 7 of [1].
%                Only method (5,0) provides a full subsidiary output; the
%                value of xPredSubsid for the others is half the
%                dimensionality of xVec, that is, the derivative terms are
%                not computed.
%
%OUTPUTS: If the function is run with two inputs as
%         [orders,isFSAL,subsidType]=RungeKNystroemGStep(order,solutionChoice)
%         Then there are just three outputs, orders, isFSAL and subsidType, 
%         described below.
%         Otherwise, all of the outputs are
%         xPredMain The entire state vector state propagated forward an
%                   interval of deltaT at the order of precision given by
%                   the input order. This is the main integration order.
%         xPredSubsid The subsidiary estimate. This is a different order
%                     than the main integration order can can be used in
%                     algorithms that adaptively adjust the stepsize.
%                     Depending on the combination of (order,
%                     solutionChoice), this is either a full state vector
%                     or just a half a state vector (i.e. the first
%                     derivatives are not always computed in the subsidiary
%                     formula depending on the algorithm chosen).
%                   g The values of the second derivatives df evaluated at
%                     various points as determined by the selected
%                     algorithm. These can sometimes be reused in
%                     interpolation routines to reduce the number of
%                     computations required. In all of
%                     the methods, g(:,1) is equal to df(xVal,t). In
%                     FSAL methods, g(:,end) is df(xPredMain,t+deltaT)
%              orders A 2X1 vector whereby order(1)=order, the input
%                     parameter, and orders(2) is the order of xPredSubsid.
%                     The value of orders(2) depends on the algorithm
%                     chosen by the combination of the order and
%                     solutionChoice inputs.
%              isFSAL Indicates whether the first evaluation of the f of
%                     the next step is equal to g(:,end). If so, this can
%                     be passed to the function on the next step rather
%                     than having to make an additional evaluation of f.
%          subsidType Indicates the nature of the output xPredSubsid.
%                     Possible values are
%                     0 xPredSubsid is a complete state vector, having the
%                       same dimensionality as xVec.
%                     1 xPredSubsid does not contain any first derivatives
%                       and thus is half the dimensionality of xVec.
%                     2 xPredSubsid only contains first derivatives and
%                       thus is half the dimensionality of xVec.
%
%The formulae are given in this file in the form of the components of their
%Butcher tables. The general formula used in discussed in all of the
%references.
%
%The formulae were referred to above using a notation such as RKNG5(4)7. In
%general, the normenclature is
%RKNGq(p)s[F]X
%where
%q is the order of the main integrating formula
%p is the order of the subsidiary formula
%s is the number of stages (derivative function evaluations)
%The presence of an F indicates that the function is first same as last
%(FSAL). That is, f evaluated at the next step is g(:,end) and thus
%g(:,end) can be passed for fCur for the next step. The absence of an F
%indicates that the function is not FSAL.
%
%REFERENCES:
%
%[1] E. Fehlberg, "Classical seventh-, sixth-, and fifth-order Runge-Kutta-
%Nyström formulas with stepsize control for general second-order
%differential equations," National Aeronautics and Space Administration,
%Marshall Space Flight Center, AL, Tech. Rep. NASA TR R-432, Oct. 1974.
%
%[2] J. R. Dormand, Numerical Methods for Differential Equations.
%Raton: CRC Press, 1996.
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

    if(nargin<5||isempty(dfCur))
        dfCur=df(xVec,t);
    end
    
    xDim=length(xVec)/2;
    x=xVec(1:xDim,1);%The position
    dx=xVec((xDim+1):end);%The first derivative 
end

switch(order)
    case 5
        switch(solutionChoice)
            case 0
                %Pg 260 of [2], Table 14.4 and Table 5.4, RKNG5(4)7FM.
                %(Has dense output..., page 265.)
                A = [0,         0,          0,          0,          0,              0,    0;
                     1/5,       0,          0,          0,          0,              0,    0;
                     3/40,      9/40,       0,          0,          0,              0,    0;
                     44/45,     -56/15,     32/9,       0,          0,              0,    0;
                     19372/6561,-25360/2187,64448/6561, -212/729,   0,              0,    0; 
                     9017/3168, -355/33,    46732/5247, 49/176,     -5103/18656,    0,    0;
                     35/384,    0,          500/1113,   125/192,    -2187/6784,     11/84,0];
                c = [0; 1/5; 3/10; 4/5; 8/9; 1; 1];
                bMain = [35/384; 0; 500/1113; 125/192; -2187/6784; 11/84; 0];
                bSubsid = [5179/57600; 0; 7571/16695; 393/640; -92097/339200; 187/2100; 1/40];

                ABar=[0,            0,      0,              0,          0,          0,  0;
                      1/50,         0,      0,              0,          0,          0,  0;
                      9/400,        9/400,  0,              0,          0,          0,  0;
                      8/225,        0,      64/225,         0,          0,          0,  0;
                      1276/59049,   16/81,  37312/295245,   4876/98415, 0,          0,  0;
                      19/132,       9/22,   -14/55,         133/660,    0,          0,  0;
                      35/384,       0,      50/159,         25/192,     -243/6784,  0,  0];

                bBarMain=[35/384;    0;      50/159; 25/192; -243/6784;  0;  0];
                bBarSubsid=[5179/57600;   0;      7571/238500;    393/3200;   -10233/339200;  0;  0];
                orders=[5;4];
                isFSAL=true;
            case 1
                %RKNG5(6)9F of [1], Table 7
                A=[0,            0,      0,      0,              0,              0,              0,              0,      0;
                      4/15,         0,      0,      0,              0,              0,              0,              0,      0;
                      1/10,         3/10,   0,      0,              0,              0,              0,              0,      0;
                      3/20,         0,      9/20,   0,              0,              0,              0,              0,      0;
                      9/40,         0,      0,      27/40,          0,              0,              0,              0,      0;
                      11/48,        0,      0,      5/8,            -5/48,          0,              0,              0,      0;
                      27112/194481, 0,      0,      56450/64827,    80000/194481,   -24544/21609,   0,              0,      0;
                      -26033/41796, 0,      0,      -236575/38313,  -14500/10449,   275936/45279,   228095/73788,   0,      0;
                      7/81,         0,      0,      0,              -250/3483,      160/351,        2401/5590,      1/10,   0];
                c=[0; 4/15; 2/5; 3/5; 9/10; 3/4; 2/7; 1; 1];
                bMain=A(end,:)';
                
                ABar=[0,           0,      0,          0,          0,              0,              0,                  0,      0;
                       8/225,       0,      0,          0,          0,              0,              0,                  0,      0;
                       1/25,        1/25,   0,          0,          0,              0,              0,                  0,      0;
                       9/160,       81/800  9/400,      0,          0,              0,              0,                  0,      0;
                       81/840       0,      729/3200,   81/1600,    0,              0,              0,                  0,      0;
                       11283/88064, 0,      3159/88064, 7275/44032, -33/688,        0,              0,                  0,      0;
                       6250/194481, 0,      0,          0,          -3400/194481,   1696/64827,     0,                  0,      0;
                       -6706/45279, 0,      0,          0,          1047925/1946997,-147544/196209, 1615873/1874886,    0,      0;
                       31/360,      0,      0,          0,          0,              64/585,         2401/7800,          -1/300, 0];

                bBarMain=ABar(end,:)';
                
                bBarSubsid=bBarMain;
                bBarSubsid(7)=0;
                bBarSubsid(8)=bBarMain(7);
                bSubsid=[];%There is no secondary formula for the
                %derivatives.
                orders=[5;6];
                isFSAL=true;
            otherwise
                error('Unknown Solution Choice');
        end
    case 6
        switch(solutionChoice)
            case 0
                %RKNG6(7)11F of [1], Table 5.
                c(1)= 0;
                c(2)= 0.10185185185185185185185185185185;
                c(3)= 0.15277777777777777777777777771778;
                c(4)= 0.22916666666666666666666666666667;
                c(5)= 0.625;
                c(6)= 0.375;
                c(7)= 0.66666666666666666666666666666667;
                c(8)= 0.16666666666666666666666666666667;
                c(9)= 0.97037314832177011063191449465592;
                c(10)=1;
                c(11)=1;

                A=zeros(11,11);
                A(2,1)=   0.10185185185185185185185185185185;
                A(3,1)=   0.38194444444444444444444444444444e-1;
                A(3,2)=   0.11458333333333333333333333333333;
                A(4,1)=   0.57291666666666666666666666666667e-1;
                A(4,2)=   0;
                A(4,3)=   0.171875;
                A(5,1)=   0.81869834710743801652892561983471;
                A(5,2)=   0;
                A(5,3)=  -0.31379132231404958677685950413223e1;
                A(5,4)=   0.29442148760330578512396694214876e1;
                A(6,1)=   0.78409090909090909090909090909091e-1;
                A(6,2)=   0;
                A(6,3)=   0;
                A(6,4)=   0.29066985645933014354066985645933;
                A(6,5)=   0.59210526315789473684210526315789e-2;
                A(7,1)=   0.89637111859334081556303778526001e-1;
                A(7,2)=   0;
                A(7,3)=   0;
                A(7,4)=   0.20414673046251993620414673046252;
                A(7,5)=   0.14243014944769330734243014944769;
                A(7,6)=   0.23045267489711934156378600823045;
                A(8,1)=   0.10180041152263374485596707818930;
                A(8,2)=   0;
                A(8,3)=   0;
                A(8,4)=   0;
                A(8,5)=  -0.31604938271604938271604938271605;
                A(8,6)=   0.14579659024103468547912992357437;
                A(8,7)=   0.23511904761904761904761904761905;
                A(9,1)=  -0.16746873130588455834971263869845e1;
                A(9,2)=   0;
                A(9,3)=   0;
                A(9,4)=  -0.78716193268293956927272248512526e1;
                A(9,5)=   0.53924880763160586316973727446928e1;
                A(9,6)=  -0.15448282105008430450113021691446e1;
                A(9,7)=  -0.32555460065369958589440586873101e1;
                A(9,8)=   0.99245659289317916591142538446550e1;
                A(10,1)= -0.34703198562189991428634815033844e1;
                A(10,2)=  0;
                A(10,2)=  0;
                A(10,4)= -0.15938792782884673788823685509672e2;
                A(10,5)=  0.11404863801986940109108344617441e2;
                A(10,6)= -0.34698562868578472962336158110987e1;
                A(10,7)= -0.74043819346894372049120079607603e1;
                A(10,8)=  0.19941980772362036997897885025319e2;
                A(10,9)= -0.63493713698019674173438857844779e-1;
                A(11,1)=  0.53710681696188942565754159955995e-1;
                A(11,2)=  0;
                A(11,3)=  0;
                A(11,4)=  0;
                A(11,5)=  0;
                A(11,6)=  0.21516847168246719445491147392647;
                A(11,7)=  0.34207435392349983829712168876606;
                A(11,8)=  0.22650729089707197573589454446475;
                A(11,9)=  0.30472907521849356793365990503756;
                A(11,10)=-0.14218987341772151898734177215190;
        
                ABar=zeros(11,11);
                ABar(2,1)=  0.51868998628257887517146776406036e-2;
                ABar(3,1)=  0.58352623456790123456790123456790e-2;
                ABar(3,2)=  0.58352623456790123456790123456790e-2;
                ABar(4,1)=  0.13129340277777777777777777777778e-1;
                ABar(4,2)=  0;
                ABar(4,3)=  0.13129340277777777777777777777778e-1;
                ABar(5,1)=  0.17755681818181818181818181818182e-1;
                ABar(5,2)=  0;
                ABar(5,3)=  0;
                ABar(5,4)=  0.17755681818181818181818181818182;
                ABar(6,1)=  0.19257489669421487603305785123967e-1;
                ABar(6,2)=  0;
                ABar(6,3)=  0.40285268291200777831793874574623e-1;
                ABar(6,4)=  0.10349608525445846020008699434537e-1;
                ABar(6,5)=  0.42013351393188854489164086687307e-3;
                ABar(7,1)=  0.50150891632373113854595336076818e-1;
                ABar(7,2)=  0;
                ABar(7,3)=  0;
                ABar(7,4)=  0.12832080200501253132832080200501;
                ABar(7,5)=  0.14277669482347844920944336149015e-1;
                ABar(7,6)=  0.29472859102488732118361747991378e-1;
                ABar(8,1)=  0.10084019204389574759945130315501e-1;
                ABar(8,2)=  0;
                ABar(8,3)=  0;
                ABar(8,4)=  0;
                ABar(8,5)= -0.20329218106995884773662551440329e-1;
                ABar(8,6)=  0.89555163629237703311777385851460e-2;
                ABar(8,7)=  0.15178571428571428571428571427142e-1;
                ABar(9,1)=  0.58788608928240748733225234592929e-1;
                ABar(9,2)=  0;
                ABar(9,3)= -0.68167418196165528792823095885550;
                ABar(9,4)=  0;
                ABar(9,5)=  0;
                ABar(9,6)=  0.39049081248066910635460201662578e-1;
                ABar(9,7)=  0.13202858101777483299456820865647;
                ABar(9,8)=  0.92261993425952482280294445343764;
                ABar(10,1)= 0.84356399126286632857568756095560e-1;
                ABar(10,2)= 0;
                ABar(10,2)=-0.13999437419511294662384532710268e1;
                ABar(10,4)= 0;
                ABar(10,5)= 0;
                ABar(10,6)= 0;
                ABar(10,7)= 0.15182217068949492865241002704702;
                ABar(10,8)= 0.16612294394931082916619461943874e1;
                ABar(10,9)= 0.25357326422396130665331704968283e-2;
                ABar(11,1)= 0.53595659842036653630856529407254e-1;
                ABar(11,2)= 0;
                ABar(11,3)= 0;
                ABar(11,4)= 0;
                ABar(11,5)= 0;
                ABar(11,6)= 0.13393181352957964677413293047098;
                ABar(11,7)= 0.11449729305181632233053748356721;
                ABar(11,8)= 0.18915603352979407729184141673934;
                ABar(11,9)= 0.79150409147660666995755819490274e-2;
                ABar(11,10)=0.90415913200723327305605786618445e-3;

                bBarMain=ABar(end,:)';
                bBarSubsid=bBarMain;
                bBarSubsid(10)=0;
                bBarSubsid(11)=bBarMain(10);
                bMain=A(end,:)';
                bSubsid=[];%There is no secondar formula for the
                %derivatives.
                orders=[6;7];
                isFSAL=true;
            otherwise
                error('Unknown Solution Choice');
        end
        
    case 7
        switch(solutionChoice)
            case 0
                %RKNG7(8)14F from [1], Table 3.
                c(1)=0;
                c(2)=0.73470804841064383606103566874183e-1;
                c(3)=0.11020620726159657540915535031127;
                c(4)=0.16530931089239486311373302546691;
                c(5)=0.5;
                c(6)=0.26628826929126164263520369439706;
                c(7)=0.63371173070873835736479630560294;
                c(8)=0.75;
                c(9)=0.5625;
                c(10)=0.125;
                c(11)=0.375;
                c(12)=0.96652805160235700834451642119099;
                c(13)=1;
                c(14)=1;
        
                A=zeros(14,14);
                A(2,1)= 0.73470804841064383606103566874183e-1;
                A(3,1)= 0.27551551815399143852288837577819e-1;
                A(3,2)= 0.82654655446197431556866512733456e-1;
                A(4,1)= 0.41327327723098715778433256366728e-1;
                A(4,2)= 0;
                A(4,3)= 0.12398198316929614733529976910018;
                A(5,1)= 0.89670558763795782658378325389875;
                A(5,2)= 0;
                A(5,3)=-0.34585915268314336799411093327799e1;
                A(5,4)= 0.30618859391934758533573260788811e1;
                A(6,1)= 0.57053369423965328229363646644660e-1;
                A(6,2)= 0;
                A(6,3)= 0;
                A(6,4)= 0.20664670695498245793175075857561;
                A(6,5)= 0.25881929123138564740892891767848e-2;
                A(7,1)= 0.22130953402732738534298054285043e-1;
                A(7,2)= 0;
                A(7,3)= 0;
                A(7,4)= 0.36666901842159380713935315648106;
                A(7,5)= 0.32075560532767078702498038199180;
                A(7,6)=-0.75843846443258975333835287154967e-1;
                A(8,1)= 0.83333333333333333333333333333333e-1;
                A(8,2)= 0;
                A(8,3)= 0;
                A(8,4)= 0;
                A(8,5)= 0;
                A(8,6)= 0.38436436964131621037911008488971;
                A(8,7)= 0.28230229702535045628755658177696;
                A(9,1)= 0.83496093750000000000000000000000e-1;
                A(9,2)= 0;
                A(9,3)= 0;
                A(9,4)= 0;
                A(9,5)= 0;
                A(9,6)= 0.38306747479397408845870063561136;
                A(9,7)= 0.13548721270602591154129936438864;
                A(9,8)=-0.39550781250000000000000000000000e-1;
        
                A(10,1)= 0.73420353223593964334705075445816e-1;
                A(10,2)= 0;
                A(10,3)= 0;
                A(10,4)= 0;
                A(10,5)= 0;
                A(10,6)= 0.98808964916022916024205710336420e-1;
                A(10,7)= 0.24153311327327749549842803451955;
                A(10,8)=-0.48707561728395061728395061728395e-1;
                A(10,9)=-0.24005486968449931412894375857339;

                A(11,1)= 0.81378441127067064904041056404207e-2;
                A(11,2)= 0;
                A(11,3)= 0;
                A(11,4)= 0;
                A(11,5)= 0;
                A(11,6)= 0;
                A(11,7)=-0.36266091174647134384031532058792;
                A(11,8)= 0.69726880597127928317272609847243e-1;
                A(11,9)= 0.37797780620763392161154341509711;
                A(11,10)=0.28181838082900278742109519000315;
        
                A(12,1)=-0.14042538922482838913280031225476e1;
                A(12,2)= 0;
                A(12,3)= 0;
                A(12,4)= 0;
                A(12,5)= 0;
                A(12,6)=-0.13555559029404957528304113342361e2;
                A(12,7)=-0.15021472824848050961721330969968e1;
                A(12,8)= 0.14767543284167949686233606841588e1;
                A(12,9)=-0.21707681965133688432577373607995e1;
                A(12,10)=0.66149759502676558681039202833030e1;
                A(12,11)=0.11507526173569321530679222376434e2;

                A(13,1)= -0.52708651815801315268176882187497e1;
                A(13,2)=  0;
                A(13,3)=  0;
                A(13,4)=  0;
                A(13,5)=  0;
                A(13,6)= -0.49965599553656833001045921529105e2;
                A(13,7)= -0.50302228928658231516135124812231e1;
                A(13,8)=  0.44548269045298760506518238622704e1;
                A(13,9)= -0.86071533124033841312406742989148e1;
                A(13,10)= 0.23840410046372287590078676456468e2;
                A(13,11)= 0.41711581466028388124069667164840e2;
                A(13,12)=-0.13297747642437995408237095558512;

                A(14,1)=  0.35099303056581883152660173681744e-1;
                A(14,2)=  0;
                A(14,3)=  0;
                A(14,4)=  0;
                A(14,5)=  0;
                A(14,6)=  0;
                A(14,7)=  0;
                A(14,8)=  0.25223475276631606400638853417712;
                A(14,9)=  0.11840033306876549234162515364336;
                A(14,10)= 0.20258133611250929893187899871888;
                A(14,11)= 0.26757025259420140796393329272621;
                A(14,12)= 0.16586384510629873791268098150965;
                A(14,13)=-0.41749822704672884309167134456960e-1;
        
        
                ABar=zeros(14,14);
                ABar(2,1)= 0.26989795819968848329994970508715e-2;
                ABar(3,1)= 0.30363520297464954371244341822304e-2;
                ABar(3,2)= 0.30363520297464954371244341822304e-2;
                ABar(4,1)= 0.68317920669296147335299769100184e-2;
                ABar(4,2)= 0;
                ABar(4,3)= 0.68317920669296147335299769100184e-2;
                ABar(5,1)=-0.10263757731977888994310872824217e-2;
                ABar(5,2)= 0;
                ABar(5,3)= 0;
                ABar(5,4)= 0.12602637577319778889943108728242;
                ABar(6,1)= 0.98909903843107417913313499064241e-2;
                ABar(6,2)= 0;
                ABar(6,3)= 0.20401758759111349514170518498571e-1;
                ABar(6,4)= 0.50265147713328703261825104735338e-2;
                ABar(6,5)= 0.13545726631277755415728360014730e-3;
                ABar(7,1)= 0.36772464695317721429741572462201e-1;
                ABar(7,2)= 0;
                ABar(7,3)= 0;
                ABar(7,4)= 0.82132294778521785827721741407693e-1;
                ABar(7,5)= 0.30087165409098963036870918119641e-1;
                ABar(7,6)= 0.51803353935993790519824105531789e-1;
                ABar(8,1)= 0.41233049088272873123221004021091e-1;
                ABar(8,2)= 0;
                ABar(8,3)= 0;
                ABar(8,4)= 0.11335100293061819105328798078376;
                ABar(8,5)= 0.56722148592237668841301677436715e-1;
                ABar(8,6)= 0.57456202064954525469376924736474e-1;
                ABar(8,7)= 0.12487597323916741512812413021961e-1;
                ABar(9,1)= 0.42146301269531250000000000000000e-1;
                ABar(9,2)= 0;
                ABar(9,3)= 0;
                ABar(9,4)= 0;
                ABar(9,5)=-0.78088073730468750000000000000000e-1;
                ABar(9,6)= 0.14104682102928772004397085536135;
                ABar(9,7)= 0.74603813736337279956029144638648e-1;
                ABar(9,8)=-0.21505737304687500000000000000000e-1;
        
                ABar(10,1)= 0.55243877171925011431184270690444e-2;
                ABar(10,2)= 0;
                ABar(10,3)= 0;
                ABar(10,4)= 0;
                ABar(10,5)= 0;
                ABar(10,6)= 0.45913375893505158838018111807029e-2;
                ABar(10,7)= 0.12009956992268139808927955623138e-1;
                ABar(10,8)=-0.24361818415637860082304526748971e-2;
                ABar(10,9)=-0.11877000457247370827617741197988e-1;

                ABar(11,1)=  0.12396099092300428073855581211091e-1;
                ABar(11,2)=  0;
                ABar(11,3)=  0;
                ABar(11,4)=  0;
                ABar(11,5)=  0;
                ABar(11,6)=  0;
                ABar(11,7)= -0.23148568834881149606828637484335e-1;
                ABar(11,8)=  0.44057716338592294670599538200368e-2;
                ABar(11,9)=  0.24164236870396086181891828927171e-1;
                ABar(11,10)= 0.52494961238325405884021273526037e-1;

                ABar(12,1)= -0.12148292337172366838692371706654;
                ABar(12,2)=  0;
                ABar(12,3)=  0;
                ABar(12,4)= -0.15948786809469047245658868763595e1;
                ABar(12,5)=  0.77089844409590354601143580630038e-1;
                ABar(12,6)=  0;
                ABar(12,7)=  0;
                ABar(12,8)=  0.98844932135442618048624328441900e-1;
                ABar(12,9)= -0.18517690177654009760124559303975;
                ABar(12,10)= 0.16665727117807342381867154630279e1;
                ABar(12,11)= 0.52611925503652552568040599022802;

                ABar(13,1)= -0.49475846764102332689603709547782;
                ABar(13,2)=  0;
                ABar(13,3)=  0;
                ABar(13,4)= -0.56513209641364305307648232070852e1;
                ABar(13,5)=  0.42750028729043677987389324310306;
                ABar(13,6)=  0;
                ABar(13,7)=  0;
                ABar(13,8)=  0.30293416726956828108567954376506;
                ABar(13,9)= -0.10280329379503342611614151091571e1;
                ABar(13,10)= 0.54254171279669182157854764162220e1;
                ABar(13,11)= 0.15340242607867031086671199895920e1;
                ABar(13,12)=-0.15763473585838266589893780962028e-1;

                ABar(14,1)=  0.35173987586306713954725819037907e-1;
                ABar(14,2)=  0;
                ABar(14,3)=  0;
                ABar(14,4)=  0;
                ABar(14,5)=  0;
                ABar(14,6)=  0;
                ABar(14,7)=  0;
                ABar(14,8)=  0.63858784354258308506882892800822e-1;
                ABar(14,9)=  0.50866724905581448754291007148603e-1;
                ABar(14,10)= 0.17703179472766752427031494226269;
                ABar(14,11)= 0.16781715613041509463911067215393;
                ABar(14,12)= 0.45385629257942440722375392950401e-2;
                ABar(14,13)= 0.71298936997666580243712730101115e-3;

                bBarMain=ABar(end,:)';
                bBarSubsid=bBarMain;
                bBarSubsid(13)=0;
                bBarSubsid(14)=bBarMain(13);

                bMain=A(end,:)';
                bSubsid=[];%There is no secondar formula for the
                %derivatives.
                orders=[7;8];
                isFSAL=true;
            otherwise
                error('Unknown Solution Choice');
        end
    otherwise
        error('Unknown order provided')
end

%Determine subsidiary type.
if(isempty(bBarSubsid))
    subsidType=2;
elseif(isempty(bSubsid))
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
    %%%%Use the selected coefficients.
    numStages=length(c);
    g=zeros(xDim,numStages);

    g(:,1)=dfCur;
    for curStage=2:numStages
        xCur=x+c(curStage)*deltaT*dx+deltaT^2*(g*ABar(curStage,:)');
        dxCur=dx+deltaT*(g*A(curStage,:)');
        xVecCur=[xCur;dxCur];

        g(:,curStage)=df(xVecCur,t+c(curStage)*deltaT);
    end

    Phi=dx+deltaT*sum(bsxfun(@times,bBarMain,g'),1)';
    dPhi=sum(bsxfun(@times,bMain,g'),1)';

    xNewMain =x+deltaT*Phi;
    dxNewMain=dx+deltaT*dPhi;

    xPredMain=[xNewMain;dxNewMain];

    Phi=dx+deltaT*sum(bsxfun(@times,bBarSubsid,g'),1)';
    xNewSubsid=x+deltaT*Phi;
    
    if(~isempty(bSubsid))
        %If the coefficients for a subsidiary second deriavtive were
        %provided.
        dPhi=sum(bsxfun(@times,bSubsid,g'),1)';
        dxNewSubsid=dx+deltaT*dPhi;
        xPredSubsid=[xNewSubsid;dxNewSubsid];
    else
        xPredSubsid=xNewSubsid;
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
