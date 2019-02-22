function E=solveKeplersEq(M,e,eqSel)
%%SOLVEKEPLERSEQ  Solve Kepler's equation for elliptical/ parabolic or
%                 hyperbolic orbits in different forms. If the eqSel
%                 parameter is not given and  e<=1, then the function
%                 solves the equation M=E-e*sin(E) for E by default. If
%                 e>1, then the function solves the equation M=e*sinh(E)-E
%                 by default. If desired, however, the user can choose
%                 which function to solve including the additional
%                 functions M=E-(1-e)*sin(E) and M=E+(e-1)*asinh(E).
%                 These equations have to be solved iteratively as no
%                 closed-form solutions exist. These equation arises when
%                 propagating orbits under ideal Keplerian dynamic models.
%
%INPUTS: M  The left-hand side of the equation, typically the mean anomaly.
%        e  The eccentricity value in the equation such that e>=0.
%     eqSel An optional parameter specifying which type of equation to
%           solve. If omitted, the equation M=E-e*sin(E) is solved for E if
%           0<=e<=1 and the equation M=e*sinh(E)-E is solved otherwise.
%           Possible values of eqSel are
%           0 Solve the equation M=E-e*sin(E)
%           1 Solve the equation M=e*sinh(E)-E (the hyperbolic equation)
%           2 Solve the equation M=E-(1-e)*sin(E)
%           3 Solve the equation M=E+(e-1)*asinh(E)
%
%OUTPUTS:    E   The solution to the selected equation.
%
%It is assumed that if eqSel=0 or eqSel=2 that 0<=e<=1.
%
%The solution to Kepler's equations is described in Chapter 4.1 of [1]. The
%implementation used is taken from [2].
%
%REFERENCES:
%[1] G. Beutler, Methods of Celestial Mechanics: Physical, Mathematical and
%    Numerical Principles. Berlin: Springer, 2005, vol. 1.
%[2] R. H. Gooding and A. W. Odell, "The hyperbolic Kepler's equation,
%    and the elliptic equations revisited," Royal Aerospace Executive,
%    Procurement Executive, Ministry of Defence, Farnborough, Hants, United
%    Kingdom, Tech. Rep. 369, Jul. 1989.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
   if(e<=1)
       eqSel=0;
   else
       eqSel=1;
   end
end

switch(eqSel)
    case 0%Solve Kepler's traditional equation for elliptical orbits
        E=solveKeplerGooding(M,e);
    case 1%Solve Kepler's equation for hyperbolic orbits.
        eL=M/e;
        g=1/e;
        g1=1-g;
        s=solveKeplerHypGooding(eL,g1);
        E=asinh(s);
    case 2
        E=solveKeplerGooding(M,1-e);
    case 3
        E=solveKeplerHypGooding(M,e);
    otherwise
        error('Invalid Value of eqSel provided.')
end
end


function E=solveKeplerGooding(M,e)
%%SOLVEKEPLERGOODING This implements Gooding's algorithm for solving
%                    Kepler's equation for elliptical orbits,
%                    M=E-(1-e1)*sin(E).
%
%INPUTS: M  The left-hand side of the equation, typically the mean anomaly.
%        e  The eccentricity value in the equation such that e<=0<=1.
%
%OUTPUTS: The value of E such that M=E-(1-e1)*sin(E).
%
%The algorithm is taken from
%R. H. Gooding and A. W. Odell, "The hyperbolic Kepler's equation,
%and the elliptic equations revisited," Royal Aerospace Executive,
%Procurement Executive, Ministry of Defence, Farnborough, Hants, United
%Kingdom, Tech. Rep. 369, Jul. 1989.
%Note that the algorithm described in the report uses 1-e instead of e. 
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

e1=1-e;
sw=0.25;

%Reduce M to be in the range -pi to pi.
EMR=mod(M,2*pi);
if(EMR<-pi)
    EMR=EMR+2*pi;
elseif(EMR>pi)
    EMR=EMR-2*pi;
end
EE=EMR;
if(EE==0)
    E=EE+M-EMR;
    return;
end

if(EE<0)
    EE=-EE;
end
e=1-e1;%The eccentricity

%EMR is range reduced EM and EE is absolute value of EMR.
%Starter by first solving the subic equation.
W=solveCubicEqGoodin(e,2*e1,3*EE);
EE=(EE^2+(pi-EE)*W)/pi;
if(EMR<0)
   EE=-EE; 
end

%Do iterations of Halley's method, each followed by Newton's method.
%Gooding used two iterations. Here, I go until convergence, which seems to
%always be in 5 or fewer iterations, but I allow to up to 30 just in case.
%At any rate, it seems to be more than the 2 iterations Gooding uses.
EEPrev=Inf;
for iter=1:30
    fDD=e*sin(EE);%Second derivative
    fDDD=e*cos(EE);%Thirs derivative
    
    if(EE^2/6+e1>sw)
        f=(EE-fDD)-EMR;
        fD=1-fDDD;%First Derivative
    else
        f=GoodingSinDiffFun(e1,EE)-EMR;
        fD=2*e*sin((1/2)*EE)^2+e1;
    end
    
    %This is the iterative update with the modifications to reduce
    %underflow problems.
    dEE=f*f/((1/2)*f*fDD-fD^2);
    if(dEE==0)
        break;
    end

    W=fD+(1/2)*dEE*(fDD+(1/3)*dEE*fDDD);
    fD=fD+dEE*(fDD+(1/2)*dEE*fDDD);
    EE =EE-(f-dEE*(fD-W))/fD;
    %Without the modfications given by Gooding to reduce underflow
    %problems, the iterative update for EE would be
    %f=f+dEE*(fE+(1/2)*dEE*(fDD+(1/3)*dEE*fDDD))
    %EE=EE+dEE-f/fD;
    
    %Check for convergence.
    if(abs(EEPrev-EE)<=eps(EE))
        break;
    end
    EEPrev=EE;
end

E=EE+M-EMR;
end


function x=solveCubicEqGoodin(a,b,c)
%SOLVECUBICEQGOODING  This solves the equation a*x^3+3*b*x-2*c=0, where
%                     a>=0 and b^3+a*c^2>=0 for the real root of x.
%                     Additionally, if a and b are both zero, then zero is
%                     generated in lieu of an indeterminate solution.
%
%INPUTS: a,b,c The coefficients in the equation a*x^3+3*b*x-2*c=0 where
%              a>=0 and b^3+a*c^2>=0.
%
%OUTPUTS: x    The one real solution to a*x^3+3*b*x-2*c=0.
%
%In the report
%R. H. Gooding and A. W. Odell, "The hyperbolic Kepler's equation,
%and the elliptic equations revisited," Royal Aerospace Executive,
%Procurement Executive, Ministry of Defence, Farnborough, Hants, United
%Kingdom, Tech. Rep. 369, Jul. 1989.
%the source of the algorithm is listed as being in 
%R. H. Gooding, "Solution of the hyperbolic Kepler's equation," Royal
%Aerospace Executive, Procurement Executive, Ministry of Defence,
%Farnborough, Hants, United Kingdom, Tech. Rep. 87042, 1987.
%However, the 1987 report appears to be unobtainable. Nonetheless, the
%solution is given in Equation 4.47 in
%G. Beutler, Methods of Celestial Mechanics: Physical, Mathematical and
%Numerical Principles. Berlin: Springer, 2005, vol. 1.
%The modification to return zero instead of infinity is consistent with how
%Gooding used the algorithm.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(a==0&&b==0||c==0)
    x=0;
else
    d=nthroot(sqrt(a)*abs(c)+sqrt(b^3+a*c^2),3)^2;
    x=2*c/(d+b+b^2/d);
end
end

function x=GoodingSinDiffFun(e,EE)
%%GOODINGSINDIFFFUN Evaluate the function EE-(1-e)*sin(EE) using Gooding's
%                   EMKEP procedure for when e and EE are close
%                   to (1,0). This is supposed to be more accurate than
%                   just directly evaluating the functions, unless EE is
%                   large as it is then supposed to worsen rounding errors.
%
%
%The algorithm is the EMKPL algorithm taken from Appendix C of
%A. W. Odell and R. H. Gooding, "Procedure for solving Kepler's
%equation," Celestial Mechanics, vol. 38, no. 4, pp. 307-334, Apr. 1986.
%modified to solve EE-(1-e)*sin(EE) instead of EE-e*sin(EE).
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

x=e*sin(EE);
EE2=-EE^2;
term=EE;
d=0;
while(1)
    d=d+2;
    term=term*EE2/(d*(d+1));
    x0=x;
    x=x-term;
    if(x==x0)
        break;
    end
end

end


function s=solveKeplerHypGooding(eL,g1)
%%SOLVEKEPLERHYP This implement's Gooding's algorithm to the hyperbolic
%                Kepler equation M=sinh(E)-(1-g1)*E for E. This is related
%                to the equation M=e*sinh(E)-E through a transformation of
%                variables.
%
%The algorithm is taken from
%R. H. Gooding and A. W. Odell, "The hyperbolic Kepler's equation,
%and the elliptic equations revisited," Royal Aerospace Executive,
%Procurement Executive, Ministry of Defence, Farnborough, Hants, United
%Kingdom, Tech. Rep. 369, Jul. 1989.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

sw=1/2;

s=eL;
if(eL~=0)
    %Starter based on Lagrange's theorem
    g=1-g1;
    cl=sqrt(1+eL^2);
    al=asinh(eL);
    w=g^2*al/cl^3;
    s=1-g/cl;
    s=eL+g*al/nthroot(s^3+w*eL*(1.5-g/0.75),3);
    
    %Iterate Halley-then Newton process.
    for iter=1:10
        s0=s^2;
        s1=s0+1;
        s2=sqrt(s1);
        s3=s1*s2;
        fDD=g*s/s3;%Second derivative
        fDDD=g*(1-2*s0)/(s1*s3);%Third derivative
        
        if(s0/6+g1>=sw)
            f=(s-g*asinh(s))-eL;
            fD=1-g/s2;%First derivative
        else
            f=GoodingHyperSinDiffFun(g1,s)-eL;
            fD=(s0/(s2+1)+g1)/s2;
        end
        dS=f*fD/((1/2)*f*fDD-fD^2);
        sTemp=s+dS;
        if(sTemp==s)
            break;
        end
        f=f+dS*(fD+(1/2)*dS*(fDD+(1/3)*dS*fDDD));
        fD=fD+dS*(fDD+(1/2)*dS*fDDD);
        s=sTemp-f/fD;
    end
end


end

function x=GoodingHyperSinDiffFun(g1,s)
%%GOODINGHYPERSINDIFFFUN Evaluate the function s-(1-g1)*asinh(s) when
%                        (g1,s) is close to (0,0) using Gooding's method.
%                        This is supposed to have a higher precision than
%                        just explicitly evaluating the function.
%
%The algorithm is the SHMKEP function taken from Appendix B of 
%R. H. Gooding and A. W. Odell, "The hyperbolic Kepler's equation,
%and the elliptic equations revisited," Royal Aerospace Executive,
%Procurement Executive, Ministry of Defence, Farnborough, Hants, United
%Kingdom, Tech. Rep. 369, Jul. 1989.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

g=1-g1;
t=s/(1+sqrt(1+s^2));
x=s*(g1+g*t^2);
term=2*g*t;
twoI1=1;
%Iterate until convergence
while(1)
    twoI1=twoI1+2;
    term=term*t^2;
    x0=x;
    x=x-term/twoI1;
    if(x==x0)
        break;
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
