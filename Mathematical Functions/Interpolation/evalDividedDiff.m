function val=evalDividedDiff(t,f,fDerivs)
%%EVALDIVIDEDDIFF Given a set of knots t and function values, f, evaluate
%               the divided difference f[t(1),...,t(end)]. The notation
%               using brackets is common when discussing divided
%               differences. If a value of t is repeated r times, then it
%               is assumed that fDerivs(t(i),r) returns the rth derivative
%               of the function. Divided differences are used in various
%               interpolation routines. 
%
%INPUTS: t NX1 or 1XN set of points (knots) at which the function values
%          are available.
%        f The function values at the points in t. This can either be an
%          NX1 or 1XN vector,or a function handle that can take the t
%          values as a parameter and return a scalar.
%  fDerivs This parameter is only needed if x values are repeated. This is
%          a function handle or vector such that fDerivs(t(i),r) returns
%          the rth derivative of the function evaluated at t(i). Repeated
%          knots and derivatives arise in hermite interpolation.
%
%OUTPUTS: val The value of the divided difference f[t(1),...,t(end)].
%
%Divded differences in terms of Newton's interpolation formula are
%discussed in Chapter 2.1.3 of [1]. The expression including allowing for
%repeated knots is given in terms of B-Splines in Equation 2.4.4.1 in
%Chapter 2.4.4 of [1].
%
%EXAMPLE:
%This example uses the values given in Table 3.1.1 of [2].
% x=[1;1.3;1.6;1.9;2.2];
% f=[0.7651977;0.6200860;0.4554022;0.2818186;0.1103623];
% val=evalDividedDiff(t,f,fDerivs)
%One will get val=0.001825102880660.
%
%REFERENCES:
%[1] J. Stoer and R. Burlisch, Introduction to Numerical Analysis, 2nd ed.
%    New York: Springer-Verlag, 1991.
%[2] R. L. Burden and J. D. Faires, Numerical Analysis, 9th ed. Boston, MA:
%    Brooks/Cole, 2011.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    fDerivs=[];
end

%Make sure that the points are sorted.
t=sort(t,'ascend');

nMax=length(t);

inLeftBranch=true(nMax,1);
idxMin=zeros(nMax,1);
idxMax=zeros(nMax,1);
leftValues=NaN(nMax,1);
goingDown=true;

idxMin(1)=1;
idxMax(1)=nMax;
%The divided difference value from the previous level down.
prevLevelRet=[];

curLevel=1;
while(curLevel>0)
    if(goingDown)
        %This is only for the case of repeated knots, indicating fixed
        %derivatives.
        if(~isempty(fDerivs)&&all(t(idxMin(curLevel))==t(idxMin(curLevel):idxMax(curLevel))))
            r=idxMax(curLevel)-idxMin(curLevel);
            prevLevelRet=fDerivs(t(idxMin(curLevel)),r)/factorial(r);
            
            goingDown=false;
            curLevel=curLevel-1;
            continue;
        end

        %If we are at the bottom level of the recursion.
        if(curLevel==nMax)
            prevLevelRet=f(idxMin(curLevel));
            goingDown=false;
            curLevel=curLevel-1;
            continue;
        else%We are not at the bottom level of the recursion.
            %We have to go down the left side of the recursion for this
            %level.
            inLeftBranch(curLevel+1)=true;
            idxMin(curLevel+1)=idxMin(curLevel)+1;
            idxMax(curLevel+1)=idxMax(curLevel);
            curLevel=curLevel+1;
            continue;
        end
    else
        if(inLeftBranch(curLevel)==false)
            %If we have both halves of the f values for the divided
            %difference at this level, then evaluate it and go up another
            %level.
            tDiff=t(idxMax(curLevel))-t(idxMin(curLevel));
            prevLevelRet=(leftValues(curLevel)-prevLevelRet)/tDiff;
            curLevel=curLevel-1;
            continue;
        else
            %If we only have the left-half of the divided difference for
            %this level, store it and go down the other half.
            leftValues(curLevel)=prevLevelRet;
            inLeftBranch(curLevel)=false;
            inLeftBranch(curLevel+1)=true;
            
            idxMin(curLevel+1)=idxMin(curLevel);
            idxMax(curLevel+1)=idxMax(curLevel)-1;
            goingDown=true;
            curLevel=curLevel+1;
            continue;
        end
    end
end

val=prevLevelRet;
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
