function [xi,w,Hm,algParams]=KorobovLatticeRules(numDim,method,params,periodizeAlg,lowerBounds,upperBounds)
%%KOROBOVLATTICERULES Obtain lattice points for multidimensional numerical
%               integration based on the Korobov lattice rules, which are
%               sometimes also just referred to as "number theoretic"
%               lattice rules. The standard points are meant for
%               integrating a function that is periodic in all dimensions
%               with the period being on the region 0-1 in all dimensions.
%               Inputs to this function can perform a transformation to
%               periodize the points (for use with non-periodic functions,
%               though periodization is not always necessary), and can
%               change the integration region to other finite or infinite
%               bounds.
%
%INPUTS: numDim The number of dimensions of the lattice points that are to
%               be generated. numDim>=2. Valid values depend on the method
%               chosen.
%        method The method of generating the lattice points that should be
%               used. Possible values are
%               0 Use the Korobov rules based on Theorem 6.3-5 in Chapter
%                 6.3 of [1]. These rules are based on minimizing a cost
%                 function H(m) over an integer m on the interval p-1 where
%                 p is a prime number. In this option, one selects a
%                 tabulated rule. The value of params is the index of the
%                 tabulated rule (increasing params increases the number of
%                 lattice points) and can vary from 1 to a maximum number
%                 depending on numDim. Possible combinations of numDim and
%                 the maximum value of params are: (3,28), (4,24), (5,18),
%                 (6,14). In all instances, the highest index of params
%                 corresponds to a lattice with 10007 points. Only
%                 numDim 3-6 are supported.
%               1 Use the Korobov rules based on Theorem 6.3-5 in Chapter
%                 6.3 of [1], but unlike in method 0, params is a prime
%                 number that is the number of lattice points to generate
%                 and this function actually performs the optimization of
%                 Theorem 6.3-5, producing results for dimensions and
%                 numbers of points not tabulated in method 0. Note that
%                 the upper limit of the sum in Equation 6.3-6 should be p
%                 and not m.
%               2 Use the Korobov rules based on Theorem 6.3-6 in Chapter
%                 6.3 of [1]. These rules are based on minimizing a cost
%                 function H*(m) that is more complicated than that used in
%                 methods 0 and 1. These formula have more points than
%                 those in methods 0 and 1. In this option, one selects a
%                 tabulated rule. The tabulated rules are taken from the
%                 tables in Chapter 6.3 of [1]. However, method 3 can often
%                 find better rules (because it does not appear that the
%                 authors of the table optimized over a). However method 3
%                 can also be extremely slow due to the joint optimization
%                 over parameters a and b. The value of params is the indx
%                 of the tabulated rule (increasing params increases the
%                 number of lattice points) and can vary from 1 to a
%                 maximum number depending on numDim. Possible combinations
%                 of numDim and the maximum value of params along with the
%                 number of lattice points in the rule for the maximum
%                 value of params are (3,6,100063), (4,6,100063),
%                 (5,6,100063), (6,6,100063), (7,8,100063), (8,6,100063),
%                 (9,6,159053), (10,6,155093). Only numDim 3-10 are
%                 supported.
%               3 Use the Korobov rules based on Theorem 6.3-6 in Chapter
%                 6.3 of [1] but unlike in method 0 params is a structure
%                 holding two prime numbers and two integer bound,
%                 params.p, params.q, params.aMin, params.aMax and an
%                 optimization is performed to determine the lattice
%                 points. There will be p*q lattice points. a is bounded by
%                 aMin and aMax for a brute-force search for an optimal
%                 value of an integer parameter a. Usually p>a>q. If aMin
%                 and aMax are omitted, then aMin=min(p,q) and aMax=
%                 max(p,q) are used. This allows for the computation of
%                 results not tabulated in [1].
%               4 In this method, params is the final output parameter of
%                 the function algParams fed back in. This mode of the
%                 function regenerates points that might have taken a very
%                 long time to optimize. algParams will contain only a few
%                 parameters that are independent of the number of
%                 dimensions of the problem, and are thus more efficient to
%                 store than saving found points.
%               5 Use the Fibonacci grid method described in Chapter 6.2 of
%                 [1] and Chapter 4.3 of [2]. This has the same format as
%                 the Korobov rules, as in Equation 6.3-4 of [1].This
%                 requires that numDim=2 and exhibits a certain optimality.
%                 In this instance, params>1 is the index of the number in
%                 the Fibonacci sequence to use. Thus, the total number of
%                 points in the lattice is FibonacciNum(params). As noted
%                 in Chapter 8.6 of [2], this grid method will be better if
%                 params is odd.
%  periodizeAlg If this parameter is omitted or an empty matrix is passed,
%               then no periodization of the points will be performed.
%               Periodization of the points can improve performance when
%               used for integration over non-periodic functions. If this
%               parameter is provided, then it is a sturcture with element
%               periodizeAlg.method and optionally periodizeAlg.param1 and
%               periodizeAlg.param2. The meaning of the parameters is the
%               same as in the function periodizeLatticePoints.
%lowerBounds,upperBounds If these p[arameters are provided (and are not
%               empty matrices), then the lattice points will be
%               transformed for integration over regions other than just
%               +/-1 in all dimensions. These are numDimX1 or 1XnumDim
%               vectors where lowerBounds are the lower integration bounds
%               and upperBounds are the upper integration bounds. If
%               transformed, the following substitutions are used:
%               Considering only the current dimensions of xi where a is
%               the lower bound and b is the upper bound:
%               * If a and b are finite:
%                 xi(curDim,:)=(b-a)*xi(curDim,:)+a, w=(b-a)*w
%               * If a=-Inf and b is finite:
%                 xi(curDim,:)=-1./xi(curDim,:)+b+1; w=w./xi(curDim,:).^2;
%                 --Note this may cause problems if any elements of xi are
%                 zero.
%               * If a is finite and b=Inf:
%                 xi(curDim,:)=1./(1-xi(curDim,:))+a-1;
%                 w=w./(1-xi(curDim,:)).^2;
%                 --Note this may cause problems if any elements of xi are
%                 1.
%               * If a=-Inf and b=Inf:
%                 xi(curDim,:)=1./(1-xi(curDim,:))-1./xi(curDim,:);
%                 w=w.*(1./(1-xi(curDim,:)).^2+1./xi(curDim,:).^2);
%                 --Note this may cause problems if any elements of xi are
%                 zero or 1.
%               * Same as any above, but the lower bound is greater than
%                 the upper bound:
%                 The bounds are flipped, but the same substitution is
%                 made, just w is multiplied by -1.
%
%OUTPUTS: xi A numDimXnumPoints matrix of the lattice points generated
%            using the selected method.
%          w A numPointsX1 set of weights associated with the lattice
%            points. If no periodization or transformation of the
%            integration region took place, then these are all just
%            1/numPoints.
%         Hm The value of the cost function used in methods 0-3. In all
%            instances, values closer to 1 are better. This is useful for
%            comparing the quality of points generated by the same
%            algorithm.
%  algParams A structure containing the parameters used to generate the
%            points. This can be passed to method 4 to reuse optimal
%            solutions produced by algorithms 1 or 3, which might take a
%            very long time to compute.
%
%EXAMPLE 1:
%Consider the examples in 4D of Section 6.6 of [1]. We start with the
%problem of integrating over the 4-dimensional cube of the function
%exp(x1*x2*x3*x4). Let us consider a few variants, varying the number of
%points used and periodizing versus not periodizing:
% %The function
% f=@(x)exp(prod(x,1));
% %Integrated over 0->1 in all dimensions, for 4 dimensions has an "exact"
% %solution of
% exactSol=1.0693976088597706235;
% %(Which comes from evaluating a hypergeometric function)
% %First, consider using tabulated formulae without periodization using a
% %low-order formula.
% numDim=4;
% method=0;
% params=1;
% [xi,w]=KorobovLatticeRules(numDim,method,params);
% est1=sum(bsxfun(@times,f(xi),w'));
% error1=abs(exactSol-est1)
% %Now a low-order formula with periodization,
% periodizeAlg=[];
% periodizeAlg.method=0;
% [xi,w]=KorobovLatticeRules(numDim,method,params,periodizeAlg);
% est2=sum(bsxfun(@times,f(xi),w'));
% error2=abs(exactSol-est2)
% %Now with the highest-order tabulated formula, without periodization.
% params=24;
% [xi,w]=KorobovLatticeRules(numDim,method,params);
% est3=sum(bsxfun(@times,f(xi),w'));
% error3=abs(exactSol-est3)
% %And finally, a high-order tabulated formula with periodization.
% [xi,w]=KorobovLatticeRules(numDim,method,params,periodizeAlg);
% est4=sum(bsxfun(@times,f(xi),w'));
% error4=abs(exactSol-est4)
%One will observe that the error decreases with the use of periodization as
%well as with the use of an increasing number of points.
%
%EXAMPLE 2:
%Method 0 is just pre-tabulated values that can be computed using method 1.
%Here, we demonstrate how  method 1 produces the same set of lattice points
%as tabulated im method 0 when given the (prime) number of points.
% numDim=4;
% method=0;
% params=3;
% xi1=KorobovLatticeRules(numDim,method,params);
% method=1;
% params=size(xi1,2);
% xi2=KorobovLatticeRules(numDim,method,params);
%One will observe that xi1 and xi2 are the same. However, this is not
%always the case. The primary reason is because multiple optimal solutions
%for a given prime number can exist and method 1 only chooses the first
%one. Also, finite precision errors can play a role.
%
%EXAMPLE 3:
%Here, we consider integration over a region that is not the unit cube. We
%consider the function sum(abs(x-1/2)) in five dimensions with bounds as
%follows:
% f=@(x)sum(abs(x-1/2),1);
% lowerBounds=[-2;-1;0;4;3];
% upperBounds=[2;1;1;6;8];
% %The exact solution against which we will compare is easily found:
% exactSol=915;
% numDim=5;
% method=0;
% params=1;%The first tabulated value uses 1069 points.
% [xi,w]=KorobovLatticeRules(numDim,method,params,[],lowerBounds,upperBounds);
% est1=sum(bsxfun(@times,f(xi),w'));
% error1=abs(exactSol-est1)
% %Now, we again consider periodizing the points.
% periodizeAlg=[];
% periodizeAlg.method=0;
% [xi,w]=KorobovLatticeRules(numDim,method,params,periodizeAlg,lowerBounds,upperBounds);
% est2=sum(bsxfun(@times,f(xi),w'));
% error2=abs(exactSol-est2)
% %Here, we consider the use of Koborov type II points instead. By
% %necessity, this uses many more points:
% method=2;
% params=1;%The first tabulated value uses 15019 points.
% [xi,w]=KorobovLatticeRules(numDim,method,params,[],lowerBounds,upperBounds);
% est3=sum(bsxfun(@times,f(xi),w'));
% error3=abs(exactSol-est3)
% %Also we consider it with the type II points and periodization.
% [xi,w]=KorobovLatticeRules(numDim,method,params,periodizeAlg,lowerBounds,upperBounds);
% est4=sum(bsxfun(@times,f(xi),w'));
% error4=abs(exactSol-est4)
% %Finally, we consider using a larger number of type II points and
% %periodization.
% params=4;%51097 points
% [xi,w]=KorobovLatticeRules(numDim,method,params,periodizeAlg,lowerBounds,upperBounds);
% est5=sum(bsxfun(@times,f(xi),w'));
% error5=abs(exactSol-est5)
%From this example, one can see that the errors again improve with
%periodization, but that the first choice of type II points performs worse
%than the type I points, despite using many more points. However, once
%periodization is used with the type II points, the results are
%significantly better.
%
%EXAMPLE 4:
%In this instance, we look at how differing search spaces improve Komorov's
%second algorithm versus the tabulated values.
% numDim=5;
% method=2;
% params=1;
% %First we get the tabulated values.
% [xi1,w1,Hm1,algParams1]=KorobovLatticeRules(numDim,method,params);
% %Now, we try to optimize over the parameters, rather than using the
% %tabulated version.
% params=[];
% params.p=algParams1.p;
% params.q=algParams1.q;
% method=3;
% %Note that without setting a bound on aMin and aMax, the following line
% %is VERY slow due to the large search space if one searches for all a
% %values from q to p. Thus, here the search space is restricted.
% params.aMin=1;
% params.aMax=100;
% [xi2,w2,Hm2,algParams2]=KorobovLatticeRules(numDim,method,params);
% %One will find that Hm1=1.003235184762100 but
% %Hm2=1.002308498901230 and algParams2 contains algParams2.method=2,
% %algParams2.numDim=5, algParams2.p=653, algParams2.q=23, algParams2.a=89
% %algParams2.b=11; The optimization method provided a better
% %result than the tabulated value. The optimized result can be used on the
% %opimization problem from example 3 as 
% f=@(x)sum(abs(x-1/2),1);
% lowerBounds=[-2;-1;0;4;3];
% upperBounds=[2;1;1;6;8];
% exactSol=915;
% numDim=5;
% method=4;
% params=algParams2;
% [xi,w]=KorobovLatticeRules(numDim,method,params,[],lowerBounds,upperBounds);
% est1=sum(bsxfun(@times,f(xi),w'));
% error1=abs(exactSol-est1)
% %One will see that the performance is better than in the previous example,
% but only by a very small amount.
%
%EXAMPLE 5:
%We consider integration in 2D.
% f=@(x)exp(prod(x,1));
% exactSol=1.3179021514544038949;
% numDim=2;
% method=5;
% params=11;%The 11th Fibonacci number is 89.
% [xi,w]=KorobovLatticeRules(numDim,method,params);
% est1=sum(bsxfun(@times,f(xi),w'));
% error1=abs(exactSol-est1)
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%[2] I. H. Sloan and S. Joe, Lattice Methods for Multiple integration.
%    Oxford, United Kingdom: Clarendon Press, 1994.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

Hm=[];
algParams=[];
switch(method)
    case 0
        num=params;
        [xi,a]=getKorobovTypeIPoints(numDim,num);
        
        %Compute Hm given the optimal a value, if desired as a return
        %value.
        if(nargout>2)
            p=size(xi,2);
            Hm=HmTypeI(numDim,p,m);
            
            algParams.method=0;
            algParams.numDim=numDim;
            algParams.p=num;
            algParams.a=a;
        end
    case 1
        primeVal=params;
        [xi,Hm,a]=findKorobovTypeIPointsViaOpt(numDim,primeVal);
        
        algParams.method=0;
        algParams.numDim=numDim;
        algParams.p=primeVal;
        algParams.a=a;
    case 2
        num=params;
        [xi,p,q,a,b]=getKorobovTypeIIPoints(numDim,num);
        
        %Compute Hm given the optimal values of p,q,a,and b, if desired as
        %a return value. This is Equation 6.3-6 in [1].
        if(nargout>2)
            Hm=HmTypeII(numDim,p,q,a,b);
        end
        
        algParams.method=2;
        algParams.numDim=numDim;
        algParams.p=p;
        algParams.q=q;
        algParams.a=a;
        algParams.b=b;
    case 3
        p=params.p;
        q=params.q;
        if(isfield(params,'aMin'))
            aMin=params.aMin;
            aMax=params.aMax;
        else
            aMin=min(p,q);
            aMax=max(p,q);
        end

        [xi,Hm,p,q,a,b]=findKorobovTypeIIPointsViaOpt(numDim,p,q,aMin,aMax);
        
        algParams.method=2;
        algParams.numDim=numDim;
        algParams.p=p;
        algParams.q=q;
        algParams.a=a;
        algParams.b=b;
    case 4
        algParams=params;
        switch(algParams.method)
            case 0%Koborov type I
                p=algParams.p;
                a=algParams.a;
                xi=points4TypeIParams(numDim,p,a);
                if(nargout>1)
                   Hm=HmTypeI(numDim,p,m);
                end
            case 2%Koborov type II
                p=algParams.p;
                q=algParams.q;
                a=algParams.a;
                b=algParams.b;
                xi=points4TypeIIParams(numDim,p,q,a,b);
                
                if(nargout>1)
                   Hm=HmTypeII(numDim,p,q,a,b); 
                end
            otherwise
                error('Invalid method in the algorithmic parameters passed')
        end
    case 5
        if(numDim~=2)
            error('A Fibonacci lattice can only be used when numDim=2')
        end
        num=params;
        xi=getFibonacciPoints(num);
    otherwise
        error('Unknown method specified')
end

%If the points should be periodized.
if(nargin>3&&~isempty(periodizeAlg))
    method=periodizeAlg.method;
    if(isfield(periodizeAlg,'param1'))
        param1=periodizeAlg.param1;
    else
        param1=[];
    end
    if(isfield(periodizeAlg,'param2'))
        param1=periodizeAlg.param2;
    else
        param2=[];
    end
    
    [xi,w]=periodizeLatticePoints(method,xi,[],param1,param2);
else%Otherwise, the xi do not change and the weights are just uniform.
    numPoints=size(xi,2);
    w=(1/numPoints)*ones(numPoints,1);
end

%If lower and upper bounds are provided.
if(nargin>4&&~isempty(lowerBounds))
    for curDim=1:numDim
        a=lowerBounds(curDim);
        b=upperBounds(curDim);
        
        if(a<b)
            flipSign=1;
        else%Make a<b and flip the sign of the result.
            temp=a;
            a=b;
            b=temp;
            flipSign=-1;
        end
        %Now, a<b.

        if(isfinite(a)&&isfinite(b))%Both bounds finite
            xi(curDim,:)=a+(b-a)*xi(curDim,:);
            w=flipSign*w*(b-a);
        elseif(~isfinite(a)&&isfinite(b))%a=-Inf
            xi(curDim,:)=-1./xi(curDim,:)+b+1;
            w=flipSign*w./xi(curDim,:).^2;
        elseif(isfinite(a)&&~isfinite(b))%b=Inf
            xi(curDim,:)=1./(1-xi(curDim,:))+a-1;
            w=flipSign*w./(1-xi(curDim,:)).^2;
        else%a=-Inf and b=Inf
            xi(curDim,:)=1./(1-xi(curDim,:))-1./xi(curDim,:);
            w=flipSign*w.*(1./(1-xi(curDim,:)).^2+1./xi(curDim,:).^2);
        end
    end
end
end

function [xi,Hm,a]=findKorobovTypeIPointsViaOpt(numDim,primeVal)
%%FINDKOROBOXTYPEIPOINTSVIAOPT Find the optimal set of of lattice points
%               based on Theorem 6.3-5 in [1], which is the first one
%               listed from Korobov. The points are for numerical
%               integration of a continuous multidimensional function that
%               is periodic in all dimensions with the first period being
%               the range 0-1, in all dimensions. The points returned are
%               wrapped to the range 0-1.            
%
%INPUTS: numDim The number of dimensions of the lattice points that are to
%               be generated. numDim>=2.
%      primeVal The number of points to generate. This must be a prime
%               number>numDim. Prime numbers can be found in Matlab using
%               the primes function.
%
%OUTPUTS: xi A numDimXprimeVal matrix of the lattice points generated using
%            Korobov's type I rule.
%         Hm The value of the cost function from Equation 6.3-6 in [1]. The
%            closer this is to 1, the better the set of points is.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isprime(primeVal)==false)
    error('To generate type I Korobov points, the number of points must be a prime number')
end

if(primeVal<numDim)
    error('To generate type I Korobov points, the number of points must be greater than the number of dimensions.')
end

p=primeVal;

%Do a brute force search over all m to find the minimum value of the cost
%function in Equation 6.3-6 in [1]. The equation has been corrected as the
%sum term should have been to p, not to m. This is Theorem 6.3-5.
minCost=Inf;
for m=1:(p-1)
    Hm=HmTypeI(numDim,p,m);
    if(Hm<minCost)
        mMin=m;
        minCost=Hm;
    end
end
Hm=minCost;
a=mMin;

xi=points4TypeIParams(numDim,p,a);
end

function Hm=HmTypeI(numDim,p,m)
%%HMTYPEI This implements the cost function in Equation 6.3-6 in [1] when
%         given all parameters.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    Hm=0;
    mPows=m.^(0:(numDim-1));
    for k=1:p
        Hm=Hm+prod((1-2*mod(k*mPows/p,1)).^2);
    end
    
    Hm=(3^numDim/p)*Hm;
end

function Hm=HmTypeII(numDim,p,q,a,m)
%%HMTYPEI This implements the cost function in Equation 6.3-7 in [1] when
%         given all parameters.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    Hm=0;
    numPoints=p*q;
    mPows=m.^(0:(numDim-1));
    aPows=a.^(0:(numDim-1));
    for k=1:numPoints
        Hm=Hm+prod((1-2*mod(k*(p*mPows+q*aPows)./numPoints,1)).^2);
    end
    Hm=3^numDim/numPoints*Hm;
end


function [xi,Hm,p,q,a,b]=findKorobovTypeIIPointsViaOpt(numDim,p,q,aMin,aMax)
%%FINDKOROBOVTYPEIIPOINTSVIAOPT Find the optimal set of of lattice points
%               based on Theorem 6.3-7 in [1], which is the second one
%               listed from Korobov. The points are for numerical
%               integration of a continuous multidimensional function that
%               is periodic in all dimensions with the first period being
%               the range 0-1, in all dimensions. The points returned are
%               wrapped to the range 0-1.  
%
%INPUTS:numDim The number of dimensions of the lattice points that are to
%              be generated. numDim>=2.
%          p,q Two prime numbers>numDim. p*q will be the number of lattice
%              points generated. Usually p>q, as this execution time of
%              this function is proportional to q.
%    aMin,aMax The value a is an integer>0, often larger than q. Varying a
%              changes the location of the optimal point (can improve Hm)
%              for a fixed p and q. aMin and aMax define the range of a to
%              search for an optimal value.
%
%OUTPUTS: xi A numDimX(p*q) matrix of the lattice points generated using
%            Korobov's type II rule.
%         Hm The value of the cost function from Equation 6.3-7 in [1]. The
%            closer this is to 1, the better the set of points is.
%
%Finding good combinations of p, q, and a affect the quality of the points.
%In most tabulated sets of points, one observes that p>a>q.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(~isprime(p)||~isprime(q))
   error('One or more of the two required prime numbers to search for a type II Korobov is not prime')
end

if(p<numDim||q<numDim)
    error('To generate type II Korobov points, the two prime numbers must be greater than the number of dimensions.')
end

%Do a brute force search over all m to find the minimum value of the cost
%function in Equation 6.3-7 in [1]. This is Theorem 6.3-6.
minCost=Inf;
for a=aMin:aMax
    for m=1:(q-1)
        Hm=HmTypeII(numDim,p,q,a,m);
        
        if(Hm<minCost)
           minCost=Hm;
           mMin=m;
           aMinVal=a;
        end
    end
end
Hm=minCost;
b=mMin;
a=aMinVal;

xi=points4TypeIIParams(numDim,p,q,a,b);
end

function [xi,a]=getKorobovTypeIPoints(numDim,num)
%%GETKOROBOVTYPEIPOINTS Obtain the lattice points based on Theorem 6.3-5 in
%               [1] for numerical integration of a continuous
%               multidimensional function that is periodic in all
%               dimensions with the first period being the range 0-1, in
%               all dimensions. The points returned are wrapped to the
%               range 0-1.
%
%INPUTS: numDim The number of dimensions of the points to generate.
%               Tabulated values are 3<=numDum<=6.
%           num For a given number of dimensions, multiple formulae are
%               tabulated. num selects which one. Higher values of num
%               produce larger sets of lattice points. In all instances,
%               the maximum value of num yields a set of 10007 lattice
%               points. num starts at 1. The maximum value of num for
%               different values of numDim are
%               numDim  the maximum num
%               3       28
%               4       24
%               5       18
%               6       14
%
%OUTPUTS: xi A numDimXnumPoints set of points for numerically integrating a
%            periodic function between 0-1 in all dimensions. Associated
%            weights would just be 1/(numPoints) (uniform).
%
%The lattice points (which all have uniform weighting -- so 1/numPoints),
%can be used for integrating continuous functions that are periodic over
%0-1 in all dimensions. For non-periodic functions, a periodization of the
%points prior to numeric integration should be performed.
%
%Solutions to the formula of Theorem 6.3-5 for specific numbers of points
%are taken from Tables 6.2 through 6.5 of [1], providing solutions for from
%three to six dimensions.
%
%Only p (the number of points) and the first non-unitary dimension of the
%generating vector theta need be tabulated (it is theta(2)), because the
%other components of theta are mod(theta(2)^idx,p). The modulo is due to
%the periodic nature of the space.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

switch(numDim)
    case 3
        %First column is p, the second columns is theta(2) from Table 6.2
        %in Chapter 6.3 of [1], pg. 202.
        optVecTable=[101,   40;
                     101,   48;
                     199    30;
                     199,   73;
                     307,   75;
                     307,   131;
                     523,   78;
                     523,   114;
                     701,   215;
                     701,   313;
                     1069,  136;
                     1069,  338;
                     1543,  355;
                     1543,  552;
                     2129,  359;
                     2129,  937;
                     3001,  276;
                     3001,  772;
                     4001,  722;
                     4001,  1934;
                     5003,  1476;
                     5003,  1949;
                     6007,  592;
                     6007,  2831;
                     8191,  739;
                     8191,  3303;
                     10007, 544;
                     10007, 3072];
    case 4
        %First column is p, the second columns is theta(2) from Table 6.3
        %in Chapter 6.3 of [1], pg. 203.
        optVecTable=[307,   42;
                     307,   95;
                     523,   178;
                     523,   238;
                     701,   82;
                     701,   265;
                     1069,  71;
                     1069,  271;
                     1543,  128;
                     1543,  663;
                     2129,  766;
                     2129,  970;
                     3001,  174;
                     3001,  1466;
                     4001,  113;
                     4001,  956;
                     5003,  792;
                     5003,  2053;
                     6007,  1351;
                     6007,  2610;
                     8191,  2488;
                     8191,  3842;
                     10007, 1206;
                     10007, 1784];
    case 5
        %First column is p, the second columns is theta(2) from Table 6.4
        %in Chapter 6.3 of [1], pg. 203.
        optVecTable=[1069,  63;
                     1069,  526;
                     1543,  58;
                     1543,  133;
                     2129,  618;
                     2129,  720;
                     3001,  408;
                     3001,  890;
                     4001,  1534;
                     4001,  1651;
                     5003,  840;
                     5003,  1352;
                     6007,  509;
                     6007,  1487;
                     8191,  1386;
                     8191,  2228;
                     10007, 198;
                     10007, 1870];
    case 6
        %First column is p, the second columns is theta(2) from Table 6.5
        %in Chapter 6.3 of [1], pg. 204.
        optVecTable=[2129,  41;
                     2129,  727;
                     3001,  233;
                     3001,  322;
                     4001,  1751;
                     4001,  1780;
                     5003,  2037;
                     5003,  2208;
                     6007,  312;
                     6007,  1521;
                     8191,  1632;
                     8191,  3699;
                     10007, 2240;
                     10007, 2399];
    otherwise
        error('Untabulated number of dimensions chosen')
end

p=optVecTable(num,1);
a=optVecTable(num,2);

xi=points4TypeIParams(numDim,p,a);
end

function [xi,p,q,a,b]=getKorobovTypeIIPoints(numDim,num)
%%GETKOROBOVTYPEIIPOINTS Obtain the lattice points based on Theorem 6.3-5
%               in [1] for numerical integration of a continuous
%               multidimensional function that is periodic in all
%               dimensions with the first period being the range 0-1, in
%               all dimensions. The points returned are wrapped to the
%               range 0-1.
%
%INPUTS: numDim The number of dimensions of the points to generate.
%               Tabulated values are 3<=numDum<=10.
%           num For a given number of dimensions, multiple formulae are
%               tabulated. num selects which one. Higher values of num
%               produce larger sets of lattice points. In all instances,
%               the maximum value of num yields a set of 10007 lattice
%               points. num starts at 1. The maximum value of num for
%               different values of numDim  and the number of points
%               generated for the maximum value are
%               numDim   the maximum num   max number of points
%               3        6                 100063
%               4        6                 100063
%               5        6                 100063
%               6        6                 100063
%               7        8                 100063
%               8        6                 100063
%               9        6                 159053
%               10       6                 155093
%
%OUTPUTS: xi A numDimXnumPoints set of points for numerically integrating a
%            periodic function between 0-1 in all dimensions. Associated
%            weights would just be 1/(numPoints) (uniform).
%
%The lattice points (which all have uniform weighting -- so 1/numPoints),
%can be used for integrating continuous functions that are periodic over
%0-1 in all dimensions. For non-periodic functions, a periodization of the
%points prior to numeric integration should be performed.
%
%Solutions to the formula of Theorem 6.3-6 for specific numbers of points
%are taken from Tables 6.6 through 6.13 of [1], providing solutions for 
%from three to ten dimensions.
%
%p, q, a and b are tabulated instead of the actual points, which can be
%recovered from p, q, a, and b, using Theorem 6.3-6. The thetas obtained
%from the tabulated values is not identical to that in the tables, because
%no effort is made to multiply the theta by a large number so that the
%first component (modulo p*q) equals 1. the results are, however, equally
%valid.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

switch(numDim)
    case 3
        %The columns are p,q,a,b from Table 6.6 in Chapter 6.3 of [1], pg.
        %204.
        optVecTable=[691,   29, 176,    20;
                     907,   31, 402,    12;
                     1259,  31, 535,    5;
                     1543,  37, 355,    14;
                     1907,  43, 275,    10;
                     2129,  47, 359,    4];
    case 4
        %The columns are p,q,a,b from Table 6.7 in Chapter 6.3 of [1], pg.
        %204.
        optVecTable=[691,   29, 320,    6;
                     907,   31, 316,    3;
                     1259,  31, 483,    9;
                     1543,  37, 128,    13;
                     1907,  43, 60,     37;
                     2129,  47, 766,    5];   
    case 5
        %The columns are p,q,a,b from Table 6.8 in Chapter 6.3 of [1], pg.
        %205.
        optVecTable=[653,   23, 193,    15;
                     691,   29, 271,    17;
                     1069,  31, 63,     17;
                     1381,  37, 480,    13;
                     1733,  41, 828,    12;
                     2129,  47, 618,    31];
    case 6
        %The columns are p,q,a,b from Table 6.9 in Chapter 6.3 of [1], pg.
        %205.
        optVecTable=[653,   23, 254,    3;
                     691,   29, 29,     18;
                     1069,  31, 63,     8;
                     1381,  37, 264,    15;
                     1733,  41, 680,    11;
                     2129,  47, 727,    20];
    case 7
        %The columns are p,q,a,b from Table 6.10 in Chapter 6.3 of [1], pg.
        %206.
        optVecTable=[653,   23, 32,     19;
                     787,   23, 173,    7;
                     829,   29, 175,    6;
                     1069,  31, 159,    16;
                     1249,  37, 430,    12;
                     1543,  37, 82,     14;
                     1733,  41, 680,    17;
                     2129,  47, 718,    30];
    case 8
        %The columns are p,q,a,b from Table 6.11 in Chapter 6.3 of [1], pg.
        %206.
        optVecTable=[829,   29, 32,     12;
                     1069,  31, 313,    17;
                     1249,  37, 351,    19;
                     1543,  37, 438,    21;
                     1733,  41, 104,    38;
                     2129,  47, 86,     20];
    case 9
        %The columns are p,q,a,b from Table 6.12 in Chapter 6.3 of [1], pg.
        %207.
        optVecTable=[1069,  31, 68,     6;
                     1249,  37, 128,    28;
                     1543,  37, 117,    11;
                     1733,  41, 459,    9;
                     2129,  47, 636,    17;
                     3001,  53, 108,    26];
    case 10
        %The columns are p,q,a,b from Table 6.13 in Chapter 6.3 of [1], pg.
        %207.
        optVecTable=[4507,  19, 1611,   9;
                     4507,  23, 1611,   3;
                     5003,  23, 431,    12;
                     4507,  29, 1611,   10;
                     5003,  29, 431,    16;
                     5003,  31, 431,    27];
    otherwise
        error('Untabulated number of dimensions chosen')
end

p=optVecTable(num,1);
q=optVecTable(num,2);
a=optVecTable(num,3);
b=optVecTable(num,4);

xi=points4TypeIIParams(numDim,p,q,a,b);
end

function xi=points4TypeIParams(numDim,p,a)
%%POINTS4TYPEIPARAMS Generate the lattice points based on Theorem 6.3-4
%               in [1] for numerical integration of a continuous
%               multidimensional function that is periodic in all
%               dimensions with the first period being the range 0-1, in
%               all dimensions. The points returned are wrapped to the
%               range 0-1. This function takes the explicit parameters of
%               Equation 6.3-6.
%
%INPUTS: numDim The number of dimensions of the lattide points to generate.
%       p,a The parameters of Equation 6.3-5 of [1]. ps is a prime number
%           and is equal to the number of points produced. a is a positive
%           integer that should minimize the cost function.
%
%OUTPUTS: xi A numDimXnumPoints set of points for numerically integrating a
%            periodic function between 0-1 in all dimensions. Associated
%            weights would just be 1/(numPoints) (uniform).
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Recovering theta from Theorem 6.3-5
theta=[1,a.^(1:(numDim-1))]';

%Creating the points from theta using the pattern of Theorem 6.3-4.
xi=zeros(numDim,p);
for k=0:(p-1)
    xi(:,k+1)= theta*(k/p);
end
%Wrapping the points to the range 0-1. This is valid due to the assumed
%periodicity of the function to be integrated.
xi=wrapRange(xi,0,1);
end

function xi=points4TypeIIParams(numDim,p,q,a,b)
%%POINTS4TYPEIIPARAMS Generate the lattice points based on Theorem 6.3-6
%               in [1] for numerical integration of a continuous
%               multidimensional function that is periodic in all
%               dimensions with the first period being the range 0-1, in
%               all dimensions. The points returned are wrapped to the
%               range 0-1. This function takes the explicit parameters of
%               Equation 6.3-7.
%
%INPUTS: numDim The number of dimensions of the lattide points to generate.
%       p,q,a,b The parameters of Equation 6.3-7 of [1]. p and q are prime
%               numbers; p*q is the number of points to generate. a and b
%               are positive integers that should minimize the cost
%               function.
%
%OUTPUTS: xi A numDimXnumPoints set of points for numerically integrating a
%            periodic function between 0-1 in all dimensions. Associated
%            weights would just be 1/(numPoints) (uniform).
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Recreating theta as in Theorem 6.3-6.
i=(1:numDim)';
theta=p*b.^(i-1)+q.*a.^(i-1);

numPoints=p*q;

%Creating the points from theta using the pattern of Theorem 6.3-4.
xi=zeros(numDim,numPoints);
for k=0:(numPoints-1)
    xi(:,k+1)= theta*(k/numPoints);
end
%Wrapping the points to the range 0-1. This is valid due to the assumed
%periodicity of the function to be integrated.
xi=wrapRange(xi,0,1);
end


function xi=getFibonacciPoints(FibPos)
%GETFIBONACCIPOINTS  Obtain the lattice points based on Fibonacci numbers
%               for numerical integration of a continuous two-dimensional
%               function that is periodic in both dimensions with the first
%               period being the range 0-1, in both dimensions. The points
%               returned are wrapped to the range 0-1.
%
%INPUTS: FibPos The number of points generated depends on the position of
%               the maximum number in the Fibonacci series used. This is
%               the inex of that maximum number. FibPos>1. The number of
%               points is equal to the Fibonacci number at the chosen
%               position.
%
%OUTPUTS: xi A numDimXnumPoints set of points for numerically integrating a
%            periodic function between 0-1 in all dimensions. Associated
%            weights would just be 1/(numPoints) (uniform).
%
%The method of generating Fibonacci points is described in Chapter 6.2 of
%[1] and Chapter 4.3 of [2]. As mentioned in Chapter 6.3 of [1], the
%Fibonacci points exhibit certain optimality properties that are difficult
%to obtain in higher numbers of dimensions.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%[2] I. H. Sloan and S. Joe, Lattice Methods for Multiple Integration.
%    Oxford, United Kingdom: Clarendon Press, 1994.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(FibPos<2)
   error('The position in the Fibonacci series used for generating 2D lattice points must be >=1') 
end

Fm=FibonacciNum(FibPos);
Fm1=FibonacciNum(FibPos-1);

numPoints=Fm;
xi=zeros(2,numPoints);
for j=0:(Fm-1)
   xi(:,j+1)=(j/Fm)*[1;Fm1];
end

%Wrap the points back into the range of 0->1.
xi=wrapRange(xi,0,1);
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
