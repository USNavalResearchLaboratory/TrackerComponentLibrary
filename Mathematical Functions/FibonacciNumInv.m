function k=FibonacciNumInv(n,choice)
%FIBONACCINUMINV Given a real number >=1, determine the k such that
%                n==FibonacciNum(k). If n given is not a Fibonacci number,
%                then the behaviour is determined by the choice option.
%
%INPUTS: n A scalar or matrix of values for which the argument of
%          FibonacciNum is desired, n>=1.
%   choice An optional parameter that determines what is returned if n
%          is not an exact Fibonacci number. If this parameter is omitted,
%          then the default is one (the next lower value).
%          0 means return the closest value.
%          1 means return the next lower value.
%          2 means return the next higher value.    
%
%OUTPUTS: k The argument of FibonacciNum that gives n, or an appropriately
%           close value selected by the choice option. k will always be
%           >=2.
%
%A non-recursive expression for Fibonacci numbers is given in [1] and can
%be inverted. However, due to finite prevision issues, for n in the mid
%70's the value returned can be off. Thus, testing needs to be done for all
%inputs in the choice option to make sure that the proper value is chosen.
%
%REFERENCES:
%[1] Chandra, Pravin and Weisstein, Eric W. "Fibonacci Number." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/FibonacciNumber.html
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<2||isempty(choice))
        choice=1; 
    end

    %The golden ratio.
    phi=(1+sqrt(5))/2;
    k=fix(log(n*sqrt(5))/log(phi)+(1/2));
    
    %Due to finite precision errors, even if n is a Fibonacci number, k
    %might be the correct value, it might be 1 less or it might be 1 too
    %much. We will adjust k to deal with the situation. Note that from
    %above, k>=2, so there are no issues with k-1 being <1.
    switch(choice)
        case 0%Return the closest value
            valsLower=FibonacciNum(k-1);
            valsCur=FibonacciNum(k);
            valsUpper=FibonacciNum(k+1);

            diffLower=abs(n-valsLower);
            diffCur=abs(n-valsCur);
            diffUpper=abs(n-valsUpper);

            sel=(diffLower<diffCur)&(diffLower<diffUpper);
            k(sel)=k(sel)-1;
            sel=(valsUpper<diffCur)&(diffUpper<diffLower);
            k(sel)=k(sel)+1;
        case 1%Return the next lower value.
            %We have to do a check, because finite precision error
            %above might make k be too large.
            valsCur=FibonacciNum(k);

            sel=valsCur>n;
            k(sel)=k(sel)-1;
        case 2%Return the next higher value.
            valsCur=FibonacciNum(k);
            sel=valsCur<n;
            k(sel)=k(sel)+1;
        otherwise
            error('Unknown choice value specified.')
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
