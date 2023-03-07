classdef Interval < matlab.mixin.CustomDisplay
    %%INTERVAL A class to implement basic interval arithmetic operations.
    %          The class is partially conformant to the IEEE Std 1788-2015
    %          for Interval arithmetic with respect to set-based intervals.
    %          The most notable standard violation is the lack of
    %          decorations. Rather than defining a single point value,
    %          Intervals define a range, or set, of possible values.  All
    %          computations rely on machine infrastructure to provide
    %          rounding accuracy at floating-point level precision.
    %          Intervals play a role in the box particle filter. One can
    %          create an manipulate arrays and matrices of Intervals in a
    %          manner similar to standard datatypes in Matlab.
    %
    %The IEEE Standard is [1]. This class does not implement different
    %decorations or flavors. NaNs are used to indicate empty intervals.
    %The use of complex numbers in intervals is not supported.
    %
    %As an example, vectors of two intervals.
    %   a=[Interval(1,2);Interval(0,0)];
    %   b=Interval(3,4);Interval(12,24)];
    %   a+b
    %   a.*b
    %   a*b
    %
    %REFERENCES:
    %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
    %    pp.1-97, June 30 2015.
    %
    %October 2015 David F. Crouse and David Karnick, Naval Research
    %Laboratory, Washington D.C.
    %(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
    properties
        lb%Lower bounds of the interval array
        ub%Upper bounds of the interval array
    end
    
    methods
        function newInt=Interval(lb,ub)
        %%INTERVAL Create a new instance of the Interval class.
        %
        %INPUTS: Three types of initializations can be performed. The first
        %        type is with no inputs, in which case the interval set if
        %        just empty.
        %        The second type of initialization has one input:
        %        val  A scalar value. The interval will be a single point.
        %             If one wishes to create an interval having the finite
        %             precision limitations of a constant, then this input
        %             method should not be used. Rather, one should use
        %             lb=val-eps(val) and ub=val+eps(val) in a two-
        %             parameter initialization.
        %        The third type of initialization  takes two inputs:
        %        lb, ub numDumX1 vectors of the lower and upper bounds of
        %               the interval in each of the dimensions.
        %
        %If something in ub or lb is a NaN, then the corresponding thing in
        %the other bound will be set to a NaN. This avoids the creation of
        %invalid Intervals.
        %
            if(nargin<1)%Create empty interval set.
                newInt.lb=[];
                newInt.ub=[];
            elseif(nargin<2)%Only a scalar is given.
                newInt.lb=lb;
                newInt.ub=lb;
            elseif(any(size(lb)~=size(ub)))
                error('Lower bound and Upper bound must be same size.')
            else%Interval bounds are given.
                %If there is a NaN in one, put it in the other.
                sel=isnan(lb) | isnan(ub);
                lb(sel)=NaN;
                ub(sel)=NaN;

                newInt.lb=lb;
                newInt.ub=ub;

                assert(all(isnan(newInt.lb(:))|newInt.lb(:)<=newInt.ub(:)));
            end
        end
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Numeric functions of intervals of Chapter 9.2 of the     %%%
%%% IEEE Std 1788-2015                                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function c=cancelMinus(a,b)
        %%CANCELMINUS Given that one knows that a=b+c (all are intervals),
        %             recover c from a and b. This is only defined if the
        %             width of a is not less than b, otherwise an empty
        %             Interval is returned. An empty interval is also
        %             returned if either a or b contains +/-Inf, because
        %             the original value cannot be recovered.
        %
        %This function is defined in Section 9.3 of [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            %Allocate space
            c=Interval.empty(size(a.lb));
            
            %Presumably, the rounding mode should be opposite the rounding
            %mode one would use with addition so as the undo any rounding
            %that may have occurred.
            setProcRoundingMode(3);%Round to +Inf
            c.lb=a.lb-b.lb;
            setProcRoundingMode(0);%Round to -Inf
            c.ub=a.ub-b.ub;

            %Next, deal with b being empty and a not empty or any of the
            %endpoints of the intervals being infinite. There is also the
            %case where the width of a is less than the width of b. In all
            %cases, we will return an interval that is the entire line.
            selNaN=a.lb==-Inf|a.ub==Inf|b.lb==-Inf|b.ub==Inf;
            %For the test of the width of a being less than the width of b,
            %we order things to minimize finite precision errors. The best
            %rounding mode to use is probably nearest, though a strict
            %proof of that would be difficult. We also throw in the issue
            %of rounding modes having made the upper bound larger than the
            %lower bound.
            setProcRoundingMode(2);%Round to nearest
            selNaN=selNaN|(isnan(a)|isnan(b))|((a.ub-b.ub)-(a.lb-b.lb)<0)|(c.lb>c.ub);
            
            s.type='()';
            s.subs{1}=find(selNaN);
            c=c.subsasgn(s,Interval(-NaN,NaN));
        end

        function c=cancelPlus(a,b)
        %%CANCELPLUS Given that one knows that a was obtained from c-b,
        %            recover the interval c from a and b.
        %
        %This function is defined in Section 9.3 of [1]as cjust
        %cancelMinus(a,-b).
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.

            c=cancelMinus(a,-b);
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Operations of Chapter 9.3 of the IEEE Std 1788-2015  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function c=intersection(a,b)
        %%INTERSECTION Pairwise intersection of two intervals. This is the
        %           region where the intervals overlap. If theere is no
        %           intersection, then an empty interval (with NaNs) is
        %           returned.
        %
        %This function is defined in Section 9.3 of [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            if(~isa(a,'Interval'))
                a=Interval(a);
            elseif(~isa(b,'Interval'))
                b=Interval(b);
            end
            
            lbNew=max(a.lb,b.lb);
            ubNew=min(a.ub,b.ub);
            noIntersect=(lbNew>ubNew | isnan(a) | isnan(b));
            lbNew(noIntersect)=NaN;
            ubNew(noIntersect)=NaN;
            c=Interval(lbNew,ubNew);
        end

        function c=convexHull(a,b)
        %%CONVEXHULL Find the pairwise convex hull of two intervals. This
        %            is the smallest interval that can contain both
        %            intervals.
        %
        %This function is defined in Section 9.3 of [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            if(~isa(a,'Interval'))
                a=Interval(a);
            elseif(~isa(b,'Interval'))
                b=Interval(b);
            end
            
            %Empty intervals (NaNs) do not contribute to the hull unless
            %both are NaNs.
            lbNew=min(a.lb,b.lb);
            ubNew=max(a.ub,b.ub);
            c=Interval(lbNew,ubNew);
        end
        
        function c=hull(a,b)
        %%HULL This is another name for the method convexHull. See comments
        %      to that method for more information.

            c=convexHull(a,b);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Numeric functions of intervals of Chapters 9.4 and       %%%
%%% and 10.5.9 of the IEEE Std 1788-2015                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function c=inf(a)
        %%INF Return the lower bound (infimum) of an Interval.
        %
        %This function is defined in Chapers 9.4 and 10.5.9 of [1].
        %       
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            c=a.lb;
        end
        
        function c=sup(a)
        %%SUP Return the upper bound (supremum) of an Interval.
        %    
        %This function is defined in Chapters  9.4 and 10.5.9  of [1].
        %       
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            c=a.ub;
        end
        
        function c=mid(a)
        %%MID Find the midpoint of an Interval. Empty intervals return a
        %     NaN.
        %
        %This function is defined in Chapters 9.4 and and 10.5.9 of [1].
        %The method used is that of [2] as making a midpoint formula
        %compliant with [1] is not trivial.
        %       
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] F. Goulard, "How do you compute the midpoint of an interval?"
        %    ACM Transaction on Mathematical Software, vol. 40, no. 2, Feb.
        %    2014.
            
            %Allocate space
            c=size(a.lb);
        
            sel0=isnan(a.lb|a.ub);
            sel1=a.lb==-Inf & b.ub==Inf;
            sel2=a.lb==-Inf &~(b.ub==Inf);
            sel3=~(a.lb==-Inf)&(b.ub==Inf);
            sel4=~(sel1|sel2|sel3);
            
            c(sel0)=NaN;
            c(sel1)=0;
            c(sel2)=-Inf;
            c(sel3)=Inf;

            setProcRoundingMode(2);%Round to nearest
            c(sel4)=0.5*(a.lb+a.ub);
            sel5=~isfinite(c)&sel4;
            c(sel5)=a.lb/2+a.ub/2;
        end
        
        function c=wid(a)
        %%WID Find the width of the interva (upper bound minus lower
        %     bound).
        %
        %This function is defined in Chapters 9.4 and and 10.5.9 of [1]. As
        %per Chapter 12.12.8, the result is rounded towards +Inf.
        %       
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
            
            if(~isa(a,'Interval'))
                c=0;
            else
                setProcRoundingMode(3);%Round to +Inf
                c=a.ub-a.lb;
            end
        end

        function c=rad(a)
        %%RAD Find the radius of the interval, element by element. This is
        %     half the difference of the upper and lower bounds.
        %
        %This function is defined in Chapters 9.4 and and 10.5.9 of [1]. As
        %per Chapter 12.12.8, the midpoint must cover both boundaries, so
        %it is assumed that rad is rounded towards +Inf.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
            
            midPoint=a.mid();

            %The midpoint must cover both boundaries
            setProcRoundingMode(3);%Round to +Inf
            c=max(a.ub-midPoint,midPoint-a.lb);
        end
        
        function c=mag(a)
        %%MAG Find the magnitude of the Interval, element-by-element. That
        %     is the maximum of the absolute value of the Interval.
        %
        %This function is defined in Chapters 9.4 and and 10.5.9 of [1].
        %       
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            c=max(abs(a.lb),abs(a.ub));
        end
        
        function c=mig(a)
        %%MIG Find the migitude of the interval. That is the minimum of
        %    the absolute value of the Interval.
        %
        %This function is defined in Chapters 9.4 and and 10.5.9 of [1].
        %       
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            c=min(abs(a.lb),abs(a.ub));
            sel=a.lb<=0&&a.ub>=0;
            c(sel)=0;
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Numeric functions of intervals of Table 9.3 in Chapter   %%%
%%% 9.5 of the IEEE Std 1788-2015                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function c=equal(a,b)
        %%EQUAL Determine equality between sets of intervals (or intervals
        %       and a point).
        %
        %This function is defined in Table 9.3 of Chapter 9.5 of [1] and is
        %the same as the equals function of Table 10.7 of [1] except the
        %results of all independent intervals are and-ed together.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
            
            c=equals(a,b);
            c=all(c(:));
        end

        function c=subset(a,b)
        %%SUBSET Determine whether a is a subset of b, taking into
        %        consideration all intervals in a set. This includes the
        %        endpoints. Compare to the contains function, which omits
        %        the endpoints and  does not and-the results of the
        %        independent comparisons together.
        %
        %This function is defined in Table 9.3 of Chapter 9.5 of [1]. In
        %Chapter 10.5.10, it is specified that an empty interval is both a
        %subset of any interval and disjoint from them. Thus, a check for
        %an empty interval makes that true here.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
            
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end

            c=(a.ub<=b.ub)&(a.lb>=b.lb) | (isnan(a.lb) |isnan(b.lb));
            c=all(c(:));
        end

        function c=interior(a,b)
        %%INTERIOR Determine whether all values in a are within b (and
        %          not on the boundary). This is the same as the
        %          containedBy function except this function considers all
        %          intervals, and-ing the results together.
        %
        %This function is defined in Table 9.3 of Chapter 9.5 of [1]. From
        %Chapter 10.5.10, comparisons of two infinite values of the same
        %sign are true and also true if a is empty. Compare to the
        %containedBy function.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
            
            c=containedBy(a,b);
            c=all(c(:));
        end
        
        function c=disjoint(a,b)
        %%DISJOINT Determine whether two sets of intervals are disjoint.
        %          That is, none of the sets intersect anywhere
        %
        %This function is defined in Table 9.3 of Chapter 9.5 of [1]. In
        %Chapter 10.5.10, it is specified that an empty interval is both a
        %subset of any interval and disjoint from them. Thus, a check for
        %an empty interval makes that true here.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
            
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end
        
            c = (a.lb>b.ub) | (a.ub<b.lb) |isnan(a.lb) |isnan(b.lb);
            c=all(c(:));
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Numeric functions of intervals of in Chapter 10.5.10 of  %%%
%%% IEEE Std 1788-2015 that are not in Table 9.3             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function c=isEmpty(a)
        %%ISEMPTY This returns true if the interval is empty and false
        %         otherwise. This works on each interval in a set
        %         separately. Empty intervals are marked using NaN values.
        %
        %This function is defined in Chapters 10.5.10 of [1].
        %       
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            c=isnan(a.lb)|isnan(a.ub);
        end
        
        function c=isnan(a)
        %%ISNAN Overloaded isnan function in Matlab. Note that an Inteval
        %       of [NaN,NaN] is used to indicate an empty Interval. This
        %       function just calls isEmpty. See the comments to isEmpty
        %       for more information.
        
            c=isEmpty(a);
        end
        
        function c=isEntire(a)
        %%ISENTIRE This returns true if an Interval goes from -Inf to Inf.
        %          This does all Intervals in a set separately.
        %
        %This function is defined in Chapters 10.5.10 of [1].
        %       
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            c=a.lb==-Inf & a.ub ==Inf;
        end
        
        function c=less(a,b)
        %%LESS Determine whether the set of Intervals a is weakly less than
        %      the set of Intervals b. This means that a.lb<=b.lb and
        %      a.ub<=b.ub for all intervals in a and b.
        %    
        %This function is defined in Chapter 10.5.10 of [1]. From Table
        %10.4, it is required that two empty intervals return true.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.

            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end
        
            c=(a.lb<=b.lb & a.ub<=b.ub) | (isnan(a.lb)&isnan(b.lb));
            c=all(c(:));
        end
        
        function c=precedes(a,b)
        %%PRECEDES Determine whether the set of Intervals a precedes the
        %       set of Intervals b. This means that a.ub<=b.lb for all
        %       intervals in a and b.
        %
        %This function is defined in Chapter 10.5.10 of [1]. From Table
        %10.4, it is defined that an empty interval preces all intervals.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.  
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end
            
            c=a.ub<=b.lb | isnan(a.lb) | isnan(b.lb);
            c=all(c(:));
        end
        
        function c=strictLess(a,b)
        %%LESS Determine whether the set of Intervals a is stricly less than
        %      the set of Intervals b. This means that a.lb<b.lb and
        %      a.ub<b.ub for all intervals in a and b.
        %    
        %This function is defined in Chapter 10.5.10 of [1]. Comparisons
        %between two infinite values of the same sign are considered true.
        %Also, comparisons between empty sets are only true if both are
        %empty sets.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end
            
            c=(((a.lb < b.lb)|(a.lb==-Inf & b.lb==-Inf)|(a.lb==Inf & b.lb==Inf))) & ((a.ub < b.ub)|(a.ub==-Inf &b.ub==-Inf)|(a.ub==Inf & b.ub==Inf))|(isnan(a.lb)&isnan(b.lb));
            c=all(c(:));
        end
        
        function c=strictPrecedes(a,b)
        %%STRICTPRECEDES Determine whether the set of Intervals a strictly
        %       precedes the set of Intervals b. This means that a.ub< b.lb
        %       for all intervals in a and b.
        %
        %This function is defined in Chapter 10.5.10 of [1]. From Table
        %10.4, it is defined that an empty interval strictly preces all
        %intervals. Compare this function to before.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end
            
            c=a.ub<b.lb | isnan(a.lb) | isnan(b.lb);
            c=all(c(:));
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Forward Mode Elementary Functions of Chapter 10.5.3  of  %%%
%%% the IEEE Std 1788-2015                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function a=neg(a)
        %%NEG The negation function (find the negative of an interval).
        %
        %The negation function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %The method for negating an interval is taken from Section III of
        %[2].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
            
            temp=a.lb;
            a.lb=-a.ub;
            a.ub=-temp;
        end
 
        function c=uminus(a)
        %%UMINUS Overloaded negation operator (the unary minus operator in
        %        Matlab). This function just calls neg. See the comments to
        %        neg for more information.
            
            c=neg(a);
        end

        function [c,undefined]=add(a,b)
        %%ADD Find the sum of two sets of Intervals. An input can also be a
        %     set of points.
        %
        %The addition function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %Basic interval addition is described in Section III of [2] and
        %proper rounding modes have been set. If an Inf - Inf or -Inf + Inf
        %occurs, then an empty interval is returned.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            
            if(~isa(b,'Interval'))
                b=Interval(b);
            end

            setProcRoundingMode(0);%Round to -Inf
            lbNew=a.lb+b.lb;
            setProcRoundingMode(3);%Round to +Inf
            ubNew=a.ub+b.ub;

            %If an Inf-Inf case occured, make sure that both upper and
            %lower bounds are set to NaN to indicate an empty interval.
            undefined=xor(isnan(lbNew),isnan(ubNew));
            lbNew(undefined)=NaN;
            ubNew(undefined)=NaN;
            
            c=Interval(lbNew,ubNew);
        end
        
        function c = plus(a,b)
        %%PLUS Overloaded addition (+) operator. This function just calls
        %      add. See the comments to add for more information.
            
            c=add(a,b);
        end

        function [c,undefined]=sub(a,b)
        %%NEG The subtraction function a-b. a and b can be set of intervals
        %     or points.
        %
        %The subtraction function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %Basic interval subtraction is described in Section III of [2] and
        %proper rounding modes have been set. If an Inf - Inf or -Inf + Inf
        %occurs, then an empty interval is returned.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            
            if(~isa(b,'Interval'))
                b=Interval(b);
            end
        
            setProcRoundingMode(0);%Round to -Inf
            lbNew=a.lb-b.ub;
            setProcRoundingMode(3);%Round to +Inf
            ubNew=a.ub-b.lb;
            
            %If an Inf-Inf case occured, make sure that both upper and
            %lower bounds are set to NaN to indicate an empty interval.
            undefined=xor(isnan(lbNew),isnan(ubNew));
            lbNew(undefined)=NaN;
            ubNew(undefined)=NaN;
            
            c=Interval(lbNew,ubNew);
        end
        
        function c=minus(a,b)
        %%MINUS Overloaded subtraction operator. This just calls the sub
        %%function. See the comments to sub for more information.
        
            c=sub(a,b);
        end

        function c=mul(a,b)
        %%MUL Multiplication of two intervals or of an interval and a
        %     scalar (element-wise when given matrices of intervals).
        %
        %The multiplication function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %The standard method of interval multilication is given in Section
        %III of [2]. However, given the use of NaN to indicate an empty
        %interval and the possibility of +/-Inf values, a few extensions
        %are necessary. If any Interval is empty (NaN), then the resulting
        %product is empty. If Inf*0 is encountered, it is ignored. For
        %example, [0,Inf]*[-Inf,0]=[-Inf,0]. Only [-Inf,Inf]*0 will result
        %in an empty Interval [NaN,NaN].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009

            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            
            if(~isa(b,'Interval'))
                b=Interval(b);
            end

            setProcRoundingMode(0);%Round to -Inf
            lbNew=min(min(min(a.lb.*b.lb,a.lb.*b.ub),a.ub.*b.lb),a.ub.*b.ub);
            setProcRoundingMode(3);%Round to +Inf
            ubNew=max(max(max(a.lb.*b.lb,a.lb.*b.ub),a.ub.*b.lb),a.ub.*b.ub);

            c=Interval(lbNew,ubNew);
        end
        
        function c=times(a,b)
        %%TIMES Overloaded element-wise multiplication operator (the .*
        %       multiplication). This function just calls mul. See the
        %       comments to mul for more information.
        
            c=mul(a,b);
        end

        function [c,d]=div(a,b)
        %%DIV  Evaluate a divided by a (element-wise), where a and/or b is
        %      an Interval. This is performed by multiplying a (with the
        %      times function) by the inverse of b (found using the recip
        %      function). The availability of two outputs reflects show the
        %      recip function handles intervals spanning zero.
        %
        %The division function is in Table 9.1 of Chapter 9.1 of [1]. The
        %ability for two-output division is discussed in Chapter 10.5.5.
        %This function is implemented using multiplication and the recip
        %function. See the somments to the recip function for more details
        %on the behaviour when dividing by 0 when considering the standard
        %in [1].
        %
        %When an interval spans 0, there should actually be two outputs,
        %because the reciprocal produces two intervals. If an interval does
        %not span 0, then the second interval will be set to the empty
        %interval [NaN,NaN]. However, if only one output of this function
        %is requested, then in such an instance -Inf to Inf will be the
        %bounds of the single interval returned by the recip function. How
        %that effects the output through multiplication with a is
        %determined by the times function.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            if(~isa(b,'Interval'))
                c=times(a,1./b);
                d=[];
            elseif(nargout<2)
                c=times(a,recip(b));
                d=[];
            else
                [b1,b2]=recip(b);
                c=times(a,b1);
                d=times(a,b2);
            end
        end
        
        function [c,d]=mulRevToPair(a,b)
        %%MULREVTOPAIR Divide a./b element-wise, obtaining two outputs.
        %          This is just another name for the div function producing
        %          two outputs. Unlike the div function, this function will
        %          not take the union of the two outputs if only one output
        %          is requested.
        %
        %The mulRevToPair function is defined in Chapter 10.5.5 of [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
            
            [c,d]=div(a,b);
        end
        
        function [c,d]=rdivide(a,b)
        %%RDIVIDE Overloaded element-wise division operation (the ./
        %         operation). This function just calls div. See the
        %         comments to div for more information.
        
            if(nargout<2)
                c=div(a,b);
            else
                [c,d]=div(a,b);
            end
        end
        
        function [c,d]=ldivide(a,b)
        %%LDIVIDE Overloaded element-wise left-hand division operator (the
        %        .\ operation). This function just calls div with the
        %        proper order of parameters.

            if(nargout<2)
                c=div(b,a);
            else
                [c,d]=div(b,a);
            end
        end

        function [c,d]=recip(a)
        %%RECIP Find the reciprocal (inverse) of an Interval (element
        %       -wise). When an interval spans 0, there should actually be
        %       two outputs, because the reciprocal produces two intervals.
        %       If an interval does not span 0, then the second interval
        %       will be set to the empty interval [NaN,NaN]. However, if
        %       only one output of this function is requested, then in such
        %       an instance -Inf to Inf will be the bounds of the single
        %       interval returned.
        %
        %The reciprocal function is in Table 9.1 of Chapter 9.1 of [1].
        %Chapter 10.5.5's discussion two-output division. Unlike the
        %example in [1], this function returns an infite (-Inf,Inf)
        %interval when dividing by 0. In [1], the example has no solution
        %when an interval not spanning zero is divided by zero, but has a
        %solution when it spans zero. What is done here is more consistent,
        %because multiplication by the inverse equals division.
        %
        %A typical implementation of the reciprocal function can be pulled
        %from the definition of divions in [2]. However, the ability to
        %obtain two intervals goes beyind that text.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] G.I. Hargreaves. Interval Analysis in MATLAB. Technical Report 416,
        %    Manchester Centre for Computational Mathematics, Dept. Math., 2002.
        %    [Online] http://www.manchester.ac.uk/mims/eprints.
            
            setProcRoundingMode(0);%Round to -Inf
            lbNew=1./a.ub;
            
            setProcRoundingMode(3);%Round to +Inf
            ubNew=1./a.lb;
            
            %Range touches 0
            %-0 case not always caught arithmetically
            lbNew(a.ub==0)=-Inf;
            ubNew(a.lb==0)=Inf;
            
            %Range crosses 0
            fullRange=(a.ub > 0)&(a.lb < 0);
            if(nargout<2)
                lbNew(fullRange)=-Inf;
                ubNew(fullRange)=Inf;
                d=[];
            else
                lbTwo=nan(size(a.lb));
                ubTwo=nan(size(a.ub));
                
                lbTwo(fullRange)=lbNew(fullRange);
                ubTwo(fullRange)=Inf;
                lbNew(fullRange)=-Inf;
                
                d=Interval(lbTwo,ubTwo);
            end

            c=Interval(lbNew,ubNew);
        end
        
        function c=sqr(a)
        %%SQR Square an interval. That is, multiply it by itself.
        %
        %The square function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %THis function just calls mul(a,a).
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            c=mul(a,a);
        end

        function a=sqrt(a)
        %%SQRT Evaluate the square root of an Interval. Intervals including
        %      negative numbers result in empty intervals being returned.
        %
        %The square root function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %The method for computing the square root of an interval is taken
        %from Section III of [2]. However, for getting the lower and upper
        %bounds of the exponent of a scalar, the sqrt values are
        %computed with different rounding modes rather than using a
        %different number of terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
            
            %Make intervals outside the valid range empty.
            sel=a.lb<0;
            a.lb(sel)=NaN;
            a.ub(sel)=NaN;
        
            setProcRoundingMode(0);%Round to -Inf
            a.lb=sqrt(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=sqrt(a.ub);
        end
        
        function d=fma(a,b,c)
        %%FMA Fused multiply add. This function evaluates a.*b+c, where one
        %     or more of the parameters is an Interval.
        %
        %The fma function is in Table 9.1 of Chapter 9.1 of [1]. If
        %implemented in assembly, this function could be potentially
        %efficient as processors often have an fma instruction for floating
        %point values. However, in Matlab, this just relies on the
        %overloaded ./ and plus operators and is just a.*b+c.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            d=a.*b+c;
        end

        function result=pown(a,b)
        %%POWN Raise an interval to a scalar integer power. The power b
        %      must be >=0.
        %
        %The pown function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %This function is implemented as described in Section III of [2].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            if(numel(b)~=1||~isnumeric(b)||~isreal(b)||fix(b)~=b||b<0)
                error('The second input must be a scalar integer >=0')
            end
        
            numEls=numel(a.lb);
            
            %Create an empty interval for the results.
            result=Interval.empty(size(a.lb));
            
            s.type='()';
            s.subs{1}=[];
            for curEl=1:numEls
                s.subs{1}=curEl;
                aCur=a.subsref(s);
                if(b==0)
                    result=result.subsasgn(s,Interval(1,1));
                elseif(aCur.lb>=0||mod(b,2)~=0)
                    setProcRoundingMode(0);%Round to -Inf
                    lowerBound=aCur.lb.^b;
                    setProcRoundingMode(3);%Round to +Inf
                    upperBound=aCur.ub.^b;
                    result=result.subsasgn(s,Interval(lowerBound,upperBound));
                elseif(aCur.ub<=0&&mod(b,2)==0)
                    setProcRoundingMode(0);%Round to -Inf
                    lowerBound=aCur.ub.^b;
                    setProcRoundingMode(3);%Round to +Inf
                    upperBound=aCur.lb.^b;
                    result=result.subsasgn(s,Interval(lowerBound,upperBound));
                else
                    setProcRoundingMode(3);%Round to +Inf
                    result=result.subsasgn(s,Interval(0,max(aCur.lb.^b,aCur.ub.^b)));
                end
            end
        end

        function c=pow(a,b)
        %%POW Raise the Interval a to a real power b (element by element)
        %     that is >=0. b can be a scalar (all Intervals raised to the
        %     same power) or a matrix having the same dimensionality as the
        %     intervals.
        %
        %The pown function is in Table 9.1 of Chapter 9.1 of [1] and is
        %defined to be exp(times(b,log(a))). Rather than going by 1 where
        %negative numbers have no value, the restrictions of the log
        %function to intervals are applied.
        %
        %This function is implemented as described in Section III of [2].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
            
            if(~isa(b,'Interval')&&(~isnumeric(b)||~isreal(b)))
                error('The second input must be a scalar integer >=0')
            end
        
            c=exp(times(b,log(a)));
        end

        function c=power(a,b)
        %POWER Overloaded eleement-wise power operation (the .^ operation).
        %      This function just calls pow. See the comments to pow for
        %      more information.

            c=pow(a,b);
        end

        function a=exp(a)
        %%EXP Evaluate the exponential of an Interval (Euler's number
        %     raised to it).
        %
        %The exponential function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %The method for computing the exponential of an interval is taken
        %from Section III of [2]. However, for getting the lower and upper
        %bounds of the exponent of a scalar, the exponent values are
        %computed with different rounding modes rather than using a
        %different number of terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009

            setProcRoundingMode(0);%Round to -Inf
            a.lb=exp(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=exp(a.ub);
        end
        
        function a=exp2(a)
        %%EXP2 Evaluate 2 raised to  an Interval.
        %
        %The exp2 function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %The method for computing the exp2 of an interval is the same as
        %that of computing exp of an interval taken from Section III of
        %[2]. However, for getting the lower and upper bounds of exp2 of a
        %scalar, the powers are computed with different rounding modes
        %rather than using a different number of terms in a series
        %expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009

            setProcRoundingMode(0);%Round to -Inf
            a.lb=2.^a.lb;
            setProcRoundingMode(3);%Round to +Inf
            a.ub=2.^a.ub;
        end

        function a=exp10(a)
        %%EXP10 Evaluate 10 raised to  an Interval.
        %
        %The exp10 function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %The method for computing the exp10 of an interval is the same as
        %that of computing exp of an interval taken from Section III of
        %[2]. However, for getting the lower and upper bounds of exp10 of a
        %scalar, the powers are computed with different rounding modes
        %rather than using a different number of terms in a series
        %expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009

            setProcRoundingMode(0);%Round to -Inf
            a.lb=10.^a.lb;
            setProcRoundingMode(3);%Round to +Inf
            a.ub=10.^a.ub;
        end

        function a=log(a)
        %%LOG Evaluate the natural logarithm of an Interval. Intervals
        %     crossing into negative territory cause empty intervals to be
        %     returned.
        %
        %The natural logarithm function is in Table 9.1 of Chapter 9.1 of
        %[1].
        %
        %The method for computing the natural logarithm of an interval is
        %done analogously to how the inverse tangent is done in Section III
        %of [2] (because the natural logarithm function is similarly
        %monotonic). However, for getting the lower and upper bounds of the
        %natural logarithm of a scalar, the log values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
            
            %Make intervals outside the valid range empty.
            sel=(a.lb<0);
            a.lb(sel)=NaN;
            a.ub(sel)=NaN;
        
            setProcRoundingMode(0);%Round to -Inf
            a.lb=log(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=log(a.ub);
        end

        function a=log2(a)
        %%LOG2 Evaluate the base-2 logarithm of an Interval. Intervals
        %     crossing into negative territory cause empty intervals to be
        %     returned.
        %
        %The base-2 logarithm function is in Table 9.1 of Chapter 9.1 of
        %[1].
        %
        %The method for computing the base-2 logarithm of an interval is
        %done analogously to how the inverse tangent is done in Section III
        %of [2] (because the base-2 logarithm is similarly monotonic).
        %However, for getting the lower and upper bounds of the base-2
        %logarithm of a scalar, the log2 values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
            
            %Make intervals outside the valid range empty.
            sel=(a.lb<0);
            a.lb(sel)=NaN;
            a.ub(sel)=NaN;
        
            setProcRoundingMode(0);%Round to -Inf
            a.lb=log2(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=log2(a.ub);
        end
        
        function a=log10(a)
        %%LOG10 Evaluate the base-10 logarithm of an Interval. Intervals
        %     crossing into negative territory cause empty intervals to be
        %     returned.
        %
        %The base-0 logarithm function is in Table 9.1 of Chapter 9.1 of
        %[1].
        %
        %The method for computing the base-10 logarithm of an interval is
        %done analogously to how the inverse tangent is done in Section III
        %of [2] (because the base-10 logarithm is similarlymonotonic).
        %However, for getting the lower and upper bounds of the base-10
        %logarithm of a scalar, the log10 values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
            
            %Make intervals outside the valid range empty.
            sel=(a.lb<0);
            a.lb(sel)=NaN;
            a.ub(sel)=NaN;
        
            setProcRoundingMode(0);%Round to -Inf
            a.lb=log10(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=log10(a.ub);
        end

        function result=sin(a)
        %%SIN Evaluate the sine of an Interval.
        %
        %The sine function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %The method for computing the sine of an interval is taken from
        %Equation 6 of [2]. However, for getting the lower and upper bounds
        %of the sine of a scalar, the sine values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009

            numEls=numel(a.lb);
            
            %Create an empty interval for the results.
            result=Interval.empty(size(a.lb));
            
            piInt=Interval.pi();
            halfPiInt=piInt./2;
            
            region1=Interval(-halfPiInt.lb,halfPiInt.lb);
            region2=Interval(halfPiInt.ub,piInt.lb);
            region3=Interval(0,piInt.lb);
            region4=Interval(-piInt.lb,0);
            
            s.type='()';
            s.subs{1}=[];
            for curEl=1:numEls
                s.subs{1}=curEl;
                aCur=a.subsref(s);
                if(aCur.subset(region1))
                    setProcRoundingMode(0);%Round to -Inf
                    sinLo=sin(aCur.lb);
                    setProcRoundingMode(3);%Round to +Inf
                    sinHi=sin(aCur.ub);
                    result=result.subsasgn(s,Interval(sinLo,sinHi));
                elseif(aCur.subset(region2))
                    setProcRoundingMode(0);%Round to -Inf
                    sinLo=sin(aCur.ub);
                    setProcRoundingMode(3);%Round to +Inf
                    sinHi=sin(aCur.lb);
                    result=result.subsasgn(s,Interval(sinLo,sinHi));
                elseif(aCur.subset(region3))
                    setProcRoundingMode(0);%Round to -Inf
                    sinLo=sin(aCur.lb);
                    sinHi=sin(aCur.ub);
                    result=result.subsasgn(s,Interval(min(sinLo,sinHi),1));
                elseif(aCur.subset(region4))
                    result=result.subsasgn(s,-sin(-aCur));
                else
                    result=result.subsasgn(s,Interval(-1,1));
                end
            end
        end
        
        function result=cos(a)
        %%COS Evaluate the cosine of an Interval.
        %
        %The cosine function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %The method for computing the cosine of an interval is taken from
        %Equation 7 of [2]. However, for getting the lower and upper bounds
        %of the cosine of a scalar, the cosine values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009

            numEls=numel(a.lb);
            
            %Create an empty interval for the results.
            result=Interval.empty(size(a.lb));
            
            piInt=Interval.pi();
            halfPiInt=piInt./2;
            
            region1=Interval(0,piInt.lb);
            region2=Interval(-piInt.lb,0);
            region3=Interval(-halfPiInt.lb,halfPiInt.ub);
            
            s.type='()';
            s.subs{1}=[];
            for curEl=1:numEls
                s.subs{1}=curEl;
                aCur=a.subsref(s);
                if(aCur.subset(region1))
                    setProcRoundingMode(0);%Round to -Inf
                    cosLo=cos(aCur.ub);
                    setProcRoundingMode(3);%Round to +Inf
                    cosHi=cos(aCur.lb);
                    result=result.subsasgn(s,Interval(cosLo,cosHi));
                elseif(aCur.subset(region2))
                    result=result.subsasgn(s,-cos(-aCur));
                elseif(aCur.subset(region3))
                    setProcRoundingMode(0);%Round to -Inf
                    cosLo=cos(aCur.lb);
                    cosHi=cos(aCur.ub);
                    result=result.subsasgn(s,Interval(min(cosLo,cosHi),1));
                else
                    result=result.subsasgn(s,Interval(-1,1));
                end
            end
        end
        
        function a=tan(a)
        %%TAN Evaluate the tangent of an Interval.
        %
        %The tangent function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %The method for computing the tangent of an interval is taken from
        %Equation 8 of [2]. However, for getting the lower and upper bounds
        %of the tangent of a scalar, the necessary sin and cosine values
        %are computed with different rounding modes rather than using a
        %different number of terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            piInt=Interval.pi();
            halfPiInt=piInt./2;
            
            sel=(a.lb>=-halfPiInt.lb) & (a.ub<=halfPiInt.lb);
            setProcRoundingMode(0);%Round to -Inf
            sinLo=sin(a.lb(sel));
            cosHi=cos(a.ub(sel));
            setProcRoundingMode(3);%Round to +Inf
            sinHi=sin(a.ub(sel));
            cosLo=cos(a.lb(sel));
            
            setProcRoundingMode(0);%Round to -Inf
            a.lb(sel)=sinLo/cosLo;
            
            setProcRoundingMode(3);%Round to +Inf
            a.ub(sel)=sinHi/cosHi;

            a.lb(~sel)=-Inf;
            a.ub(~sel)=Inf;
        end

        function a=asin(a)
        %%ASIN Evaluate the inverse sine of an Interval. Intervals
        %      extending beyond +/-1 cause empty intervals to be returned.
        %
        %The inverse sine function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %The method for computing the inverse sine of an interval is done
        %analogously to how the inverse tangent is done in Section III of
        %[2] (because the inverse sine function is similarly monotonic).
        %However, for getting the lower and upper bounds of the inverse
        %sine of a scalar, the asin values are computed with different
        %rounding modes rather than using a different number of terms in a
        %series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
            
            %Make intervals outside the valid range empty.
            sel=(a.ub<-1)|(a.lb>1);
            a.lb(sel)=NaN;
            a.ub(sel)=NaN;
            
            setProcRoundingMode(0);%Round to -Inf
            a.lb=asin(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=asin(a.ub);
        end
        
        function a=acos(a)
        %%ACOS Evaluate the inverse cosine of an Interval. Intervals
        %      extending beyond +/-1 cause empty intervals to be returned.
        %
        %The inverse cosine function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %The method for computing the inverse cosine of an interval is done
        %analogously to how the inverse tangent is done in Section III of
        %[2] (because the inverse cosine function is similarly monotonic).
        %However, for getting the lower and upper bounds of the inverse
        %cosine of a scalar, the acos values are computed with different
        %rounding modes rather than using a different number of terms in a
        %series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
            
            %Make intervals outside the valid range empty.
            sel=(a.ub<-1)|(a.lb>1);
            a.lb(sel)=NaN;
            a.ub(sel)=NaN;
            
            setProcRoundingMode(0);%Round to -Inf
            temp=a.lb;
            a.lb=acos(a.ub);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=acos(temp);
        end

        function a=atan(a)
        %%ATAN Evaluate the inverse tangent of an Interval.
        %
        %The inverse tangent function is in Table 9.1 of Chapter 9.1 of
        %[1].
        %
        %The method for computing the inverse tangent of an interval is
        %taken from Section III of [2]. However, for getting the lower and
        %upper bounds of the inverse tangent of a scalar, the atan values
        %are computed with different rounding modes rather than using a
        %different number of terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009

            setProcRoundingMode(0);%Round to -Inf
            a.lb=atan(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=atan(a.ub);
        end
        
        function [result,result2]=atan2(yVals,xVals)
        %%ATAN2 Evaluate the four-quadrant inverse tangent of an interval.
        %       This ranges from -pi to pi with the break on the negative x
        %       axis. The function can be called with one output, in which
        %       case an interval going across quadrants 1 and 4 returns and
        %       Interval from -pi to pi. When two outputs are requested,
        %       then such an interval is split into two halves. When an
        %       interval in quadrant 1 ends at y=0, then it is not
        %       considered crossing. When an interval in quadrant 3 ends at
        %       y=0, then it is also not considered crossing.
        %
        %The inverse tangent function is in Table 9.1 of Chapter 9.1 of
        %[1].
        %
        %The implementation of this function is actually quite tricky,
        %because one must pay attention to which quadrants the input
        %interval crosses. The inverse tangent values go from -pi to pi.
        %
        %The four quadrants are
        %    y
        %    |
        %  1 | 2
        %    |
        % ___0___ x
        %    |
        %  4 | 3
        %    |
        %The value of the angle increases counterclockwise starting from
        %the x axis.
        %
        %Example: One of each type of type of quadrant plus one of each
        %spanning quadrant plus the origin.
%         xLb=zeros(10,1);
%         xUb=zeros(10,1);
%         yLb=zeros(10,1);
%         yUb=zeros(10,1);
%         %Quadrant 1
%         xLb(1)=-1;
%         xUb(1)=-0.5;
%         yLb(1)=0.5;
%         yUb(1)=1;
%         %Quadrant 2
%         xLb(2)=0.5;
%         xUb(2)=1;
%         yLb(2)=0.5;
%         yUb(2)=1;
%         %Quadrant 3
%         xLb(3)=0.5;
%         xUb(3)=1;
%         yLb(3)=-1;
%         yUb(3)=-0.5;
%         %Quadrant 4
%         xLb(4)=-1;
%         xUb(4)=-0.5;
%         yLb(4)=-1;
%         yUb(4)=-0.5;
%         %Quadrants 1->2
%         xLb(5)=-1;
%         xUb(5)=1;
%         yLb(5)=0.5;
%         yUb(5)=1;
%         %Quadrants 3->4
%         xLb(6)=-1;
%         xUb(6)=1;
%         yLb(6)=-1;
%         yUb(6)=-0.5;
%         %Quadrants 2->3
%         xLb(7)=0.5;
%         xUb(7)=1;
%         yLb(7)=-1;
%         yUb(7)=1;
%         %Quadrants 1->4
%         xLb(8)=-1;
%         xUb(8)=-0.5;
%         yLb(8)=-1;
%         yUb(8)=1;
%         %All quadrants
%         xLb(9)=-1;
%         xUb(9)=1;
%         yLb(9)=-1;
%         yUb(9)=1;
%         %The origin
%         xLb(10)=0;
%         xUb(10)=0;
%         yLb(10)=0;
%         yUb(10)=0;
%         xVals=Interval(xLb,xUb);
%         yVals=Interval(yLb,yUb);
%         [s1,s2]=atan2(yVals,xVals)
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.

             numEls=numel(xVals.lb);
            
            %Create an empty interval for the results.
            result=Interval.empty(size(xVals.lb));
            result2=Interval.empty(size(xVals.lb));
            
            s.type='()';
            s.subs{1}=[];
            for curEl=1:numEls
                s.subs{1}=curEl;
                x=xVals.subsref(s);
                y=yVals.subsref(s);
                
                if(y.lb>=0&&x.lb<=0&&x.ub<=0)
                %Quadrant 1 only.
                    setProcRoundingMode(0);%Round to -Inf
                    minVal=atan2(y.ub,x.ub);
                    setProcRoundingMode(3);%Round to +Inf
                    maxVal=atan2(y.lb,x.lb);
                elseif(y.lb>=0 &&x.lb>0)
                %Quadrant 2 only
                    setProcRoundingMode(0);%Round to -Inf
                    minVal=atan2(y.lb,x.ub);
                    setProcRoundingMode(3);%Round to +Inf
                    maxVal=atan2(y.ub,x.lb);
                elseif(y.lb>=0)
                %Spanning quadrants 1 and 2.
                    setProcRoundingMode(0);%Round to -Inf
                    minVal=atan2(y.lb,x.ub);
                    setProcRoundingMode(3);%Round to +Inf
                    maxVal=atan2(y.lb,x.lb);
                elseif(y.ub<0&&x.lb<=0&&x.ub<=0)
                %Quadrant 4 only
                    setProcRoundingMode(0);%Round to -Inf
                    minVal=atan2(y.ub,x.lb);
                    setProcRoundingMode(3);%Round to +Inf
                    maxVal=atan2(y.lb,x.lb);
                elseif(y.ub<0&&x.lb>=0&&x.ub>=0)
                %Quadrant 3 only.
                    setProcRoundingMode(0);%Round to -Inf
                    minVal=atan2(y.lb,x.lb);
                    setProcRoundingMode(3);%Round to +Inf
                    maxVal=atan2(y.ub,x.ub);
                elseif(y.ub<0)
                %Spanning quadrants 3 and 4
                    setProcRoundingMode(0);%Round to -Inf
                    minVal=atan2(y.ub,x.lb);
                    setProcRoundingMode(3);%Round to +Inf
                    maxVal=atan2(y.ub,x.ub);
                elseif(x.lb>=0)
                %Spanning quadrants 2 and 3.
                    setProcRoundingMode(0);%Round to -Inf
                    minVal=atan2(y.lb,x.lb);
                    setProcRoundingMode(3);%Round to +Inf
                    maxVal=atan2(y.ub,x.lb);
                elseif(x.ub<=0&&nargout>1||(y.ub==0&&x.ub==0))
                %Spanning quadrants 1 and 4 if we want two outputs. We also
                %get here if we are in the third quadrant and y.ub=0 and
                %x.ub=0.
                    piInt=Interval.pi();

                    setProcRoundingMode(0);%Round to -Inf
                    minVal=atan2(y.ub,x.ub);
                    maxVal=piInt.ub;
                    
                    minVal2=-piInt.ub;
                    setProcRoundingMode(3);%Round to +Inf
                    maxVal2=atan2(y.lb,x.ub);
                    
                    if(y.ub==0&&x.ub==0)
                        %If the Interval only touches the axis and does not
                        %actually cross, then there is no top half.
                        minVal=minVal2;
                        maxVal=maxVal2;
                    else
                        result2=result2.subsasgn(s,Interval(minVal2,maxVal2));
                    end
                else
                 %Spanning all four quadrants or just quadrants 1 and 4 but
                 %we do not want two outputs.
                    piInt=Interval.pi();
                    setProcRoundingMode(0);%Round to -Inf
                    minVal=-piInt.lb;
                    setProcRoundingMode(3);%Round to +Inf
                    maxVal=piInt.ub;
                end
                
                result=result.subsasgn(s,Interval(minVal,maxVal));
            end
            
        end

        function a=sinh(a)
        %%SINH Evaluate the hyperbolic sine of an Interval.
        %
        %The hyperbolic sine function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %The method for computing the hyperbolic sine of an interval is
        %done analogously to how the inverse tangent is done in Section III
        %of [2] (because the hyperbolic sine function is similarly
        %monotonic). However, for getting the lower and upper bounds of the
        %hyperbolic sine of a scalar, the sinh values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
            
            setProcRoundingMode(0);%Round to -Inf
            a.lb=sinh(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=sinh(a.ub);
        end
        
        function a=cosh(a)
        %%COSH Evaluate the hyperbolic cosine of an Interval.
        %
        %The hyperbolic cosine function is in Table 9.1 of Chapter 9.1 of
        %[1].
        %
        %The cosh function is symmetric about the origin. Thus, if the
        %interval contains the origin, then the minimum is zero. If the
        %interval does not contain the origin, then cosh increases
        %monotonically and the method to use is similar to that of atan in
        %[2].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
                        
            %Case 1 spans the origin.
            sel1=a.lb<=0&a.ub>=0;
            
            %The origin is the lower bound.
            temp=a.lb(sel1);
            a.lb(sel1)=0;
            setProcRoundingMode(3);%Round to +Inf
            a.ub(sel1)=max(cosh(temp),a.ub(sel1));
            
            %Case 2 means both are on one side or the other of the origin.
            a(~sel1)=abs(a(~sel1));
            setProcRoundingMode(0);%Round to -Inf
            a.lb(~sel1)=cosh(a.lb(~sel1));
            setProcRoundingMode(3);%Round to +Inf
            a.ub(~sel1)=cosh(a.ub(~sel1));
        end
        
        function a=tanh(a)
        %%TANH Evaluate the hyperbolic tangent of an Interval.
        %
        %The hyperbolic tangent function is in Table 9.1 of Chapter 9.1 of
        %[1].
        %
        %The method for computing the hyperbolic tangent of an interval is
        %done analogously to how the inverse tangent is done in Section III
        %of [2] (because the hyperbolic tangent function is similarly
        %monotonic). However, for getting the lower and upper bounds of the
        %hyperbolic tangent of a scalar, the tanh values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
            
            setProcRoundingMode(0);%Round to -Inf
            a.lb=tanh(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=tanh(a.ub);
        end

        function a=asinh(a)
        %%ASINH Evaluate the inverse hyperbolic sine of an Interval.
        %
        %The inverse hyperbolic sine function is in Table 9.1 of Chapter
        %9.1 of [1].
        %
        %The method for computing the inverse hyperbolic sine of an
        %interval is done analogously to how the inverse tangent is done in
        %Section III of [2] (because the inverse hyperbolic sine function
        %is similarly monotonic). However, for getting the lower and upper
        %bounds of the inverse hyperbolic sine of a scalar, the ashinh 
        %values are computed with different rounding modes rather than
        %using a different number of terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
            
            setProcRoundingMode(0);%Round to -Inf
            a.lb=asinh(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=asinh(a.ub);
        end
        
        function a=acosh(a)
        %%ACOSH  Evaluate the inverse hyperbolic cosine of an Interval.
        %        Intervals going below 1 cause empty intervals to be
        %        returned.
        %
        %The inverse hyperbolic cosine function is in Table 9.1 of Chapter
        %9.1 of [1].
        %
        %The method for computing the inverse hyperbolic cosine of an
        %interval is done analogously to how the inverse tangent is done in
        %Section III of [2] (because the inverse hyperbolic cosine function
        %is similarly monotonic). However, for getting the lower and upper
        %bounds of the inverse hyperbolic cosine of a scalar, the acosh
        %values are computed with different rounding modes rather than
        %using a different number of terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            %Make intervals outside the valid range empty.
            sel=a.ub<1;
            a.lb(sel)=NaN;
            a.ub(sel)=NaN;
            
            setProcRoundingMode(0);%Round to -Inf
            a.lb=acosh(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=acosh(a.ub);
        end
        
        function a=atanh(a)
        %%ATANH Evaluate the inverse hyperbolic tangent of an Interval.
        %       Intervals going outside the range of [-1,1] cause empty
        %       intervals to be returned.
        %
        %The inverse hyperbolic tangent function is in Table 9.1 of Chapter
        %9.1 of [1].
        %
        %The method for computing the inverse hyperbolic tangent of an
        %interval is done analogously to how the inverse tangent is done in
        %Section III of [2] (because the inverse hyperbolic tangent
        %function is similarly monotonic). However, for getting the lower
        %and upper bounds of the inverse hyperbolic tangent of a scalar,
        %the atanh values are computed with different rounding modes rather
        %than using a different number of terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            %Make intervals outside the valid range empty.
            sel=(a.ub<-1)|(a.lb>1);
            a.lb(sel)=NaN;
            a.ub(sel)=NaN;
            
            setProcRoundingMode(0);%Round to -Inf
            a.lb=atanh(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=atanh(a.ub);
        end

        function b=sign(a)
        %%SIGN Compute the signum function for an interval. This consists
        %      of an Interval that is either a point or a span. If an 
        %      entire interval is negative, then -1 is returned. If it is
        %      all positive, then +1 is returned. If the interval is just
        %      the point 0, then 0 is returned. Otherwise, the return is an
        %      Interval that spans the appropriate regions.
        %
        %The sign function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            b=Interval(sign(a.lb),sign(a.ub));
        end
        
        function a=ceil(a)
        %%CEIL Compute the ceiling function for an Interval. This only
        %      affects the endpoints of the interval; non-integer values
        %      are rounded to the next highest integer.
        %
        %The ceil function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.  
            
           a=Interval(ceil(a.lb),ceil(a.ub)); 
        end
        
        function a=floor(a)
        %%FLOOR Compute the floor function for an Interval. This only
        %      affects the endpoints of the interval; non-integer values
        %      are rounded to the next lowest integer.
        %
        %The ceil function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.  
            
           a=Interval(floor(a.lb),floor(a.ub)); 
        end
        
        function a=trunc(a)
        %%TRUNC Truncate an Interval. This only affects the endpoints of
        %      the interval; non-integer value are truncated to integers.
        %      This function rounds positive and negative numbers towards
        %      zero.            
        %
        %The trunc function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            a=Interval(fix(a.lb),fix(a.ub)); 
        end
        
        function a=fix(a)
        %%FIX Overloaded fix function in Matlab. This function just calls
        %     trunc. See the comments to trunc for more details.
        
            a=trunc(a);
        end

        function a=roundTiesToEven(a)
        %%ROUNDTIESTOEVEN Round an interval with ties being rounded to even
        %      integers. This only affects the endpoints of the interval;
        %      non-integer value are rounded to integers.
        %
        %The roundTiesToEven function is in Table 9.1 of Chapter 9.1 of
        %[1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.

            a=Interval(roundTies2EvenOdd(a.lb,0),roundTies2EvenOdd(a.ub,0));
        end
        
        function a=roundTiesToAway(a)
        %%ROUNDTIESTOAWAY Round an interval with ties being rounded away
        %      from zero. This only affects the endpoints of the interval;
        %      non-integer value are rounded to integers.
        %
        %The roundTiesToAway function is in Table 9.1 of Chapter 9.1 of
        %[1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
            
            a=Interval(round(a.lb),round(a.ub)); 
        end

        function a=round(a,N)
        %%ROUND Overloaded round function in Matlab. This rounds an
        %       interval with ties being rounded away from zero.
        %
        %INPUTS: a The interval(s) to round.
        %        N An optional parameter specifying the number of digits
        %          around the decimal point to which the number should be
        %          rounded. 0 means round to an integer (the default if
        %          omitted). A Positive number means round to N digits
        %          right of the decimal point and a negative number to N
        %          digits left of the decimal point.
        
            if(nargin<2)
                a=Interval(round(a.lb),round(a.ub)); 
            else
                a=Interval(round(a.lb,N),round(a.ub,N)); 
            end
        end

        function a=abs(a)
        %%ABS Evaluate the absolute value of an Interval. This makes all
        %     points in the interval >=0.
        %
        %The absolute value function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %The method for computing the tangent of an interval is taken from
        %Section 3 of [2]. 
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            sel=a.lb*a.ub>0;
            
            a.lb(sel)=min(abs(a.lb(sel)),abs(a.ub(sel)));
            a.ub=max(abs(a.lb),abs(a.ub));
            a.lb(~sel)=0;
        end
        
        function a=min(a,b,dim)
        %%MIN Find the minimum elements in an array of intervals or when
        %     comparing two intervals.
        %
        %The funcation can be called in three ways, the same as the min
        %function in Matlab:
        %1)
        %c=min(a)
        %In this case the min function is applied separetely to the upper
        %and lower bounds of a. Thus, if a holds a vector of intervals,
        %then the minimum upper and lower bounds of the vector are taken.
        %If a is a matrix of intervals, then the minimum across a column is
        %found. For multidimensional arrays, the minimum is across the
        %first non-singleton dimension.
        %2)
        %c=min(a,b)
        %In this case, the minimum is taken for the upper and lower bounds
        %for all intervals in a and b.
        %3)
        %c=min(a,[],dims) This is the same as min(a), except dims is an
        %integer specifying the dimension across which the minimum should
        %be found.
        %
        %The minimum function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.

            if(~isa(a,'Interval'))
               a=Interval(a); 
            end
            
            if(nargin==1)
               a.lb=min(a.lb);
               a.ub=min(a.ub);
            elseif(nargin==2)
                if(~isa(b,'Interval'))
                   b=Interval(b); 
                end
                a.lb=min(a.lb,b.lb);
                a.ub=min(a.ub,b.ub);
            else%nargin=3
                if(~builtin('isempty',b))
                    error('The second argument must be an empty matrix when using three input arguments')
                end
                
                a.lb=min(a.lb,[],dim);
                a.ub=min(a.ub,[],dim);
            end
            %If the multiplication ceated a NaN for one bound, make sure
            %there is a NaN for the whole bound (an empty interval).
            sel=isnan(a.lb)|isnan(a.ub);
            a.lb(sel)=NaN;
            a.ub(sel)=NaN;
        end
        
        function a=max(a,b,dim)
        %%MAX Find the maximum elements in an array of intervals or when
        %     comparing two intervals.
        %
        %The funcation can be called in three ways, the same as the max
        %function in Matlab:
        %1)
        %c=max(a)
        %In this case the max function is applied separetely to the upper
        %and lower bounds of a. Thus, if a holds a vector of intervals,
        %then the maximum upper and lower bounds of the vector are taken.
        %If a is a matrix of intervals, then the maximum across a column is
        %found. For multidimensional arrays, the maximum is across the
        %first non-singleton dimension.
        %2)
        %c=max(a,b)
        %In this case, the maximum is taken for the upper and lower bounds
        %for all intervals in a and b.
        %3)
        %c=max(a,[],dims) This is the same as max(a), except dims is an
        %integer specifying the dimension across which the maximum should
        %be found.
        %
        %The maximum function is in Table 9.1 of Chapter 9.1 of [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.

            if(~isa(a,'Interval'))
               a=Interval(a); 
            end
            
            if(nargin==1)
               a.lb=max(a.lb);
               a.ub=max(a.ub);
            elseif(nargin==2)
                if(~isa(b,'Interval'))
                   b=Interval(b); 
                end
                a.lb=max(a.lb,b.lb);
                a.ub=max(a.ub,b.ub);
            else%nargin=3
                if(~builtin('isempty',b))
                    error('The second argument must be an empty matrix when using three input arguments')
                end
                
                a.lb=max(a.lb,[],dim);
                a.ub=max(a.ub,[],dim);
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Forward Mode Elementary Functions of Chapter 10.6.1  of  %%%
%%% the IEEE Std 1788-2015                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function a=expm1(a)
        %%EXPM1 Evaluate the function exp(x)-1 without a loss of precision
        %       near x=0.
        %
        %This function is defined Table 10.5 of Chapter 10.6.1 of [1].
        %
        %Becuase the function is monotonically increasing, we can get
        %bounds on it in a similar manner to how the inverse tangent is
        %handled in Secion III of [2]. However, for getting the lower and
        %upper bounds of expm1 of a scalar, the values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            setProcRoundingMode(0);%Round to -Inf
            a.lb=expm1(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=expm1(a.ub);
        end
        
        function a=exp2m1(a)
        %%EXP2M1 Evaluate the function 2^x-1 without a loss of precision
        %       near x=0.
        %
        %This function is defined Table 10.5 of Chapter 10.6.1 of [1].
        %
        %Becuase the function is monotonically increasing, we can get
        %bounds on it in a similar manner to how the inverse tangent is
        %handled in Secion III of [2]. However, for getting the lower and
        %upper bounds of exp2m1 of a scalar, the values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            setProcRoundingMode(0);%Round to -Inf
            a.lb=exp2m1(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=exp2m1(a.ub);
        end

        function a=exp10m1(a)
        %%EXP10M1 Evaluate the function 10^x-1 without a loss of precision
        %       near x=0.
        %
        %This function is defined Table 10.5 of Chapter 10.6.1 of [1].
        %
        %Becuase the function is monotonically increasing, we can get
        %bounds on it in a similar manner to how the inverse tangent is
        %handled in Secion III of [2]. However, for getting the lower and
        %upper bounds of exp2m1 of a scalar, the values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            setProcRoundingMode(0);%Round to -Inf
            a.lb=exp10m1(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=exp10m1(a.ub);
        end

        function a=logp1(a)
        %%LOGP1 Evaluate the function log(x+1) without a loss of precision
        %       near x=0. Values of a passing below -1 cause an empty
        %       Interval to be returned.
        %
        %This function is defined Table 10.5 of Chapter 10.6.1 of [1].
        %
        %Becuase the function is monotonically increasing, we can get
        %bounds on it in a similar manner to how the inverse tangent is
        %handled in Secion III of [2]. However, for getting the lower and
        %upper bounds of log1p of a scalar, the values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            %Make intervals outside the valid range empty.
            sel=(a.lb<0);
            a.lb(sel)=NaN;
            a.ub(sel)=NaN;
        
            setProcRoundingMode(0);%Round to -Inf
            a.lb=log1p(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=log1p(a.ub);
        end
        
        function a=log1p(a)
        %%LOG1P Overloaded Matlab function logp1 for evaluating log(1+p).
        %       This just calls log1p.
        
            a=logp1(a);
        end
    
        function a=log2p1(a)
        %%LOG2P1 Evaluate the function log2(x+1) without a loss of
        %        precision near x=0. Values of a passing below -1 cause an
        %        empty Interval to be returned.
        %
        %This function is defined Table 10.5 of Chapter 10.6.1 of [1].
        %
        %Becuase the function is monotonically increasing, we can get
        %bounds on it in a similar manner to how the inverse tangent is
        %handled in Secion III of [2]. However, for getting the lower and
        %upper bounds of log2p1 of a scalar, the values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            %Make intervals outside the valid range empty.
            sel=(a.lb<0);
            a.lb(sel)=NaN;
            a.ub(sel)=NaN;
        
            setProcRoundingMode(0);%Round to -Inf
            a.lb=log2p1(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=log2p1(a.ub);
        end
        
        function a=log10p1(a)
        %%LOG10P1 Evaluate the function log10(x+1) without a loss of
        %        precision near x=0. Values of a passing below -1 cause an
        %        empty Interval to be returned.
        %
        %This function is defined Table 10.5 of Chapter 10.6.1 of [1].
        %
        %Becuase the function is monotonically increasing, we can get
        %bounds on it in a similar manner to how the inverse tangent is
        %handled in Secion III of [2]. However, for getting the lower and
        %upper bounds of log10p1 of a scalar, the values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            %Make intervals outside the valid range empty.
            sel=(a.lb<0);
            a.lb(sel)=NaN;
            a.ub(sel)=NaN;
        
            setProcRoundingMode(0);%Round to -Inf
            a.lb=log10p1(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=log10p1(a.ub);
        end
        
        function a=hypot(x,y)
        %%HYPOT Evaluate the hypotenuse function sqrt(x^2+y^2) for all
        %       corresponding intervals in x and y. This also overloads the
        %       built-in Matlab function.
        %
        %This function is defined Table 10.5 of Chapter 10.6.1 of [1].
        %
        %Since multiplication and square roots are already defined in the
        %Interval class, nothing special needs to be done for this
        %implementation.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            a=sqrt(x.*x+y.*y);
        end
        
        function b=sinPi(a)
        %%SINPI Evaluate sin(pi*a). This takes into account the fact that
        %       with finite precision, the value of pi is a range, not a
        %       point. Thus, this just calls the sin function after finding
        %       the biggest range of inputs given the range of pi.
        %
        %This function is defined Table 10.5 of Chapter 10.6.1 of [1].
        %
        %Since multiplication and square roots are already defined in the
        %Interval class, nothing special needs to be done for this
        %implementation.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            piInt=Interval.pi();
            
            %For a positive lower bound, we want to use the lowest value of
            %pi to get the lowest lower bound on the input to sin(pi*a).
            sel=a.lb>0;
            setProcRoundingMode(0);%Round to -Inf
            a.lb(sel)=piInt.lb*a.lb(sel);
            %For a negatvce lower bound, we want to use the highest value
            %of pi to get the lowest lower bound on the input to sin(pi*a).
            a.lb(~sel)=piInt.ub*a.lb(~sel);
            
            %For a positive upper bound, we want to use the highest value
            %of pi to get the highest upper bound on the input to
            %sin(pi*a).
            sel=a.ub>0;
            setProcRoundingMode(3);%Round to +Inf
            a.ub(sel)=piInt.ub.*a.ub(sel);
            
            %For a negative upper bound, we want to use the lowest value
            %of pi to get the highest upper bound on the input to
            %sin(pi*a).
            a.ub(~sel)=piInt.lb.*a.ub(~sel);
            
            b=sin(a);
        end
        
        function b=cosPi(a)
        %%COSPI Evaluate cos(pi*a). This takes into account the fact that
        %       with finite precision, the value of pi is a range, not a
        %       point. Thus, this just calls the cos function after finding
        %       the biggest range of inputs given the range of pi.
        %
        %This function is defined Table 10.5 of Chapter 10.6.1 of [1].
        %
        %Since multiplication and square roots are already defined in the
        %Interval class, nothing special needs to be done for this
        %implementation.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            piInt=Interval.pi();
            
            %For a positive lower bound, we want to use the lowest value of
            %pi to get the lowest lower bound on the input to cos(pi*a).
            sel=a.lb>0;
            setProcRoundingMode(0);%Round to -Inf
            a.lb(sel)=piInt.lb*a.lb(sel);
            %For a negatvce lower bound, we want to use the highest value
            %of pi to get the lowest lower bound on the input to cos(pi*a).
            a.lb(~sel)=piInt.ub*a.lb(~sel);
            
            %For a positive upper bound, we want to use the highest value
            %of pi to get the highest upper bound on the input to
            %cos(pi*a).
            sel=a.ub>0;
            setProcRoundingMode(3);%Round to +Inf
            a.ub(sel)=piInt.ub.*a.ub(sel);
            
            %For a negative upper bound, we want to use the lowest value
            %of pi to get the highest upper bound on the input to
            %cos(pi*a).
            a.ub(~sel)=piInt.lb.*a.ub(~sel);
            
            b=cos(a);
        end
        
        function b=tanPi(a)
        %%TANPI Evaluate tan(pi*a). This takes into account the fact that
        %       with finite precision, the value of pi is a range, not a
        %       point. Thus, this just calls the tan function after finding
        %       the biggest range of inputs given the range of pi.
        %
        %This function is defined Table 10.5 of Chapter 10.6.1 of [1].
        %
        %Since multiplication and square roots are already defined in the
        %Interval class, nothing special needs to be done for this
        %implementation.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            piInt=Interval.pi();
            
            %For a positive lower bound, we want to use the lowest value of
            %pi to get the lowest lower bound on the input to tan(pi*a).
            sel=a.lb>0;
            setProcRoundingMode(0);%Round to -Inf
            a.lb(sel)=piInt.lb*a.lb(sel);
            %For a negatvce lower bound, we want to use the highest value
            %of pi to get the lowest lower bound on the input to tan(pi*a).
            a.lb(~sel)=piInt.ub*a.lb(~sel);
            
            %For a positive upper bound, we want to use the highest value
            %of pi to get the highest upper bound on the input to
            %tan(pi*a).
            sel=a.ub>0;
            setProcRoundingMode(3);%Round to +Inf
            a.ub(sel)=piInt.ub.*a.ub(sel);
            
            %For a negative upper bound, we want to use the lowest value
            %of pi to get the highest upper bound on the input to
            %tan(pi*a).
            a.ub(~sel)=piInt.lb.*a.ub(~sel);
            
            b=tan(a);
        end
        
        function a=asinPi(a)
        %%ASINPI Evaluate asin(a)/pi. This takes into account the fact that
        %       with finite precision, the value of pi is a range, not a
        %       point.
        %
        %This function is defined Table 10.5 of Chapter 10.6.1 of [1].
        %
        %Since multiplication and square roots are already defined in the
        %Interval class, nothing special needs to be done for this
        %implementation.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            a=asin(a);

            piInt=Interval.pi();
            
            %For a positive lower bound, we want to use the largest value of
            %pi to get the lowest lower bound on the input to asin(a)/pi.
            sel=a.lb>0;
            setProcRoundingMode(0);%Round to -Inf
            a.lb(sel)=a.lb(sel)/piInt.ub;
            %For a negatice lower bound, we want to use the lowest value
            %of pi to get the lowest lower bound on the input to 
            %asin(a)/pi.
            a.lb(~sel)=a.lb(~sel)/piInt.lb;
            
            %For a positive upper bound, we want to use the lowest value
            %of pi to get the highest upper bound on the input to
            %asin(a)/pi.
            sel=a.ub>0;
            setProcRoundingMode(3);%Round to +Inf
            a.ub(sel)=a.ub(sel)/piInt.lb;
            
            %For a negative upper bound, we want to use the highest value
            %of pi to get the highest upper bound on the input to
            %asin(a)/pi.
            a.ub(~sel)=a.ub(~sel)/piInt.ub;
        end
        
        function a=acosPi(a)
        %%ACOSPI Evaluate acos(a)/pi. This takes into account the fact that
        %       with finite precision, the value of pi is a range, not a
        %       point.
        %
        %This function is defined Table 10.5 of Chapter 10.6.1 of [1].
        %
        %Since multiplication and square roots are already defined in the
        %Interval class, nothing special needs to be done for this
        %implementation.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            a=acos(a);

            piInt=Interval.pi();
            
            %For a positive lower bound, we want to use the largest value of
            %pi to get the lowest lower bound on the input to acos(a)/pi.
            sel=a.lb>0;
            setProcRoundingMode(0);%Round to -Inf
            a.lb(sel)=a.lb(sel)/piInt.ub;
            %For a negatice lower bound, we want to use the lowest value
            %of pi to get the lowest lower bound on the input to 
            %acos(a)/pi.
            a.lb(~sel)=a.lb(~sel)/piInt.lb;
            
            %For a positive upper bound, we want to use the lowest value
            %of pi to get the highest upper bound on the input to
            %acos(a)/pi.
            sel=a.ub>0;
            setProcRoundingMode(3);%Round to +Inf
            a.ub(sel)=a.ub(sel)/piInt.lb;
            
            %For a negative upper bound, we want to use the highest value
            %of pi to get the highest upper bound on the input to
            %acos(a)/pi.
            a.ub(~sel)=a.ub(~sel)/piInt.ub;
        end
        
        function a=atanPi(a)
        %%ATANPI Evaluate atan(a)/pi. This takes into account the fact that
        %       with finite precision, the value of pi is a range, not a
        %       point.
        %
        %This function is defined Table 10.5 of Chapter 10.6.1 of [1].
        %
        %Since multiplication and square roots are already defined in the
        %Interval class, nothing special needs to be done for this
        %implementation.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            a=atan(a);

            piInt=Interval.pi();
            
            %For a positive lower bound, we want to use the largest value of
            %pi to get the lowest lower bound on the input to atan(a)/pi.
            sel=a.lb>0;
            setProcRoundingMode(0);%Round to -Inf
            a.lb(sel)=a.lb(sel)/piInt.ub;
            %For a negatice lower bound, we want to use the lowest value
            %of pi to get the lowest lower bound on the input to 
            %atan(a)/pi.
            a.lb(~sel)=a.lb(~sel)/piInt.lb;
            
            %For a positive upper bound, we want to use the lowest value
            %of pi to get the highest upper bound on the input to
            %atan(a)/pi.
            sel=a.ub>0;
            setProcRoundingMode(3);%Round to +Inf
            a.ub(sel)=a.ub(sel)/piInt.lb;
            
            %For a negative upper bound, we want to use the highest value
            %of pi to get the highest upper bound on the input to
            %atan(a)/pi.
            a.ub(~sel)=a.ub(~sel)/piInt.ub;
        end
        
        function a=atan2Pi(a,b)
        %%ATAN2PI Evaluate atan(a,b)/pi. This takes into account the fact
        %       that with finite precision, the value of pi is a range, not
        %       a point.
        %
        %This function is defined Table 10.5 of Chapter 10.6.1 of [1].
        %
        %Since multiplication and square roots are already defined in the
        %Interval class, nothing special needs to be done for this
        %implementation.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            a=atan2(a,b);

            piInt=Interval.pi();
            
            %For a positive lower bound, we want to use the largest value of
            %pi to get the lowest lower bound on the input to atan(a)/pi.
            sel=a.lb>0;
            setProcRoundingMode(0);%Round to -Inf
            a.lb(sel)=a.lb(sel)/piInt.ub;
            %For a negatice lower bound, we want to use the lowest value
            %of pi to get the lowest lower bound on the input to 
            %atan(a)/pi.
            a.lb(~sel)=a.lb(~sel)/piInt.lb;
            
            %For a positive upper bound, we want to use the lowest value
            %of pi to get the highest upper bound on the input to
            %atan(a)/pi.
            sel=a.ub>0;
            setProcRoundingMode(3);%Round to +Inf
            a.ub(sel)=a.ub(sel)/piInt.lb;
            
            %For a negative upper bound, we want to use the highest value
            %of pi to get the highest upper bound on the input to
            %atan(a)/pi.
            a.ub(~sel)=a.ub(~sel)/piInt.ub;
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Slope Functions of Chapter 10.6.2 of the IEEE Std        %%%
%%% 1788-2015                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function a=expSlope1(a)
        %%EXPSLOPE1 Evaluate the function (exp(x)-1)/x. This function plays
        %           a role in slope algorithms for implementing improved
        %           range enclosures. Care is taken at the x=0 singularity,
        %           which equals 1 by using the expSlope1 function (written
        %           for non-intervals).
        %
        %This function is defined Table 10.6 of Chapter 10.6.2 of [1].
        %
        %Becuase the function is monotonically increasing, we can get
        %bounds on it in a similar manner to how the inverse tangent is
        %handled in Secion III of [2]. Proper rounding for the operations
        %has been added. However, for getting the lower and upper bounds of
        %the expSlope1 function of a scalar, the values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            setProcRoundingMode(0);%Round to -Inf
            a.lb=expSlope1(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=expSlope1(a.ub);
        end

        function a=expSlope2(a)
        %%EXPSLOPE2 Evaluate the function 2*(exp(x)-1-x)/x^2. This function
        %           plays a role in slope algorithms for implementing
        %           improved range enclosures. Care is taken at the x=0
        %           singularity, which equals 1 by using the expSlope2
        %           function (written for non-intervals).
        %
        %This function is defined Table 10.6 of Chapter 10.6.2 of [1].
        %
        %Becuase the function is monotonically increasing, we can get
        %bounds on it in a similar manner to how the inverse tangent is
        %handled in Secion III of [2]. Proper rounding for the operations
        %has been added. However, for getting the lower and upper bounds of
        %the expSlope2 function of a scalar, the values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            setProcRoundingMode(0);%Round to -Inf
            a.lb=expSlope2(a.lb);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=expSlope2(a.ub);
        end
        
        function a=logSlope1(a)
        %%LOGSLOPE1 Evaluate the function -2*(log(1+x)-x)/x^2. This function
        %           plays a role in slope algorithms for implementing
        %           improved range enclosures. Care is taken at the x=0
        %           singularity, which equals 1 by using the logSlope1
        %           function (written for non-intervals). Intervals with
        %           bounds <-1 cause an empty interval to be returned.
        %
        %This function is defined Table 10.6 of Chapter 10.6.2 of [1].
        %
        %Becuase the function is monotonically decreasing, we can get
        %bounds on it in a similar manner to how the inverse tangent is
        %handled in Secion III of [2]. Proper rounding for the operations
        %has been added. However, for getting the lower and upper bounds of
        %the logSlope1 function of a scalar, the values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            %Make intervals outside the valid range empty.
            sel=(a.lb<-1);
            a.lb(sel)=NaN;
            a.ub(sel)=NaN;
        
            setProcRoundingMode(0);%Round to -Inf
            temp=a.lb;
            a.lb=logSlope1(a.ub);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=logSlope1(temp);
        end
        
        function a=logSlope2(a)
        %%LOGSLOPE2 Evaluate the function 3*(log(1+x)-x+x^2/2)/x^3. This
        %           function plays a role in slope algorithms for
        %           implementing improved range enclosures. Care is taken
        %           at the x=0 singularity, which equals 1 by using the
        %           logSlope2 function (written for non-intervals).
        %           Intervals with bounds <-1 cause an empty interval to be
        %           returned.
        %
        %This function is defined Table 10.6 of Chapter 10.6.2 of [1].
        %
        %Becuase the function is monotonically decreasing, we can get
        %bounds on it in a similar manner to how the inverse tangent is
        %handled in Secion III of [2]. Proper rounding for the operations
        %has been added. However, for getting the lower and upper bounds of
        %the logSlope2 function of a scalar, the values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            %Make intervals outside the valid range empty.
            sel=(a.lb<-1);
            a.lb(sel)=NaN;
            a.ub(sel)=NaN;
        
            setProcRoundingMode(0);%Round to -Inf
            temp=a.lb;
            a.lb=logSlope2(a.ub);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=logSlope2(temp);
        end        
        
        function a=sinSlope3(a)
        %%SINSLOPE3 Evaluate the function -6*(sin(x)-x)/x^3. This
        %           function plays a role in slope algorithms for
        %           implementing improved range enclosures. Care is taken
        %           at the x=0 singularity, which equals 1 by using the
        %           sinSlope3 function (written for non-intervals).
        %
        %This function is defined Table 10.6 of Chapter 10.6.2 of [1].
        %
        %Becuase the function is symmetric about 0 and monotonically
        %decreasing with respect to the absolute value of the input, we can
        %get bounds on it in a similar manner to how the inverse tangent is
        %handled in Secion III of [2]. Proper rounding for the operations
        %has been added. However, for getting the lower and upper bounds of
        %the sinSlope3 function of a scalar, the values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
            
            minVal=min(abs(a.ub),abs(a.lb));
            maxVal=max(abs(a.ub),abs(a.lb));
            
            setProcRoundingMode(0);%Round to -Inf
            a.lb=sinSlope3(maxVal);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=sinSlope3(minVal);
        end   

        function a=asinSlope3(a)
        %%ASINSLOPE3 Evaluate the function 6*(asin(x)-x)/x^3. This
        %           function plays a role in slope algorithms for
        %           implementing improved range enclosures. Care is taken
        %           at the x=0 singularity, which equals 1 by using the
        %           asinSlope3 function (written for non-intervals).
        %           Intervals leaving the -1 to 1 valid range of inputs
        %           cause empty intervals to be returned.
        %
        %This function is defined Table 10.6 of Chapter 10.6.2 of [1].
        %
        %Becuase the function is symmetric about 0 and monotonically
        %increasing with respect to the absolute value of the input, we can
        %get bounds on it in a similar manner to how the inverse tangent is
        %handled in Secion III of [2]. Proper rounding for the operations
        %has been added. However, for getting the lower and upper bounds of
        %the asinSlope3 function of a scalar, the values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
            
        %Make intervals outside the valid range empty.
            sel=a.lb<-1 | a.ub> 1;
            a.lb(sel)=NaN;
            a.ub(sel)=NaN;
        
            minVal=min(abs(a.ub),abs(a.lb));
            maxVal=max(abs(a.ub),abs(a.lb));
            
            setProcRoundingMode(0);%Round to -Inf
            a.lb=asinSlope3(minVal);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=asinSlope3(maxVal);
        end
        
        function a=atanSlope3(a)
        %%ATANSLOPE3 Evaluate the function -3*(atan(x)-x)/x^3. This
        %           function plays a role in slope algorithms for
        %           implementing improved range enclosures. Care is taken
        %           at the x=0 singularity, which equals 1 by using the
        %           atanSlope3 function (written for non-intervals).
        %
        %This function is defined Table 10.6 of Chapter 10.6.2 of [1].
        %
        %Becuase the function is symmetric about 0 and monotonically
        %decreasing with respect to the absolute value of the input, we can
        %get bounds on it in a similar manner to how the inverse tangent is
        %handled in Secion III of [2]. Proper rounding for the operations
        %has been added. However, for getting the lower and upper bounds of
        %the atanSlope3 function of a scalar, the values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
            
            minVal=min(abs(a.ub),abs(a.lb));
            maxVal=max(abs(a.ub),abs(a.lb));
            
            setProcRoundingMode(0);%Round to -Inf
            a.lb=atanSlope3(maxVal);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=atanSlope3(minVal);
        end
        
        function a=coshSlope2(a)
        %%COSHSLOPE2 Evaluate the function 2*(cosh(x)-1)/x^2. This
        %            function plays a role in slope algorithms for
        %            implementing improved range enclosures. Care is taken
        %            at the x=0 singularity, which equals 1 by using the
        %            coshSlope2 function (written for non-intervals).
        %
        %This function is defined Table 10.6 of Chapter 10.6.2 of [1].
        %
        %Becuase the function is symmetric about 0 and monotonically
        %increasing with respect to the absolute value of the input, we can
        %get bounds on it in a similar manner to how the inverse tangent is
        %handled in Secion III of [2]. Proper rounding for the operations
        %has been added. However, for getting the lower and upper bounds of
        %the coshSlope2 function of a scalar, the values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
            
            minVal=min(abs(a.ub),abs(a.lb));
            maxVal=max(abs(a.ub),abs(a.lb));
            
            setProcRoundingMode(0);%Round to -Inf
            a.lb=coshSlope2(minVal);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=coshSlope2(maxVal);
        end
        
        function a=sinhSlope3(a)
        %%SINHSLOPE2 Evaluate the function 6*(sinh(x)-x)/x^3. This
        %            function plays a role in slope algorithms for
        %            implementing improved range enclosures. Care is taken
        %            at the x=0 singularity, which equals 1 by using the
        %            sinhSlope2 function (written for non-intervals).
        %
        %This function is defined Table 10.6 of Chapter 10.6.2 of [1].
        %
        %This is two times one of the slope function in Table 10.6 of
        %Section 10.6.2 of [1]. As the standard says the value at 0 is 1,
        %but the value of the function in the standard at 0 is 1/2, it is
        %assumed that there is a typo in the standard. Hence this function
        %is twice that listed in [1].
        %
        %Becuase the function is symmetric about 0 and monotonically
        %increasing with respect to the absolute value of the input, we can
        %get bounds on it in a similar manner to how the inverse tangent is
        %handled in Secion III of [2]. Proper rounding for the operations
        %has been added. However, for getting the lower and upper bounds of
        %the sinhSlope3 function of a scalar, the values are computed with
        %different rounding modes rather than using a different number of
        %terms in a series expansion.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        %[2] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
            
            minVal=min(abs(a.ub),abs(a.lb));
            maxVal=max(abs(a.ub),abs(a.lb));
            
            setProcRoundingMode(0);%Round to -Inf
            a.lb=sinhSlope3(minVal);
            setProcRoundingMode(3);%Round to +Inf
            a.ub=sinhSlope3(maxVal);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Boolean functions of intervals from Chapter 10.6.3  of   %%%
%%% the IEEE Std 1788-2015                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function c=isSingleton(a)
        %%ISSINGLETON This function returns true if the interval represents
        %             a single point. Empty intervals are not considered
        %             singleton
        %
        %This function is defined Chapter 10.6.3 of [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            c=(a.lb==a.ub);
        end
        
        
        function c=isMember(a,val)
        %%ISMEMBER Determine whether a value is within the interval. This
        %          return a boolean value. For sets of intervals, this doe
        %          not AND-together all of the results.
        %
        %This function is defined Chapter 10.6.3 of [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
            
            c=a.lb<=val & a.ub>=val;
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The 16 state of overlapping interval situations of Table %%%
%%% 10.7 of Chapter 10.6.4 of the IEEE Std 1788-2015         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function c=bothEmpty(a,b)
        %%BOTHEMPTY Given Intervals a and b, this function returns true if
        %           both intervals are empty. An empty Interval is one
        %           whose upper and lower bounds are both NaNs. On the
        %           other hand, if the upper and lower bound arrays are
        %           empty matrices, this function returns an empty matrix.
        %
        %This function is defined in Table 10.7 of Chapter 10.6.4 [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
            
            c=isnan(a.lb)&isnan(a.ub)&isnan(b.lb)&isnan(b.ub);
        end
        
        function c=firstEmpty(a,b)
        %%FIRSTEMPTY Given Intervals a and b, this function returns true if
        %           Interval a is empty and Interval b is not empty. An
        %           empty Interval is one whose upper and lower bounds are
        %           both NaNs. On the other hand, if the upper and lower
        %           bound arrays are empty matrices, this function returns
        %           an empty matrix.
        %
        %This function is defined in Table 10.7 of Chapter 10.6.4 [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
            
            c=(isnan(a.lb)&isnan(a.ub))&~(isnan(b.lb)&isnan(b.ub));
        end
        
        function c=secondEmpty(a,b)
        %%SECONDEMPTY Given Intervals a and b, this function returns true
        %           if Interval a is not empty and Interval b is empty. An
        %           empty Interval is one whose upper and lower bounds are
        %           both NaNs. On the other hand, if the upper and lower
        %           bound arrays are empty matrices, this function returns
        %           an empty matrix.
        %
        %This function is defined in Table 10.7 of Chapter 10.6.4 [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
            
            c=~(isnan(a.lb)&isnan(a.ub))&(isnan(b.lb)&isnan(b.ub));
        end
        
        function c=before(a,b)
        %%BEFORE Returns c(i)=true if Interval a comes before Interval b.
        %        This means that at no point do they overlap and a.ub<b.lb.
        %        Input arguments can also be points.
        %
        %This function is defined in Table 10.7 of Chapter 10.6.4 [1]. To
        %be consistent with the definition of handling empty intervals in
        %the precedes function in Table 10.5, this function returns true if
        %any input is an empty interval.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end
            c=a.ub<b.lb| isnan(a.lb) | isnan(b.lb);
        end
        
        function c=lt(a,b)
        %%LT Overloaded strictly less than < comparison This is the same as
        %    the before function. See the before function for more details.
        
            c=before(a,b);
        end
        
        function c=meets(a,b)
        %%MEETS  Returns c(i)=true if Interval a meets Interval b. This
        %        means that the upper bound of a equals the lower
        %        bound of b, and that neither a not b is a single point.
        %        Input arguments can also be points.
        %
        %This function is defined in Table 10.7 of Chapter 10.6.4 [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end
            c=(a.lb<a.ub)&(a.ub==b.lb)&(b.lb<b.ub);
        end

        function c=overlaps(a,b)
        %%OVERLAPS Returns c(i)=true if Interval a overlaps Interval b.
        %        This means that a.lb<b.lb, b.lb <a.ub and a.ub<b.ub. Thus,
        %        this is a specific type of overlapping. Input arguments
        %        can also be points.
        %
        %This function is defined in Table 10.7 of Chapter 10.6.4 [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end
        
            c=(a.lb<b.lb)&(b.lb<a.ub)&(a.ub<b.ub);
        end
        
        function c=starts(a,b)
        %%STARTS Returns c(i)=true if Interval a starts Interval b. That is
        %        if all of a is contained in b and if the lower bound of a
        %        equals the lower bound of b. The upper bound of a must be
        %        strictly less than the upper bound of b. Input arguments
        %        can also be points.
        %
        %This function is defined in Table 10.7 of Chapter 10.6.4 [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end
        
            c = (a.lb==b.lb) & (a.ub<=b.ub);
        end
        
        function c=containedBy(a,b)
        %%CONTAINEDBY Returns c(i)=true if Interval a is entirely contained
        %             by Interval b without overlap at the endpoints of b.
        %             Inputs can also be points.
        %
        %This function is defined in Table 10.7 of Chapter 10.6.4 [1].
        %Compare to the interior function. In order to be consistent with
        %the comparisons with empty intervals in Table 10.4 (with respect
        %to the interior function), two empty intervals return true and if
        %a is empty and b is not, then this returns true. Also, comparisons
        %between infinite values of the same sign are true.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015
                
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end
            
            c = ((a.lb>b.lb)|(a.lb==-Inf & b.lb==-Inf)|(a.lb==Inf & b.lb==Inf))&((a.ub<b.ub)|(a.ub==-Inf &b.ub==-Inf)|(a.ub==Inf & b.ub==Inf))|isnan(a.lb);
        end
        
        function c=finishes(a,b)
        %%FINISHES Returns c(i)=true if Interval a finishes Interval b.
        %          This means that all of interval a is contained in b
        %          without overlapping at the lower bound, and the upper
        %          bound of a coincides with the upper bound of b. inputs
        %          can also be points.
        %
        %This function is defined in Table 10.7 of Chapter 10.6.4 [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end
            c=(a.ub==b.ub) & (a.lb>b.lb);
        end
        
        function c=equals(a,b)
        %%EQUALS Returns c(i)=true if Interval a is equal to Interval b.
        %        This means that the bounds of both intervals are equal.
        %        Inputs can also be points.
        %
        %This function is defined in Table 10.7 of Chapter 10.6.4 [1].
        %Table 10.4 of Chapter 10.4.10 specifies that two empty intervals
        %must be considered equal, so that is the case here.
        %
        %Compare to the equal function.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end
            
            c=(a.ub==b.ub)&(a.lb==b.lb) |(isnan(a.lb)&isnan(a.ub));
        end
        
        function c=eq(a,b)
        %%EQ Overloaded == operator in Matlab. This is the same as the
        %    equals method. See the comments to that method for more
        %    information.
        
            c=equals(a,b);
        end
        
        function c=finishedBy(a,b)
        %%FINISHEDBY Returns c(i)=true if Interval a is finished by
        %          Interval b. This means that all of interval b is
        %          contained in a and the upper bound of b coincides with
        %          the upper bound of a without overlapping the lower
        %          bound. Inputs can also be points.
        %
        %This function is defined in 10.7 of Chapter 10.6.4 [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end
            
            c=(a.ub==b.ub)&(a.lb<b.lb);
        end

        function c=contains(a,b)
        %%CONTAINS Returns c(i)=true if Interval a contains Interval b.
        %          This means that all of interval b is within a and does
        %          not overlap with the upper or lower bounds of a.
        %          Inputs can also be points.
        %
        %This function is defined in 10.7 of Chapter 10.6.4 [1]. In order
        %to be consistent with the comparisons with empty intervals in
        %Table 10.4 (with respect to the interior function), two empty
        %intervals return true and if a is not empty and b is empty, then
        %this returns true. Also, comparisons between two infinite
        %quantities of the same sign are true.
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end

            c=((a.lb<b.lb)|(a.lb==-Inf & b.lb==-Inf)|(a.lb==Inf & b.lb==Inf))&((b.ub<a.ub)|(a.ub==-Inf &b.ub==-Inf)|(a.ub==Inf & b.ub==Inf))|isnan(a.lb);
        end
        
        function c=startedBy(a,b)
        %%STARTEDBY Returns c(i)=true if Interval a is started by Interval
        %           b. This means that the lower bound of b equals the
        %           lower bound of a and the upper bound of b is less than
        %           the upper bound of a. Inputs can also be points.
        %
        %This function is defined in 10.7 of Chapter 10.6.4 [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
   
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end
            
            c=(a.lb==b.lb)&(b.ub<a.ub);
        end
        
        function c=overlappedBy(a,b)
        %%OVERLAPPEDBY Returns c(i)=true if Interval a is overlapped by
        %          Interval b. This means that the lower bound of b is less
        %          than the lower bound of a, the upper bound of b is
        %          greater than the lower bound of a and the upper bound of
        %          b is less than the upper bound of a. OInputs can also be
        %          points.
        %
        %This function is defined in 10.7 of Chapter 10.6.4 [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end
            
            c=(b.lb<a.lb)&(a.lb<b.ub)&(b.ub<a.ub);
        end        

        function c=metBy(a,b)
        %%METBY Returns c(i)=true if Interval a is met by Interval b. This
        %       means that neither a nor b can be a single point and the
        %       upper bound of b equals the lower bound of a. Inputs can
        %       also be points.
        %
        %This function is defined in 10.7 of Chapter 10.6.4 [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end
            
            c=(b.lb<b.ub)&(b.ub==a.lb)&(a.lb<a.ub);
        end
        
        function c=after(a,b)
        %%AFTER Returns c(i)=true if Interval a comes after Interval b.
        %       This means that the upper bound of b is less than the lower
        %       bound of a. Inputs can also be points.
        %
        %This function is defined in 10.7 of Chapter 10.6.4 [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            if(~isa(a,'Interval'))
                a=Interval(a);
            end
            if(~isa(b,'Interval'))
                b=Interval(b);
            end

            c=b.ub<a.lb;
        end
        
        function c=gt(a,b)
        %%GT Overloaded strictly greater than > comparison This is the same as
        %    the after function. See the after function for more details.
        
            c=after(a,b);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Misc Overloaded Matlab functions                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function C=mtimes(A,B)
        %%MTIMES Overloaded matrix multiplication operator. Multiply A*B. A
        %        and/ or B can be Intervals or doubles.
            
            C=matMultiply(A,B,'Interval.zeros');
        end

        function c=isfinite(a)
        %%ISFINITE Overloading the isfinite operator in Matlab. This will
        %          be true if either the upper or lower bounds are Int or
        %          if the interval is empty, meaning contains NaNs.

            c=isfinite(a.lb)&isfinite(a.ub);
        end
        
        function c=subsref(a,s)
        %%SUBSREF Overloaded subscript reference operation. When Interval a
        %         contains multiple objects, this function lets one use
        %         a(c:d) to access objects from c to d.
        %	
            switch s(1).type
                case '()'
                    lbNew=builtin('subsref',a.lb,s);
                    ubNew=builtin('subsref',a.ub,s);
                    c=Interval(lbNew,ubNew);
                otherwise
                    c=builtin('subsref',a,s);
            end
        end
        
        function a=subsasgn(a,s,b)
        %%SUBSASGN Overloaded subscript assignment operation.
            if(isempty(a))
                a=Interval();
            end
            switch(s(1).type)
                case '()'
                    if(isa(b,'Interval'))
                        a.lb=builtin('subsasgn',a.lb,s,b.lb);
                        a.ub=builtin('subsasgn',a.ub,s,b.ub);
                    else
                        a.lb=builtin('subsasgn',a.lb,s,b);
                        a.ub=builtin('subsasgn',a.ub,s,b);
                    end
                otherwise
                    builtin('subsasgn',a,s.b);
            end
        end
        
        function c=end(a,k,~)
        %%END Overloaded last index operation.
            c=size(a,k);
        end

        function c=sum(a,varargin)
        %%SUM Overloaded sum over a particular dimension. If any
        %Interval is empty (NaN), then the resulting sum is empty.
            setProcRoundingMode(0);%Round to -Inf
            lbNew=sum(a.lb,varargin{:});
            setProcRoundingMode(3);%Round to +Inf
            ubNew=sum(a.ub,varargin{:});
            c=Interval(lbNew,ubNew);
        end
        
        function c=prod(a,dim)
        %%PROD Overloaded product over a particular dimension. If any
        %Interval is empty (NaN), then the resulting product is empty.
            sz=size(a);
            if(nargin<2)
                [~,dim]=find(sz>1,1);
                if(isempty(dim))
                    dim=1;
                end
            end
            
            %Define full indices for dynamic indexing
            %Overloaded subsref not automatically called from within class
            s.type='()';
            s.subs=cell(1,length(sz));
            for n=1:length(sz)
                s.subs{n}=1:sz(n);
            end
            
            c=Interval(1,1);
            for m=1:sz(dim)
                s.subs{dim}=m;
                c=c.*subsref(a,s);
            end
        end
        
        function c=horzcat(varargin)
        %%HORZCAT Overloaded horizontal concatenation operation.
        
            lbNew=cell2mat(cellfun(@(x) Interval.infimumS(x),varargin,'UniformOutput',false));
            ubNew=cell2mat(cellfun(@(x) Interval.supremumS(x),varargin,'UniformOutput',false));
            c=Interval(lbNew,ubNew);
        end
        
        function c=vertcat(varargin)
        %%VERTCART Overloaded vertical concatenation operation.
        
            lbNew=cell2mat(cellfun(@(x) Interval.infimumS(x),varargin.','UniformOutput',false));
            ubNew=cell2mat(cellfun(@(x) Interval.supremumS(x),varargin.','UniformOutput',false));
            c=Interval(lbNew,ubNew);
        end
        
        function c=cat(dim,varargin)
        %%CAT Overloaded concatenation operation.
        
            lb_varargin=cellfun(@(x) Interval.infimumS(x),varargin,'UniformOutput',false);
            lbNew=cat(dim,lb_varargin{:});
            ub_varargin=cellfun(@(x) Interval.supremumS(x),varargin,'UniformOutput',false);
            ubNew=cat(dim,ub_varargin{:});
            c=Interval(lbNew,ubNew);
        end
        
        function c=transpose(a)
        %%TRANSPOSE Overloaded transpose operation(the .' operation, not
        %           the ' operation). This transposes the matrices holding
        %           the upper and lower bounds of the intervals.
            
            c=Interval(a.lb.',a.ub.');
        end
        
        function c=ctranspose(a)
        %%CTRANSPOSE Overloaded complex conjugate transpose operation (the
        %           ' operation, not the .' operation). This transposes the
        %           matrices holding the upper and lower bounds of the
        %           intervals. Intervals must be real, so this is the same
        %           as the transpose operator.
            
            c=Interval(a.lb',a.ub');
        end
        
        function varargout=size(a,dim)
        %%SIZE Overloaded size operation. This returns the size of the
        %      lower-bound matrix of the intervals, which is the same as
        %      the size of the upper bound matrix of the intervals.
        
            d=size(a.lb);
            if(nargout>1)
                varargout=num2cell(d);
            elseif(nargin>1)
                varargout={d(dim)};
            else
                varargout={d};
            end
        end
        
        function L=length(a)
        %%LENGTH Overloaded length operation. This returns the length
        %of the lower-bound matrix of the intervals, which is the same as
        %the size of the upper bound matrix of the intervals.
        
            L=length(a.lb);
        end
        
        function n=numel(a)
        %%NUMEL Overloaded element numel operation. This is the number of
        %intervals in a set. This just returns the length of the lower-
        %bound matrix of the intervals, which is the same as the length of
        %the upper bound matrix of the intervals. For versions of Matlab
        %prior to to 2015b, numel is called from subsref and subsasgn and
        %must be set to 1.
            
            v=ver;
            version=str2double(v(1).Version);
            %R2015b, numel no longer called for overloaded indexing
            if(version>=9)
                n=numel(a.lb);
            else
            %Otherwise must be set to 1 to not break overloaded subsref
            %indexing
                n=1;
            end
        end
        
        function n=numArgumentsFromSubscript(~,~,~)
        %%NUMARGUMENTSFROMSUBSCRIPT Overloaded method to allow proper
        %indexing from overloaded subsref and subsasgn. Requires Matlab
        %R2015b or later to work.
        
            n=1;
        end
        
        function plot(a,b)
        %%PLOT Overloaded plot function. This function plots intervals as
        %      blue lines or green boxes. If one input is provided, the
        %      function plots the intervals in a as vertical lines
        %      equispaced. If a and b are Intervals, then boxes are plotted
        %      in x-y cooridnates.
        %
        %Example:
        % a=[Interval(1,2);Interval(3,4);Interval(5,6)];
        % b=[Interval(1,2);Interval(4,5);Interval(1,3)];
        % plot(a,b)
        %Three boxes will be shown with one taller than the others.

            
            %Must use prod(size(y)) and not numel(y) because numel will
            %return 1 for MATLAB versions earlier than 2015b
            
            if(nargin==1)
                y=a;
                numEl=prod(size(y));
                x=1:numEl;
                xWd=zeros(size(x));
                yWd=Interval.widthS(y);
            else
                x=a;
                y=b;
                numEl=prod(size(x));
                xWd=Interval.widthS(a);
                yWd=Interval.widthS(y);
            end
            
            xInf=Interval.infimumS(x);
            xSup=Interval.supremumS(x);
            yInf=Interval.infimumS(y);
            ySup=Interval.supremumS(y);
            
            for n=1:numEl
                if(xWd(n)==0 && yWd(n)==0)
                    plot(xInf(n),yInf(n),'.b','linewidth',2)
                elseif(xWd(n)==0 || yWd(n)==0)
                    plot([xInf(n),xSup(n)],[yInf(n),ySup(n)],'linewidth',2,'color','b')
                else
                    fill([xInf(n),xSup(n),xSup(n),xInf(n)],[yInf(n),yInf(n),ySup(n),ySup(n)],'g')
                end
                if n==1
                    tf=ishold;
                    if(~tf)
                        hold on
                    end
                end
            end
            if(~tf)
                hold off
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Misc Other functions                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [empty,c]=emptyIntersect(a,b)
        %%EMPTYINTERSECT Determine whether the intersection of two
        %       intervals (or an interval and a point) is empty and return
        %       the intersection region (or an empty interval). Empty
        %       intervals are denoted by having NaN for their bounds.
        %
        %OUTPUTS: empty A boolean value indicating whether the intersection
        %               region of the two intervals is empty.
        %             c The intersection interval of the two intervals.
        
            c=intersection(a,b);
            empty=isnan(c);
        end

        function [c,d]=relComp(aInt,bInt)
        %%RELCOMP Find the relative complement of Interval bInt in Interval
        %         aInt. That is, the part of Interval aInt that is not in
        %         bInt (set difference). There are two outputs as Interval
        %         bInt can be contained within Interval aInt.
        %
        %As an example, we can consider many possible scenarios
%         b=Interval(2,4);
%         aInt=Interval.empty(7,1);%Allocate space
%         aInt(1)=Interval(0,1);%a is before b
%         aInt(2)=Interval(5,6);%a is after b
%         aInt(3)=Interval(0,3);%Second half of a intersects b.
%         aInt(4)=Interval(3,5);%First half of a intersects b.
%         aInt(5)=Interval(0,5);%b is completely within a.
%         aInt(6)=Interval(2.5,3.5);%a is cmpletely within b.
%         bInt=repmat(b,size(aInt));
%         [c,d]=relComp(aInt,bInt)
        
            %Allocate space
            c=Interval.empty(size(aInt.lb));
            d=Interval.empty(size(aInt.lb));
        
            numEls=numel(aInt.lb);
            
            s.type='()';
            s.subs{1}=[];
            for curEl=1:numEls
                s.subs{1}=curEl;
                a=aInt.subsref(s);
                b=bInt.subsref(s);
                
                if(a.ub<b.lb || b.ub<a.lb)
                %The regions do not intersect at all.
                    minVal=a.lb;
                    maxVal=a.ub;
                elseif(a.lb<b.lb && a.ub<=b.ub)
                %Second half of a intersects with b.
                    minVal=a.lb;
                    maxVal=b.lb;
                elseif(a.lb>b.lb &&a.lb <b.ub && a.ub>b.ub)
                %First half of a intersects with b
                    minVal=b.ub;
                    maxVal=a.ub;
                elseif(a.lb<b.lb && b.ub <a.ub)
                %b is completely within a. This results in two intervals.
                    minVal=a.lb;
                    maxVal=b.lb;
                    d=d.subsasgn(s,Interval(b.ub,a.ub));
                else
                %a is completely within b. Return the empty set.
                	minVal=NaN;
                    maxVal=NaN;
                end
                
                c=c.subsasgn(s,Interval(minVal,maxVal));
            end
        end
        
        function c=mince(a,n)
        %%MINCE Given a column vector a of intervals that are n equally
        %       spaced subintervals of each element of a. Each row of the
        %       output corresponds to the set of subintervals for that row
        %       on the input. Subintervals are returned in ascending
        %       order. Overlap due to the finite precision (and rounding)
        %       of the upper and lower bounds is possible. n is the number
        %       of subintervals to create. If n is omitted or an empty
        %       matrix is passed, n=100 is used.

            if(nargin<2||isempty(n))
                n=100;
            end
            %Allocate space for the result.
            c=Interval.empty(length(a),n);
    
            numEls=numel(a.lb);
            s.type='()';
            s.subs{1}=[];
            sC.type='()';
            sC.subs=[];
            for curEl=1:numEls
                s.subs{1}=curEl;
                aCur=a.subsref(s);

                setProcRoundingMode(0);%Round to -Inf
                lbNew=linspace(aCur.lb,aCur.ub,n+1);
            
                setProcRoundingMode(3);%Round to +Inf
                ubNew=linspace(aCur.lb,aCur.ub,n+1);

                sC.subs={curEl,1:n};
                c=c.subsasgn(sC,Interval(lbNew(1:end-1),ubNew(2:end)));
            end
        end
    end
    
    methods (Static)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Constant functions of Chapter 10.5.2 of the IEEE Std     %%%
%%% 1788-2015                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function a=empty(varargin)
        %%EMPTY Returns an empty Interval or a matrix of empty Intervals
        %       with the given dimensions.
        %
        %INPUTS: This function is called like the zeros function. The input
        %        can either be a collection of a,b,c... specifying the
        %        size of each dimensions, or it can be 
        %        N A row vector specifying the size of each of the
        %          dimensions of the matrix of empty intervals that is
        %          created. If a scalar is passed, then an NXN matrix is
        %          generated. If omitted, a 1X1 empty interval is created.
        %
        %OUTPUTS: a A matrix of empty intervals having the given
        %           dimensions. Empty intervals have NaN values as their
        %           dimensions.
        %
        %This function is defined in Chapter 10.5.2 of [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            if(isempty(varargin))
                varargin{1}=1;
            end
            
            a=Interval(NaN(varargin{:}),NaN(varargin{:}));
        end
        
        function a=entire(varargin)
        %%ENTIRE Returns an Interval spanning from -Inf to Inf or a matrix
        %        of empty Intervals from -Inf to Inf. with the given
        %        dimensions.
        %
        %INPUTS: This function is called like the zeros function. The input
        %        can either be a collection of a,b,c... specifying the
        %        size of each dimensions, or it can be 
        %        N A row vector specifying the size of each of the
        %          dimensions of the matrix of entire intervals that is
        %          created. If a scalar is passed, then an NXN matrix is
        %          generated. If omitted, a 1X1 entire interval is created.
        %
        %OUTPUTS: a A matrix of entire intervals having the given
        %           dimensions. 
        %
        %This function is entire in Chapter 10.5.2 of [1].
        %
        %REFERENCES:
        %[1] IEEE Standard for Interval Arithmetic," in IEEE Std 1788-2015,
        %    pp.1-97, June 30 2015.
        
            if(isempty(varargin))
                varargin{1}=1;
            end
            
            a=Interval(-Inf(varargin{:}),Inf(varargin{:}));
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Other static functions                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function c=zeros(varargin)
        %%ZEROS Create a matrix of Intervals where both the upper and lower
        %       bounds of each interval equals zero.
        %
        %INPUTS: This function is called like the zeros function built into
        %        Matlab. The input can either be a collection of a,b,c...
        %        specifying the size of each dimensions, or it can be 
        %        N A row vector specifying the size of each of the
        %          dimensions of the matrix of zero intervals that is
        %          created. If a scalar is passed, then an NXN matrix is
        %          generated. If omitted, a 1X1 zero interval is created.

            if(isempty(varargin))
                varargin{1}=1;
            end

            zeroMat=zeros(varargin{:});
            c=Interval(zeroMat,zeroMat);
        end

        function c=infimumS(a)
        %%INFIMUMS Return the lower bound (infimum) of an Interval or a
        %     constant. Though there is a non-static infimum function, this
        %     one is a convenience function that can be used with Intervals
        %     or non-Interval data types.
        
            if(~isa(a,'Interval'))
                c=a;
            else
                c=a.lb;
            end
        end
        
        function c=supremumS(a)
        %%SUPREMUMS Return the upper bound (supremum) of an Interval or a
        %     constant. Though there is a non-static supremum function,
        %     this one is a convenience function that can be used with Intervals
        %     or non-Interval data types.
        
            if(~isa(a,'Interval'))
                c=a;
            else
                c=a.ub;
            end
        end

        function c=widthS(a)
        %%WIDTHS     Find the width of the interva (upper bound minus lower
        %            bound). Though there is a non-static width function,
        %            this one is a convenience function that can be used
        %            with Intervals or non-Interval data types.
        
            if(~isa(a,'Interval'))
                c=a;
            else
                c=wid(a);
            end
        end

        function c=pi()
        %%PI Obtain an interval for the constant pi. This is an interval to
        %    take into account the finite precision limitations of the
        %    processor.
        %
        %IN Section II of 1, it is suggested that one obtain bounds on pi
        %by using the atan function with different rounding modes. Here,
        %that is done, but one is just using 4*atan(1) as opposed to the
        %more complicated method in [1].
        %
        %REFERENCES
        %[1] M. Daumas, D. Lester, C. Muñoz, "Verified Real Number
        %    Calculations: A Library for Interval Arithmetic", IEEE
        %    Transactions on Computers, vol.58, no. 2, pp. 226-237,
        %    Feb. 2009
        
            setProcRoundingMode(0);%Round to -Inf
            lb=4*atan(1);%The lower bound.

            setProcRoundingMode(3);%Round to +Inf
            ub=4*atan(1);
            c=Interval(lb,ub);
        end
    end
    
    methods (Access=protected)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Overloaded methods to enable custom display of Interval objects %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function displayScalarObject(a)
            %%DISPLAYSCALAROBJECT Overloaded display operation for scalar
            %Intervals
            fprintf('[%s,%s]\n\n',num2str(a.lb),num2str(a.ub))
        end
        
        function displayNonScalarObject(a)
        %%DISPLAYSCALAROBJECT Overloaded display operation for non-scalar
        %Intervals
            dispCellLb=num2cell(a.lb);
            dispCellUb=num2cell(a.ub);
            disp(cellfun(@(x,y) sprintf('[%s,%s]',num2str(x),num2str(y)),dispCellLb,dispCellUb,'UniformOutput',false))
        end
        
        function displayEmptyObject(a)
        %%DISPLAYEMPTYOBJECT Overloaded display operation for empty
        %Intervals
            dimstr=matlab.mixin.CustomDisplay.convertDimensionsToString(a);
            className = matlab.mixin.CustomDisplay.getClassNameForHeader(a);
            emptyHeader = [dimstr,' empty ',className];
            header = sprintf('%s\n',emptyHeader);
            disp(header)
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
