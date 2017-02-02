function [xPredMain,xPredSubsid,k,orders]=RosenbrockStep(xCur,tCur,f,deltaT,fCur,J,dfdt,order)
%%ROSENBROCKSTEP Perform a single step of a Rosenbrock method whose main
%                integrating order is given by order. All implemented
%                Rosenbrock routines are embedded routines, meaning that a
%                secondary result of a different order, either one order
%                higher or one lower, can also be returned along with the
%                order of the second output. The second output can be used
%                in an algorithm that has an adaptive set size, such as in
%                the function RosenbrockAdaptiveOverRange. Rosenbrock
%                methods are diagonal implicit Runge Kutta methods that are
%                often well-suited for use with stiff ordinary differential
%                equations. That is, the algorithm is for integrating
%                dx/dt=f(x,t) given initial conditions (x0,t0) when the
%                form of f makes it such that standard Runge-Kutta methods
%                will fail/ will be required to take extremely small
%                stepsizes. However, unlike Runge-Kutta methods, Rosenbrock
%                methods require that derivative information be provided.
%                To integrate over a range of values adaptively choosing
%                the stepsize, the function RosenbrockAdaptiveOverRange
%                should be used.
%
%INPUTS: The function can either be called with one input as
%        [orders,isFSAL]=RosenbrockStep(order,solutionChoice)
%        to just get both of the orders of a particular method whose main
%        order is given by order (see below for details), or it can be run
%        as
%        [xPred,xPredSubsid,dxdt,orders,isFSAL]=RosenbrockStep(xCur,curT,f,fCur,deltaT,Jacobian,dfdt,order)
%        to actually do a full step and provide the orders. All of the
%        inputs in the full step are:
%           xCur    The value of the (scalar or vector) state over which
%                   integration is being performed.
%           tCur    The time at which xCur is taken.
%           f       f(x,t) returns the derivative of x with respect to time
%                   taken at time t.
%           fCur    The value f(xCur,tCur). This is explicitely provided
%                   rather than just using the function f, because when
%                   performing multiple steps, sometimes the value is
%                   computed while computing xPredSubsid and can thus be
%                   reused on the next step. If an empty matrix is
%                   provided, then this function will just compute fCur
%                   from f.
%           deltaT  The size of the single (time) step over which the 
%                   integration is performed.
%               J   The symmetric derivative matrix of f with respect to
%                   the components of xCur. For example, the function
%                   numjac could be used to numerically approximate the
%                   Jacobian.
%             dfdt  The derivative of f with respect to t.
%            order  The main integration order of the Rosenbrock method.
%                   If this parameter is omitted, then the default order of
%                   2 is used. Orders can range from 2 to 3. Given order 2,
%                   the Rosenbrock formula of [2] is used. Given order 3,
%                   the Rosenbrock formula of [1] is used.
%
%OUTPUTS: If the function is run with one input as
%         [orders,isFSAL]=RosenbrockStep(order)
%         Then there are just two outputs, orders and isFSAL, described
%         below.
%         Otherwise, if the stepsize is such that a certain matrix based on
%         the jacobian and stepsize is not full rank, then all empty
%         matrices (except for order) are returned. In such an instance,
%         change the step size and try again. Otherwise, all of the
%         outputs are
%           xPredMain The state propagated forward an interval of deltaT
%                     at the order of precision given by the input order.
%                     This is the main integration order.
%         xPredSubsid The subsidiary estimate. This is a different order
%                     than the main integration order can can be used in
%                     algorithms that adaptively adjust the stepsize.
%                   k The values of the derivatives f evaluated at various
%                     points as determined by the selected algorithm. These
%                     can sometimes be reused in interpolation routines to
%                     reduce the number of computations required. In all of
%                     the methods, k(:,1) is equal to f(xVal,curT). In
%                     FSAL methods, k(:,end) is
%                     f(xPredMain,curT+deltaT)
%              orders A 2X1 vector whereby order(1)=order, the input
%                     parameter, and orders(2) is the order of xPredSubsid.
%                     The value of orders(2) depends on the algorithm
%                     chosen by the combination of the order and
%                     solutionChoice inputs.
%              isFSAL Indicates whether the first evaluation of the f of
%                     the next step is equal to k(end,:)'. If so, this can
%                     be passed to the function on the next step rather
%                     than having to make an additional evaluation of f.
%
%Rosenbrock methods are diagonal implicit Runge-Kutta methods that perform
%only a single iteration of Newton's method each step to solve the implicit
%equations. The Rosenbrock formula used are the 2(3) formula in [2] used in
%Matlab's ode23s function, and the 3(4) formula in Section 5 of [1].
%
%REFERENCES:
%[1] L. F. Shampine, "Implementation of Rosenbrock methods," ACM Transactions
%on Mathematical Software, vol. 8, no. 2, pp. 93-113, Jun. 1982.
%[2] L. F. Shampine and M. W. Reichelt, "The Matlab ODE suite," Journal on
%Scientific Computing, vol. 18, no. 1, pp. 1-22, Jan. 1997.
%
%February 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin==1)
        %If the function was just called as
        %orders=RosenbrockStep(order);
        order=xCur;
        if(order==2)
            xPredMain=[2;3];%orders
            xPredSubsid=true;%isFSAL
        else
            xPredMain=[3;4];%orders
            xPredSubsid=false;%isFSAL
        end
        return;
    end
    
    if(nargin<8)
        order=2;
    end
    
    if(isempty(fCur))
       fCur=f(xCur,tCur); 
    end

    xDim=size(xCur,1);
    switch(order)
        case 2
            %Use the 2(3) formula of [2]
            orders=[2;3];

            d=1/(2+sqrt(2));
            e32=6+sqrt(2);

            %The coefficient due to Wolfbrandt.
            W=eye(xDim,xDim)-deltaT*d*J;

            if(rank(W)==xDim)
                F0=fCur;
                k1=W\(F0+deltaT*d*dfdt);
                F1=f(xCur+(1/2)*deltaT*k1,tCur+0.5*deltaT);
                k2=W\(F1-k1)+k1;
                %The second-order estimate.
                xPredMain=xCur+deltaT*k2;

                %This is F0 on the next estimate since this is a FSAL  
                %(first same as last) method.
                F2=f(xPredMain,tCur+deltaT);

                %The third order corrector estimate.
                k3=W\(F2-e32*(k2-F1)-2*(k1-F0)+deltaT*d*dfdt);
                xPredSubsid=xPredMain-(deltaT/6)*(k1-2*k2+k3);
                k=zeros(xDim,3);
                k(:,1)=F0;
                k(:,2)=F1;
                k(:,3)=F2;
            else
                xPredMain=[];
                xPredSubsid=[];
                k=[];
            end
        case 3
            orders=[3;4];
            %Use the 3(4) formula from [1]

            W=eye(xDim,xDim)-(1/2)*deltaT*J;

            if(rank(W)==xDim)
                %If the matrix is not singular, then the step can be tried.
                F0=fCur;
                k1=W\(F0+deltaT*(1/2)*dfdt);
                F1=f(xCur+deltaT*k1,tCur+deltaT);
                k2=W\(F1-(3/2)*deltaT*dfdt-4*k1);
                F2=f(xCur+(24/25)*deltaT*k1+(3/25)*deltaT*k2,tCur+(3/5)*deltaT);
                k3=W\(F2+(121/50)*deltaT*dfdt+(186/25)*k1+(6/5)*k2);
                k4=W\(F2+(29/250)*deltaT*dfdt-(56/125)*k1-(27/125)*k2-(1/5)*k3);

                %The third order estimate
                xPredMain=xCur+deltaT*((98/108)*k1+(11/72)*k2+(25/216)*k3);
                %Thr fourth-order corrector.
                xPredSubsid=xCur+deltaT*((19/18)*k1+(1/4)*k2+(25/216)*k3+(125/216)*k4);
                k=zeros(xDim,3);
                k(:,1)=F0;
                k(:,2)=F1;
                k(:,3)=F2;
            else%If the matrix is singular, then the step size should be
                %adjusted.
                xPredMain=[];
                xPredSubsid=[];
                k=[];
            end
        otherwise
            error('Invalid order given')
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
