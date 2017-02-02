function xSim=WeakStochRungeKStep(method,xCur,curT,a,D,deltaT,W,V,dadt,dadx,dDdx,d2adx,d2Ddx)
%%WEAKSTOCHRUNGEKSTEP Simulate a nonlinear continuous-time random process
%                     specified by the Langevin equation forward in time by
%                     a step-size of deltaT using a weak Taylor scheme.
%
%INPUTS: method  An integer defining which weak Taylor method to use.
%                1  The Euler-Maruyama method
%                2  The explicit order 2.0 weak Taylor scheme with additive
%                   noise. This requires the diffusion function to be
%                   additive.
%                3  The explicit order 2.0 weak Taylor scheme
%                4  The simplified order 2.0 weak Taylor scheme. This
%                   requires the diffusion function to be autonomous.
%        xCur    The initial target state at time curT.
%        curT    The time of the initial state xCur.
%        a       The drift function in the continuous-time stochastic
%                dynamic model. It takes the state and a time variable
%                as its arguments, a(x,t).
%        D       The diffusion function in the continuous-time stochastic
%                dynamic model. It takes the state and a time variable as
%                its arguments, D(x,t). If the process noise must be
%                additive, D(x,t) should not depend on x. If the process
%                noise must be autonomous, D(x,t) should not depend on t.
%       deltaT   The size of the single step over which the strong
%                stochastic Runge-Kutta integration is performed.
%         W      A matrix of independent random variables which satisfies
%                the moment conditions in Eq. 14.2.3 of Kloeden. If not
%                supplied, Gaussian random variables are used.
%         V      A matrix of independent two-point distributed random
%                variables satisfying Eq. 14.2.8-10 in Kloeden. If not
%                supplied, these are calculated.
%     dadt       A xDimX1 matrix of the derivative of the drift function
%                with respect to time at state xCur a   nd time curT. This can
%                also be a function that takes the state and a time
%                variable as its arguments, dadt(x,t).
%     dDdt       A xDimXdColDim matrix of the derivative of the diffusion
%                function with respect to time at state xCur and time curT.
%                This can also be a function that takes the state and a
%                time variable as its arguments, dDdt(x,t). If the process
%                noise must be additive, dDdt(x,t) should not depend on x.
%                If the process noise must be autonomous, dDdt(x,t) should
%                not depend on t.
%     dadx       A xDimXxDim matrix of the derivative of the drift function
%                with respect to the state at state xCur and time curT.
%                This can also be a function that takes the state and a
%                time variable as its arguments, dadx(x,t).
%     dDdx       A xDimXxDimXdColDim matrix of the derivative of the
%                diffusion function with respect to state at state xCur and
%                time curT. This can also be a function that takes the
%                state and a time variable as its arguments, dDdx(x,t). If
%                the process noise must be additive, dDdx(x,t) should not
%                depend on x. If the process noise must be autonomous,
%                dDdx(x,t) should not depend on t. If not provided, this is
%                assumed to be zero.
%    d2adx       A xDimXxDimXxDim matrix of the second derivative of the
%                drift function with respect to the state at state xCur and
%                time curT. The value at point (m,k,l) represents
%                d2a(m)/dx(k)dx(l). This can also be a function that takes
%                the state and a time variable as its arguments,
%                d2adx(x,t). If not provided, this is assumed to be zero.
%    d2Ddx       A xDimXxDimXxDimXdColDim matrix of the second derivative
%                of the diffusion function with respect to state at state
%                xCur and time curT. The value at point (m,k,l,j)
%                represents d2D(m,j)/dx(k)dx(l). This can also be a
%                function that takes the state and a time variable as its
%                arguments, d2Ddx(x,t). If the process noise must be
%                additive, d2Ddx(x,t) should not depend on x. If the
%                process noise must be autonomous, d2Ddx(x,t) should not
%                depend on t. If not provided, this is assumed to be zero.
%OUTPUTS: xSim   The simulated target state at time curT+deltaT.
%
%The algorithm for the Euler-Marayama method is described in Chapter 10.2
%of [1].
%
%The algorithms for both the additive and non-additive explicit order 2.0
%Taylor scheme are described in Section 15.1 of [1]..
%
%The algorithm for the simplified order 2.0 Taylor scheme is described in
%Section 14.2 of [1]..
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%
%April 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

aCur=a(xCur,curT);
DCur=D(xCur,curT);
dColDim=size(DCur,2);
xDim=size(aCur,1);

%Generate the random component following Eq. 11.2.4,11.2.8-10 in Kloeden
if(nargin<7 || isempty(W))
    W=sqrt(deltaT)*randn(dColDim);
end
if(nargin<8 || isempty(V))
    V=-deltaT*eye(dColDim);
    for j=1:dColDim
        for j2=(j+1):dColDim
            V(j,j2)=sign(randn)*deltaT;
            V(j2,j)=-V(j,j2);
        end
    end
end

switch method
    case 1 % Euler-Marayama
        xSim=xCur+aCur*deltaT+DCur*W(:,1);
    case 2 % Explicit additive 2.0 - Eq. 15.1.4
        y=xCur+aCur*deltaT+sum(DCur*W,2);
        xSim=xCur+(deltaT/2)*(a(y,curT)+aCur)+sum(DCur*W,2);
    case 3 % Explicit 2.0 - Eq. 15.1.3
        y=xCur+aCur*deltaT+sum(DCur*W,2);
        f1Val=xCur+(deltaT/2)*(a(y,curT)+aCur);
        
        Rp=bsxfun(@plus,xCur+aCur*deltaT,DCur*sqrt(deltaT));
        Rm=bsxfun(@minus,xCur+aCur*deltaT,DCur*sqrt(deltaT));
        Up=bsxfun(@plus,xCur,DCur*sqrt(deltaT));
        Um=bsxfun(@minus,xCur,DCur*sqrt(deltaT));
        
        DRp=zeros(size(DCur));
        DRm=zeros(size(DCur));
        DUpm=zeros(xDim,dColDim,dColDim);
        for j=1:dColDim
            Dp=D(Rp(:,j),curT);
            Dm=D(Rm(:,j),curT);
            DRp(:,j)=Dp(:,j);
            DRm(:,j)=Dm(:,j);
            DUpm(:,:,j)=D(Up(:,j),curT)-D(Um(:,j),curT);
        end
        
        f2Val=(DRp+DRm+((4-2)*DCur))*W;
        f3Val=(DRp-DRm)*(W.^2-deltaT);
        for j=1:dColDim
            for r=1:dColDim
                if r==j
                    continue;
                end
                f2Val(:,j)=f2Val(:,j)+DUpm(:,:,r)*W(:,j);
                f3Val(:,j)=f3Val(:,j)+DUpm(:,:,r)*(W(:,j).*W(:,r)+V(r,j));
            end
        end
        
        xSim=f1Val+(1/4)*sum(f2Val,2)+(1/(4*sqrt(deltaT)))*sum(f3Val,2);
    case 4 % Order 2.0
        if nargin<13
            d2Ddx=zeros(xDim,xDim,xDim,dColDim);
        end
        if nargin<12
            d2adx=zeros(xDim,xDim,xDim);
        end
        
        if(isa(dadt,'function_handle'))
            dadt=dadt(xCur,curT);
        end
        if(isa(dadx,'function_handle'))
            dadx=dadx(xCur,curT);
        end
        if(isa(dDdx,'function_handle'))
            dDdx=dDdx(xCur,curT);
        end
        if(isa(d2adx,'function_handle'))
            d2adx=d2adx(xCur,curT);
        end
        if(isa(d2Ddx,'function_handle'))
            d2Ddx=d2Ddx(xCur,curT);
        end
        
        %L operators
        L0a=taylorOperator(aCur,DCur,0,dadx,dadt,d2adx);
        L0D=taylorOperator(aCur,DCur,0,dDdx,[],d2Ddx);
        Lja=taylorOperator(aCur,DCur,[],dadx);
        LjD=taylorOperator(aCur,DCur,[],dDdx);
        
        %Simplified order 2.0 weak Taylor scheme(Eq. 14.2.7 in Kloeden)
        f1Val=xCur+aCur*deltaT+(deltaT^2/2)*L0a;
        f2Val=DCur + (deltaT/2)*(L0D-Lja);
        f3Val=0;
        for j1=1:dColDim
            f3Val=f3Val+(1/2)*LjD(:,:,j1)*(bsxfun(@plus,bsxfun(@times,W(:,j1),W),V(:,j1)));
        end
        
        xSim=f1Val+sum(f2Val*W,2)+sum(f3Val,2);
    otherwise
        error('Not a valid method')
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
