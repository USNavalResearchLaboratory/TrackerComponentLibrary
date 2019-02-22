function val=superfactorial(n)
%%SUPERFACTORIAL Evaluate the superfactorial of a nonnegative real integer.
%                This is factorial(1)*factorial(2)*...factorial(n).
%
%INPUTS: n A scalar or matrix of real integers >=0.
%
%OUTPUTS: val A matrix having the same dimensions as n holding the
%             superfactorial of the values in n.
%
%The superfactorial is defined in [1]. For values 18 and below, the value
%is taken from a tabel. For higher values, the value is computed from the
%definition starting from the value at 18, maintaining the current
%factorial value at each step. Values 27 and higher overflow with double
%precision arithmetic.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Superfactorial." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/Superfactorial.html
%
%September 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%A table of values for 0<=n<=18
valTable=[1;
          1;
          2;
          12;
          288;
          34560;
          24883200;
          125411328000;
          5056584744960000;
          1834933472251084800000;
          6658606584104736522240000000;
          265790267296391946810949632000000000;
          127313963299399416749559771247411200000000000;
          792786697595796795607377086400871488552960000000000000;
          69113789582492712943486800506462734562847413501952000000000000000;
          90378331112371142262979521568630736335023247731599748366336000000000000000000;
          1890966832292234727042877370627225068196418587883634153182519380410368000000000000000000000;
          672593129192865130334217631473916658864122332882577979675277211683839238972899328000000000000000000000000;
          4306192564997715382115598640379294845786123319603755168023536027873932927153136831171640950784000000000000000000000000000];

if(any(n(:)<0)||any(n(:)~=fix(n(:)))||any(~isreal(n(:))))
    error('n must be a positive integer.') 
end
val=zeros(size(n));

numVals=numel(n);

for curEl=1:numVals
    nCur=n(curEl);
    
    if(nCur<=18)
        val(curEl)=valTable(nCur+1);
    else
        factVal=6402373705728000;%factorial(18)
        val(curEl)=valTable(18+1);
        
        for curN=19:nCur
            %Quit if there was an overflow.
            if(~isfinite(val(curEl)))
                break;
            end
            factVal=factVal*curN;
            val(curEl)=val(curEl)*factVal;
        end
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
