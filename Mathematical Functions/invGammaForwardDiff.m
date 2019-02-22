function g=invGammaForwardDiff(z,epsVal)
%%INVGAMMAFORWARDDIFF This function evaluate the forward difference for two
%         inverse gamma functions: (1/gamma(z+epsVal)-1/gamma(z))/epsVal .
%         This in done in such a manner that the result is numerically
%         stable for epsVal>=0, whereas simply evaluating
%         (1/gamma(z+epsVal)-1/gamma(z))/epsVal would fail as epsVal
%         approaches zero, whereas this approaches the asymptotic limit of
%         -(polygamma(z)/gamma(z)).
%
%INPUTS:      z A real scalar value with z>0.
%        epsVal A positive real value with epsVal>=0.
%
%OUTPUTS: g The value (1/gamma(z+epsVal)-1/gamma(z))/epsVal
%
%This function implements the algorithm for evaluating such finite
%differences in Chapter III or [1]. The algorithm adjusts the value of z to
%be in a certain range (create an additional multiplicative coefficient for
%the final solution in the process) and then utilizes Chebyshev
%interpolation to solve the problem in that range. The interpolation
%coefficients are not explicitly tabulated in [1], but were provided by
%the author of [1].
%
%REFERENCES:
%[1] R. C. Forrey, "Computing the hypergeometric function," Journal of
%    Computational Physics, vol. 137, no. 1, pp. 79-100, Oct. 1997.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

x=z-1;
sumVal=0;
coeff=1;
negX=0;
while(1)
    if (x>1)
        sumVal=sumVal+coeff*gamma(x+epsVal);
        coeff=coeff*x;
        x=x-1;
    elseif(x<0)
        x1=x+epsVal+2;
        coeff1=1;
        while(x1<1)
          coeff1=x1*coeff1;
          x1=x1+1;
        end
        sumVal=sumVal+coeff*coeff1/gamma(x1);
        coeff=coeff*(x+1);
        x=x+1;
        negX=1;
    else
        break;
    end
end

if((x>=0.25)&&(x<=0.75))
    c1=getChebyshevCoeffs(1);
    numCoeffs=length(c1);
    T1=zeros(numCoeffs,1);
    F1=zeros(numCoeffs,1);

    T1(0+1)=1;
    T1(1+1)=2*(x+epsVal)-1;
    F1(0+1)=0;
    F1(1+1)=2;
    temp=c1(1+1)*F1(1+1);
    for i=2:(numCoeffs-1)
        T1(i+1)=(4*(x+epsVal)-2)*T1(i-1+1)-T1(i-2+1);
        F1(i+1)=4*T1(i-1+1)+(4*x-2)*F1(i-1+1)-F1(i-2+1);
        temp=temp+c1(i+1)*F1(i+1);
    end
elseif ((x>=0)&&(x<0.25))
    c1=getChebyshevCoeffs(2);
    numCoeffs=length(c1);
    T1=zeros(numCoeffs,1);
    F1=zeros(numCoeffs,1);

    T1(0+1)=1;
    T1(1+1)=2*(x+epsVal);
    F1(0+1)=0;
    F1(1+1)=2;
    temp=c1(1+1)*F1(1+1);
    for i=2:(numCoeffs-1)
        T1(i+1)=4*(x+epsVal)*T1(i-1+1)-T1(i-2+1);
        F1(i+1)=4*T1(i-1+1)+4*x*F1(i-1+1)-F1(i-2+1);
        temp=temp+c1(i+1)*F1(i+1);
    end
elseif((x>0.75)&&(x<=1))
    c1=getChebyshevCoeffs(3);
    numCoeffs=length(c1);
    T1=zeros(numCoeffs,1);
    F1=zeros(numCoeffs,1);

    T1(0+1)=1;
    T1(1+1)=2*(x+epsVal)-2;
    F1(0+1)=0;
    F1(1+1)=2;
    temp=c1(1+1)*F1(1+1);
    for i=2:(numCoeffs-1)
        T1(i+1)=(4*(x+epsVal)-4)*T1(i-1+1)-T1(i-2+1);
        F1(i+1)=4*T1(i-1+1)+(4*x-4)*F1(i-1+1)-F1(i-2+1);
        temp=temp+c1(i+1)*F1(i+1);
    end
end

if(negX==0)
    x1=z;
    coeff1=1;
    while(x1<1)
        coeff1=x1*coeff1;
        x1=x1+1;
    end
    x2=z+epsVal;
    coeff2=1;
    while(x2<1)
        coeff2=x2*coeff2;
        x2=x2+1;
    end
    temp=sumVal+coeff*temp;
    g=-temp*coeff1*coeff2/gamma(x1)/gamma(x2);
else
    x1=x+1;
    coeff1=1;
    while(x1<1)
        coeff1=x1*coeff1;
        x1=x1+1;
    end
    x2=x+1+epsVal;
    coeff2=1;
    while(x2<1)
        coeff2=x2*coeff2;
        x2=x2+1;
    end
    coeff=-coeff*coeff1*coeff2/gamma(x1)/gamma(x2);
    g=sumVal+coeff*temp;
end
end

function c=getChebyshevCoeffs(coeffNum)
switch(coeffNum)
    case 1
        %The CI interpolating coefficients of Section III
        c=[0.94178559779549466571096003120435196;
           0.44153813248410067571913157711414607e-2;
           0.56850436815993633786326645888162378e-1;
          -0.42198353964185605010125001866024699e-2;
           0.13268081812124602205840067963889683e-2;
          -0.18930245297988804325239470239464680e-3;
           0.36069253274412452565780822094442805e-4;
          -0.60567619044608642184855483216922771e-5;
           0.10558295463022833447318234541645507e-5;
          -0.18119673655423840482918555144273961e-6;
           0.31177249647153222777902517006137963e-7;
          -0.53542196390196871408740949118221475e-8;
           0.91932755198595889468877475468573503e-9;
          -0.15779412802883397617671187106425584e-9;
           0.27079806229349545432695717700017206e-10;
          -0.46468186538257301439531283506784063e-11;
           0.79733501920074196555512936759234830e-12;
          -0.13680782098309160264738694164685656e-12;
           0.23473194865638006534799539031857605e-13;
          -0.40274326149490669507857892267787757e-14;
           0.69100517473721009958174457696435176e-15;
          -0.11855845002219929396593062972684083e-15;
           0.20341485424963760969383490105975402e-16;
          -0.34900543417173691101844936408331408e-17;
           0.59879938564842634972645168624438135e-18;
          -0.10273780578716378747008169519685451e-18;
           0.17627028160574041125936108594612916e-19;
          -0.30243206536626379817809691872233988e-20;
           0.51889146600668142375785699199940389e-21;
          -0.89027708392150216484577040964212789e-22;
           0.15274740724470977041487116294681806e-22;
          -0.26207312865170684216151526387496724e-23;
           0.44964644619824783627762340991300087e-24;
          -0.77147147879836211531329396406348717e-25;
           0.13236365808260955301316348853544449e-25;
          -0.22709797413377406198008958539204735e-26;
           0.38966913277073699893252807432563276e-27;
          -0.66795989154793901466615113245736539e-28;
           0.11456694360946249087722449327564468e-28;
          -0.20956088513945987438866120550893160e-29;
           0.34345153487326051089311279207743562e-30;
          -0.74448389617685196161619686887550341e-31];
    case 2
        %The CII coefficients of Section III. The first element is not
        %used.
        c= [0.11528686913857579339872890819003657e1;
           -0.39836641427188668813550502856567435;
            0.16381491849746834445969671065563396;
           -0.41349972584595838242416447164595642e-1;
            0.11739888104509743948748485834561229e-1;
           -0.31509159742825717845846783104528302e-2;
            0.85084809366682540330028115184077086e-3;
           -0.22845443192182297253614554810213881e-3;
            0.61296656896858907270916323759970391e-4;
           -0.16433766723011959082591541534833589e-4;
            0.44046701847148520660258125028242579e-5;
           -0.11803851479587223345492859134791582e-5;
            0.31630339312403588488305625683201151e-6;
           -0.84755796666686117564957022251013564e-7;
            0.22710572677209079780536954678987573e-7;
           -0.60853209609268373214751556259951644e-8;
            0.16305620921375867864482570008163625e-8;
           -0.43690846345047718022878883179027790e-9;
            0.11706935476739890379554689241357534e-9;
           -0.31368649843198552351255033209421610e-10;
            0.84052057618382692960217222664957228e-11;
           -0.22521682699590609081199019088965996e-11;
            0.60346669123807723976181127096882828e-12;
           -0.16169841538137032176079290114309245e-12;
            0.43326960175123609635570088625382667e-13;
           -0.11609424034675431553315176322024985e-13;
            0.31107358004300087572452155428660087e-14;
           -0.83351914632193111475558815401948979e-15;
            0.22334078222557889355389486422061460e-15;
           -0.59843982246058550382747881611851515e-16;
            0.16035146716190080240936859943115090e-16;
           -0.42966046133076898235808019603294715e-17;
            0.11512717363557431988678458870224873e-17;
           -0.30848233202835882015258583966299712e-18;
            0.82657591746540727258216017499064442e-19;
           -0.22148034956862123422799663231945171e-19;
            0.59345480806145642339133686333296721e-20;
           -0.15901573656881585725893714030807897e-20;
            0.42608138203898096080539369435375448e-21;
           -0.11416816226321087557458906349840213e-21;
            0.30591266842950015571055286508657438e-22;
           -0.81969053674548061989664444282339330e-23;
            0.21963543471485197662543467891802004e-23;
           -0.58851140572211577956963471197095354e-24;
            0.15769121438531798083082131134888596e-24;
           -0.42253211944581570323425035302537635e-25;
            0.11321706791574145306428072576766804e-25;
           -0.30335842761477973373797446515125892e-26;
            0.81281383350578045680446098123885346e-27;
           -0.21782407988772728568103833180457024e-27;
            0.58395544064782062129754390403734767e-28;
           -0.15729062977489325257494410942884130e-28;
            0.42390612257722955199550993363196147e-29;
           -0.11242203351086692027388616387423238e-29;
            0.27892280419588143241883200553486195e-30;
           -0.75766427928255356179910217971637866e-31];
    case 3
        %The CIII coefficents of Section III.
        c= [0.10532770878177862619534128247576828e1;
            0.21902166104535936497306369004840667;
            0.53885821783347712865216341722976574e-1;
            0.25387290658986838596948519579519148e-2;
            0.61466596479014144199820446583715941e-3;
           -0.32319247384294465724865638122474435e-5;
            0.60054921157267140200751871810266970e-5;
           -0.41824428090189489334617924547407754e-6;
            0.74607235650174366232051332482639985e-7;
           -0.84349526185192483560074198183789434e-8;
            0.11322169721817117406057072389666464e-8;
           -0.14175349900034682206860980369914924e-9;
            0.18156967683771854495445069753509525e-10;
           -0.23052163748763990586386231147733255e-11;
            0.29327030584105892891631030300077869e-12;
           -0.37268590170679729030689484336505900e-13;
            0.47360432581610222494078892575939043e-14;
           -0.60172423075766780010690060490450222e-15;
            0.76443979970650480527157880770622904e-16;
           -0.97108892590783757664936380167684001e-17;
            0.12335488659810502174628042595177563e-17;
           -0.15668997427423797214874298423999374e-18;
            0.19902969432180950952170748993213290e-19;
           -0.25280701093316992983208535829903356e-20;
            0.32111217127088658654008440525466587e-21;
           -0.40787027055654288157193053732139852e-22;
            0.51806681115442807351458062924762066e-23;
           -0.65803415226414646040514708695329147e-24;
            0.83581632724068042390791744946381128e-25;
           -0.10616267321620223331012310816058461e-25;
            0.13484159784261929973156667845986312e-26;
           -0.17130640476670792317750095910458264e-27;
            0.21720215147689502411187819143753676e-28;
           -0.27633054946463729557612727034555572e-29;
            0.26664265210535867308016959008022142e-30];
    otherwise
        error('Unknown coefficient number given')
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
