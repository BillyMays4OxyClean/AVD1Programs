%% Wattegedara, Prabhasha Dinesh B A
clear
close all
clc

%% **INPUT**


M0 = 0.8;
alt = 15240;     % meters
p0 = 7.495E-4;     % atm
t0 = 270.65;    % Kelvin

hPR = 42800;
R = 287;
gamma = 1.4;
cp = 1.004;
a0 = sqrt(gamma*R*t0); % m/s
v0 = M0*a0;
Treq = 26690;      % Newtons
range = 7960000;   % meters


for alpha = 1:0.5:5
    for piC = 5:1:30
        % FREE STREAM (Station 0)


        tauR = 1+((gamma-1)/2)*M0^2;
        piR = (1+((gamma-1)/2)*M0^2)^(gamma/(gamma-1));


        Tt0 = t0*tauR;
        pto = p0*piR;



        %% DIFFUSER INLET (Station 2)


        tauD = 1;
        piD = 0.96;


        Tt2 = Tt0*tauD;
        pt2 = pto*piD;


        %% COMPRESSOR EXIT (Station 3)

        %piC = 5:1:30
        tauC = (piC).^((gamma-1)./gamma);


        Tt3 = Tt2.*tauC;
        pt3 = pt2.*piC;



        %% FAN EXIT (Station 13)


        piF = 1:0.2:6;   % Fan pressure Ratio
        tauF = (piF).^((gamma-1)./gamma);


        Tt13 = Tt2.*tauF;
        pt13 = pt2.*piF;


        %% COMBUSTOR EXIT (Station 4)

        Tt4 = 1500;
        piB = 1;
        tauB = Tt4/Tt3;
        taulambda = Tt4./t0;


        pt4 = pt3.*piB;     % atm


        % Fuel-Air Ratio


        f = ((cp*t0)/hPR)*(taulambda-(tauR*tauC));




        %% COMPRESSOR-TURBINE COUPLING (Station 5)


        %alpha = 1:0.5:13.5
        tauT = 1-(tauR./taulambda).*((tauC-1)+alpha.*(tauF-1));
        piT = (tauT).^(gamma./(gamma-1));


        Tt5 = Tt4.*tauT;
        pt5 = pt4.*piT;




        %% PRIMARY NOZZLE EXIT (Station 9)


        piNP = 1;
        tauNP=1;
        Tt9 = Tt5.*tauNP;
        pt9 = pt5.*piNP;
        NPRp = piR.*piC.*piT;


        %M9 = sqrt( (2./(gamma-1)) .* ((NPRp).^((gamma-1)./gamma) - 1) );
        M9 = (5.*((NPRp).^(2/7))).^(1/2);
        Tt9_T9 = 1+((gamma-1)./2).*M9.^2;
        T9 = Tt9./Tt9_T9;
        a9 = sqrt(gamma.*R.*T9);
        v9 = M9.*a9;
        v9_a0 = v9./a0;


        %% SECONDARY NOZZLE EXIT (Station 19)


        piNS = 1;
        tauNS=1;
        Tt19 = Tt13.*tauNS;
        pt19 = pt13.*piNS;
        NPRs = piR.*piF;


        M19 = sqrt((2./(gamma-1)).*((NPRs).^((gamma-1)./gamma)-1));
        Tt19_T19 = 1+((gamma-1)./2).*M19.^2;
        T19 = Tt19./Tt19_T19;
        a19 = sqrt(gamma.*R.*T19);
        v19 = M19.*a19;
        v19_a0 = v19./a0;


        %% ENGINE PERFORMANCE


        %alpha = 1:0.5:13.5;
        Sp_Thrst = a0.*(((1./(1+alpha)).*((v9./a0)-M0))+((alpha./(1+alpha)).*((v19./a0)-M0)));
        SFC = (f./((Sp_Thrst).*(1+alpha))).*1000;
        nTH = 1-1./(tauR.*tauC);


        nP = (Sp_Thrst.*(1+alpha).*v0)./(0.5.*(1+f).*(v9.^2)+0.5.*alpha.*(v19.^2)-0.5.*(1+alpha).*(v0.^2));
        nO = nTH.*nP;


        mdot = Treq./Sp_Thrst;
        mdotf = mdot.*f;


        %Data1 = table(SFC(:),Sp_Thrst(:),nO(:),nTH(:),nP(:),piC(:),Tt4(:),alpha(:),piF(:),mdot(:),mdotf(:));


        Area = mdot.*sqrt(t0)/(((p0*101325)*M0*sqrt(gamma/R)*((1/(1+((gamma-1)*M0^2/2)))^((gamma+1)/(2*(gamma-1))))));
        Diameter = sqrt(4*Area/pi());
    end
    plot(5:1:30,Sp_Thrst)
    xlabel('\pi_C')
    ylabel('')
    hold on
end