%% Quickie for Miguel (hope it gives you a hardon)
clear
close all
clc

%% Inputs

ATMe = zeros(length(0:5:50),5); % Stored in the format given by Justin Rhinestone
card = 1; % we playin cards cuh
T = zeros(11,1);
P = zeros(11,1);
rho = zeros(11,1);
mu = zeros(11,1);
nu = zeros(11,1);

for alt = 0:1524:15240 % Altitude in thousands of feet.
    if alt == 0 % Since Justin's code has an aneurysm when you put in 0 meters for sea level, I just put in 1 meter lol.
        [~,ATMe(card,:)] = StandardATM(1);
        [~, ~, ~, T(card,1), P(card,1), rho(card,1), ~, ~, mu(card,1), nu(card,1), ~, ~, ~] = atmo(alt/1000,1524,2);
    else
        [~,ATMe(card,:)] = StandardATM(alt);
        [~, ~, ~, T(card,1), P(card,1), rho(card,1), ~, ~, mu(card,1), nu(card,1), ~, ~, ~] = atmo(alt/1000,1524,2);
    end
    card = card + 1; % Cardinal number representation since alt doesn't keep track of the number iteration we are at. (2 apples, etc)
end

Altidude = ATMe(:,1);
Temperature = ATMe(:,2);
Pressure = ATMe(:,3);
Density = ATMe(:,4);
a = ATMe(:,5);
ATM_Table = table(Altidude,Temperature,Pressure,Density,a); % Turn it into a pretty table

AoA = 5;

[Altitude,RS,Mach,V,Weight,q,n,Tmax,t,delV] = PablosWholeAssCode();

for i=1:length(Altidude)
    CD(i,1) = CallDisWholeAssFunction(V(i)/a(i),rho(i),mu(i),V(i),AoA);
end

Draggin = q.*CD*1316;

plot(0:5:50,Draggin,'-k','LineWidth',2)
xlabel('Altitude, kft')
ylabel('Drag, lbf')

DragginDeezNutsAcrossYourFace = table(Draggin);
plot(Altidude,Draggin,'-k','LineWidth',2)
title("Drag vs. Altitude ")
xlabel('Altitude, feet')
ylabel('Drag, lbf')
xtickangle(45)
ax = gca;
ax.XAxis.Exponent = 0;
writetable(DragginDeezNutsAcrossYourFace,'MiguelsPeen.xlsx','Sheet','DragDatalbf')

ExportTable = table(Altitude,RS,Mach,V,Weight,q,n,Tmax,t,delV);
writetable(ExportTable,'MiguelsPeen.xlsx','Sheet','PablosVelocityShit')
clc
% winopen('MiguelsPeen.xlsx')

function CD = CallDisWholeAssFunction(Mach,rho,mu,V,AoA)
    %% Created by Luke Patterson
    %% for the purpose of calculating total Vehicle drag during tandem flight

    %% Flight Condition
    % WK2 Cruises at 15200 meters at Mach 0.55
    % The flight condition class definition will contain and calculate all
    % values relating to atmospheric properties and airspeed

    unitSystem = 'FPS';

%     Cruise = FlightCondition;
%     Cruise.Altitude = 15200;
%     Cruise.Units = 'SI';
%     Cruise.Name = 'Cruise';
%     Cruise.AoA = 0;
%     Cruise = Cruise.SetSpeed('Mach',0.55);
%     Cruise = Cruise.ConvertUnits(unitSystem);
% 
%     rho = Cruise.Atmosphere.rho;
%     a = Cruise.Atmosphere.a;
%     M = Cruise.AirSpeed.Mach;
%     mu = Cruise.Atmosphere.mu;
%     V = Cruise.AirSpeed.V;

    %% Geometry Definitions

    MainWing = WingGeometry;
    MainWing.Name = 'Main Wing';
    MainWing.b = 141.0833;
    MainWing.Rc = 10.3937;
    MainWing.RootAirfoil = NACA4WingSection('2421',MainWing.Rc,128);
    MainWing.RootAirfoil.Parent = MainWing;
    MainWing.Tc = 5.1968;
    MainWing.TipAirfoil = NACA4WingSection('2309',MainWing.Tc,128);
    MainWing.TipAirfoil.Parent = MainWing;
    MainWing.Taper = 0.5;
    MainWing.AR = 15.1365;
    MainWing.MAC = 8.084;
    MainWing.Sweep = 0;
    MainWing.Swet = 2141.3;
    MainWing.Se = 1041.3;
    MainWing.S = 1316.882;
    MainWing.V = 2573;
    MainWing.Units = 'FPS';
    MainWing = MainWing.ConvertUnits(unitSystem);

    HT = WingGeometry;
    HT = HT.ImportFromCell(readcell('DataDump.xlsx','Sheet','HT'));
    HT = HT.ConvertUnits(unitSystem);

    VT = WingGeometry;
    VT = VT.ImportFromCell(readcell('DataDump.xlsx','Sheet','VT'));
    VT = VT.ConvertUnits(unitSystem);

    Fuselage = FuselageGeometry;
    Fuselage = Fuselage.ImportFromCell(readcell('DataDump.xlsx','Sheet','Fuselage'));
    Fuselage.Units = 'FPS';
    Fuselage = Fuselage.ConvertUnits(unitSystem);


    SS2 = WingGeometry;
    SS2.S = 47.63;
    SS2.Units = 'SI';
    SS2 = SS2.ConvertUnits(unitSystem);

    %% Drag Coefficient Calculations for sub-sonic wings and bodies
    
    CLw = @(AoA) 6.4247*(deg2rad(AoA)+0.0363);

    CDMainWing = WingDrag(MainWing,Mach,rho,mu,V,CLw,AoA);

    CDFuselage = BodyDrag(Fuselage,rho,V,mu);
    
    CLHT = @(AoA) FWCLa(Mach,VT.AR,VT.Sweep)*deg2rad(AoA);
    CDHT = WingDrag(HT,Mach,rho,mu,V,CLHT,AoA);

    CLVT = @(AoA) FWCLa(Mach,VT.AR,VT.Sweep)*deg2rad(AoA);
    CDVT = WingDrag(VT,Mach,rho,mu,V,CLVT,AoA);
    
    CDSS2 = 0.0314;
    Sreff = MainWing.S;
    
%     CD = CDMainWing*Sreff + 2*CDFuselage*Fuselage.Swet + CDHT*HT.S/Sreff + CDVT*VT.S/Sreff;
    CD = CDMainWing + 2*CDFuselage*Fuselage.Sb/Sreff + CDSS2*SS2.S/Sreff + CDHT*HT.S/Sreff + CDVT*VT.S/Sreff;

%     D = Cruise.AirSpeed.q*(CDMainWing*MainWing.S + 2*CDFuselage*Fuselage.Sb + CDSS2*SS2.S);
% 
%     L = Cruise.AirSpeed.q*CLw(AoA)*MainWing.S;


% 
%     if strcmp(unitSystem,'SI')
%         unit = 'Newtons';
%     elseif strcmp(unitSystem,'FPS')
%         unit = 'Pounds';
%     end
% 
%     fprintf('WhiteKnight Two drag during M = %.3f @ AoA = %d is %.2f %s\n',M,AoA,D,unit)
%     fprintf('WhiteKnight Two lift during M = %.3f @ AoA = %d is %.2f %s\n',M,AoA,L,unit)


    %% Data Export

%     FC = Cruise.Data2Cell();
%     MW = MainWing.Data2Cell();
%     F = Fuselage.Data2Cell();

%     writecell(FC,'MiguelsPeen.xlsx','Sheet','AerodynamicResults')
%     writecell(MW,'MiguelsPeen.xlsx','Sheet','MainWing')
%     writecell(F,'MiguelsPeen.xlsx','Sheet','Fuselage')

    %% Drag Functions
    function CDWing = WingDrag(Wing,Mach,rho,mu,V,CLw,AoA)
        tc = Wing.RootAirfoil.t;
        xc = Wing.RootAirfoil.p;
        Swet = Wing.Swet;
        Sref = Wing.S;
        sweep = Wing.Sweep;
        cmac = Wing.MAC;
        AR = Wing.AR;
        e = 0.7;
        Cf = frictionCoefficient(rho,V,cmac,mu);
        L = thiccParameter(xc);
        R = LiftingSurfaceCorrelationFactor(sweep,Mach);
        CD0Wing = Cf*(1 + L*tc + 100*tc^4 ) * R * Swet / Sref;
        CDLWing = CLw(AoA)^2/(pi*AR*e);
        CDWing = CD0Wing + CDLWing;
    end

    function R = LiftingSurfaceCorrelationFactor(sweep,Mach)
        % Until the process is developed from DATC0M, 1.146 will be used. This
        % assumes Mach = 0.55 and sweep is 0 degrees.
        peTransonic = @(M) abs(0.6-M)/M;
        peSubSonic = @(M) abs(0.25-M)/M;

        if peTransonic(Mach) > peSubSonic(Mach)
            R = 1.065;
        else
            R = 1.146;
        end
    end

    function L = thiccParameter(xc)
        if xc >= 0.3
            L  = 1.2;
        elseif xc <= 0.3
            L = 2.0;
        end
    end

    function CD0Body = BodyDrag(body,rho,V,mu)
        Cf = frictionCoefficient(rho,V,body.MaxDiameter,mu);
        CDf = Cf*(1 + 60/(body.FinenessRatio)^3 + 0.0025*(body.FinenessRatio) ) * body.Swet/body.Sb;
        CDb = 0.029 * (body.BaseDiameter/body.MaxDiameter)^3/sqrt(CDf);
        CD0Body = CDf + CDb;
    end

    function Cf = frictionCoefficient(rho,V,cmac,mu)
        Re = rho*V*cmac/mu;
        if Re < 5*10^5
            Cf = 1.328/sqrt(Re);
        elseif Re > 5*10^5
            Cf = 0.455/((log10(Re))^2.58);
        end
    end
end


function [Altitude,RS,Mach,V,Weight,q,n,Tmax,t,delV] = PablosWholeAssCode()
    %%Inputs everyhting is in english uunits 
    %Inputs/Initialize Variables
    g0=32.1741; %ft/s^2
    T0=518.67; %R
    rho0=0.0023769; %slug/ft^3
    k=1.4;
    Re=6378E3;
    Re=Re*3.28;
    R=1717;
    i=0;
    EheightM=[50:50:400].*10^3; %%% dont need to do 
    Vmax=100000;
    S= 1315;  %%%  
    W0=64992.265;  %%% Inital wight
    Th=7809.9*4; % lbf
    Isp=344; %%%

    PM=0;

    %%%
    %Data Henry 


    %Optimal Trajectory for a given Eheight calculation
    h2=[0,5,10,15,20,25,30,35,40,45,50]*10^3;
    can=0;

    for H=1:length(h2)  %H1=0:.1:350 for H=1:length(h2)-40
        h=h2(1,H);          %h=H1*10^3;    %h=h2(1,H);  %H=H+1;
        for Mach=0.1:0.001:0.8
            can=can+1;
            M=Mach;
            if h<36089
                theta=(1-h/145442);
                T=theta*T0;
                sigma=(1-h/145442)^(4.2561);
                rho=sigma*rho0;
            elseif h<=65617
                theta=(0.751865);
                T=theta*T0;
                sigma=0.297076*exp(-(h-36089)/20806);
                rho=sigma*rho0;
            elseif h<=104987
                theta=(0.682457+h/945374);
                T=theta*T0;
                sigma=(0.978261+h/659515)^(-35.16320);
                rho=sigma*rho0;  
            elseif h<=154199
                theta=0.482561+h/337634;
                T=theta*T0;
                sigma=(0.857003+h/190115)^(-13.20114);
                rho=sigma*rho0;         
            elseif h<=167232
                theta=0.939268;
                T=theta*T0;
                sigma=0.00116533*exp(-(h-154199)/25992);
                rho=sigma*rho0;
            elseif h<=232940
                theta=1.434843-h/337634;
                T=theta*T0;
                sigma=(0.798990-h/606330)^(11.20114);
                rho=sigma*rho0;
            elseif h<=278.386
                theta=1.237723-h/472687;
                T=theta*T0;
                sigma=(0.900194-h/649922)^(16.08160);
                rho=sigma*rho0;
            end
            a=sqrt(k*R*T);
            V=M*a;
            q=0.5*rho*V^2;

           % Aerodynamics data 
           Cdo = 0.266;

           % rocket propulsion
           H1=h/10^3;
    %        Isp_rocket=interIsp(h); % might be H1=h/10^3. Ambiguous in Chudoba's work
           hblah = ceil(H/3.28);
           [trust,mdot] = Thrust(hblah,Mach,'Imperial');
           ST = trust*4/mdot;
           if H==1
             W=W0;
             n=1;
             delV=0;
             t=0;
           else
             %HB=HB+1;  
             delV=V-Data2(H-1,4);
             %g=g0*(Re/(Re+h))^2;
             W=(Data2(H-1,5))/exp(delV/(g0*ST));
             Ps=(V*(Th-D))/Data2(H-1,5);
             t=((h-Data2(H-1,1)*10^3)/Ps)+Data2(H-1,9);
             %W=W0;
             %W=((Data2(H-1,5)))/(exp(delV/(g0*Isp)));
             %W=(Th/Isp_rocket)*((Data2(H-1,5))/(V*(Th-D)))*(h-Data2(H-1,1)*10^3);
             %W=2*Data2(H-1,5)-Data2(H-1,5)*exp(-(delV)/(Isp_rocket*g0));
             %W(i,j)=W(i,j-1)/(exp(V(i,j)-V(i,j-1)))/(32.2*Isp_rocket(j));
             %W=Data2(H-1,5)/(exp(delV))/(g0*Isp_rocket);
           end

           D=Cdo*q*S;
            RS=(ST*V*(Th-D))/(Th*W);      %this eqn with w0 works best
            %RS=(g0*Isp_rocket*(D-Th)*V)/(Th);
            %RS=1/RS;
            %mdot=Th/(Isp*g0);
            %RS=(V*(Th-D))/mdot;
            n=(Th-D-sind(10))/W;       %this eqn with W0 gives better ans       % OR 
    %         n=((Th-D)/W)-sind(10);           % OR 
            %n=((Th-D)/W);
           %{
            if H==1
                n=1;
            end
            %}
            Tmax=T+0.2*T*M^2;
            AllData(can,:)=[RS M V W q n Tmax t delV h];
            if Mach < 0.8
                if q <=300 && n>0 && n<4 && delV>=0 Tmax<1500; %&& H1 <= 50 %300 -1000 -4
                    %Change q to <=300
                    PM=1+PM;
                    MidData(PM,:)=[RS M V W q n Tmax t delV h];
                    %{
                    elseif q <=300 && Tmax<1000 && n>0 && n<4 && delV>0 && H1 > 50 && M>Data2(H-1,3)
                    PM=1+PM;
                    MidData(PM,:)=[RS M V W q n Tmax t delV h];
                    [Value Loc]=max(MidData(:,1)); %add stop
                    RSbest=MidData(Loc,1);
                    Mbest=MidData(Loc,2);
                    Vbest=MidData(Loc,3);
                    Wbest=MidData(Loc,4);
                    qbest=MidData(Loc,5);
                    nbest=MidData(Loc,6);
                    Tmaxbest=MidData(Loc,7);
                    tbest=MidData(Loc,8);
                    delVbest=MidData(Loc,9);
                    Data2(H,:)=[h/10^3 RSbest Mbest Vbest Wbest qbest nbest Tmaxbest tbest delVbest];
                    vpa(Data2)
                    PM=0;
                    can=0;
                    clear MidData AllData
                    break
                    %}
                else
                    continue
                end
            else   %add stop
                if q <=300 && n>0 && n<4 && delV>=0 Tmax<1500;
                   %Change q to <=300
                   PM=1+PM;
                   MidData(PM,:)=[RS M V W q n Tmax t delV h];
                   vpa(MidData);
                   [Value Loc]=max(MidData(:,1)); %add stop
                   RSbest=MidData(Loc,1);
                   Mbest=MidData(Loc,2);
                   Vbest=MidData(Loc,3);
                   Wbest=MidData(Loc,4);
                   qbest=MidData(Loc,5);
                   nbest=MidData(Loc,6);
                   Tmaxbest=MidData(Loc,7);
                   tbest=MidData(Loc,8);
                   delVbest=MidData(Loc,9);
                   Data2(H,:)=[h RSbest Mbest Vbest Wbest qbest nbest Tmaxbest tbest delVbest];
                   PM=0;
                end
                vpa(MidData);
                [Value Loc]=max(MidData(:,1)); %add stop
                RSbest=MidData(Loc,1);
                Mbest=MidData(Loc,2);
                Vbest=MidData(Loc,3);
                Wbest=MidData(Loc,4);
                qbest=MidData(Loc,5);
                nbest=MidData(Loc,6);
                Tmaxbest=MidData(Loc,7);
                tbest=MidData(Loc,8);
                delVbest=MidData(Loc,9);
                Data2(H,:)=[h/10^3 RSbest Mbest Vbest Wbest qbest nbest Tmaxbest tbest delVbest];
                vpa(Data2);
                PM=0;
                can=0;
                clear MidData AllData
                % if H==11;
                % delV=Data2(11,4)-Data2(1,4);
                % W=W/exp(delV/(g0*Isp));
                % end
            end
        end
    end

    figure(1)
    plot(Data2(:,3),Data2(:,1),'-b','LineWidth',2)
    ylabel('Altitude, ft*10^3')
    xlabel('Mach')


    figure(2)
    plot(Data2(:,4),Data2(:,1),'-b','LineWidth',2)
    ylabel('Altitude, ft*10^3')
    xlabel('Velocity, ft/s')

    figure(3)
    RC= Data2(:,4)*sind(10)
    plot(RC,Data2(:,1),'-b','LineWidth',2)
    ylabel('Altitude, ft*10^3')
    xlabel('Rate of Climb, ft/s')


    figure(4)
    plot(Data2(:,5),Data2(:,1),'-b','LineWidth',2)
    ylabel('Altitude, ft*10^3')
    xlabel('Weight, lbf')

    Altitude = Data2(:,1);
    RS=Data2(:,2);
    Mach=Data2(:,3);
    V=Data2(:,4);
    Weight=Data2(:,5);
    q=Data2(:,6);
    n=Data2(:,7);
    Tmax=Data2(:,8);
    t=Data2(:,9);
    delV=Data2(:,10);


%     ExportTable = table(Altitude,RS,Mach,V,Weight,q,n,Tmax,t,delV);
%     writetable(ExportTable,'MiguelsPeen.xlsx','Sheet','PablosVelocityShit')


    function Ispout = interIsp(h)
        RH=[0,921,2500,5000,10000,15000,20000,25000,30000,35000,36089,40000,...
            45000,50000,55000,60000,65000,65617,70000,75000,80000,85000,90000,...
            91000,96000,101000,104987,111000,121000,131000,140000,146000,150000,...
            154199,160000,170000,179200,180000,190000,200000,200131,225000,250000,...
            275000,300000,310000,320000,330000,340000,350000,375000,400000];

        RIsp_LNG_LOX=[127.2,132.7,141.8,155.6,180.4,202.1,220.9,237.4,251.6,264,...
            266.5,274.6,283.4,290.9,297.1,302.5,307.1,307.6,311,314.5,317.5,...
            320.2,322.7,323.1,325.3,327.3,328.7,330.6,333.3,335.7,337.6,338.8,...
            339.5,340.3,341.2,342.7,343.9,344,345.2,346.4,346.4,348.8,351.4,...
            353.5,355.1,355.4,356.1,356.3,356.6,356.6,356.6,356.6];

         Ispout = interp1(RH,RIsp_LNG_LOX,h,'spline','extrap');    
    end
end

function CLa = FWCLa(M0,AR,delta)
    B = sqrt(1 - M0^2);
    CLa = 2*pi*AR/(2+sqrt(4+AR^2*B^2*(1+(tand(delta)^2/B^2)) ));
end