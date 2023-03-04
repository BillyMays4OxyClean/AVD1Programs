clear
close all
clc

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

%% Geometry Definitions
unitSystem = 'FPS';

MainWing = WingGeometry;
MainWing = MainWing.ImportFromCell(readcell('DataDump.xlsx','Sheet','MainWing'));
MainWing = MainWing.ConvertUnits(unitSystem);

HT = WingGeometry;
HT = HT.ImportFromCell(readcell('DataDump.xlsx','Sheet','HT'));
HT = HT.ConvertUnits(unitSystem);

VT = WingGeometry;
VT = VT.ImportFromCell(readcell('DataDump.xlsx','Sheet','VT'));
VT = VT.ConvertUnits(unitSystem);

Fuselage = FuselageGeometry;
Fuselage = Fuselage.ImportFromCell(readcell('DataDump.xlsx','Sheet','Fuselage'));
Fuselage = Fuselage.ConvertUnits(unitSystem);

SS2 = WingGeometry;
SS2.S = 47.63;
SS2.AR = 1.46;
SS2.Sweep = 45;
SS2.Rc = 28.2;
SS2.RootAirfoil = NACA4WingSection('0006',SS2.Rc,128);
SS2.Units = 'SI';
SS2 = SS2.ConvertUnits(unitSystem);

%Optimal Trajectory for a given Eheight calculation
h2=[0:5:50]*10^3;
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
        
        

        Cruise = FlightCondition;
        Cruise.Altitude = h; % Altitude in feet
        Cruise.Units = unitSystem;
        Cruise.Name = 'Climb';
        Cruise.AoA = 7.5;

        Cruise = Cruise.SetSpeed('Mach',M);
        Cruise = Cruise.ConvertUnits(unitSystem);


        %% Drag Coefficient Calculations for sub-sonic wings and bodies
        
        [CDMW,Dmw] = MainWing.ReturnDrag(Cruise);

        [CDHT,Dht] = HT.ReturnDrag(Cruise);

        [CDVT,Dvt] = VT.ReturnDrag(Cruise);

        [CDfuse,Dfuse] = Fuselage.ReturnDrag(Cruise);

        CD0SS2 = 0.0314;

        [CLSS2,LSS2] = SS2.ReturnLift(Cruise);

        CDLSS2 = CLSS2^2/(pi*SS2.AR*SS2.RootAirfoil.e);

        CDSS2 = CD0SS2 + CDLSS2;

        DSS2 = Cruise.AirSpeed.q*CDSS2*SS2.S;

        CD =  CDMW + 2*CDfuse + 2*CDHT + 2*CDVT + CDSS2;

        D = Dmw + 2*Dfuse + 2*Dht + 2*Dvt + DSS2;

        [Clmw,Lmw] = MainWing.ReturnLift(Cruise);

        [Clht,Lht] = HT.ReturnLift(Cruise);

        [Clvt,Lvt] = VT.ReturnLift(Cruise);

        CL = Clmw + CLSS2;

        L = Lmw + LSS2;
        
        
       % Aerodynamics data 
       Cdo = 0.266;
    
       % rocket propulsion
       H1=h/10^3;
%        Isp_rocket=interIsp(h); % might be H1=h/10^3. Ambiguous in Chudoba's work
       hblah = ceil(H/3.28);
       [trust,mdot] = Thrust(hblah,Mach,'Imperial');
       ST = trust/mdot;
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
        n=L/W;       %this eqn with W0 gives better ans       % OR 
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
RC= Data2(:,4)*sind(10);
plot(RC,Data2(:,1),'-b','LineWidth',2)
ylabel('Altitude, ft*10^3')
xlabel('Rate of Climb, ft/s')


figure(4)
plot(Data2(:,5),Data2(:,1),'-b','LineWidth',2)
ylabel('Altitude, ft*10^3')
xlabel('Weight, lbf')
xtickangle(35)
ax = gca;
ax.XAxis.Exponent = 0;

Altitude = Data2(:,1);
RS=Data2(:,2);
Mach=Data2(:,3);
Velocity=Data2(:,4);
Weight=Data2(:,5);
q=Data2(:,6);
n=Data2(:,7);
Tmax=Data2(:,8);
t=Data2(:,9);
delV=Data2(:,10);


ExportTable = table(Altitude,RS,Mach,Velocity,Weight,q,n,Tmax,t,delV);
writetable(ExportTable,'VelocityProfileExport.xlsx')

save("PablosVelocityProfile.mat","Altitude","Velocity","Mach","q")

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