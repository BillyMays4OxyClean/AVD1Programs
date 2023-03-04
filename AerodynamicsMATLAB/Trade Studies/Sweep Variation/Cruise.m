%% Created by Luke Patterson for the purpose of conducting a trade study during cruise.
clear
close all
clc

%% Unit Definition
unitSystem = 'FPS';

%% Geometry Definitions
MainWing = WingGeometry;
MainWing = MainWing.ImportFromCell(readcell('WK2Geometry.xlsx','Sheet','MainWing'));
MainWing = MainWing.ConvertUnits(unitSystem);

HT = WingGeometry;
HT = HT.ImportFromCell(readcell('WK2Geometry.xlsx','Sheet','HT'));
HT = HT.ConvertUnits(unitSystem);

VT = WingGeometry;
VT = VT.ImportFromCell(readcell('WK2Geometry.xlsx','Sheet','VT'));
VT = VT.ConvertUnits(unitSystem);

Fuselage = FuselageGeometry;
Fuselage = Fuselage.ImportFromCell(readcell('WK2Geometry.xlsx','Sheet','Fuselage'));
Fuselage = Fuselage.ConvertUnits(unitSystem);

SS2 = WingGeometry;
SS2.S = 47.63;
SS2.AR = 1.46;
SS2.Sweep = 45;
SS2.Rc = 28.2;
SS2.RootAirfoil = NACA4WingSection('0006',SS2.Rc,128);
SS2.Units = 'SI';
SS2 = SS2.ConvertUnits(unitSystem);
SS2ExportCell = SS2.Data2Cell();

%% Flight Condition
% WK2 Cruises at 15200 meters at Mach 0.55
% The flight condition class definition will contain and calculate all
% values relating to atmospheric properties and airspeed

AoA_i = 0;
AoA_inc = 1;
AoA_f = 15;
AoA = AoA_i:AoA_inc:AoA_f;

Cruz = FlightCondition;
Cruz.Altitude = 50000; % Altitude in feet
Cruz.Units = 'FPS';
Cruz.Name = 'Cruise';

Cruz = Cruz.SetSpeed('Mach',0.55);
Cruz = Cruz.ConvertUnits(unitSystem);

%% Aerodynamic Calculations
supress = 0;
n = 1;

CDMW = zeros(length(AoA),1);
CD0MW = zeros(length(AoA),1);

CD0 = zeros(length(AoA),1);
CD = zeros(length(AoA),1);

Clmw = zeros(length(AoA),1);
CMy = zeros(length(AoA),1);
NP = zeros(length(AoA),1);

CL = zeros(length(AoA),1);

D = zeros(length(AoA),1);
L = zeros(length(AoA),1);

Cm = zeros(length(AoA),1);
xac = zeros(length(AoA),1);
b = 1;
for Sweep = 0:15:60
    MainWing.Sweep = Sweep;
    n = 1;
    uaR(b) = Sweep;
    for AoA = AoA_i:AoA_inc:AoA_f

        Cruz.AoA = AoA;

        [CDMW(n,b),CD0MW(n,b),CDLMW,Dmw] = MainWing.ReturnDrag(Cruz);

        [CDHT,CD0HT,CDLHT,Dht] = HT.ReturnDrag(Cruz);

        [CDVT,CD0VT,CDLVT,Dvt] = VT.ReturnDrag(Cruz);

        [CDfuse,CD0fuse,Dfuse] = Fuselage.ReturnDrag(Cruz);

        CD0SS2 = 0.0314;

        [CLSS2,~,LSS2] = SS2.ReturnLift(Cruz);

        CDLSS2 = CLSS2^2/(pi*SS2.AR*SS2.RootAirfoil.e);

        CDSS2 = CD0SS2 + CDLSS2;

        DSS2 = Cruz.AirSpeed.q*CDSS2*SS2.S;

        CD0(n,b) = CD0MW(n,b) + 2*CD0fuse + 2*CD0HT + 2*CD0VT + CD0SS2;

        CD(n,b) =  CDMW(n,b) + 2*CDfuse + 2*CDVT + 2*CDHT + CDSS2; % + CDHT; ;

        D(n,b) = Dmw + 2*Dfuse + 2*Dht + 2*Dvt + DSS2;

        [Clmw(n,b),~,Lmw] = MainWing.ReturnLift(Cruz);

        [Clht,~,Lht] = HT.ReturnLift(Cruz);

        [Clvt,~,Lvt] = VT.ReturnLift(Cruz);

        CL(n,b) = Clmw(n,b) + CLSS2;

        L(n,b) = Lmw + LSS2;

        [Cm(n,b),xac(n,b)] = Cm0(Clmw(n,b),MainWing);

        if strcmp(unitSystem,'SI')
            unit = 'N';
        elseif strcmp(unitSystem,'FPS')
            unit = 'lbf';
        end

        if supress == 1
            fprintf('WhiteKnight Two drag during M = %.3f @ AoA = %.1f is %.2f %s\n',Cruz.AirSpeed.Mach,Cruz.AoA,D(n),unit)
            fprintf('WhiteKnight Two lift during M = %.3f @ AoA = %.1f is %.2f %s\n',Cruz.AirSpeed.Mach,Cruz.AoA,L(n),unit)
        end
        n = n + 1;
    end

    Ao = AoA_i:AoA_inc:AoA_f;
    b = b + 1;
end

figure()
grid on
hold on
Cdim = size(CD);
dragon = string;
for n=1:Cdim(2)
    plot(Ao,CL(:,n),'LineWidth',2)
    dragon(n) = num2str(uaR(n));
end
title('Total Vehicle C_L vs. AoA')
xlabel('Angle of Attack, degrees')
ylabel('Section lift coefficient, C_L')
l=legend(dragon,'Location','SouthEast');
title(l,'Wing Sweep, degrees')


figure()
grid on
hold on
for n=1:Cdim(2)
    plot(Ao,L(:,n)./D(:,n),'LineWidth',2)
    dragon(n) = num2str(uaR(n));
end
title('L/D vs. AoA')
xlabel('Angle of Attack, degrees')
ylabel('L/D')
l=legend(dragon,'Location','SouthEast');
title(l,'Wing Sweep, degrees')

%{
a1 = plot(CD(:,n),CL(:,n),'-Vk','LineWidth',2);

grid on
title({'Drag Polar';'';''})
xlabel('C_D')
ylabel('C_L')

xtickangle(45)

ax = gca;
ax.XAxis.Exponent = 0;

new_axis = axes('Position',ax.Position,'XAxisLocation','top','Color','none');
a2 = line(L./D,CL,'Parent',new_axis,'Color','b','LineWidth',2,'Marker','o');
xlabel('L/D')
set(new_axis,{'xcolor'},{'b'})
%}

function [Cm0, xac] = Cm0(CL,Wing)
    xac = Wing.AR/6 * (1+2*Wing.Taper)/(1+Wing.Taper)*tand(Wing.Sweep);
    Cm0 = Wing.RootAirfoil.Cmac - xac*CL;
end