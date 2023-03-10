%% Created by Luke Patterson for the analysis of the Boeing 747 during cruise
clear
close all
clc

%% Unit Definition
unitSystem = 'FPS';

%% Geometry Definitions for the Boeing 747
MainWing = WingGeometry;
MainWing.S = 550.36;
MainWing.b = 195.66929;
MainWing.Rc = 47.24409;
MainWing.Tc = 13.32021;

bac = GenericAirfoil;
bac = bac.ImportFromDatFile('bacxxx.txt');
bac.PlotFoil();
figure()

MainWing.RootAirfoil = bac;
MainWing.TipAirfoil = bac;
MainWing.Taper = 0.2819;
MainWing.AR = 6.9607;
MainWing.MAC = 33.4491;
MainWing.meanThicc = 0.1981;
MainWing.Sweep = 38.8596;
MainWing.Swet = 11332.944;
MainWing.Se = 5399.1;
MainWing.Kw = 2.0604;
MainWing.V = 1950.2;
MainWing.Units = 'FPS';
MainWing = MainWing.ConvertUnits(unitSystem);
exp2Excel = MainWing.Data2Cell();
writecell(exp2Excel,'B747Geometry.xlsx','Sheet','MainWing');

Fuselage = FuselageGeometry;
Fuselage.SlendernessRatio = 0.0479;
Fuselage.FinenessRatio = 10.5548;
Fuselage.e = 1.5072;
Fuselage.MaxDiameter = 21.333;
Fuselage.BaseDiameter = 7.5;
Fuselage.Sb = pi/4 * Fuselage.MaxDiameter^2;
Fuselage.L = 225.166;
Fuselage.W = 21.333;
Fuselage.H = 32.1522;
Fuselage.Swet = 14787.8356;
Fuselage.V = 65231;
Fuselage.Units = 'FPS';
Fuselage = Fuselage.ConvertUnits(unitSystem);
exp2Excel = MainWing.Data2Cell();
writecell(exp2Excel,'B747Geometry.xlsx','Sheet','Fuselage');

HT = WingGeometry;
HT.b = 72.73;
HT.Rc = 31.56168;
HT.Tc = 8.39895;
HT.RootAirfoil = NACA4WingSection('0012',HT.Rc,128);
HT.TipAirfoil = NACA4WingSection('0012',HT.Tc,128);
HT.AR = 3.6389;
HT.Sweep = 34.2025;
HT.MAC = 22.218;
HT.Swet = 3110.3;
HT.Se = 1493.4;
HT.S = 1453;
HT.Units = 'FPS';
HT = HT.ConvertUnits(unitSystem);
exp2Excel = HT.Data2Cell();
writecell(exp2Excel,'B747Geometry.xlsx','Sheet','HT');

VT = WingGeometry;
VT.S = 1051.957;
VT.b = 38.05774;
VT.Rc = 42.6509;
VT.Tc = 12.63123;
VT.RootAirfoil = NACA4WingSection('0012',VT.Rc,128);
VT.TipAirfoil = NACA4WingSection('0012',VT.Tc,128);
VT.Taper = 0.2961;
VT.AR = 1.3768;
VT.Sweep = 45;
VT.meanThicc = 0.2375;
VT.MAC = 30.358;
VT.Swet = 2228.8;
VT.Se = 1051.957;
VT.Units = 'FPS';
VT = VT.ConvertUnits(unitSystem);
exp2Excel = MainWing.Data2Cell();
writecell(exp2Excel,'B747Geometry.xlsx','Sheet','VT');

%% Flight Condition
AoA_f = 15;
AoA = 0:AoA_f;

begin = FlightCondition;
begin.Altitude = 35105; % Altitude in feet, condition for cruise is 35kft.
begin.Units = 'FPS';
begin = begin.SetSpeed('V',850.7); % Takeoff velocity in ft/s
begin = begin.ConvertUnits(unitSystem);
begin.Name = 'Takeoff';

%% Aerodynamic Calculations
supress = 1;
n = 1;

V = zeros(length(AoA),1);

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

for A = 1:length(AoA)

    begin.AoA = AoA(A);
    
    [CDMW(n,1),CD0MW(n,1),CDLMW,Dmw] = MainWing.ReturnDrag(begin);
    [CDHT,CD0HT,CDLHT,Dht] = HT.ReturnDrag(begin);
    [CDVT,CD0VT,CDLVT,Dvt] = VT.ReturnDrag(begin);
    [CDfuse,CD0fuse,CDLfuse,Dfuse] = Fuselage.ReturnDrag(begin);
    
    CD0(n,1) = CD0MW(n,1) + CD0fuse + CD0HT + CDVT ;
    CD(n,1) =  CDMW(n,1) + CDfuse + CDHT + CDVT;
    D(n,1) = Dmw + Dfuse + Dht + Dvt;

    [Clmw(n,1),~,Lmw] = MainWing.ReturnLift(begin);
    [Clht,~,Lht] = HT.ReturnLift(begin);
    [Clvt,~,Lvt] = VT.ReturnLift(begin);
    CL(n,1) = Clmw(n,1);
    L(n,1) = Lmw;

    [Cm(n,1),xac(n,1)] = Cm0(Clmw(n,1),MainWing);
    
    if strcmp(unitSystem,'SI')
        unit = 'N';
    elseif strcmp(unitSystem,'FPS')
        unit = 'lbf';
    end
    
    if supress == 1
        fprintf('Boeing 747 drag during M = %.3f @ AoA = %.1f is %.2f %s\n',begin.AirSpeed.Mach,begin.AoA,D(n),unit)
        fprintf('Boeing 747 lift during M = %.3f @ AoA = %.1f is %.2f %s\n',begin.AirSpeed.Mach,begin.AoA,L(n),unit)
    end
    n = n + 1;
end

plot(AoA,CL,'-b','LineWidth',2)
grid on
title('C_L vs. AoA')
xlabel('Angle of Attack, degrees')
ylabel('Section lift coefficient, C_L')
figure()

plot(AoA,CD,'-k','LineWidth',2)
grid on
title('C_D vs. AoA')
xlabel('Angle of Attack, degrees')
ylabel('Section drag coefficient, C_D')

figure()
plot(AoA,L,'-b','LineWidth',2)
ax = gca;
ax.YAxis.Exponent = 0;
grid on
title('Total Vehicle Lift vs. AoA')
xlabel('Angle of Attack, degrees')
ylabel(['Tandem Vehicle Lift, ',unit])

figure()
plot(AoA,D,'-k','LineWidth',2)
ax = gca;
ax.YAxis.Exponent = 0;
grid on
title('Total Vehicle Drag vs. AoA')
xlabel('Angle of Attack, degrees')
ylabel(['Tandem Vehicle Drag, ',unit])

figure()
plot(CD,CL,'-b','LineWidth',2)
xlabel('C_D')
ylabel('C_L')
title('C_L vs. C_D')

figure()
plot(AoA,L./D,'-b','LineWidth',2)
xlabel('Angle of Attack, degrees')
ylabel('L/D')
title('L/D vs. AoA')
Drag = D;
Lift = L;

figure()
plot(AoA,Cm,'--k','LineWidth',2)
xlabel('Angle of Attack, degrees')
ylabel('Moment Coefficient')
title('C_M vs. AoA')

figure()
a1 = plot(CD,CL,'-Vk','LineWidth',2);

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



function [Cm0, xac] = Cm0(CL,Wing)
    xac = Wing.AR/6 * (1+2*Wing.Taper)/(1+Wing.Taper)*tand(Wing.Sweep);
    Cm0 = Wing.RootAirfoil.Cmac - (xac*CL);
end