%% Created by Luke Patterson for Senior Design Drag Calculations
clear
close all
clc

%% Flight Condition
A = 15200;

%% CD Input from VSPAero
[m,e] = StandardATM(A);
rho = e(1,4);
a = e(1,5);
M = 0.55;
V = M*a;
Sref = 1316/(3.28)^2; % m^2

CD =  [0.001720	0.000930	0.000340	0.000010	-0.000020	-0.000130	-0.000050	-0.000030	-0.000010	-0.000010	-0.000010	0.000000	0.000000	0.000000	0.000000	0.001720	0.000930	0.000340	0.000010	-0.000020	-0.000130	-0.000050	-0.000030	-0.000010	-0.000010	-0.000010	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000010	0.000010	0.000010	0.000000	0.000000	0.000010	0.000010	0.000010	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000];
%CD = [0.001270	0.000660	0.000230	-0.000010	-0.000020	-0.000080	-0.000030	-0.000010	-0.000010	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.001270	0.000660	0.000230	-0.000010	-0.000020	-0.000080	-0.000030	-0.000010	-0.000010	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000010	0.000000	0.000000	0.000000	0.000000	0.000010	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000];

CDTot = 0.022;

D = 1/2*rho*V^2*CDTot*Sref; % Drag in Newtons
%D = 3.28/2.2*D; % Drag in lbf

CD0 = 0.1;
SS2Sref = 47.63*(3.28)^2;
Dss2 = 1/2*rho*V^2*CD0*SS2Sref;
%Dss2 = 3.28/2.2*Dss2;
fprintf('WhiteKnight Two drag during cruise (M = %.3f) is %.2f Newtons\n',M,D)
fprintf('SpaceShip Two drag during cruise (M = %.3f) is %.2f Pounds\n',M,Dss2)