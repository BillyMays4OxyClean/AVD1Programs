clear
close all
clc


V=[1:1:200];
rho=1.225;
W=29500;
%Cdo=0.0061;
Cdo=0.226
e=0.82;
Sw=122.165;
AR=15.1365;
Cdol=0;
Q=0.108

%Landing Weight
WF=9381


%Take-off
PR=(1/2)*rho*V.^3*Sw*Cdo+((W*9.81)^2./((1/2)*rho.*V*Sw))*(1/(pi*AR*e));
ThruR=(((1/2)*(rho).*(V).^2.*(Cdo))/(((W*9.81)/Sw)))+Cdol+(((W*9.81)/Sw)./((1/2)*pi*e*AR*rho.*V.^2)).*(W*(9.81));

%Climb
W10=287150.5/9.81;
W20=283912.19/9.81;
W30=280291.34/9.81;
%W30=100000.00
W40=275011.3/9.81;
W50=267440.43/9.81;
rho10=0.9046;
rho20=0.6527;
rho30=0.4583;
rho40=0.3016;
rho50=0.1865;
Cdooo=0.266;

%PR3=(1/2)*rho30*V.^3*Sw*Cdo+((W30*9.81)^2./((1/2)*rho30.*V*Sw))*(1/(pi*AR*e));

PR10=(((Cdooo)*(rho10).*(V).^3)./(2*(W10*(9.81)/Sw)))+Cdol.*V+((2*(W10*9.81)/Sw)./(pi*e*AR*rho10.*V))*(W10*(9.81));
PR20=(((Cdooo)*(rho20).*(V).^3)./(2*(W20*(9.81)/Sw)))+Cdol.*V+((2*(W20*9.81)/Sw)./(pi*e*AR*rho20.*V))*(W20*(9.81));
PR30=(((Cdooo)*(rho30).*(V).^3)./(2*(W30*(9.81)/Sw)))+Cdol.*V+((2*(W30*9.81)/Sw)./(pi*e*AR*rho30.*V))*(W30*(9.81));
PR40=(((Cdooo)*(rho40).*(V).^3)./(2*(W40*(9.81)/Sw)))+Cdol.*V+((2*(W40*9.81)/Sw)./(pi*e*AR*rho40.*V))*(W40*(9.81));
PR50=(((Cdooo)*(rho50).*(V).^3)./(2*(W50*(9.81)/Sw)))+Cdol.*V+((2*(W50*9.81)/Sw)./(pi*e*AR*rho50.*V))*(W50*(9.81));

TR10=(((1/2)*(rho10).*(V).^2.*(Cdooo))/((W10*(9.81)/Sw)))+Cdol+(((W10*9.81)/Sw)./((1/2)*pi*e*AR*rho10.*V.^2)).*(W10*(9.81));
TR20=(((1/2)*(rho20).*(V).^2.*(Cdooo))/((W20*(9.81)/Sw)))+Cdol+(((W20*9.81)/Sw)./((1/2)*pi*e*AR*rho20.*V.^2)).*(W20*(9.81));
TR30=(((1/2)*(rho30).*(V).^2.*(Cdooo))/((W30*(9.81)/Sw)))+Cdol+(((W30*9.81)/Sw)./((1/2)*pi*e*AR*rho30.*V.^2)).*(W30*(9.81));
TR40=(((1/2)*(rho40).*(V).^2.*(Cdooo))/((W40*(9.81)/Sw)))+Cdol+(((W40*9.81)/Sw)./((1/2)*pi*e*AR*rho40.*V.^2)).*(W40*(9.81));
TR50=(((1/2)*(rho50).*(V).^2.*(Cdooo))/((W50*(9.81)/Sw)))+Cdol+(((W50*9.81)/Sw)./((1/2)*pi*e*AR*rho50.*V.^2)).*(W50*(9.81));

%Fuel consumption Q dot the totalof fuel consumed per unit time
%ovee the power avaible.

%Trying to find the chnage in fuel consumption


qpp10=Q./PR10;
qpp20=Q./PR20;
qpp30=Q./PR30;
qpp40=Q./PR40;
qpp50=Q./PR50;


qpt10=Q./TR10;
qpt20=Q./TR20;
qpt30=Q./TR30;
qpt40=Q./TR40;
qpt50=Q./TR50;


SFC0=2.76E-04
2.63E-04
SFC10=2.44E-04
2.20E-04
SFC20=2.04E-04
1.84E-04
SFC30=1.51E-04
1.39E-04
SFC40=1.23E-04
1.22E-04
SCF50=1.10E-04



FuelC10=SFC10*TR10
FuelC20=SFC20*TR20
FuelC30=SFC30*TR30
FuelC40=SFC30*TR40
FuelC50=SFC30*TR50







%%Cruise
% Cdoo=0.023;
Cdoo=0.226

% rho2=1.056;
rho2=0.1865;
% V2=183.92;
V2=162;
% W2=[92027.61:196.2:99483.21];

W2=[99483.21:-196.2:92027.61];

L=137824.9192

n=L./W2;
V4=[188.65:-1.32:137.2];
phi=acos(1./n);
Radius= V4./(n.*9.81.*sin(phi));




TurnRate=(9.81.*sqrt(n.^2-1))./V4;


W3=99483.21;

WC=99483.21;
WCC=95755.41;
WCCC=92027.61;

%Power Avaiable
PA=(((Cdoo)*(rho2)*(V2)^3)./(2.*(W2.*(9.81)/Sw)))+Cdol*V2+((2.*(W2.*9.81)/Sw)/(pi*e*AR*rho2*V2)).*(W2.*(9.81));
%Thrust Avaiable
TR=(((1/2)*(rho2)*(V2)^2*(Cdoo))./((W2.*(9.81)/Sw)))+Cdol+(((W2.*9.81)/Sw)/((1/2)*pi*e*AR*rho2*V2^2)).*(W2.*(9.81));






TRC=(((1/2)*(rho2)*(V).^2*(Cdoo))./((WC*(9.81)/Sw)))+Cdol+(((WC*9.81)/Sw)./((1/2)*pi*e*AR*rho2*V.^2))*(WC*(9.81));
TRCC=(((1/2)*(rho2)*(V).^2*(Cdoo))./((WCC*(9.81)/Sw)))+Cdol+(((WCC*9.81)/Sw)./((1/2)*pi*e*AR*rho2*V.^2))*(WCC*(9.81));
TRCCC=(((1/2)*(rho2)*(V).^2*(Cdoo))./((WCCC*(9.81)/Sw)))+Cdol+(((WCCC*9.81)/Sw)./((1/2)*pi*e*AR*rho2*V.^2))*(WCC*(9.81));
qdot=1.86*10^-5;

FC=qdot*TRC;
FCC=qdot*TRCC;
FCCC=qdot*TRCCC;

%Thrust and Power SPecific Fuel Consumption
V3=[20:1:100];
PAC=(((Cdoo)*(rho2)*(V3).^3)/(2*(W3*(9.81)/Sw)))+Cdol.*V3+((2*(W3*9.81)/Sw)./(pi*e*AR*rho2.*V3)).*(W3*(9.81));
qp=Q./PAC;

TRC=(((1/2)*(rho2)*(V3).^2*(Cdoo))/((W3*(9.81)/Sw)))+Cdol+(((W3*9.81)/Sw)./((1/2)*pi*e*AR*rho2*V3.^2))*(W3*(9.81));
qp2=Q./TRC;




%Decent




%Landing 
PRF=(1/2)*rho*V.^3*Sw*Cdo+((WF*9.81)^2./((1/2)*rho.*V*Sw))*(1/(pi*AR*e));
ThruRF=(((1/2)*(rho).*(V).^2.*(Cdo))/((WF*(9.81)/Sw)))+Cdol+(((WF*9.81)/Sw)./((1/2)*pi*e*AR*rho.*V.^2)).*(WF*(9.81));




%rate of Climb



%%







%Take-off Power Required 
figure(1)
plot(V,PR)
xlim([50 200])
xlabel('Airspeed(m/s)')
ylabel('Power Required (Watts)')
title ('Power Req at Take-off')

%Take-off Thrust required
figure(2)
plot(V,ThruR)
xlabel('Airspeed(m/s)')
ylabel('Thrust Required (N)')
title ('Thrust Required at Take-off')
xlim([50 110])

%Climb Power Required
figure(3)
plot(V,PR10) 
hold on 
plot(V,PR20)
hold on 
plot(V,PR30)
hold on 
plot(V,PR40)
hold on 
plot(V,PR50)
hold off
xlim([40 150])
legend('10,000 ft', '20,000ft','30,000ft','40,000ft','50,000ft');
xlabel('Airspeed(m/s)')
ylabel('Power Required(Watts)')
title ('Power Required at Climb')


%Climb Thrust Required
figure(4)
plot(V,TR10) 
hold on 
plot(V,TR20)
hold on 
plot(V,TR30)
hold on 
plot(V,TR40)
hold on 
plot(V,TR50)
hold off
% xlim([0 80])
xlim([40 150])
legend('10,000 ft', '20,000ft','30,000ft','40,000ft','50,000ft');
xlabel('Airspeed(m/s)')
ylabel('Thrust Avalaible(N)')
title ('Thrust Required at Climb')

%Cruis both avaible and require power
figure(5)
plot(W2,PA)
xlabel('Weight(N)')
ylabel('Power Avaiable (Watts)')
title ('Power Avaiable at Cruise')

%Cruis both avaible and require thrust
figure(6)
plot(W2,TR)
xlabel('Weight(N)')
ylabel('Thrust Avaiable (N)')
title ('Thrust Avaiable at Cruise')

figure(7)
plot(V3,qp)
xlabel('Airspeed (m/s)')
ylabel('Power-Specific Fuel Consumption (kg/(s*N))')
title ('Power-Specific Fuel Consumptionat Cruise')


figure(8)
plot(V3,qp2)
xlabel('Airspeed (m/s)')
ylabel('Thrust-Specific Fuel Consumption (kg/(s*N))')
title ('Thrust-Specific Fuel Consumption Cruise')


figure(9)
plot(V,qpp10) 
hold on 
plot(V,qpp20)
hold on 
plot(V,qpp30)
hold on 
plot(V,qpp40)
hold on 
plot(V,qpp50)
hold off
xlabel('Airspeed (m/s)')
ylabel('Power-Specific Fuel Consumption (kg/(s*N))')
title ('Power-Specific Fuel Consumption Climb')
xlim([40 150])
legend('10,000 ft', '20,000ft','30,000ft','40,000ft','50,000ft','Location', 'Best');


figure(10)
plot(V,qpt10) 
hold on 
plot(V,qpt20)
hold on 
plot(V,qpt30)
hold on 
plot(V,qpt40)
hold on 
plot(V,qpt50)
hold off
xlabel('Airspeed (m/s)')
ylabel('Thrust-Specific Fuel Consumption (kg/(s*N))')
title ('Thrust-Specific Fuel Consumptionat Climb')
xlim([20 200])
legend('10,000 ft', '20,000ft','30,000ft','40,000ft','50,000ft','Location', 'Best');

figure(11)
plot(V,FC)
hold on 
plot(V,FCC)
hold on 
plot(V,FCCC)
hold off
xlabel('Airspeed (m/s)')
ylabel('Fuel Consumption (kg/s)/N')
title ('Fuel Consumption at Climb')
legend('10,000 ft', '20,000ft','30,000ft','40,000ft','50,000ft','Location', 'Best');
xlim([100 200])



figure(12)
plot(V,PRF)
xlabel('Airspeed (m/s)')
ylabel('Power Required (Watts)')
title ('Power Required at Landing')
xlim([20 150])

figure(13)
plot(V,ThruRF)
xlabel('Airspeed (m/s)')
ylabel('Thrust Required (N)')
title ('Thrust Required at Landing')
xlim([20 150])


% figure(14)
% plot(V,PR)
% xlim([50 200])
% xlabel('Airspeed(m/s)')
% ylabel('Power Required (Watts)')
% title ('Power Req at Take-off')
% yyaxis right
% plot(V,ThruR)
% ylabel('Thrust Required (N)')



figure(14)
plot(W2,n)
xlabel('Weight (N)')
ylabel('Load Factor')
title ('Load factor during Cruise')


figure(15)
plot(V4,Radius)
xlabel('Velocity of Aircraft (m/s)');
ylabel('Turning Radius (m)');
%title('
xlim([138 190])


figure(16)
plot(V4,TurnRate)
xlabel('Velocity of Aircraft (m/s)');
ylabel('Turn Rate (degrees/sec)');
%title('
xlim([138 190])


