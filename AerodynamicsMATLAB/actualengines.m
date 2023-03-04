clc
clear all
close all

[ATMOm,ATMOe] = StandardATM(15240,'true');

Altitude=ATMOm(:,1);
P=ATMOm(:,3);
T=ATMOm(:,2);
p=ATMOm(:,4);
a0=ATMOm(:,5);

% T=217;
% p=1.125;
% a0=344;

%AE3007A1 Engine%
%inputs%
T0=216.7;
y=1.4;
cp=1004;
hpr=42800000;
Tt4=1800;
prc=24;
prf=1.8;
M0=.215:.05:.55;
a=5.3;
gc=1;
A=2.199; 

%Equations%
R=(y-1)/y*cp;
tr=1+(y-1)/2*M0.^2;
tv=Tt4/T0;
tc=prc^((y-1)/y);
tf=prf^((y-1)/y);
v9a0=sqrt(2/(y-1)*(tv-tr*(tc-1+a*(tf-1))-tv./(tr*tc)));
v9=v9a0.*a0;
v19a0=sqrt(2/(y-1)*(tr*tf-1));
v19=v19a0.*a0;
V=M0.*a0;
mair=p.*V.*A;
mc=mair/(1+a);
mf=mair-mc;


%thrust%
F=mc/gc.*(v9-V)+mf/gc.*(v19-V)/10;
Fmo=(a0./gc).*1./(1+a).*(v9./a0-M0+a*(v19./a0-M0));
f=cp*T0/hpr*(tv-tr*tc);
S=f./((1+a).*(Fmo));

%efficiency%
Nt=1-1./(tr*tc);
Np=(v9./V-1+a.*(v19./V-1))./((v9.^2)./V.^2-1+a.*((v19.^2)./V.^2-1));
No=Nt.*Np;

%Power%
tt=1-tr/tv*(tc-1+a*(tf-1));
Ve=(F+mair.*V)./(mair+mf);
W=1/2.*((mair+mf).*Ve.^2-mair.*V.^2);

%ALF502L Engine%
%inputs%
T0=216.7;
y=1.4;
cp=1004;
hpr=42800000;
Tt4=1800;
prc2=13.6;
prf2=1.74;
%M0=.215;
a2=5.2;
gc=1;
A2=1.27;

%Equations%
R=(y-1)/y*cp;
tr2=1+(y-1)/2*M0.^2;
tv2=Tt4/T0;
tc2=prc2^((y-1)/y);
tf2=prf2^((y-1)/y);
v9a02=sqrt(2/(y-1)*(tv2-tr2*(tc2-1+a2*(tf2-1))-tv2./(tr2*tc2)));
v92=v9a02.*a0;
v19a02=sqrt(2/(y-1)*(tr2*tf2-1));
v192=v19a02.*a0;
V2=M0.*a0;
mair2=p.*V2.*A2;
mc2=mair2/(1+a2);
mf2=mair2-mc2;


%thrust%
F2=mc2/gc.*(v92-V2)+mf2/gc.*(v192-V2)/10;
Fmo2=(a0./gc).*1./(1+a2).*(v92./a0-M0+a*(v192./a0-M0));
f2=cp*T0/hpr*(tv-tr*tc);
S2=f2./((1+a2).*(Fmo2));

%efficiency%
Nt2=1-1./(tr2*tc2);
Np2=(v92./V2-1+a2.*(v192./V2-1))./((v92.^2)./V2.^2-1+a.*((v192.^2)./V2.^2-1));
No2=Nt2.*Np2;

%Power%
tt2=1-tr2/tv2*(tc2-1+a2*(tf2-1));
Ve2=(F2+mair2.*V2)./(mair2+mf2);
W2=1/2.*((mair2+mf2).*Ve2.^2-mair2.*V2.^2);

% 

%PW308A%

%inputs%
T0=216.7;
y=1.4;
cp=1004;
hpr=42800000; %70000
Tt4=1800; %2400
prc3=18;
prf3=1.7;
%M0=.215;
a3=4.1;
gc=1;
A3=1.299;

%Equations%
R3=(y-1)/y*cp;
tr3=1+(y-1)/2*M0.^2;
tv3=Tt4/T0;
tc3=prc3^((y-1)/y);
tf3=prf3^((y-1)/y);
v9a03=sqrt(2/(y-1)*(tv3-tr3*(tc3-1+a3*(tf3-1))-tv3./(tr3*tc3)));
v93=v9a03.*a0;
v19a03=sqrt(2/(y-1).*(tr3*tf3-1));
v193=v19a03.*a0;
V3=M0.*a0;
mair3=p.*V3*A3;
mc3=mair3./(1+a3);
mf3=mair3-mc3;


%thrust%
F3=mc3/gc.*(v93-V3)+mf3/gc.*(v193-V3)/10;
Fmo3=(a0/gc).*1./(1+a3).*(v93./a0-M0+a3.*(v193./a0-M0));
f3=cp*T0/hpr*(tv3-tr3*tc3);
S3=f3./((1+a3).*(Fmo3));

%efficiency%
Nt3=1-1./(tr3*tc3);
Np3=2.*M0.*(v93./a0-M0+a3.*(v193./a0-M0))./(v93.^2./a0.^2-M0.^2+a3.*(v193.^2./a0.^2-M0.^2));
No3=Nt3.*Np3;

Ve3=(F3+mair3.*V3)./(mair3+mf3);
W3=1/2.*((mair3+mf3).*Ve3.^2-mair3.*V3.^2);

%F33=[transpose(F3(1:138,1)),transpose(F3(139:276,2)),transpose(F3(277:414,3)),transpose(F3(415:552,4)),transpose(F3(553:690,5)),transpose(F3(691:828,6)),transpose(F3(829:966,7)),transpose(F3(967:1104,8)),transpose(F3(1105:1242,9)),transpose(F3(1243:1380,10)),transpose(F3(1381:1520,11))];
% plot(Altitude,W3)
% title('PW308A Power at .55 Mach')
% xlabel('Altitude, meters')
% ylabel('Power, Watt')


%PW305B Engine%
%inputs%
T0=216.7;
y=1.4;
cp=1004;
hpr4=42800000; %30000
Tt4=1800; %1000
prc4=15.5;
prf4=1.8;
%M0=.215;
a4=4.3;
gc=1;
A4=.927; 

%Equations%
R=(y-1)/y*cp;
tr4=1+(y-1)/2*M0.^2;
tv4=Tt4/T0;
tc4=prc4^((y-1)/y);
tf4=prf4^((y-1)/y);
v9a04=sqrt(2/(y-1)*(tv4-tr4*(tc4-1+a4*(tf4-1))-tv4./(tr4*tc4)));
v94=v9a04.*a0;
v19a04=sqrt(2/(y-1)*(tr4*tf4-1));
v194=v19a04.*a0;
V4=M0.*a0;
mair4=p.*V.*A4;
mc4=mair4/(1+a4);
mf4=mair4-mc4;


%thrust%
F4=mc4/gc.*(v94-V4)+mf4/gc.*(v194-V4)/10;
Fmo4=(a0./gc).*1./(1+a4).*(v94./a0-M0+a4*(v194./a0-M0));
f4=cp*T0/hpr4*(tv4-tr4*tc4);
S4=f4./((1+a4).*(Fmo4));

%efficiency%
Nt4=1-1./(tr4*tc4);
Np4=(v94./V4-1+a4.*(v194./V4-1))./((v94.^2)./V4.^2-1+a4.*((v194.^2)./V4.^2-1));
No4=Nt4.*Np4;

%Power%
tt4=1-tr4/tv4*(tc4-1+a4*(tf4-1));
Ve4=(F4+mair4.*V4)./(mair4+mf4);
W4=1/2.*((mair4+mf4).*Ve4.^2-mair4.*V4.^2);
% %table%
% M0=M0';
% F=F';
% Fmo=Fmo';
% S=S';
% Nt=Nt';
% Np=Np';
% No=No';
% table(M0,F,Fmo,S,Nt,Np,No)

plot(Altitude,W3,'LineWidth',2)
title('Watts varying by Altitude at Mach .215')
xtickangle(30)
% legend('AE3007A1','ALF502L','PW308A','PW305B','Location','Southeast')
xlabel('Altitude, meters')
ylabel('Power, Watts')

surf(M0,Altitude,V3)