%% Testing
clear
close all
clc

nseg = 2^5;

NACA = NACA4WingSection('2412',1,nseg);
NACA.ExportToDat('NACA2412');
% NACA.PlotFoil();
n=1;
for x = 0:1/(nseg-2):1
    dydxn(n,1) = NACA.MeanSlope(x);
    n = n + 1;
end

importFoil = GenericAirfoil;
% importFoil.xc = NACA.xc;
% importFoil.yc = NACA.yc;
% importFoil.t = NACA.t;
% importFoil.xu = NACA.xu;
% importFoil.yu = NACA.yu;
% importFoil.xl = NACA.xl;
% importFoil.yl = NACA.yl;
importFoil = importFoil.ImportFromDatFile('NACA2412.dat');
hold on
% importFoil.PlotFoil()
dycdx = importFoil.CamberLineSlope();

pe = abs(dydxn-dycdx)./dydxn * 100;

pe_s = abs(NACA.yc-importFoil.yc)./NACA.yc * 100;

table(dydxn,dycdx,pe)
table(NACA.yc,importFoil.yc,pe_s)
paoa = abs(NACA.a0l-importFoil.a0l)/NACA.a0l * 100

plot(NACA.xc,NACA.yc,'-b',importFoil.xu,importFoil.yc,'-r')
figure()
plot(NACA.xu(1:(length(NACA.xu)-1),1),dydxn(1:length(dydxn),1),'-b',importFoil.xu(1:(length(importFoil.xu)-1),1),dycdx,'-r')