clear
close all
clc

PR = @(Rho,V,W) Rho*V*W;

pr = PR(0.008,200,10000)