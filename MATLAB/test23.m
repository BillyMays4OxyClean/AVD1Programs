clear
close all
clc

alt = 15200;

[ATMOm,ATMOe] = StandardATM(alt);

rho = ATMOm(:,4);