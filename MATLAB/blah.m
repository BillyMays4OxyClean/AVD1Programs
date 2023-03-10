clear
close all
clc

t = 0.5*10^-2;
h = 15*10^-2;
b = 20*10^-2;
Vz = 100*1000;

A1 = 5*10^-2;
A4 = A1;
A2 = 7.5*10^-2;
A3 = A2;

% Iy = t*(2*h)^3/12 + 2*(b*t^3/12 + b*t*h^2)
% 
% q1 = Vz*t*b*h/Iy
% 
% q2 = q1 + Vz*t*h^2/(2*Iy)

q1 = -Vz*A1/(2*h*(A1+A2))

q2 = -Vz/(2*h)

q3 = - Vz*A1/(2*h*(A1 + A2))