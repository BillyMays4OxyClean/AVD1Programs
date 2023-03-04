%take-off
mach1=0.222;
altitude1= 0;

%Before drop
mach2=0.55;
altitude2= 50000;

slope=(altitude2-altitude1)/(mach2-mach1);

% for example
y1 = altitude1; % y coordinate of base point
x1 = mach1;% x coordinate of base point
m = slope; % slope
% define some x values between zero and 10
x = mach1:0.01:mach2
% evalualte the point slope formulae
y = m*(x-x1) + y1;
% plot the fitted line
plot(x,y,'ob',x,y)