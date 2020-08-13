% Copyright Ashish Tewari (c) 2006
global dtr;
global mu;
global omega;
global S;
global c;
global rm;
global b;
global CLmax;
global sweep;
global f8;

dtr = pi/180;
mu = 3.986004e14;
omega = 2*pi/(23*3600+56*60+4.0905);
S =  56.5;
c = 3.9926;
rm = 6378140;
b = 14.1512;
CLmax = 1.5;
sweep = 50*pi/180;
f8 = fopen('data8.mat', 'a');

long = 70*dtr;
lat = 20*dtr;
radint = rm+200;
velint = 69;
fpaint = 40*dtr;
chiint = 270*dtr;
m = 17350;

options = odeset('RelTol',1e-4);
init = [long; lat; radint; velint; fpaint; chiint; m];
[t, o] = ode45('accelclimb',[0, 80], init,options);
fclose('all');
