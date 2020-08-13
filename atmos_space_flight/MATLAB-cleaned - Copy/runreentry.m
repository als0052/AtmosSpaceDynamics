% Copyright Ashish Tewari (c) 2006
global dtr;
global mu;
global S;
global c;
global m;
global rm;
global omega;
global Gamma;
global f8;

dtr = pi/180;
mu = 3.986004e14;
S = 4;
c = 0.5;
m = 350;
rm = 6378140;
omega = 2*pi/(23*3600+56*60+4.0905);
Gamma = 1.41;
f8 = fopen('data8.mat', 'a');

long = -10*dtr;
lat = -79.8489182889*dtr;
rad = 6579.89967e3;
vel = 7589.30433867;
fpa = 0.54681217*dtr;
chi = 99.955734*dtr;
options = odeset('RelTol', 1e-8);
orbinit = [long; lat; rad; vel; fpa; chi];
[t, o] = ode45('reentry',[0, 1750],orbinit,options);
fclose('all');
