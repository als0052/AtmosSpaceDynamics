% Copyright Ashish Tewari (c) 2006
global dtr;
global mu;
global omega;
global S;
global c;
global rm;
global CD0;
global K;
global b;
global CLmax;
global CLG;
global fT0;
global tsfc0;
global mur;


dtr = pi/180;
mu = 3.986004e14;
omega = 2*pi/(23*3600+56*60+4.0905);
S = 223.0815;
c = 5.42;
rm = 6378140;
CD0 = 0.02;
K = 0.055;
b = 41.2;
CLmax = 1.6;
CLG = 0.1;
fT0 = 211128.17;
tsfc0 = 0.4;
mur = 0.03;

long = 0*dtr;
lat = 51.45*dtr;
rad = rm+3;
vel = 0.001;
fpa = 0;
chi = 270*dtr;
m = 84890.909;
init = [long; lat; rad; vel; fpa; chi; m];
[t, o] = ode45('takeoff',[0, 60], init);
fclose('all');
