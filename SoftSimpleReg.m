winstyle = 'docked';
% winstyle = 'normal';

set(0,'DefaultFigureWindowStyle',winstyle) %sets default figure properties
set(0,'defaultaxesfontsize',18)
set(0,'defaultaxesfontname','Times New Roman')
% set(0,'defaultfigurecolor',[1 1 1])

% clear VARIABLES;
clear
global spatialFactor;
global c_eps_0 c_mu_0 c_c c_eta_0
global simulationStopTimes;
global AsymForcing
global dels
global SurfHxLeft SurfHyLeft SurfEzLeft SurfHxRight SurfHyRight SurfEzRight



dels = 0.75;
spatialFactor = 1;

c_c = 299792458;                  % speed of light
c_eps_0 = 8.8542149e-12;          % vacuum permittivity
c_mu_0 = 1.2566370614e-6;         % vacuum permeability
c_eta_0 = sqrt(c_mu_0/c_eps_0);


tSim = 200e-15 %total time the simulation is set to run for
f = 230e12;  %frequency of 2300 THz
lambda = c_c/f;

xMax{1} = 20e-6; %simulation domain (20 microns)
nx{1} = 200; % x grid points
ny{1} = 0.75*nx{1}; % y grid points


Reg.n = 1;

mu{1} = ones(nx{1},ny{1})*c_mu_0;

epi{1} = ones(nx{1},ny{1})*c_eps_0;
%epi{1}(125:150,55:95)= c_eps_0*11.3; %this is the permitity being set to ne different in certain region

%epi{1}(75:100,55:95)= c_eps_0*11.3; %modified 
%epi{1}(130:170,55:95)= c_eps_0*11.3; %modified
%epi{1}(65:75,55:95)= c_eps_0*11.3; %modified
%epi{1}(35:55,55:95)= c_eps_0*11.3; %modified

epi{1}(125:150,55:95)= c_eps_0*11.3;
epi{1}(75:100,20:95)= c_eps_0*11.3; %modified 
epi{1}(100:170,20:95)= c_eps_0*11.3; %modified
epi{1}(170:75,20:95)= c_eps_0*11.3; %modified
epi{1}(75:55,20:95)= c_eps_0*11.3; %modified





sigma{1} = zeros(nx{1},ny{1});
sigmaH{1} = zeros(nx{1},ny{1});

dx = xMax{1}/nx{1};
dt = 0.25*dx/c_c;
nSteps = round(tSim/dt*2);
yMax = ny{1}*dx;
nsteps_lamda = lambda/dx

movie = 1;
Plot.off = 0;
Plot.pl = 0;
Plot.ori = '13';
Plot.N = 100;
Plot.MaxEz = 1.1;
Plot.MaxH = Plot.MaxEz/c_eta_0;
Plot.pv = [0 0 90];
Plot.reglim = [0 xMax{1} 0 yMax];


bc{1}.NumS = 1;
bc{1}.s(1).xpos = nx{1}/(4) + 1;
bc{1}.s(1).type = 'ss';
bc{1}.s(1).fct = @PlaneWaveBC;

bc{1}.NumS = 2;
bc{1}.s(2).xpos = nx{1}/(4) + 6;  %modified
bc{1}.s(2).type = 'ss';
bc{1}.s(2).fct = @PlaneWaveBC;
%mag = -1/c_eta_0;
mag = 1;
phi = 0;
omega = f*2*pi;
betap = 0;
t0 = 30e-15;
st = 15e-15;
s = 0;
y0 = yMax/2; % A plane wave is specified as a source at the mid-point of y-axis
sty = 1.5*lambda;
bc{1}.s(1).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'};  %modified
bc{1}.s(2).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'};

Plot.y0 = round(y0/dx);

bc{1}.xm.type = 'a'; %These are absorbing boundary conditions
bc{1}.xp.type = 'a';
bc{1}.ym.type = 'a';
bc{1}.yp.type = 'a';

pml.width = 20 * spatialFactor; %sets up perfectly matched layer for absorbing outgoing waves
pml.m = 3.5;

Reg.n  = 1;
Reg.xoff{1} = 0;
Reg.yoff{1} = 0;

RunYeeReg






