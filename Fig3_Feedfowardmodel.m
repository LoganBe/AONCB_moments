clear; clc; close all;
rng(6895)
set(0,'DefaultAxesFontSize',14)

% DESCRIPTION HERE

%% Time Info and constants
dt = 0.1; %Time step (ms)
tmax = 2e3; %Max time (ms)
ttrim = 100; %Time to remove to be asymptotic (ms) 
t = 0:dt:tmax+ttrim; %time vector
tt = t(1:end-ttrim/dt-1)*1e-3; %Time Vector Trimed in in Seconds
trimidx = binsearch(t,ttrim); %trim time idx

%Time Constants
tauleak = 15; % Leak Time Constant 
tausyn = 2*dt; % Synaptic Time Constant 

%Biophysical Paramters (Normalized to Leak = 0)
V0 = 0; % Rest voltage
Ve_r = 60; %Excitatory reversal potential
Vi_r = -10; %Inhibitory reversal potential

% FF weights
We = 0.015; Wi = 0.06;

% Input correlations (large for show)
rhoe = 0.1; rhoi = 0.1;

% Number of neurons 
Ne = 3; Ni = 0;
N = Ni+Ne;

% Set Information into Structure
theta.tauleak = tauleak;
theta.tausyn = tausyn;
theta.trim = trimidx;
theta.ereversal = Ve_r;
theta.ireversal = Vi_r;
theta.lreversal = 0;

theta.corrinfo.corr = "corr";
theta.corrinfo.corridx = 0;
theta.corrinfo.ae1 = 0; theta.corrinfo.ae2 = (1/rhoe)-theta.corrinfo.ae1-1; % Convert corr to beta (dist paramter)
theta.corrinfo.ai1 = 0; theta.corrinfo.ai2 = (1/rhoi)-theta.corrinfo.ai1-1; % Convert corr to beta (dist parameter)
beta_e = theta.corrinfo.ae2; beta_i = theta.corrinfo.ai2;

%% Presynaptic Inputs
% Excitatory Inputs
Ke = 1000; re = 10; %Hz 
theta.K.ee = Ke; theta.r.ee = re; 

% Inhibitory Inputs
Ki = 250; ri = 10; %Hz
theta.K.ei = Ki; theta.r.ei = ri;

% Make inputs using Copula for 3 populations 
d = 2*N;  corrmat = zeros(d);
corrmat = corrmat + diag(diag(ones(d)));
corrmat(corrmat == 0) = 0.6; % Arbitary large cross neuron corr [note 0.6 ~= rho = 0.6]

% MAKE INPUTS
[ys,yraw] = makecopinput(theta,t,dt,d/2,d/2,corrmat);
ye1 = ys{1}; ye2 = ys{2}; ye3 = ys{3}; % Excitatory Inputs
yi1 = ys{4}; yi2 = ys{5}; yi3 = ys{6}; % Inhibitory Inputs

% Filter inputs with box car of length tausyn ms
yfilt = ones(floor(tausyn/dt),1); 
yefsamp1 = filter(yfilt,1,ye1); yefsamp2 = filter(yfilt,1,ye2); yefsamp3 = filter(yfilt,1,ye3); 
yifsamp1 = filter(yfilt,1,yi1); yifsamp2 = filter(yfilt,1,yi2); yifsamp3 = filter(yfilt,1,yi3); 

% Weighted input
wb_fe1 = yefsamp1.*We/tausyn; wb_fe2 = yefsamp2*We/tausyn; wb_fe3 = yefsamp3*We/tausyn;
wb_fi1 = yifsamp1.*Wi/tausyn; wb_fi2 = yifsamp2*Wi/tausyn; wb_fi3 = yifsamp3*Wi/tausyn;

% Put input together to make raster
inplot_e = [yraw{1}(:,1:350), yraw{2}(:,1:350), yraw{3}(:,1:350)];
idx = randperm(size(inplot_e,2)); inplot_e = inplot_e(:,idx);
inplot_i = [yraw{4}(:,1:100), yraw{5}(:,1:100), yraw{6}(:,1:100)];
idx = randperm(size(inplot_i,2)); inplot_i = inplot_i(:,idx);

%SIM V
input.feedfowarde = wb_fe1;
input.feedfowardi = wb_fi1;
V1 = cifv(t,Ne,Ni,input,theta);

input.feedfowarde = wb_fe2;
input.feedfowardi = wb_fi2;
V2 = cifv(t,Ne,Ni,input,theta);

input.feedfowarde = wb_fe3;
input.feedfowardi = wb_fi3;
V3 = cifv(t,Ne,Ni,input,theta);

%% Plot Results
figure(1); clf; hold on
raster_plot(inplot_i,dt*1e-3,'b');
raster_plot(inplot_e,dt*1e-3,'r',100*3);
xlim([0,2]); ylim([0,350*3])

figure(2); clf; 
plot(tt,V1,'k','linewidth',2); box off
xlabel('Time (s)'); ylabel('Voltage (mV)')

figure(3); clf; 
plot(tt,V2,'g','linewidth',2); box off
xlabel('Time (s)'); ylabel('Voltage (mV)')

figure(4); clf; 
plot(tt,V3,'m','linewidth',2); box off
xlabel('Time (s)'); ylabel('Voltage (mV)')


