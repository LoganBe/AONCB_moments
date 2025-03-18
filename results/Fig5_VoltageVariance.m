clear; clc; close all;
rng(6895)
set(0,'DefaultAxesFontSize',14)

% Figure 5
% Uses Monte Carlo simulation to show the effects of input correlation on
% subthreshold variability
% First uses large weights for asynchronous/within-pool correlations/across-pool correlations
% Second uses small weights for asynchronous/within-pool correlations/across-pool correlations

%% Time Info and constants
dt = 0.1; %Time step (ms)
tmax = 1e3; %Max time (ms)
ttrim = 200; %Time to remove so in asymtotic state (ms)
t = 0:dt:tmax+ttrim; %time vector (ms)
tt = t(1:end-ttrim/dt-1)*1e-3; %Time Vector Trimed (s)
trimidx = binsearch(t,ttrim); % idx of ttrim

%Time Constants
tauleak = 15; % Leak constant (ms)
tausyn = 2*dt; % Syanptic time constant (ms)

%Biophysical Paramters (Normalized to 0)
V0 = 0; % Rest voltage
Vth = 15; %Threshold voltage 
Ve_r = 60; %Excitatory reversal potential
Vi_r = -10; %Inhibitory reversal potential

% Input Paramters 
N = 1;

% Save parameters into a structure to load into model
theta.tauleak = tauleak;
theta.tausyn = tausyn;
theta.trim = trimidx;
theta.ereversal = Ve_r;
theta.ireversal = Vi_r;
theta.lreversal = 0;

%% Example 1: Uncorr - Large Weights
% Input parameters - Small weight
Ke = 100;  We = 0.15; 
Ki = 25; Wi = 0.6;  

theta.corrinfo.corr = "uncorr";
rbank = [1,10,30]; % Rates to test
V = zeros(length(tt),length(rbank));
for i = 1:length(rbank)
    theta.K.ee = Ke; theta.r.ee = rbank(i);
    theta.K.ei = Ki; theta.r.ei = rbank(i);

    [ye, yi, yre, yri] = makeffinput(theta,N,0,t,dt);
    wb_fe = ye.*We./tausyn; % Weighted Excitatory input
    wb_fi = yi.*Wi./tausyn; % Weighted Inhibitory input

    input.feedfowarde = wb_fe;
    input.feedfowardi = wb_fi;
    V(:,i) = cifv(t, N, 0, input, theta);
end

% Theoretical Grid using small weight approximation
theta.Wf = [We,Wi];
rbank = 0:0.1:50;
theta.r.ee = rbank; theta.r.ei = rbank; 
theta.rhoWe = 0; theta.rhoWi = 0; theta.rhoWei = 0;
[meanVt,varVt,skewVt] = swstats(theta);

% Plot Results
figure(1); clf;
plot(tt,V,'linewidth',2)
box off
ylim([-4,30])
xlabel('Time (s)'); ylabel('Voltage (mV)')

figure(2); clf; hold on
imagesc(rbank,rbank,varVt')
plot(rbank,rbank,'k','linewidth',2)
xlim([0,50]); ylim([0,50])
xlabel('r_i (Hz)'); ylabel('r_e (Hz)')
box off
colorbar
set(gca,'Ydir','normal')
%% Example 2: Uncorr - Moderate Weights
% Input parameters - Small weight
Ke = 1000; Ki = 250; 
We = 0.015; Wi = 0.06;  

theta.corrinfo.corr = "uncorr";
rbank = [1,10,30]; % Rates to test
V = zeros(length(tt),length(rbank));
for i = 1:length(rbank)
    theta.K.ee = Ke; theta.r.ee = rbank(i);
    theta.K.ei = Ki; theta.r.ei = rbank(i);

    [ye, yi, yre, yri] = makeffinput(theta,N,0,t,dt);
    wb_fe = ye.*We./tausyn; % Weighted Excitatory input
    wb_fi = yi.*Wi./tausyn; % Weighted Inhibitory input

    input.feedfowarde = wb_fe;
    input.feedfowardi = wb_fi;
    V(:,i) = cifv(t, N, 0, input, theta);
end

% Theoretical Grid
theta.Wf = [We,Wi];
rbank = 0:0.1:50;
theta.r.ee = rbank; theta.r.ei = rbank; 
theta.rhoWe = 0; theta.rhoWi = 0; theta.rhoWei = 0;
[~,varVt] = swstats(theta);

% Plot Results
figure(3); clf;
plot(tt,V,'linewidth',2)
box off
ylim([-4,30])
xlabel('Time (s)'); ylabel('Voltage (mV)')

figure(4); clf; hold on
imagesc(rbank,rbank,varVt')
plot(rbank,rbank,'k','linewidth',2)
xlim([0,50]); ylim([0,50])
xlabel('r_i (Hz)'); ylabel('r_e (Hz)')
box off
colorbar
set(gca,'Ydir','normal')

%% Example 3: Corr W - Large Weights
% Input parameters - Small weight
Ke = 100; Ki = 25; 
We = 0.15; Wi = 0.6;  

theta.corrinfo.corr = "corr";
rho_e = 0.03; rho_i = 0.03; 
theta.corrinfo.corridx = 0; 
rho_ei = 0; 
theta.corrinfo.ae1 = 0; theta.corrinfo.ae2 = (1/rho_e)-theta.corrinfo.ae1-1;
theta.corrinfo.ai1 = 0; theta.corrinfo.ai2 = (1/rho_i)-theta.corrinfo.ai1-1;

rbank = [1,10,30]; % Rates to test
V = zeros(length(tt),length(rbank));
for i = 1:length(rbank)
    theta.K.ee = Ke; theta.r.ee = rbank(i);
    theta.K.ei = Ki; theta.r.ei = rbank(i);

    [ye, yi, yre, yri] = makeffinput(theta,N,0,t,dt);
    wb_fe = ye.*We./tausyn; % Weighted Excitatory input
    wb_fi = yi.*Wi./tausyn; % Weighted Inhibitory input

    input.feedfowarde = wb_fe;
    input.feedfowardi = wb_fi;
    V(:,i) = cifv(t, N, 0, input, theta);
end

% Theoretical Grid
theta.Wf = [We,Wi];
rbank = 0:0.1:50;
theta.r.ee = rbank; theta.r.ei = rbank; 
theta.rhoWe = rho_e; theta.rhoWi = rho_i; theta.rhoWei = 0;
[~,varVt] = swstats(theta);

% Plot Results
figure(5); clf;
plot(tt,V,'linewidth',2)
box off
ylim([-4,30])
xlabel('Time (s)'); ylabel('Voltage (mV)')

figure(6); clf; hold on
imagesc(rbank,rbank,varVt')
plot(rbank,rbank,'k','linewidth',2)
xlim([0,50]); ylim([0,50])
xlabel('r_i (Hz)'); ylabel('r_e (Hz)')
box off
colorbar
set(gca,'Ydir','normal')

%% Example 4: Corr W - Moderate Weights
% Input parameters - Small weight
Ke = 1000; Ki = 250; 
We = 0.015; Wi = 0.06;  

theta.corrinfo.corr = "corr";
rho_e = 0.03; rho_i = 0.03; 
theta.corrinfo.corridx = 0; 
rho_ei = 0; 
theta.corrinfo.ae1 = 0; theta.corrinfo.ae2 = (1/rho_e)-theta.corrinfo.ae1-1;
theta.corrinfo.ai1 = 0; theta.corrinfo.ai2 = (1/rho_i)-theta.corrinfo.ai1-1;

rbank = [1,10,30]; % Rates to test
V = zeros(length(tt),length(rbank));
for i = 1:length(rbank)
    theta.K.ee = Ke; theta.r.ee = rbank(i);
    theta.K.ei = Ki; theta.r.ei = rbank(i);

    [ye, yi, yre, yri] = makeffinput(theta,N,0,t,dt);
    wb_fe = ye.*We./tausyn; % Weighted Excitatory input
    wb_fi = yi.*Wi./tausyn; % Weighted Inhibitory input

    input.feedfowarde = wb_fe;
    input.feedfowardi = wb_fi;
    V(:,i) = cifv(t, N, 0, input, theta);
end

% Theoretical Grid
theta.Wf = [We,Wi];
rbank = 0:0.1:50;
theta.r.ee = rbank; theta.r.ei = rbank; 
theta.rhoWe = rho_e; theta.rhoWi = rho_i; theta.rhoWei = rho_ei;
[~,varVt] = swstats(theta);

% Plot Results
figure(7); clf;
plot(tt,V,'linewidth',2)
box off
ylim([-4,30])
xlabel('Time (s)'); ylabel('Voltage (mV)')

figure(8); clf; hold on
imagesc(rbank,rbank,varVt')
plot(rbank,rbank,'k','linewidth',2)
xlim([0,50]); ylim([0,50])
xlabel('r_i (Hz)'); ylabel('r_e (Hz)')
box off
colorbar
set(gca,'Ydir','normal')

%% Example 5: Corr X - Large Weights
% Input parameters - Small weight
Ke = 100; Ki = 25; 
We = 0.15; Wi = 0.6;  

theta.corrinfo.corr = "corr";
rho_e = 0.03; rho_i = 0.03; 
theta.corrinfo.corridx = 1; 
rho_ei = 0.03; 
theta.corrinfo.ae1 = 0; theta.corrinfo.ae2 = (1/rho_e)-theta.corrinfo.ae1-1;
theta.corrinfo.ai1 = 0; theta.corrinfo.ai2 = (1/rho_i)-theta.corrinfo.ai1-1;

rbank = [1,10,30]; % Rates to test
V = zeros(length(tt),length(rbank));
for i = 1:length(rbank)
    theta.K.ee = Ke; theta.r.ee = rbank(i);
    theta.K.ei = Ki; theta.r.ei = rbank(i);

    [ye, yi, yre, yri] = makeffinput(theta,N,0,t,dt);
    wb_fe = ye.*We./tausyn; % Weighted Excitatory input
    wb_fi = yi.*Wi./tausyn; % Weighted Inhibitory input

    input.feedfowarde = wb_fe;
    input.feedfowardi = wb_fi;
    V(:,i) = cifv(t, N, 0, input, theta);
end

% Theoretical Grid
theta.Wf = [We,Wi];
rbank = 0:0.1:50;
theta.r.ee = rbank; theta.r.ei = rbank; 
theta.rhoWe = rho_e; theta.rhoWi = rho_i; theta.rhoWei = rho_ei;
[~,varVt] = swstats(theta);

% Plot Results
figure(9); clf;
plot(tt,V,'linewidth',2)
box off
ylim([-4,30])
xlabel('Time (s)'); ylabel('Voltage (mV)')

figure(10); clf; hold on
imagesc(rbank,rbank,varVt')
plot(rbank,rbank,'k','linewidth',2)
xlim([0,50]); ylim([0,50])
xlabel('r_i (Hz)'); ylabel('r_e (Hz)')
box off
colorbar
set(gca,'Ydir','normal')

%% Example 6: Corr X - Moderate Weights
% Input parameters - Small weight
Ke = 1000; Ki = 250; 
We = 0.015; Wi = 0.06;  

theta.corrinfo.corr = "corr";
rho_e = 0.03; rho_i = 0.03; 
theta.corrinfo.corridx = 1; 
rho_ei = 0.03; 
theta.corrinfo.ae1 = 0; theta.corrinfo.ae2 = (1/rho_e)-theta.corrinfo.ae1-1;
theta.corrinfo.ai1 = 0; theta.corrinfo.ai2 = (1/rho_i)-theta.corrinfo.ai1-1;

rbank = [1,10,30]; % Rates to test
V = zeros(length(tt),length(rbank));
for i = 1:length(rbank)
    theta.K.ee = Ke; theta.r.ee = rbank(i);
    theta.K.ei = Ki; theta.r.ei = rbank(i);

    [ye, yi, yre, yri] = makeffinput(theta,N,0,t,dt);
    wb_fe = ye.*We./tausyn; % Weighted Excitatory input
    wb_fi = yi.*Wi./tausyn; % Weighted Inhibitory input

    input.feedfowarde = wb_fe;
    input.feedfowardi = wb_fi;
    V(:,i) = cifv(t, N, 0, input, theta);
end

% Theoretical Grid
theta.Wf = [We,Wi];
rbank = 0:0.1:50;
theta.r.ee = rbank; theta.r.ei = rbank; 
theta.rhoWe = rho_e; theta.rhoWi = rho_i; theta.rhoWei = rho_ei;
[~,varVt] = swstats(theta);

% Plot Results
figure(11); clf;
plot(tt,V,'linewidth',2)
box off
ylim([-4,30])
xlabel('Time (s)'); ylabel('Voltage (mV)')

figure(12); clf; hold on
imagesc(rbank,rbank,varVt')
plot(rbank,rbank,'k','linewidth',2)
xlim([0,50]); ylim([0,50])
xlabel('r_i (Hz)'); ylabel('r_e (Hz)')
box off
colorbar
set(gca,'Ydir','normal')




