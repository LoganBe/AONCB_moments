clear; clc; close all
rng(6895)
set(0,'DefaultAxesFontSize',14)

% DESCRIPTION HERE
%% Time Info and constants
dt = 0.1; %Time step (ms)
tmax = 3e3; %Max time (ms)
ttrim = 100; %Time to remove (ms)
t = 0:dt:tmax+ttrim; %time vector
tt = t(1:end-ttrim/dt-1)*1e-3; %Time Vector Trimed in in Seconds
trimidx = binsearch(t,ttrim); %trim time idx

%Time Constants
tauleak = 15; %Leak Time constant 
tausyn = 2*dt; % Synaptic time constant

%Biophysical Paramters (Normalized to 0)
V0 = 0; % Rest voltage
Vth = 15; %Threshold voltage 
Ve = 60; %Excitatory reversal potential
Vi = -10; %Inhibitory reversal potential

% Number of neurons/trials
N = 100; 

% FF weights
We = 0.15; Wi = 0.6;
W_f = [We, Wi];

% Smallweight functions to estimate stats
meanVtf = @(ke,ki,re,ri) (ke.*re.*We.*Ve + ki.*ri.*Wi.*Vi)./(1+ke.*re.*We + ki.*ri.*Wi); % MEAN
varVtf = @(ke,ki,re,ri,mu) (ke.*re.*We.^2.*(Ve-mu).^2 + ki.*ri.*Wi.^2.*(Vi-mu).^2)./(2*tauleak.*(1+ke.*re.*We+ki.*ri.*Wi)); % VAR
q = @(ke,ki,re,ri,mu) ke.*re.*We.^2.*(Ve-mu).^2./(ke.*re.*We.^2.*(Ve-mu).^2 + ki.*ri.*Wi.^2.*(Vi-mu).^2); % Q
corrVtf = @(ke,ki,re,ri,mu,fe,fi) fe.*q(ke,ki,re,ri,mu) + fi.*(1-q(ke,ki,re,ri,mu)); % CORR(Q)

%% Excitatory| K_e vs F_e
Ke = 200; re = 10; % excitatory
ki_ = 25; ri = 10; f_i = 0.4; % inhibitory

rre = re*1e-3; rri = ri*1e-3; % rate in spikes/ms

% Grid over Ke and Fe
kbank = 1:Ke; fbank = linspace(0,1,1000);
[KK,FF] = meshgrid(kbank,fbank);
meanVgrid = meanVtf(KK,ki_,rre,rri);
corrVgrid = corrVtf(KK,ki_,rre,rri,meanVgrid,FF,f_i);

figure(1); clf;
imagesc(0:0.1:1,kbank,corrVgrid')
colorbar; ylim([1,200]); xlim([0,1]); clim([0,1])
set(gca,'Ydir','normal')
xlabel('F_e'); ylabel('K_e')

% Cross Section - Constant Ke
Ke = 100; fe = (0:Ke)./Ke;
meanVf = meanVtf(Ke,ki_,rre,rri);
varVf = varVtf(Ke,ki_,rre,rri,meanVf);
corrVf = corrVtf(Ke,ki_,rre,rri,meanVf,fe,f_i);

figure(2); clf; hold on
plot(fe*100,meanVf.*ones(length(fe),1),'k','linewidth',2); box off;
xlabel('F_e'); ylabel('\mu')

figure(3); clf; hold on
plot(fe*100,varVf.*ones(length(fe),1),'k','linewidth',2); box off;
xlabel('F_e'); ylabel('\sigma^2')

figure(4); clf; hold on
plot(fe*100,corrVf,'k','linewidth',2); box off;
ylim([0,1])
xlabel('F_e'); ylabel('\rho')

% Cross Section - Constant Fe
Ke = 0:200; fe = 0.75;
meanVf = meanVtf(Ke,ki_,rre,rri);
varVf = varVtf(Ke,ki_,rre,rri,meanVf);
corrVf = corrVtf(Ke,ki_,rre,rri,meanVf,fe,f_i);

figure(5); clf; hold on
plot(Ke,meanVf.*ones(length(fe),1),'k','linewidth',2); box off;
xlabel('K_e'); ylabel('\mu')

figure(6); clf; hold on
plot(Ke,varVf.*ones(length(fe),1),'k','linewidth',2); box off;
xlabel('K_e'); ylabel('\sigma^2')

figure(7); clf; hold on
plot(Ke,corrVf,'k','linewidth',2); box off;
ylim([0,1])
xlabel('K_e'); ylabel('\rho')

%% Inhibitory | K_i vs F_i
Ki = 50; ri = 10; % Inhibitory
ke_ = 100; re = 10; f_e = 0.5; % Excitatory

rre = re*1e-3; rri = ri*1e-3; % Rate in spikes/ms

% Grid over Ki and Fi
kbank = 1:Ki; fbank = linspace(0,1,1000);
[KK,FF] = meshgrid(kbank,fbank);
meanVgrid = meanVtf(ke_,KK,rre,rri);
corrVgrid = corrVtf(ke_,KK,rre,rri,meanVgrid,f_e,FF);

figure(8); clf;
imagesc(0:0.1:1,kbank,corrVgrid')
colorbar; ylim([1,50]); xlim([0,1]); clim([0,1])
set(gca,'Ydir','normal')
xlabel('F_i'); ylabel('K_i')

% Cross Section over Fi
Ki = 25; fi = (0:Ki)./Ki;
meanVf = meanVtf(ke_,Ki,rre,rri);
varVf = varVtf(ke_,Ki,rre,rri,meanVf);
corrVf = corrVtf(ke_,Ki,rre,rri,meanVf,f_e,fi);

figure(9); clf; hold on
plot(fi*100,meanVf.*ones(length(fi),1),'k','linewidth',2); box off;
xlabel('F_i'); ylabel('\mu')

figure(10); clf; hold on
plot(fi*100,varVf.*ones(length(fi),1),'k','linewidth',2); box off;
xlabel('F_i'); ylabel('\sigma^2')

figure(11); clf; hold on
plot(fi*100,corrVf,'k','linewidth',2); box off;
ylim([0,1])
xlabel('F_i'); ylabel('\rho')

% Cross Section over Ki
Ki = 0:50; fi = 0.75;
meanVf = meanVtf(ke_,Ki,rre,rri);
varVf = varVtf(ke_,Ki,rre,rri,meanVf);
corrVf = corrVtf(ke_,Ki,rre,rri,meanVf,f_e,fi);

figure(12); clf; hold on
plot(Ki,meanVf.*ones(length(fi),1),'k','linewidth',2); box off;
xlabel('K_i'); ylabel('\mu')

figure(13); clf; hold on
plot(Ki,varVf.*ones(length(fi),1),'k','linewidth',2); box off;
xlabel('K_i'); ylabel('\sigma^2')

figure(14); clf; hold on
plot(Ki,corrVf,'k','linewidth',2); box off;
ylim([0,1])
xlabel('K_i'); ylabel('\rho')

%% Change in E and I rate 
Ke = 100; f_e = 0.65; % Excitatory
Ki = 25; f_i = 0.25;  % Inhibitory
rbank = 0:0.1:50; % Input rate (joint over E and I)

[R1,R2] = meshgrid(rbank,rbank);
R1 = R1*1e-3; R2 = R2*1e-3; % Rate in (spikes/ms)
meanVgrid = meanVtf(Ke,Ki,R1,R2); % Mean
corrVgrid = corrVtf(Ke,Ki,R1,R2,meanVgrid,f_e,f_i); % Corr

% Grid over Re and Ri
figure(13); clf;
imagesc(rbank,rbank,corrVgrid')
colorbar; clim([0,1]); xlim([0,50]); ylim([0,50])
set(gca,'Ydir','normal')
xlabel('r_i'); ylabel('r_e')

% Cross section over joint R
rr = rbank*1e-3;
meanVf = meanVtf(Ke,Ki,rr,rr);
varVf = varVtf(Ke,Ki,rr,rr,meanVf);
corrVf = corrVtf(Ke,Ki,rr,rr,meanVf,f_e,f_i);

figure(13); clf; hold on
plot(rr*1e3,meanVf,'k','linewidth',2); box off;
xlabel('rate (Hz)'); ylabel('\mu')

figure(14); clf; hold on
plot(rr*1e3,varVf,'k','linewidth',2); box off;
xlabel('rate (Hz)'); ylabel('\sigma^2')

figure(15); clf; hold on
plot(rr*1e3,corrVf,'k','linewidth',2); box off;
ylim([0,1])
xlabel('rate (Hz)'); ylabel('\rho')
