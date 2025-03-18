clear; clc; close all
rng(6895)
set(0,'DefaultAxesFontSize',14)

% Figure 8
% Examine the effects of input correlations on voltage skewness. 
% First shows example voltage traces for aysnch/synch case
% Second uses full moment to show theory and simulation agreement
% Third shows the effects of drive on skewness
% Lastly shows the effects of input correlations on skewness
%% Time Info and constants
dt = 0.1; %Time step (ms)
tmax = 5e3; %Max time (ms)
ttrim = 100; %Time to remove (ms)
t = 0:dt:tmax+ttrim; %time vector (ms)
tt = t(1:end-ttrim/dt-1)*1e-3; %Time Vector Trimed in in Seconds
trimidx = binsearch(t,ttrim); %trim time idx

%Time Constants
tauleak = 15;  %Leak time constant 
tausyn = 2*dt; %Synaptic Time constant 

%Biophysical Paramters (Normalized to 0)
V0 = 0; % Rest voltage
Vth = 15; %Threshold voltage 
Ve_r = 60; %Excitatory reversal potential
Vi_r = -10; %Inhibitory reversal potential

% Number of neurons
N = 1;

% FF weights
We = 0.015; Wi = 0.06;
W_f = [We, Wi];
theta.Wf = W_f;

% Set Information into Structure
theta.tauleak = tauleak;
theta.tausyn = tausyn;
theta.trim = trimidx;
theta.ereversal = Ve_r;
theta.ireversal = Vi_r;
theta.lreversal = 0;

%% Example Voltage Traces
Ke = 1000; re = 10; % Excitatory input
Ki = 250; ri = 10;  % Inhibitory input

theta.K.ee = Ke; theta.r.ee = re;
theta.K.ei = Ki; theta.r.ei = ri;

% Uncorr Input
theta.corrinfo.corr = "uncorr";
[yef,yif] = makeffinput(theta,N,0,t,dt); 
wb_fe = yef.*We./tausyn; wb_fi = yif.*Wi./tausyn; %Weighted input

input.feedfowarde = wb_fe; input.feedfowardi = wb_fi;
Vuc = cifv(t,N,0,input,theta);

% Corr Input
rhoe = 0.03; rhoi = 0.03; rhoei = 0.03;
beta_e = 1./rhoe - 1; beta_i = 1./rhoi - 1;
theta.corrinfo.corr = "corr";
theta.corrinfo.corridx = 1;
theta.corrinfo.ae1 = 0; theta.corrinfo.ae2 = beta_e;
theta.corrinfo.ai1 = 0; theta.corrinfo.ai2 = beta_i;
theta.rhoWe = rhoe; theta.rhoWi = rhoi; theta.rhoWei = rhoei;

[yef,yif] = makeffinput(theta,N,0,t,dt); 
wb_fe = yef.*We/tausyn; wb_fi = yif.*Wi/tausyn; %Weighted input

input.feedfowarde = wb_fe; input.feedfowardi = wb_fi;
Vcorr = cifv(t,N,0,input,theta);

figure(1); clf;
subplot(1,4,[1,3.2]); hold on
plot(tt,Vuc)
plot(tt,Vcorr); 
legend({'Asynch','Synch'},'box','off')
ylim([0, 20]); xlim([0,1])
xlabel('Time (s)'); ylabel('Voltage (mV)')

subplot(1,4,[4,4]); hold on
histogram(Vuc(:),'binwidth',0.2,'normalization','pdf','orientation','horizontal')
histogram(Vcorr(:),'binwidth',0.2,'normalization','pdf','orientation','horizontal')
ylim([0, 20]); 
set(gca,'Yticklabels',[])

%% Calc moments for uncorr input data
theta.K.ee = 1000; theta.r.ee = 10; 
theta.K.ei = 250; theta.r.ei = 10;

theta.corrinfo.corr = "uncorr";
rhoe = 0; rhoi = 0; rhoei = 0;
theta.rhoWe = rhoe; theta.rhoWi = rhoi; theta.rhoWei = rhoei;
beta_e = 1./rhoe - 1; beta_i = 1./rhoi - 1;
theta.corrinfo.ae1 = 0; theta.corrinfo.ae2 = beta_e;
theta.corrinfo.ai1 = 0; theta.corrinfo.ai2 = beta_i;
theta.corrinfo.corridx = 0; 

T = 1e6; trim = 100;

% Calculate Moments
kk = 6; % number of moms
rp = 1; cenmomuc_ = zeros(kk,rp);
for z = 1:rp
    theta.corrinfo.corr = "uncorr";
    [~,~,nbstat] = cifvNB(theta,T,ttrim,kk); % Calculate stats using No Bin method
    mus = [1;nbstat.moms];
    % Convert to central moms
    for i = 2:kk
        for j = 1:i+1
            temp = nchoosek(i,j-1)*(-1)^(i-(j-1))*mus(j)*mus(2).^(i-(j-1));
            cenmomuc_(i,z) = cenmomuc_(i,z) + temp;
        end
    end
    [momzuc,tempm] = genmom(theta,kk);
end
cenmomluc_ = cenmomuc_;

% pact data 
datauc = zeros(kk,2);
for i = 1:kk
    datauc(i,:) = [momzuc(i),mean(cenmomluc_(i,:))];
end

%% Calc moments for corr data
theta.corrinfo.corr = "corr";
rhoe = 0.03; rhoi = 0.03; rhoei = 0.03;
theta.rhoWe = rhoe; theta.rhoWi = rhoi; theta.rhoWei = rhoei;
theta.corrinfo.ae1 = 0; theta.corrinfo.ae2 = 1/rhoe-1;
theta.corrinfo.ai1 = 0; theta.corrinfo.ai2 = 1/rhoi-1;
theta.corrinfo.corridx = 1;

% Calculate Moments
kk = 6; 
rp = 1; cenmomc_ = zeros(kk,rp);
for z = 1:rp    
    [~,~,nbstat] = cifvNB(theta,T,ttrim,kk); % Calculate stats using No Bin method
    mus = [1;nbstat.moms];
    for i = 2:kk
        for j = 1:i+1
            temp = nchoosek(i,j-1)*(-1)^(i-(j-1))*mus(j)*mus(2).^(i-(j-1));
            cenmomc_(i,z) = cenmomc_(i,z) + temp;
        end
    end
    momzc = genmom(theta,kk); % Theory
end
cenmomlc_ = cenmomc_;

% pact data 
datac = zeros(kk,2);
for i = 1:kk
    datac(i,:) = [momzc(i),mean(cenmomlc_(i,:))];
end
%% Plot results
% Combind both uc and c
data = [datauc,datac];

xx = 1:kk;
figure(2); clf; hold on
bar(data(2:end,:)); box off
set(gca,'YScale','log')
ylim([1e-3,1e4])
legend({'Theory UC','Sim UC','Theory C','Sim C'},'box','off')

%% Eff of skew for E, I and EI
rbank = 0:0.1:20;
skewVs = zeros(length(rbank),3);
for i = 1:length(rbank)
    % E only
    theta.K.ee = 1000; theta.K.ei = 0;
    theta.r.ee = rbank(i); theta.r.ei = 0;
    theta.rhoWe = 0.03; theta.rhoWi = 0; theta.rhoWei = 0;
    [~,~,skewVs(i,1)] = swstats(theta);
    
    % I only
    theta.K.ee = 0; theta.K.ei = 250;
    theta.r.ee = 0; theta.r.ei = rbank(i);
    theta.rhoWe = 0; theta.rhoWi = 0.03; theta.rhoWei = 0;
    [~,~,skewVs(i,2)] = swstats(theta);
    
    % E and I
    theta.K.ee = 1000; theta.K.ei = 250;
    theta.r.ee = rbank(i); theta.r.ei = rbank(i);
    
    theta.rhoWe = 0.03; theta.rhoWi = 0.03; theta.rhoWei = 0;
    [~,~,skewVs(i,3)] = swstats(theta);
end

figure(3); clf;
plot(rbank,skewVs,'linewidth',2)
legend({'e only','i only','e + i'},'box','off')
xlabel('r (Hz)'); ylabel('Skew')
box off

%% Eff of skew for difference in corr
rbank = 0:0.1:20;
skewVs = zeros(length(rbank),3);
for i = 1:length(rbank)
    % uncorr
    theta.K.ee = 1000; theta.K.ei = 250;
    theta.r.ee = rbank(i); theta.r.ei = rbank(i);
    theta.rhoWe = 0; theta.rhoWi = 0; theta.rhoWei = 0;
    [~,~,skewVs(i,1)] = swstats(theta);
    
    % Within corr
    theta.K.ee = 1000; theta.K.ei = 250;
    theta.r.ee = rbank(i); theta.r.ei = rbank(i);
    theta.rhoWe = 0.03; theta.rhoWi = 0.03; theta.rhoWei = 0;
    [~,~,skewVs(i,2)] = swstats(theta);
    
    % X corr
    theta.K.ee = 1000; theta.K.ei = 250;
    theta.r.ee = rbank(i); theta.r.ei = rbank(i);
    
    theta.rhoWe = 0.03; theta.rhoWi = 0.03; theta.rhoWei = 0.03;
    [~,~,skewVs(i,3)] = swstats(theta);
end

figure(4); clf;
plot(rbank,skewVs,'linewidth',2)
legend({'Uncorr','Within corr','Cross Corr'},'box','off')
xlabel('r (Hz)'); ylabel('Skew')
box off

%% Grid of K and Rho at low rate (1 Hz)
kbank = 100:1500; rhobank = linspace(0,0.04,150); 
theta.K.ee = 1000; theta.r.ee = 1; theta.rhoWe = 0.03; 
theta.K.ei = 250; theta.r.ei = 1; theta.rhoWi = 0.03;
theta.rhoWei = 0.03;
theta.Wf = W_f;

skewVkrho = zeros(length(rhobank),length(kbank));
for i = 1:length(rhobank)
    theta.rhoWe = rhobank(i); theta.rhoWi = rhobank(i); theta.rhoWei = rhobank(i);
    for j = 1:length(kbank)
        theta.K.ee = kbank(j); theta.K.ei = kbank(j)/4;
        [~,vartemp,skewVkrho(i,j)] = swstats(theta);
    end
end

figure(5); clf;
imagesc(rhobank,kbank,skewVkrho')
yticks(100:200:1500)
xlabel('\rho'); ylabel('K')
set(gca,'Ydir','normal')
box off; colorbar

%% Grid of rate and rho
rbank = 5:0.1:60; 
theta.K.ee = 1000; theta.K.ei = 250;
rhobank = linspace(0,0.04,100);
skewVkrho = zeros(length(rhobank),length(rbank));
for i = 1:length(rhobank)
    theta.rhoWe = rhobank(i); theta.rhoWi = rhobank(i); 
    theta.rhoWei = rhobank(i);
    for j = 1:length(rbank)
        theta.r.ee = rbank(j); theta.r.ei = rbank(j);
        [~,vartemp,skewVkrho(i,j)] = swstats(theta);
    end
end

figure(6); clf;
imagesc(rhobank,rbank,skewVkrho')
set(gca,'Ydir','normal')
xlabel('r_i (Hz)'); ylabel('r_e (Hz)')
box off; colorbar
