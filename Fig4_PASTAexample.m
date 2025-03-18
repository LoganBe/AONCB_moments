clear; clc; close all
rng(6895)
set(0,'DefaultAxesFontSize',14)

% DESCRIPTION HERE

%% Time Info and constants
dt = 0.1; %Time step (ms)
tmax = 5e3; %Max time (ms)
ttrim = 100; %Time to remove (ms)
t = 0:dt:tmax+ttrim; %time vector
tt = t(1:end-ttrim/dt-1)*1e-3; %Time Vector Trimed in in Seconds
trimidx = binsearch(t,ttrim); %trim time idx

%Time Constants
tauleak = 15; %Leak Time constant 
tausyn = 2*dt;

%Biophysical Paramters (Normalized to 0)
V0 = 0; % Rest voltage
Vth = 15; %Threshold voltage 
Ve_r = 60; %Excitatory reversal potential
Vi_r = -10; %Inhibitory reversal potential

% Number Neurons/trials
N = 2;

% FF Weights
We = 0.015; Wi = 0.06;

% Input correlation
rho_e = 0.03; rho_i = 0.03;

% Set Information into Structure
theta.tauleak = tauleak;
theta.tausyn = tausyn;
theta.trim = trimidx;
theta.ereversal = Ve_r;
theta.ireversal = Vi_r;
theta.lreversal = 0;

theta.corrinfo.corr = "corr";
theta.corrinfo.ae1 = 0; theta.corrinfo.ae2 = (1/rho_e)-theta.corrinfo.ae1-1; % Convert corr to beta (dist paramter)
theta.corrinfo.ai1 = 0; theta.corrinfo.ai2 = (1/rho_i)-theta.corrinfo.ai1-1; % Convert corr to beta (dist paramter)
theta.corrinfo.corridx = 1;

%% Presynaptic Input Info
Ke = 1000; re = 10; % E (Hz)
theta.K.ee = Ke; theta.r.ee = re; 

Ki = 250; ri = 10; % I (Hz)
theta.K.ei = Ki; theta.r.ei = ri;

% Make copula input
d = 2*N; corrmat = ones(d);

% MAKE INPUT
[ys,yraw] = makecopinput(theta,t,dt,d/2,d/2,corrmat);
ye1 = ys{1}; ye2 = ys{2}; % E input
yi1 = ys{3}; yi2 = ys{4}; % I input

% boxcar filtered input
yfilt = ones(floor(tausyn/dt),1);
yefsamp1 = filter(yfilt,1,ye1); yefsamp2 = filter(yfilt,1,ye2); 
yifsamp1 = filter(yfilt,1,yi1); yifsamp2 = filter(yfilt,1,yi2); 

% Weigthed input
wb_fe1 = yefsamp1.*We/tausyn; wb_fe2 = yefsamp2.*We/tausyn; 
wb_fi1 = yifsamp1.*Wi/tausyn; wb_fi2 = yifsamp2.*Wi/tausyn; 

% Simulate voltage
input.feedfowarde = wb_fe1;
input.feedfowardi = wb_fi1;
V1 = cifv(t,1,0,input,theta);

input.feedfowarde = wb_fe2;
input.feedfowardi = wb_fi2;
V2 = cifv(t,1,0,input,theta);

% Time the input
ye1_ = ye1(trimidx+1:end); ye2_ = ye2(trimidx+1:end); 
yi1_ = yi1(trimidx+1:end); yi2_ = yi2(trimidx+1:end); 

%% Plot Inputs
% Inputs for V1 
figure(1); clf; 
subplot(2,1,1); 
plot(tt,ye1_,'r'); box off;
xlim([0,1]); xticks([])
subplot(2,1,2)
plot(tt,yi1_,'b'); box off
xlim([0,1]);

% Inputs for V2
figure(2); clf; 
subplot(2,1,1); 
plot(tt,ye2_,'r'); box off;
xlim([0,1]); xticks([])
subplot(2,1,2)
plot(tt,yi2_,'b'); box off
xlim([0,1]);

% Voltages for V1 and V2
figure(3); clf;
ax1 = subplot(2,1,1); hold on
plot(tt,V1,'k');
yticks([]); xticks([]);
xlim([0,1]);
ax2 = subplot(2,1,2); hold on
plot(tt,V2,'k');
yticks([]); xlim([0,1])
xlabel('Time (s)')
%% Calc and plot densisites
% Find all points where V1 has input
idx1 = [find(ye1_ > 0)-1; find(yi1_ > 0)-1];
V1sp = V1(idx1); V2nsp = V2(idx1);

% Same for V2
idx2 = [find(ye2_ > 0)-1; find(yi2_ > 0)-1];
V1nsp = V1(idx2); V2sp = V2(idx2);

% Calc cov and make mvnpdf
mu1 = mean(V1nsp); mu2 = mean(V2sp);
mus = [mu1,mu2];
covv = cov(V1nsp,V2sp);
[X1,X2] = meshgrid(0:0.1:12,0:0.1:12);
X = [X1(:) X2(:)];
y = mvnpdf(X,mus,covv);
y = reshape(y,size(X1,1),size(X1,2));

figure(4); clf;
contour(X1,X2,y,'linewidth',2)
xlabel('V1'); ylabel('V2')
box off
%% Hypothetical example to show voltages
% Random ISIs spaced out
isie1 = [23.9234, 53.0867,84,12,500]; esamp1 = ones(size(isie1));
isii1 = 43.1903; isamp1 = ones(size(isii1));

% Spike Times and counts to discretize
ste = cumsum(isie1)'; sti = cumsum(isii1)'; 
ye1 = histcounts(ste,t); ye1 = ye1'; 
yi1 = histcounts(sti,t); yi1 = yi1';

% Filter input (inst input)
wb_fe1 = ye1*We/dt; wb_fi1 = yi1*Wi/dt;

% 2nd example of inputs
isie2 = [38.2312,19,121,500]; esamp2 = ones(size(isie2));
isii2 = 142; isamp2 = ones(size(isii2));

% Spike times and counts to discretize
ste = cumsum(isie2)'; sti = cumsum(isii2)'; 
ye2 = histcounts(ste,t); ye2 = ye2'; 
yi2 = histcounts(sti,t); yi2 = yi2';

% Filter input (inst input)
wb_fe2 = ye2*We/dt; wb_fi2 = yi2*Wi/dt;

% Adjusted time
tp = t(1:end-1);

% Sim Voltage
V1 = zeros(length(tp),1); V2 = zeros(length(tp),1);
for i = 2:length(tp)
    drive_fe = Ve_r-V1(i-1,:); drive_fi = Vi_r-V1(i-1,:);
    g_e = wb_fe1(i-1,:); g_i = wb_fi1(i-1,:);
    network1 = g_e.*drive_fe + g_i.*drive_fi; 
    leak = -V1(i-1,:); 
    V1(i,:) = V1(i-1,:) + (1./tauleak).*dt*(leak + network1); 
    
    drive_fe = Ve_r-V2(i-1,:);  drive_fi = Vi_r-V2(i-1,:);
    g_e = wb_fe2(i-1,:); g_i = wb_fi2(i-1,:);
    network2 = g_e.*drive_fe + g_i.*drive_fi;
    leak = -V2(i-1,:); 
    V2(i,:) = V2(i-1,:) + (1./tauleak).*dt*(leak + network2);
end

% Find peaks of V1 and V2
tt_ = t(1:end-1);
[~,idxV1] = findpeaks(V1);
[~,idxV2] = findpeaks(-V2);
mp = idxV1(1:end-1) + ceil(diff(idxV1)/2);

% Plot Results
figure(5); clf; 
subplot(2,1,1); hold on
plot(tt_,V1,'k','linewidth',3)
xline(tt_(idxV1(2))-dt*1e-3,'k--')
xline(tt_(mp(2))-dt*1e-3,'k--')
xline(tt_(idxV2(2))-dt*1e-3,'k--')
yticks([]); xticks([]);
xlim([0,250])
ylim([-0.05,0.1])
ylabel('Voltage')

subplot(2,1,2); hold on
plot(tt_,V2,'k','linewidth',3)
xline(tt_(idxV1(2))-dt*1e-3,'k--')
xline(tt_(mp(2))-dt*1e-3,'k--')
xline(tt_(idxV2(2))-dt*1e-3,'k--')
yticks([]); xticks([]);
xlim([0,250])
ylim([-0.05,0.1])
xlabel('Time'); ylabel('Voltage')