clear; clc; close all
rng(6895)

set(0,'DefaultAxesFontSize',14)

% Figures 7 and 9
% Examine the effects of cross neuron correlated inputs on voltage correlations. 
% First holds inhibitory paramters constant
% Second holds excitatory paramters constant
% Lastly changes in both excitatory and inhibitory firing rate
% Fig 7 - only within-pool correlations
% Fig 9 - cross pool correlations

% nomenclature:
% For correlations _w = within same neuron, _x = across pair of neurons
% For example, rhoe_w is excitatory<-> excitatory within same neuron
% while rhoe_x is excitatoy <-> excitatory across neurons

%% Time Info and constants
dt = 0.1; %Time step (ms)
tmax = 3e3; %Max time (ms)
ttrim = 100; %Time to remove (ms)
t = 0:dt:tmax+ttrim; %time vector
tt = t(1:end-ttrim/dt-1)*1e-3; %Time Vector Trimed in in Seconds
trimidx = binsearch(t,ttrim); %trim time idx

%Time Constants
tauleak = 15; %Leak Time constant 
tausyn = 2*dt; % Synaptic Time Constant

%Biophysical Paramters (Normalized to 0)
V0 = 0; % Rest voltage
Vth = 15; %Threshold voltage 
Ve = 60; %Excitatory reversal potential
Vi = -10; %Inhibitory reversal potential

% Number of neurons
N = 1;

% Feedfoward weights
We = 0.015; Wi = 0.06;
W_f = [We, Wi];

% Input correlations (no cross corr)
rhoe = 0.03; rhoi = 0.03; rhoei = 0;

% define theory stats 
meanVt_f = @(ke,ki,re,ri) ...
           (ke.*re.*We.*Ve + ki.*ri.*Wi.*Vi)./(1 + ke.*re.*We + ki.*ri.*Wi); % Mean

varVt_f = @(ke,ki,re,ri,mu,rhoe_w,rhoi_w,rhoei_w) ...
           (ke.*re.*We.^2.*(Ve-mu).^2.*(1+rhoe_w.*(ke-1)) + ki.*ri.*Wi.^2.*(Vi-mu).^2.*(1+rhoi_w.*(ki-1)) ...
            - 2.*rhoei_w.*ke.*ki.*sqrt(re.*ri).*We.*Wi.*(Ve-mu).*(mu-Vi)) ...
            ./(2.*tauleak.*(1+ke.*re.*We + ki.*ri.*Wi)); % Var
   
corrVt_f = @(ke,ki,re,ri,mu,rhoe_w,rhoi_w,rhoei_w,rhoe_x,rhoi_x,rhoei_x) ...
          (ke.^2.*re.*We.^2.*rhoe_x.*(Ve-mu).^2 + ki.^2.*ri.*Wi.^2.*rhoi_x.*(Vi-mu).^2 ...
            + 2.*rhoei_x.*ke.*ki.*sqrt(re.*ri).*We.*Wi.*(Ve-mu).*(Vi-mu)) ...
        ./(ke.*re.*We.^2.*(1+rhoe_w.*(ke-1)).*(Ve-mu).^2 + ki.*ri.*Wi.^2.*(1+rhoi_w.*(ki-1)).*(Vi-mu).^2 ...
                - 2.*rhoei_w.*sqrt(re.*ri).*ke.*We.*ki.*Wi.*(Ve - mu).*(mu - Vi)); % Corr

% Same as above but with larger weights 
WeL = We*10; WiL = Wi*10;
% Large Weight
meanVt_fL = @(ke,ki,re,ri) ...
           (ke.*re.*WeL.*Ve + ki.*ri.*WiL.*Vi)./(1 + ke.*re.*WeL + ki.*ri.*WiL);

corrVt_fL = @(ke,ki,re,ri,mu,rhoe_w,rhoi_w,rhoei_w,rhoe_x,rhoi_x,rhoei_x) ...
          (ke.^2.*re.*WeL.^2.*rhoe_x.*(Ve-mu).^2 + ki.^2.*ri.*WiL.^2.*rhoi_x.*(Vi-mu).^2 ...
            + 2.*rhoei_x.*ke.*ki.*sqrt(re.*ri).*WeL.*WiL.*(Ve-mu).*(Vi-mu)) ...
        ./(ke.*re.*WeL.^2.*(1+rhoe_w.*(ke-1)).*(Ve-mu).^2 + ki.*ri.*WiL.^2.*(1+rhoi_w.*(ki-1)).*(Vi-mu).^2 ...
                - 2.*rhoei_w.*sqrt(re.*ri).*ke.*WeL.*ki.*WiL.*(Ve - mu).*(mu - Vi));
%% Excitatory | Ke vs rho_e
Ke = 1600; re = 10; rhoe_w = 0.04; % Excitatory
Ki = 250; ri = 10; rhoi_w = 0.04; rhoi_x = 0.03; % Inhibitory 
rhoei_w = 0; rhoei_x = 0; % Cross Neuron Corr

rre = re*1e-3; rri = ri*1e-3; % Rate in Spikes/ms

% Grid of Ke vs rhoe_x
kbank = 0:Ke; rho = linspace(0,rhoe_w,1000);
[K,RHO] = meshgrid(kbank,rho);

meanVt = meanVt_f(K,Ki,rre,rri); 
rhoVt = corrVt_f(K,Ki,rre,rri,meanVt,rhoe_w,rhoi_w,rhoei_w,RHO,rhoi_x,rhoei_x);

figure(1); clf;
imagesc(rho,kbank,rhoVt')
set(gca,'Ydir','normal'); clim([0,1])
colorbar; box off;
xlabel('\rho_{ex}'); ylabel('K_e')

% Cross Section with Constant K_e
ke_ = 1000; rhoe_x = linspace(0,rhoe_w,1000);
meanVt = meanVt_f(ke_,Ki,rre,rri);
varVt = varVt_f(ke_,Ki,rre,rri,meanVt,rhoe_w,rhoi_w,rhoei_w);
corrVt = corrVt_f(ke_,Ki,rre,rri,meanVt,rhoe_w,rhoi_w,rhoei_w,rhoe_x,rhoi_x,rhoei_x);

figure(2); clf; hold on
plot(rhoe_x,meanVt.*ones(length(rhoe_x),1),'k','linewidth',2); box off;
xlabel('\rho_{ex}'); ylabel('\mu')

figure(3); clf; hold on
plot(rhoe_x,varVt.*ones(length(rhoe_x),1),'k','linewidth',2); box off;
ylim([6,9]); 
xlabel('\rho_{ex}'); ylabel('\sigma^2')

figure(4); clf; hold on
plot(rhoe_x,corrVt,'k','linewidth',2); box off;
ylim([0,1])
xlabel('\rho_{ex}'); ylabel('\rho')

% Cross Section with constant rhoe_x
ke_ = 0:Ke; rhoe_x = 0.03;
meanVt = meanVt_f(ke_,Ki,rre,rri);
varVt = varVt_f(ke_,Ki,rre,rri,meanVt,rhoe_w,rhoi_w,rhoei_w);
corrVt = corrVt_f(ke_,Ki,rre,rri,meanVt,rhoe_w,rhoi_w,rhoei_w,rhoe_x,rhoi_x,rhoei_x);

figure(5); clf; hold on
plot(ke_,meanVt.*ones(length(rhoe_x),1),'k','linewidth',2); box off;
xlim([0,Ke])
xlabel('K_e'); ylabel('\mu')

figure(6); clf; hold on
plot(ke_,varVt.*ones(length(rhoe_x),1),'k','linewidth',2); box off;
xlim([0,Ke])
xlabel('K_e'); ylabel('\sigma^2')

figure(7); clf; hold on
plot(ke_,corrVt,'k','linewidth',2); box off;
ylim([0,1]); xlim([0,Ke])
xlabel('K_e'); ylabel('\rho')

%% Inhibitory | Ki vs rhoi
Ke = 1000; re = 10; rhoe_w = 0.04; rhoe_x = 0.03; % excitatory
Ki = 400; ri = 10; rhoi_w = 0.04;  % inhibitory
rhoei_w = 0; rhoei_x = 0;

rre = re*1e-3; rri = ri*1e-3; % Rate in spikes/ms

% Grid of Ki vs rhoi
kbank = 0:Ki; rho = linspace(0,rhoi_w,1000);
[K,RHO] = meshgrid(kbank,rho);

meanVt = meanVt_f(Ke,K,rre,rri); 
rhoVt = corrVt_f(Ke,K,rre,rri,meanVt,rhoe_w,rhoi_w,rhoei_w,rhoe_x,RHO,rhoei_x);

figure(8); clf;
imagesc(rho,kbank,rhoVt')
set(gca,'Ydir','normal'); clim([0,1])
colorbar; box off;
xlabel('\rho_{ix}'); ylabel('K_i')

% Cross Section with constant Ki
ki_ = 250; rhoi_x = linspace(0,rhoi_w,1000);
meanVt = meanVt_f(Ke,ki_,rre,rri);
varVt = varVt_f(Ke,ki_,rre,rri,meanVt,rhoe_w,rhoi_w,rhoei_w);
corrVt = corrVt_f(Ke,ki_,rre,rri,meanVt,rhoe_w,rhoi_w,rhoei_w,rhoe_x,rhoi_x,rhoei_x);

figure(9); clf; hold on
plot(rhoi_x,meanVt.*ones(length(rhoi_x),1),'k','linewidth',2); box off;
xlabel('\rho_{ix}'); ylabel('\mu')

figure(10); clf; hold on
plot(rhoi_x,varVt.*ones(length(rhoi_x),1),'k','linewidth',2); box off;
xlabel('\rho_{ix}'); ylabel('\sigma^2')

figure(11); clf; hold on
plot(rhoi_x,corrVt,'k','linewidth',2); box off;
ylim([0,1])
xlabel('\rho_{ix}'); ylabel('\rho')

% Cross section with constant rho_i
ki_ = 0:Ki; rhoi_x = 0.03;
meanVt = meanVt_f(Ke,ki_,rre,rri);
varVt = varVt_f(Ke,ki_,rre,rri,meanVt,rhoe_w,rhoi_w,rhoei_w);
corrVt = corrVt_f(Ke,ki_,rre,rri,meanVt,rhoe_w,rhoi_w,rhoei_w,rhoe_x,rhoi_x,rhoei_x);

figure(12); clf; hold on
plot(ki_,meanVt.*ones(length(rhoi_x),1),'k','linewidth',2); box off;
xlim([0,Ki])
xlabel('K_i'); ylabel('\mu')

figure(13); clf; hold on
plot(ki_,varVt.*ones(length(rhoi_x),1),'k','linewidth',2); box off;
xlim([0,Ki])
xlabel('K_i'); ylabel('\sigma^2')

figure(14); clf; hold on
plot(ki_,corrVt,'k','linewidth',2); box off;
ylim([0,1]); xlim([0,Ki])
xlabel('K_i'); ylabel('\rho')

%% Change in E and I rate 
ke = 1000; ki = 250;
rhoe_w = 0.02; rhoi_w = 0.02; rhoei_w = 0.02; 
rhoe_x = 0.013; rhoi_x = 0.013; rhoei_x = 0.013; 

% Grid over re and ri
rbank = 0:0.1:50; 
[Re,Ri] = meshgrid(rbank,rbank); Re = Re*1e-3; Ri = Ri*1e-3;
meanVt = meanVt_f(ke,ki,Re,Ri); 
rhoVt = corrVt_f(ke,ki,Re,Ri,meanVt,rhoe_w,rhoi_w,rhoei_w,rhoe_x,rhoi_x,rhoei_x);

figure(8); clf;
imagesc(rbank,rbank,rhoVt')
set(gca,'Ydir','normal'); clim([0,1])
colorbar; box off;
xlabel('r_i'); ylabel('r_e')

% Cross section over joint Re and Ri
meanVt = meanVt_f(ke,ki,rbank*1e-3,rbank*1e-3);
varVt = varVt_f(ke,ki,rbank*1e-3,rbank*1e-3,meanVt,rhoe_w,rhoi_w,rhoei_w);
corrVt = corrVt_f(ke,ki,rbank*1e-3,rbank*1e-3,meanVt,rhoe_w,rhoi_w,rhoei_w,rhoe_x,rhoi_x,rhoei_x);

figure(15); clf; hold on
plot(rbank,meanVt,'k','linewidth',2); box off;
xlabel('r (Hz)'); ylabel('\mu')

figure(16); clf; hold on
plot(rbank,varVt,'k','linewidth',2); box off;
xlabel('r (Hz)'); ylabel('\sigma^2')

figure(17); clf; hold on
plot(rbank,corrVt,'k','linewidth',2); box off;
ylim([0,1])
xlabel('r (Hz)'); ylabel('\rho')

%% FIGURE 9 -- Effects of CROSS-POOL E vs CROSS-POOL I

% Moderate Weight
Ke = 1000; re = 10; rhoe_w = 0.04; % Excitatory
Ki = 250; ri = 10; rhoi_w = 0.04;  % Inhibitory
rhoei_w = 0; rhoei_x = 0; 
rre = re*1e-3; rri = ri*1e-3; 

rho_e = linspace(0,rhoe_w,1000); rho_i = linspace(0,rhoi_w,1000);
[RHOE,RHOI] = meshgrid(rho_e,rho_i);

meanVt = meanVt_f(Ke,Ki,rre,rri); 
rhoVt = corrVt_f(Ke,Ki,rre,rri,meanVt,rhoe_w,rhoi_w,rhoei_w,RHOE,RHOI,rhoei_x);

figure(18); clf;
imagesc(rho_i,rho_e,rhoVt')
set(gca,'Ydir','normal'); clim([0,1])
colorbar;box off;
xlabel('\rho_{ix}'); ylabel('\rho_{ex}')

% Sliver for constant rhoi_x and rhoe_x
rho_i_x = 0.03; 
meanVt = meanVt_f(Ke,Ki,rre,rri);
corrVt_swE = corrVt_f(Ke,Ki,rre,rri,meanVt,rhoe_w,rhoi_w,rhoei_w,rho_e,rho_i_x,rhoei_x);

rho_e_x = 0.03;
meanVt = meanVt_f(Ke,Ki,rre,rri);
corrVt_swI = corrVt_f(Ke,Ki,rre,rri,meanVt,rhoe_w,rhoi_w,rhoei_w,rho_e_x,rho_i,rhoei_x);

% Large Weight
Ke = 100; re = 10; rhoe_w = 0.04;
Ki = 25; ri = 10; rhoi_w = 0.04; 
rhoei_w = 0; rhoei_x = 0;

rre = re*1e-3; rri = ri*1e-3; 

rho_e = linspace(0,rhoe_w,1000); rho_i = linspace(0,rhoi_w,1000);
[RHOE,RHOI] = meshgrid(rho_e,rho_i);

meanVt = meanVt_fL(Ke,Ki,rre,rri); 
rhoVt = corrVt_fL(Ke,Ki,rre,rri,meanVt,rhoe_w,rhoi_w,rhoei_w,RHOE,RHOI,rhoei_x);

figure(19); clf;
imagesc(rho_i,rho_e,rhoVt')
set(gca,'Ydir','normal'); clim([0,1])
colorbar; box off;
xlabel('\rho_{ix}'); ylabel('\rho_{ex}')

% Sliver for constant rhoi_x and rhoe_x
rho_i_x = 0.03; 
meanVt = meanVt_fL(Ke,Ki,rre,rri);
corrVt_lwE = corrVt_fL(Ke,Ki,rre,rri,meanVt,rhoe_w,rhoi_w,rhoei_w,rho_e,rho_i_x,rhoei_x);

rho_e_x = 0.03;
meanVt = meanVt_fL(Ke,Ki,rre,rri);
corrVt_lwI = corrVt_fL(Ke,Ki,rre,rri,meanVt,rhoe_w,rhoi_w,rhoei_w,rho_e_x,rho_i,rhoei_x);

figure(20); clf; hold on
plot(rho_e,corrVt_swE)
plot(rho_e,corrVt_lwE)
legend({'Small Weight','Large Weight'},'Box','off')
ylim([0,1])
xlabel('\rho_{ex}'); ylabel('\rho')

figure(21); clf; hold on
plot(rho_i,corrVt_swI)
plot(rho_i,corrVt_lwI)
legend({'Small Weight','Large Weight'},'Box','off')
ylim([0,1])
xlabel('\rho_{ix}'); ylabel('\rho')

%% FIGURE 9 -- Effects of CROSS-POOL vs CROSS-NEURON
% Moderate Weight
Ke = 1000; re = 10; rhoe_w = 0.04; % Excitatory
Ki = 250; ri = 10; rhoi_w = 0.04;  % Inhibitory
rre = re*1e-3; rri = ri*1e-3; 

meanVt = meanVt_f(Ke,Ki,rre,rri);

% Grid for rho_w with rho_x
rho_x_ind = linspace(0,rhoe_w,500); 
rhoVt = zeros(length(rho_x_ind),length(rho_x_ind));
for i = 1:length(rho_x_ind)
    rho_x_x = linspace(0,rho_x_ind(i),500);
    rhoVt(i,:) = corrVt_f(Ke,Ki,rre,rri,meanVt,rhoe_w,rhoi_w,rho_x_x,rho_x_ind(i),rho_x_ind(i),rho_x_x);
end

figure(22); clf;
imagesc(linspace(0,1,length(rho_x_ind)),rho_x_ind,rhoVt)
set(gca,'Ydir','normal'); clim([0,1])
colorbar; box off;
xlabel('\rho_x'); ylabel('\rho_w')

% Cross section with constant rho_x and rho_w
rho_x_x = rho_x_ind.*0.6; 
meanVt = meanVt_f(Ke,Ki,rre,rri);
corrVt_swW = corrVt_f(Ke,Ki,rre,rri,meanVt,rhoe_w,rhoi_w,rho_x_x,rho_x_ind,rho_x_ind,rho_x_x);

rho_x_ind = 0.03; rho_x_x = linspace(0,rho_x_ind,500);
meanVt = meanVt_f(Ke,Ki,rre,rri);
corrVt_swX = corrVt_f(Ke,Ki,rre,rri,meanVt,rhoe_w,rhoi_w,rho_x_x,rho_x_ind,rho_x_ind,rho_x_x);

% Large Weight
Ke = 100; re = 10; 
Ki = 25; ri = 10; 
rre = re*1e-3; rri = ri*1e-3; 

meanVt = meanVt_fL(Ke,Ki,rre,rri);

rho_x_ind = linspace(0,rhoe_w,500); 
rhoVt = zeros(length(rho_x_ind),length(rho_x_ind));
for i = 1:length(rho_x_ind)
    rho_x_x = linspace(0,rho_x_ind(i),500);
    rhoVt(i,:) = corrVt_fL(Ke,Ki,rre,rri,meanVt,rhoe_w,rhoi_w,rho_x_x,rho_x_ind(i),rho_x_ind(i),rho_x_x);
end

figure(23); clf;
imagesc(linspace(0,1,length(rho_x_ind)),rho_x_ind,rhoVt)
set(gca,'Ydir','normal'); clim([0,1])
colorbar; box off;
xlabel('\rho_x'); ylabel('\rho_w')

rho_x_x = rho_x_ind.*0.6; 
meanVt = meanVt_fL(Ke,Ki,rre,rri);
corrVt_lwW = corrVt_fL(Ke,Ki,rre,rri,meanVt,rhoe_w,rhoi_w,rho_x_x,rho_x_ind,rho_x_ind,rho_x_x);

rho_x_ind = 0.03; rho_x_x = linspace(0,rho_x_ind,500);
meanVt = meanVt_f(Ke,Ki,rre,rri);
corrVt_lwX = corrVt_fL(Ke,Ki,rre,rri,meanVt,rhoe_w,rhoi_w,rho_x_x,rho_x_ind,rho_x_ind,rho_x_x);

rho_x_ind = linspace(0,rhoe_w,500); 
figure(24); clf; hold on
plot(rho_x_ind,corrVt_swW)
plot(rho_x_ind,corrVt_lwW)
legend({'Small Weight','Large Weight'},'Box','off')
ylim([0,1])
xlabel('\rho_w'); ylabel('\rho')

rho_x_x = linspace(0,0.03,500);
figure(25); clf; hold on
plot(rho_x_x,corrVt_swX)
plot(rho_x_x,corrVt_lwX)
legend({'Small Weight','Large Weight'},'Box','off')
ylim([0,1])
xlabel('\rho_x'); ylabel('\rho')
