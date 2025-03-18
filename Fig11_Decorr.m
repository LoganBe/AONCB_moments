clear; clc; close all;
rng(6895)
set(0,'DefaultAxesFontSize',14)
% DESCRIPTION HERE

%% Set basic parameters
dt = 0.1; %Time step (ms)
tmax = 5e3; %Max time (ms)
ttrim = 100; %Time to remove (ms)
t = 0:dt:tmax+ttrim; %time vector (ms)
tt = t(1:end-ttrim/dt-1)*1e-3;
trimidx = binsearch(t,ttrim); %trim time idx

%Time Constants
tauleak = 15; %Leak constant
tausyn = 2*dt;

%Biophysical Paramters (Normalized to 0)
V0 = 0; % Rest voltage
Vth = 15; %Threshold voltage 
Ve = 60; %Excitatory reversal potential
Vi = -10; %Inhibitory reversal potential

% FF weights
We = 0.015; %Excitatory Weight
Wi = 0.06; %Inhibitory Weight
Wf = [We,Wi];

% Number of neurons
N = 5;

% Set Information into Structure
theta.tauleak = tauleak;
theta.tausyn = tausyn;
theta.trim = trimidx;
theta.ereversal = Ve;
theta.ireversal = Vi;
theta.lreversal = 0;
theta.Wf = Wf;

%% Decorrelate voltage by having indep input
% Parameters for correlated inputs
Ke = 1000; Ki = 250; % Number Presynaptic Synapses
re = 5; ri = 10; % Rate (Hz)
Re = re*1e-3; Ri = ri*1e-3; % Rate in /ms

rhoe = 0.03; rhoi = 0.03; rhoei = 0; % Input corrs 
beta_e = 1./rhoe - 1; beta_i = 1./rhoi - 1; % Betabin paramters

% Save into struct
theta.K.ee = Ke; theta.r.ee = re; 
theta.K.ei = Ki; theta.r.ei = ri;
theta.corrinfo.ae1 = 0; theta.corrinfo.ae2 = beta_e;
theta.corrinfo.ai1 = 0; theta.corrinfo.ai2 = beta_i;

% Transfer function between desired level of correlation and copula correlation matrix. 
% Was previously fit and parameters are located and loaded in 'raparm'
d = 4; 
modelfunc = @(b,x) b(1)*log((x-b(2))/b(3)); load('raparm'); 
rho_x = 0.023;

alpha = modelfunc(raparm,rho_x);
if alpha > 1; alpha = 1; end
if alpha < 0; alpha = 0; end

corrmat = zeros(d);
corrmat(eye(d,'logical')) = 1;
corrmat(1,2) = alpha; %E1 <-> E2 
corrmat(3,4) = alpha; %I1 <-> I2 
corrmat = corrmat + corrmat' - diag(diag(corrmat));

% Make correlated input via copula
[ys,yraw] = makecopinput(theta,t,dt,2,2,corrmat,N);
ye1 = ys{1}; ye2 = ys{2};
yi1 = ys{3}; yi2 = ys{4};

% Caluclate correlation of input spiking activity
ccin_ = zeros(1,3); ccx_ = zeros(1,3);
for i = 1:N 
    c1 = mean(nonzeros(triu(corr(yraw{i,1},'rows','pairwise'),1)),'all','omitnan'); % E1 -> E1
    c2 = mean(nonzeros(triu(corr(yraw{i,2},'rows','pairwise'),1)),'all','omitnan'); % E2 -> E2
    ccin_(i,1) = mean([c1,c2]);
    
    c1 = mean(nonzeros(triu(corr(yraw{i,3},'rows','pairwise'),1)),'all','omitnan'); % I1 -> I1
    c2 = mean(nonzeros(triu(corr(yraw{i,4},'rows','pairwise'),1)),'all','omitnan'); % I2 -> I2
    ccin_(i,2) = mean([c1,c2]);
    
    c1 = mean(nonzeros(triu(corr(yraw{i,1},yraw{i,3},'rows','pairwise'),1)),'all','omitnan'); % E1 -> I1
    c2 = mean(nonzeros(triu(corr(yraw{i,2},yraw{i,4},'rows','pairwise'),1)),'all','omitnan'); % E2 -> I2
    ccin_(i,3) = mean([c1,c2]);

    ccx_(i,1) = mean(nonzeros(triu(corr(yraw{i,1},yraw{i,2},'rows','pairwise'),1)),'all','omitnan'); % E1 -> E2
    ccx_(i,2) = mean(nonzeros(triu(corr(yraw{i,3},yraw{i,4},'rows','pairwise'),1)),'all','omitnan'); % I1 -> I2

    c1 = mean(nonzeros(triu(corr(yraw{i,1},yraw{i,4},'rows','pairwise'),1)),'all','omitnan'); % E1 -> I1
    c2 = mean(nonzeros(triu(corr(yraw{i,2},yraw{i,3},'rows','pairwise'),1)),'all','omitnan'); % E2 -> I2
    ccx_(i,3) = mean([c1,c2]);
    
end
cc_in = zeros(1,6); cc_in(:,1:3) = mean(ccin_); 
cc_x = zeros(1,4); cc_x(:,1:3) = mean(ccx_); 

% Filtered and weighted sum of correlated inputs
yfilt = ones(floor(tausyn/dt),1);
yef = filter(yfilt,1,ye1); yif = filter(yfilt,1,yi1); 
wbfe_c1 = yef.*We/tausyn; wbfi_c1 = yif.*Wi/tausyn;
yef = filter(yfilt,1,ye2); yif = filter(yfilt,1,yi2); 
wbfe_c2 = yef.*We/tausyn; wbfi_c2 = yif.*Wi/tausyn;

% Paramters for uncorrelated inputs
Ke_uc = 100; Ki_uc = 0;
re_uc = 5; ri_uc = 0;
We_uc = 10*We;

bank = 0:0.2:10; % Excitatory input rates (Hz)
cc_in = repmat(cc_in,length(bank),1); cc_x = repmat(cc_x,length(bank),1);
varV = zeros(length(bank),3);
corrV = zeros(length(bank),2);
for rr = 1:length(bank)
    theta.r.ee = bank(rr); 
    Re_uc = theta.r.ee*1e-3; 
      
    % Simulation
    % Make input
    theta.corrinfo.corr = "uncorr";
    theta.K.ee = Ke_uc; 
    theta.K.ei = 0; theta.r.ei = 0;
    [yef,~,yef_raw1] = makeffinput(theta,N,0,t,dt); 
    wbfe_uc1 = yef.*We_uc/tausyn;
    [yef,~,yef_raw2] = makeffinput(theta,N,0,t,dt); 
    wbfe_uc2 = yef.*We_uc/tausyn;

    % Sim voltage
    input.feedfowarde = wbfe_c1+wbfe_uc1; input.feedfowardi = wbfi_c1;
    V1 = cifv(t,N,0,input,theta);
    varV(rr,1) = mean(var(V1,[],1));

    input.feedfowarde = wbfe_c2+wbfe_uc2; input.feedfowardi = wbfi_c2;
    V2 = cifv(t,N,0,input,theta);
    varV(rr,2) = mean(var(V2,[],1));

    corrV(rr,1) = mean(diag(corr(V1,V2)),'all');

    % Estimate exact corr of 'uncorr'
    if Re_uc > 0
        ccin_ = zeros(N,3);
        ccx_ = zeros(N,1);
        for i = 1:N
            yuncorr_e1 = zeros(length(t),Ke_uc);
            yuncorr_e2 = zeros(length(t),Ke_uc);
            for tt = 1:length(t)
                yuncorr_e1(tt,randperm(Ke_uc,yef_raw1(tt,i))) = 1;
                yuncorr_e2(tt,randperm(Ke_uc,yef_raw2(tt,i))) = 1;
            end
            ycorr_e1 = yraw{i,1};
            ycorr_i1 = yraw{i,3};
    
            ccin_(i,1) = mean(nonzeros(triu(corr(yuncorr_e1,'rows','pairwise'),1)),'all','omitnan');
            ccin_(i,2) = mean(nonzeros(triu(corr(yuncorr_e1,ycorr_e1,'rows','pairwise'),1)),'all','omitnan');
            ccin_(i,3) = mean(nonzeros(triu(corr(yuncorr_e1,ycorr_i1,'rows','pairwise'),1)),'all','omitnan');
            
            ccx_(i,1) = mean(nonzeros(triu(corr(yuncorr_e1,yuncorr_e2,'rows','pairwise'),1)),'all','omitnan');
        end
        cc_in(rr,4:end) = mean(ccin_);
        cc_x(rr,4) = mean(ccx_);
    end
   
    % Theory using Smallweight
    muT = (Ke*Re*We*Ve + Ki*Ri*Wi*Vi + Ke_uc*Re_uc*We_uc*Ve)./(1 + Ke*Re*We + Ki*Ri*Wi + Ke_uc*Re_uc*We_uc);
    varT = (Ke*Re*We^2*(1 + cc_in(rr,1)*(Ke-1))*(Ve-muT)^2 + Ke_uc*Re_uc*We_uc^2*(1 + cc_in(rr,4)*(Ke_uc-1))*(Ve-muT)^2 + Ki*Ri*Wi^2*(1+cc_in(rr,2)*(Ki-1))*(Vi-muT)^2 ... 
            -cc_in(rr,3).*Ke.*We.*Ki.*Wi.*sqrt(Re*Ri).*(Ve-muT).*(muT-Vi)-cc_in(rr,5).*Ke.*We*Ke_uc.*We_uc.*sqrt(Re*Re_uc).*(Ve-muT).*(muT-Ve)...
            -cc_in(rr,6).*Ke_uc.*We_uc*Ki.*Wi.*sqrt(Re_uc*Ri).*(Ve-muT).*(muT-Vi))...
                /(2*tauleak*(1+Ke*Re*We + Ki*Ri*Wi + Ke_uc*Re_uc*We_uc));
    varV(rr,3) = varT;

    covV = (cc_x(rr,1).*Ke.*Ke.*sqrt(Re.*Re).*We.*We.*(Ve-muT).*(Ve-muT)...
                 +cc_x(rr,2).*Ki.*Ki.*sqrt(Ri.*Ri).*Wi.*Wi.*(Vi-muT).*(Vi-muT) ...
                 +2*cc_x(rr,3).*Ke.*Ki.*sqrt(Re.*Ri).*We.*Wi.*(Ve-muT.*(Vi-muT) ...
                 +cc_x(rr,4).*Ke_uc.^2.*Re_uc.*We_uc.^2.*(Ve-muT.^2)))...
                    ./(2*tauleak.*(1 + We.*Ke.*Re + Wi.*Ki.*Ri + We_uc.*Ke_uc.*Re_uc));
    corrV(rr,2) = covV./varV(rr,3);
    pcdone(rr,length(bank))
end
%% PLot Results
figure(1); clf; hold on
plot(bank,varV(:,1:2)','linewidth',2)
plot(bank,varV(:,3),'k--','linewidth',2)

figure(2); clf; hold on
plot(bank,corrV(:,1),'linewidth',2)
plot(bank,corrV(:,2),'k--','linewidth',2)
ylim([0,1])