clear; clc; close all;
rng(6895)
set(0,'DefaultAxesFontSize',14)
% DESCRIPTION HERE

%% Time Info and constants
dt = 0.1; %Time step (ms)
tmax = 1e3; %Max time (ms)
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
Ve_r = 60; %Excitatory reversal potential
Vi_r = -10; %Inhibitory reversal potential

% FF weights
Wee = 0.015; %Excitatory Weight
Wei = 0.06; %Inhibitory Weight
Wie = 0; Wii = 0;
Wf = [Wee,Wei];

% Number of neurons
N = 1;

% Set Information into Structure
theta.tauleak = tauleak;
theta.tausyn = tausyn;
theta.trim = trimidx;
theta.ereversal = Ve_r;
theta.ireversal = Vi_r;
theta.lreversal = 0;
theta.Wf = Wf;

%% Presynaptic inputs
Ke = 1000; re = 10; % Excitatory
Ki = 0; ri = 0; % Inhibitory (E only)

theta.corrinfo.corr = "corr";
theta.corrinfo.corridx = 0;
theta.corrinfo.ai1 = 0; theta.corrinfo.ai2 = 0;

theta.K.ee = Ke; theta.K.ei = Ki;
theta.r.ee = re; theta.r.ei = ri;
%% Skewness for Jitter for varying rate
sig = 50; % STD for jitter (ms)

rbank = linspace(5,25,25);
VexAll_ = []; VexAllJ = [];
VexInd_ = zeros(length(rbank),N*length(tt)); VexIndJ = zeros(length(rbank),N*length(tt));
for rr = 1:length(rbank)
    theta.r.ee = rbank(rr);
    
    % Make input for instant corr example
    rhoe = 0.03;
    theta.corrinfo.ae1 = 0; theta.corrinfo.ae2 = (1/rhoe)-1;
    y_fee = makeffinput(theta,N,0,t,dt);
    wb_fe = y_fee*Wee/tausyn; wb_fi = zeros(size(wb_fe));

    input.feedfowarde = wb_fe; input.feedfowardi = wb_fi;
    V_ = cifv(t, N, 0, input, theta);
    VexAll_ = [VexAll_,V_(:,1)];
    VexInd_(rr,:) = V_(:);
    
    % Make input for large corr and add jitter
    rhoe = 0.2; 
    theta.corrinfo.ae1 = 0; theta.corrinfo.ae2 = (1/rhoe)-1;
    [~,~,ye] = makeffinput(theta,N,0,t,dt);
    
    % Apply jitter through Normal Dist on spike counts
    yetemp = zeros(size(ye)); 
    for nn = 1
        spidx = find(ye(:,nn)>1); % spk idx
        for i = 1:length(spidx) 
            stemp = spidx(i); 
            yetemp(stemp,nn) = 0; % Remove spikes 
            tidx = int32((randn(ye(stemp,nn),1)*sig + t(stemp))/dt); % Find jittered times
            tidx(tidx <= 0) = []; tidx(tidx > length(t)) = []; % If jitter is neg or too long, remove spike
            unnum = unique(tidx); nc = histc(tidx,unnum); % Counts within bins
            yetemp(unnum,nn) = yetemp(unnum,nn) + nc; % Add it to spike train
        end
    end
    % Weighted Filtered Input
    yfilt = ones(floor(tausyn/dt),1); 
    y_fee = filter(yfilt,1,yetemp);
    wb_fe = y_fee*Wee/tausyn; wb_fi = zeros(size(wb_fe));
    input.feedfowarde = wb_fe; input.feedfowardi = wb_fi;

    [V, y] = cifv(t, N, 0, input, theta);
    VexAllJ = [VexAllJ,V(:,1)];
    VexIndJ(rr,:) = V(:);
end

%% Make Figures
% Color bank (default colors)
cl = [0.8141484044598232, 0.2196847366397539, 0.3048058439061899;...
      0.9330257593233372, 0.3913110342176086, 0.27197231833910035;...
      0.9817762399077278, 0.6073817762399076, 0.3457900807381776;...
      0.3600153787004998, 0.7161860822760476, 0.6655132641291811;...
      0.21299500192233756, 0.5114186851211072, 0.730795847750865];

% Plot voltage over time 
figure(1); clf; hold on;
tbank = 0:1e3:10*1e3;
for i = 1:length(rbank)
    ts = tbank(i); te = tbank(i+1);
    tp = ts:dt:te-dt;
    plot(tp*1e-3,VexAll_(:,i)-mean(VexAll_(:,i)),'color',cl(i,:))
    plot(tp*1e-3,VexAllJ(:,i)-mean(VexAllJ(:,i)),'color',cl(i,:),'linewidth',2)
    if i < length(rbank)
        plot([tp(end),te]*1e-3,[VexAll_(end,i)-mean(VexAll_(:,i)),VexAll_(1,i+1)-mean(VexAll_(:,i+1))],'color',cl(i,:))
        plot([tp(end),te]*1e-3,[VexAllJ(end,i)-mean(VexAllJ(:,i)),VexAllJ(1,i+1)-mean(VexAllJ(:,i+1))],'color',cl(i,:))
    end
end
xlabel('Time (s)'); ylabel('Voltage'); 

% Plot example histograms
idx1 = binsearch(rbank,5); idx2 = binsearch(rbank,15); idx3 = binsearch(rbank,25);
idxs = [idx1,idx2,idx3];
for i = 1:length(idxs)
    figure(i+1); clf; hold on
    [h,x] = histcounts(VexInd_(idxs(i),:),'binwidth',0.2,'normalization','pdf');
    stairs(x(1:end-1),h,'color',cl(idxs(i),:),'linewidth',2)
    
    [h,x] = histcounts(VexIndJ(idxs(i),:),'binwidth',0.2,'normalization','pdf');
    stairs(x(1:end-1),h,'color',cl(idxs(i),:),'linewidth',4)
    xlim([0,30])
    xlabel('Voltage'); ylabel('Density'); 

end

%% Correlation between 2 neurons with jitter (Uses copula)
rbank = linspace(5,25,5);
% A key feature is we need to adjust the cross-neuron corr slightly to
% adjust for lower rates. This is given in abank.
abank = linspace(0.0268,0.025,length(rbank)); 

% Correlation window size in ms and std in ms
Dt = 25; sig = 50; 

rp = 25; % Number of populations
d = 4; % Number of pops (E1, E2, I1, I2) - will only use E1 and E2 here

% Transfer function between desired level of correlation and copula correlation matrix. 
% Was previously fit and parameters are located and loaded in 'raparm'
modelfunc = @(b,x) b(1)*log((x-b(2))/b(3));
load('raparm');
V1plot_ = []; V2plot_ = [];
V1plot = []; V2plot = [];

% Log info for theory
thetaTh = theta; 
thetaTh.rhoWi = 0; thetaTh.rhoWei = 0;
rhom = zeros(2,4); 
thetaTh.rhoX.ee = rhom(:,1); thetaTh.rhoX.ei = rhom(:,2); 
thetaTh.rhoX.ie = rhom(:,3); thetaTh.rhoX.ii = rhom(:,4);
thetaTh(2) = thetaTh;

corrout = zeros(length(rbank),2); 
VcellJ = cell(length(rbank),2); VcellI = cell(length(rbank),2);
V1plot = zeros(length(tt),length(rbank),2);
V2plot = zeros(length(tt),length(rbank),2);
for rr = 1:length(rbank)
    re = rbank(rr); theta.r.ee = re; 
    
    % First do jittered input with inst corr of 0.25
    rhoe = 0.25; 
    theta.corrinfo.ae1 = 0; theta.corrinfo.ae2 = (1/rhoe)-1;
    theta.corrinfo.ai1 = 0; theta.corrinfo.ai2 = (1/rhoe)-1;

    %MAKE INPUT using Copula
    alpha = modelfunc(raparm,abank(rr));
    if alpha > 1; alpha = 1; end
    if alpha < 0; alpha = 0; end

    % Correlation matrix
    corrmat = zeros(d);
    corrmat(1,1) = 1; %E1 -> E1 
    corrmat(2,2) = 1; %E2 -> E2
    corrmat(1,2) = alpha; %E1 <-> E2 
    corrmat(2,1) = alpha; %E2 <-> E1
    
    % Get spike counts
    ys = makecopinput(theta,t,dt,d/2,d/2,corrmat,rp);
    ye1 = ys{1}; ye2 = ys{2}; 
    
    % Apply Jiter
    ye1temp = zeros(size(ye1)); ye2temp = zeros(size(ye2));
    cc_ = zeros(N,4,4);
    for nn = 1:rp
        spidx = find(ye1(:,nn)>1); %Find when i have an event
        for i = 1:length(spidx) %For each event
            stemp = spidx(i); %Take that event
            ye1temp(stemp,nn) = 0; %Remove that events
            tidx = int32((randn(ye1(stemp,nn),1)*sig + t(stemp))/dt);
            tidx(tidx <= 0) = []; tidx(tidx > length(t)) = [];
            unnum = unique(tidx); nc = histc(tidx,unnum);
            ye1temp(unnum,nn) = ye1temp(unnum,nn) + nc;
        end

        spidx = find(ye2(:,nn)>1); %Find when i have an event
        for i = 1:length(spidx) %For each event
            stemp = spidx(i); %Take that event
            ye2temp(stemp,nn) = 0; %Remove that events
            tidx = int32((randn(ye2(stemp,nn),1)*sig + t(stemp))/dt);
            tidx(tidx <= 0) = []; tidx(tidx > length(t)) = [];
            unnum = unique(tidx); nc = histc(tidx,unnum);
            ye2temp(unnum,nn) = ye2temp(unnum,nn) + nc;
        end

        % Get spike trains
        yspkse1 = zeros(length(t),theta.K.ee); yspkse2 = zeros(length(t),theta.K.ee);
        for j = 1:length(t)
            yspkse1(j,randperm(theta.K.ee,ye1temp(j,nn))) = 1;
            yspkse2(j,randperm(theta.K.ee,ye2temp(j,nn))) = 1;
        end

        % Filter data to win size of DT
        bin_per_newbin = floor(Dt/dt);
        num_nb = int32(floor(size(yspkse1,1) / bin_per_newbin));
        newsize = num_nb*bin_per_newbin;
        y1rsh = reshape(yspkse1(1:newsize,:),bin_per_newbin,num_nb,size(yspkse1,2));
        outpute1 = squeeze(sum(y1rsh,1));
        y2rsh = reshape(yspkse2(1:newsize,:),bin_per_newbin,num_nb,size(yspkse2,2));
        outpute2 = squeeze(sum(y2rsh,1));

        % Calculate spike train correlations
        cc_(nn,1,1) = mean(corr(outpute1,'rows','pairwise'),'all','omitnan'); 
        cc_(nn,2,2) = mean(corr(outpute2,'rows','pairwise'),'all','omitnan');
        cc_(nn,1,2) = mean(corr(outpute1,outpute2,'rows','pairwise'),'all','omitnan');
    end
    % Take average spike train correlations
    cc = squeeze(mean(cc_,1));
    
    % Filtered weighted input
    yfilt = ones(floor(tausyn/dt),1);
    y_fee1 = filter(yfilt,1,ye1temp); y_fee2 = filter(yfilt,1,ye2temp);
    wb_fe1 = y_fee1*Wee/tausyn; wb_fi1 = zeros(size(wb_fe1));
    wb_fe2 = y_fee2*Wee/tausyn; wb_fi2 = zeros(size(wb_fe2));

    % Sim V1 and V2
    input.feedfowarde = wb_fe1; input.feedfowardi = wb_fi1;
    V1J = cifv(t,rp,0,input,theta); VcellJ{rr,1} = V1J;
    input.feedfowarde = wb_fe2; input.feedfowardi = wb_fi2;
    V2J = cifv(t,rp,0,input,theta); VcellJ{rr,2} = V2J;
    
    % Save V to plots
    V1plot(:,rr,1) = V1J(:,1)-mean(V1J(:));
    V2plot(:,rr,1) = V2J(:,1)-mean(V2J(:));

    % Calc Jitter Corr
    corrout(rr,1) = mean(diag(corr(V1J,V2J)),'all');

    % Instant corr is given by corr using Jitter so they are equal
    rhoe = mean([cc(1,1),cc(2,2)]);
    theta.corrinfo.ae1 = 0; theta.corrinfo.ae2 = (1/rhoe)-1;
    theta.corrinfo.ai1 = 0; theta.corrinfo.ai2 = (1/rhoe)-1;

    % Theoretical calculations
    thetaTh(1).r.ee = re;
    thetaTh(1).rhoWe = rhoe; 
    rhom(:,1) = cc(1,2).*ones(2,1); thetaTh(1).rhoX.ee = rhom(:,1);
    thetaTh(2) = thetaTh(1);
    [~,~,~,~,~,corrV_] = swstats(thetaTh);
    corrout(rr,2) = corrV_(1,2);

    % Make the voltage to show example
    ys = makecopinput(theta,t,dt,d/2,d/2,corrmat,rp);
    ye1 = ys{1}; ye2 = ys{2}; 

    yfilt = ones(floor(tausyn/dt),1);
    y_fee1 = filter(yfilt,1,ye1); y_fee2 = filter(yfilt,1,ye2);
    wb_fe1 = y_fee1*Wee/tausyn; wb_fi1 = zeros(size(wb_fe1));
    wb_fe2 = y_fee2*Wee/tausyn; wb_fi2 = zeros(size(wb_fe2));

    input.feedfowarde = wb_fe1; input.feedfowardi = wb_fi1;
    V1I = cifv(t,rp,0,input,theta); VcellI{rr,1} = V1I;
    
    input.feedfowarde = wb_fe2; input.feedfowardi = wb_fi2;
    V2I = cifv(t,rp,0,input,theta); VcellI{rr,1} = V2I;

    V1plot(:,rr,2) = V1I(:,1)-mean(V1I(:));
    V2plot(:,rr,2) = V2I(:,1)-mean(V2I(:));
end

%% Plot Results
figure(20); clf;
plot(rbank,corrout)
legend({'Jitter','Instant'},'box','off')
ylim([0,1]); xlim([rbank(1),rbank(end)])
xticks(2:4:26)
xlabel('Rate (Hz)'); ylabel('\rho_V'); 
box off

cl = [0.21299500192233756, 0.5114186851211072, 0.730795847750865;...
      0.3600153787004998, 0.7161860822760476, 0.6655132641291811;...
      0.9817762399077278, 0.6073817762399076, 0.3457900807381776;...
      0.9330257593233372, 0.3913110342176086, 0.27197231833910035;...
      0.8141484044598232, 0.2196847366397539, 0.3048058439061899];

Vplot1J = squeeze(V1plot(:,:,1)); Vplot2J = squeeze(V2plot(:,:,1));
Vplot1I = squeeze(V1plot(:,:,2)); Vplot2I = squeeze(V2plot(:,:,2));

figure(21); clf; hold on;
tbank = 0:1e3:10*1e3;
for i = 1:length(rbank)
    ts = tbank(i); te = tbank(i+1);
    tp = ts:dt:te-dt;
    plot(tp*1e-3,Vplot1I(:,i),'color',cl(i,:))
    plot(tp*1e-3,Vplot1J(:,i),'color',cl(i,:),'linewidth',2)
    if i < length(rbank)
        plot([tp(end),te]*1e-3,[Vplot1I(end,i),Vplot1I(1,i+1)],'color',cl(i,:))
        plot([tp(end),te]*1e-3,[Vplot1J(end,i),Vplot1J(1,i+1)],'color',cl(i,:))
    end
end
xlabel('Time (s)'); ylabel('Voltage'); 

figure(22); clf; hold on;
tbank = 0:1e3:10*1e3;
for i = 1:length(rbank)
    ts = tbank(i); te = tbank(i+1);
    tp = ts:dt:te-dt;
    plot(tp*1e-3,Vplot2I(:,i),'color',cl(i,:))
    plot(tp*1e-3,Vplot2J(:,i),'color',cl(i,:),'linewidth',2)
    if i < length(rbank)
        plot([tp(end),te]*1e-3,[Vplot2I(end,i),Vplot2I(1,i+1)],'color',cl(i,:))
        plot([tp(end),te]*1e-3,[Vplot2J(end,i),Vplot2J(1,i+1)],'color',cl(i,:))
    end
end
xlabel('Time (s)'); ylabel('Voltage'); 

% Make Contour Plots
% Example idxs (5 Hs, 15 Hz, 25 Hz)
idxs = binsearch(rbank,[5,15,25]); 
for i = 1:length(idxs)
    idx = idxs(i);
    [X1, Y1] = meshgrid(linspace(min(Vcell{idx,1}(:)), max(Vcell{idx,1}(:)), 50), linspace(min(Vcell{idx,2}(:)), max(Vcell{idx,2}(:)), 50));
    Z1 = ksdensity([Vcell{idx,1}(:), Vcell{idx,2}(:)], [X1(:), Y1(:)]);
    Z1 = reshape(Z1, size(X1));
    maxDensity1 = min(Z1(:));

    [XJ, YJ] = meshgrid(linspace(min(VcellJ{idx,1}(:)), max(VcellJ{idx,1}(:)), 50), linspace(min(VcellJ{idx,2}(:)), max(VcellJ{idx,2}(:)), 50));
    ZJ = ksdensity([VcellJ{idx,1}(:), VcellJ{idx,2}(:)], [XJ(:), YJ(:)], 'Bandwidth', 0.4);
    ZJ = reshape(ZJ, size(XJ));

    figure(i); clf; hold on
    contour(X1, Y1, Z1, [0,0.01],'LineColor',cl(idx,:));  % 20 contour levels
    contour(XJ, YJ, ZJ, [0,0.02],'linewidth',2,'LineColor',cl(idx,:));  % 20 contour levels
    xlim([0,30]); ylim([0,30])
    xlabel('V1'); ylabel('V2'); 
end


