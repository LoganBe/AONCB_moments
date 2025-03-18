function [y_fe, y_fi, y_e, y_i] = makeffinput(inputinfo, Ne, Ni, t, dt)

% makeffinput(inputinfo,Ne,Ni,t,dt,tausyn) Create external spiking input.
%                                          Input is summed poisson (or compound poisson) process
%                                          over multiple input cells to same cell type. 
%                                          Input is then filtered for synaptic integration
%
% Input: 
%   Ne -- scalar number of excitatory cells
%   Ni -- scalar number of inhibitory cells
%   t -- time vector [1xn] (ms)
%   dt -- time step (ms)
%   inputinfo -- strcture with input paramters
%           K.ee = number of E->E synapses
%           r.ee = rate of E->E neurons (Hz)
%           K.ie = number of E->I synapses
%           r.ie = rate of E->I neurons (Hz)
%           K.ei = number of I->E synapses
%           r.ei = rate of I->E neurons (Hz)
%           K.ii = number of I->I synapses
%           r.ii = rate of I->I rate
%           filt = spike train filter type: 'box' or 'exp' (default)
%           tausyn = synaptic time constant (ms)
%           syndel = synaptic time delay (bin): default = 0
%           brate = rate type: 'homogenous' (default), 'inhomogenous'
%           corrinfo.corr = string if input is correlated: "corr" or "uncorr" -> (default)
%           corrinfo.dist = string of prior dist on compound poiss: "betabin" (default)
%           corrinfo.ae1 = excitatory correlation alpha parameter
%           corrinfo.ae2 = excitatory correlation beta parameter
%           corrinfo.ai1 = inhibitory correlation alpha parameter
%           corrinfo.ai2 = inhibitory correlation beta parameter
%           corrinfo.corridx = crosscorrelation [0 = no cross, 1 = rhoe = rhoi = rhoei]
%
% Output:
%   ye_fe -- summed spike trains for the [Nee; Nie] excitatory inputs filtered for exponential
%            decay of size [Ne+Ni,length(t)]
%   ye_fi -- summed spike trains for the [Nei; Nii] inhibitory inputs filtered for exponential
%            decay of size [Ne+Ni,length(t)]
%   ye_e -- Raw Excitatory Spike Trains
%   ye_i -- Raw Inhibitory Spike Trains
%
%

%Check Conditions
if ~isfield(inputinfo,'brate'); inputinfo.brate = 'homogenous'; end %Rate type
if ~isfield(inputinfo,'corr'); inputinfo.corr = 'uncorr'; end %Correlation or not
if ~isfield(inputinfo,'filt'); inputinfo.filt = 'box'; end %Synaptic filter type
if ~isfield(inputinfo,'syndel'); syndel = 0; else; syndel = inputinfo.syndel; end %Synaptic delay
if ~isfield(inputinfo,'delay'); delay = 0; else; delay = inputinfo.delay; end %Synaptic delay
if ~isfield(inputinfo.corrinfo,'dist'); dist = "betabin"; else; dist = inputinfo.corrinfo.dist; end %Prior on compound poisson
if ~isfield(inputinfo.K,'ii'); inputinfo.K.ii = 0; inputinfo.r.ii = 0; end %Default I -> I input
if ~isfield(inputinfo.K,'ie'); inputinfo.K.ie = 0; inputinfo.r.ie = 0; end %Default E -> I input

%Load In Paramters
Nt = Ne+Ni; 
tausyn = inputinfo.tausyn;
%Number of cells
Kee = inputinfo.K.ee; %E->E
Kie = inputinfo.K.ie; %E->I
Kei = inputinfo.K.ei; %I->E
Kii = inputinfo.K.ii; %I->I
K = [Kee, Kie, Kei, Kii];

%Rates for poisson process. If rate = 0, change flag
rate_ee = inputinfo.r.ee; rates{1} = rate_ee; %E->E 
rate_ie = inputinfo.r.ie; rates{2} = rate_ie; %E->I
rate_ei = inputinfo.r.ei; rates{3} = rate_ei; %I->E
rate_ii = inputinfo.r.ii; rates{4} = rate_ii; %I->I
%rates = [rate_ee, rate_ie, rate_ei, rate_ii];

%Set rate to 0 if no synaptic inputs
for i = 1:length(K)
    if K(i) == 0; rates{i} = zeros(size(rates{i})); end
end

%Set adjustable rates
bee = rates{1}; bie = rates{2};
bei = rates{3}; bii = rates{4};

%If no inputs -> set all outputs to 0 and exit
if all(rates{1}(:)==0) && all(rates{2}(:) == 0) && all(rates{3}(:) == 0) && all(rates{4}(:) == 0)
    y_fe = zeros(length(t),Nt);
    y_fi = zeros(length(t),Nt);
    y_e = zeros(length(t),Nt);
    y_i = zeros(length(t),Nt);
    return
end

%Make input (Uncorr = Poisson, Corr = Compound Poisson)
switch inputinfo.corrinfo.corr
    case "uncorr"
        y_ee = makepoissinput(Ne,Kee,bee,t,dt); %E->E
        y_ie = makepoissinput(Ni,Kie,bie,t,dt); %E->I
        y_ei = makepoissinput(Ne,Kei,bei,t,dt); %I->E
        y_ii = makepoissinput(Ni,Kii,bii,t,dt); %I->I
    case "corr"
        [y_ee,y_ei] = makecorrinput(Ne,Kee,Kei,bee,bei,t,dt,inputinfo); %Excitatory
        [y_ie,y_ii] = makecorrinput(Ni,Kie,Kii,bie,bii,t,dt,inputinfo); %Inhibitory
end
%y_ei = [zeros(delay,Nt);y_ei(1:end-delay,:)];

% Filter Synaptic Inputs 
if strcmp(inputinfo.filt,'exp') %Exp filter
    %**SHOULD PROB ESTIMATE FIRST POINT THROUGH A GIVEN DISTIRBUTION**
    y_fee = zeros(length(t)+1,Ne); y_fei = zeros(length(t)+1,Ne); 
    y_fie = zeros(length(t)+1,Ni); y_fii = zeros(length(t)+1,Ni); 
    for tt = 2:length(t)+1
        if tt > syndel %Only add filter if we are at time > delay
            y_fee(tt,:) = y_fee(tt-1,:) - y_fee(tt-1,:)/tausyn*dt + y_ee(tt-1-syndel,:);
            y_fei(tt,:) = y_fei(tt-1,:) - y_fei(tt-1,:)/tausyn*dt + y_ei(tt-1-syndel,:);
            y_fie(tt,:) = y_fie(tt-1,:) - y_fie(tt-1,:)/tausyn*dt + y_ie(tt-1-syndel,:);
            y_fii(tt,:) = y_fii(tt-1,:) - y_fii(tt-1,:)/tausyn*dt + y_ii(tt-1-syndel,:);
        end
    end 
    y_fee(1,:) = []; y_fei(1,:) = []; 
    y_fie(1,:) = []; y_fii(1,:) = []; 

elseif strcmp(inputinfo.filt,'box') %Box filter
    yfilt = ones(floor(tausyn/dt),1);
    y_fee = filter(yfilt,1,y_ee);
    y_fei = filter(yfilt,1,y_ei);
    y_fie = filter(yfilt,1,y_ie);
    y_fii = filter(yfilt,1,y_ii);
elseif strcmp(inputinfo.filt,'inst')
    y_fee = y_ee/dt;
    y_fei = y_ei/dt;
    y_fie = y_ie/dt;
    y_fii = y_ii/dt;
end
%Wrap up final output
y_fe = [y_fee, y_fie]; 
y_fi = [y_fei, y_fii];

y_e = [y_ee, y_ie];
y_i = [y_ei, y_ii];



