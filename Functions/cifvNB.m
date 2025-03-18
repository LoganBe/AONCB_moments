function [V,input,stats,Vpt] = cifvNB(theta,T,ttrim,n)

%CIFVNB: Instant form of the cifv model (no time bins)
%[V, input, stats] = cifNB(theta,T,ttrim,n) Computes the the voltage of the
%                    cifv neuron whenever there is an input then computes the statistical voltage moments. 
%
% Input to model is genearted through ISIs given by an exponetials
% distribution with rate (lambda). Voltage is then updated given the jump
%                       J = (1-exp(-(We+Wi)/tau))*((WeVe + WiVi)/(We+Wi)-V)
% Voltage than decays to V = V0(1-exp(-(T1-T0)/tau)) at the next time of
% input
%
% Input:
%       theta -- structure with model paratmers given by:
%                theta.K.ee = #E->E Synapses. Scalar
%                theta.r.ee = Excitatory Firing rate (Hz). Scalar
%                theta.ereversal = Excitatory reversal potential (mV). Scalar 
%                theta.K.ei = #I->E Synapses. Scalar
%                theta.r.ei = Inhibitory Firing rate (Hz). Scalar
%                theta.ireversal = Inhibitory reversal potential (mV). Scalar
%                theta.Wf = Synaptic weights for [E->E,I->E; E->I, I->I]. [2x2] mat
%                theta.tauleak = Leak time constant (ms). Scalar
%                theta.Ie = External Current (A). Scalar (If not present, assumes = 0)
%                theta.corrinfo.ae1 = alpha (excitatory) paramter for Beta(). Scalar
%                theta.corrinfo.ae2 = beta (excitatory) paramter for Beta(). Scalar
%                theta.corrinfo.ai1 = alpha (inhibitory) paramter for Beta(). Scalar
%                theta.corrinfo.ai2 = beta (inhibitory) paramter for Beta(). Scalar
%                theta.corrinfo.corridx = arho (Determins if there is cross corr. Assumes = 0). Scalar
%                theta.corrinfo.corr = Correlation indicator. String. "uncorr" or "corr"
%       T -- Max time of input
%       ttrim -- Time to trim off (asymototic state)
%       n -- number of moments to calc. Default = 3
%
% Output:
%   V -- Voltages (mV)
%   input -- excitatory and inhibitory isis
%   stats -- structure with .meanV, .varV, .skewV, and .moms/.cenmoms
%

% Default number for moments to calc
if nargin < 4; n = 3; end 

% Set paramters
Ke = theta.K.ee; Ki = theta.K.ei;
re = theta.r.ee; ri = theta.r.ei;
tauleak = theta.tauleak;
We = theta.Wf(1,1); Wi = theta.Wf(1,2);
Ve = theta.ereversal; Vi = theta.ireversal;
beta_e = theta.corrinfo.ae2; beta_i = theta.corrinfo.ai2;
if ~isfield(theta,'I'); I = 0; else; I = theta.I; end
assert(I == 0,'Warning! Need to adjust NB method with injected current I')

if Ke==0 && Ki == 0
    V = [];
    input = [];
    stats.meanV = 0; stats.varV = 0; stats.skewV = 0;
    Vpt = [];
    return
end
    

% Skip everything in no input
if re == 0 && ri == 0
    V = 0; input = []; stats.meanV = 0; stats.varV = 0; stats.skewV = 0;
    stats.moms = zeros(n,1); 
else
    % Make the input (isis)
    if strcmp(theta.corrinfo.corr,"corr")
        if theta.corrinfo.corridx == 0
            be = re*beta_e*(psi(beta_e+Ke)-psi(beta_e))*1e-3; %adjusted E rate
            isie = exprnd(1/be,3*ceil(T*Ke*re*1e-3),1); % Get some ISIs
            esamp = betabinomrnd(Ke,length(isie),[0,beta_e]); esamp = esamp{1}; % Get jump sizes by betabin dist
            if size(esamp,1) ~= length(isie); esamp = esamp'; end % Make sure right size

            bi = ri*beta_i*(psi(beta_i+Ki)-psi(beta_i))*1e-3; %adjusted I rate
            isii = exprnd(1/bi,3*ceil(T*Ki*ri*1e-3),1); % Get some ISIs
            isamp = betabinomrnd(Ki,length(isii),[0,beta_i]); isamp = isamp{1}; % Get jump sizes by betabin dist
            if size(isamp,1) ~= length(isii); isamp = isamp'; end % Make sure right size
            
            input.e = isie; input.i = isii;
        else
            b = re*beta_e*(psi(beta_e+Ke+Ki)-psi(beta_e))*1e-3; % Adjusted rate
            %p = theta.pfit; b = polyval(p,ri);
            isi = exprnd(1/b,3*ceil(T*(Ke+Ki)*re*1e-3),1); % Get some ISIs
            [esamp,isamp] = betabinombivrnd(Ke,Ki,length(isi),[0,beta_e]); esamp = esamp{1}; isamp = isamp{1}; % Get jump sizes by joint betabin
            if size(esamp,1) ~= length(isi); esamp = esamp'; end %Make sure right size 
            if size(isamp,1) ~= length(isi); isamp = isamp'; end %Make sure right size
            
            input.e = isi; input.i = isi;
        end
    else
        isie = exprnd(1/(Ke*re*1e-3),3*ceil(T*Ke*re*1e-3),1); % Get some ISIs
        esamp = ones(length(isie),1); % Jump is always 1 here (no corr)
        isii = exprnd(1/(Ki*ri*1e-3),3*ceil(T*Ki*ri*1e-3),1); 
        isamp = ones(length(isii),1); % Jump is always 1 here (no corr)
        input.e = isie; input.i = isii;
    end

    % Need to combine isis to one vector
    if theta.corrinfo.corridx == 0
        ste = cumsum(isie); % Get E spike times
        rm = find(ste > T); ste(rm) = []; isie(rm) = []; esamp(rm) = []; %Remove long times
        sti = cumsum(isii); % Get I spike times
        rm = find(sti > T); sti(rm) = []; isii(rm) = []; isamp(rm) = []; %Remove long times

        % Combine samples
        stf = [ste;sti]; 
        [st,idx] = sort(stf); % Sort to order the spike times
        esampf = zeros(size(stf)); isampf = zeros(size(stf)); 
        esampf(idx <= length(ste)) = esamp; isampf(idx > length(ste)) = isamp; % Set input to match the ordered Spike times

        % Get isis for joint inputs
        isi = diff(st); 
        if isempty(isie); isi = [isii(1);isi];
        elseif isempty(isii); isi = [isie(1);isi];
        else; isi = [min(isie(1),isii(1));isi];
        end
    else
        st = cumsum(isi);
        esampf = esamp; isampf = isamp;
    end

    isi(st > T*1e3) = []; % Remove if data is too long
    isi_ = isi(2:end); % Start at 2nd isi

    % Sim V
    mom_ = zeros(n,length(isi_));
    V = zeros(length(isi),1); % V at each spike
    Vleak = I;
    Vpt = zeros(length(isi_),3);
    for i = 1:length(isi_)
        V(i) = (1-exp(-(esampf(i)*We + isampf(i)*Wi)/tauleak)).*((esampf(i)*We.*Ve + isampf(i)*Wi.*Vi)/(esampf(i)*We + isampf(i)*Wi) - Vleak) + Vleak ; % Jump
        Vleak = (V(i)-I)*exp(-isi_(i)/tauleak) + I; % Calculate V after exp decay till tauleak
        %If no current 
        if I == 0
            for j = 1:n
                mom_(j,i) = V(i).^j*tauleak/j*(1-exp(-j*isi_(i)/tauleak)); % Moments
            end
        else
            mom_(1,i) = (V(i)-I)*tauleak*(1-exp(-isi_(i)/tauleak)) + I*isi_(i);
            if n == 2
                mom_(2,i) = (V(i)-I)^2*tauleak/2*(1-exp(-2*isi_(i)/tauleak)) + 2*tauleak*I*(V(i)-I)*(1-exp(-isi_(i)/tauleak)) + I^2*isi_(i);
            end
        end
        % NEED DO FOR SKEW
        
        Vpt(i,1) = Vleak;
        Vpt(i,2) = V(i)*exp(-isi_(i)/(3*tauleak));
        Vpt(i,3) = V(i);
    end

    % Trim the data to get in asymtotic state
    isitidx = find(st(2:end) > ttrim,1,'first');
    mom_(:,1:isitidx-1) = [];
    st_ = cumsum(isi);
    st_ = st_(isitidx+1:end)-ttrim;
    V(1:isitidx-1) = [];

    mom = sum(mom_,2)/st_(end);
    meanV = mom(1);
    varV = mom(2)-meanV.^2;
    cenmom = zeros(3,1);
    cenmom(2) = varV;
    cenmom(3) = mom(3) - 3*mom(2).*meanV + 2*meanV.^3;
    skewV = cenmom(3)./varV.^(3/2);
    
    % Pack output
    stats.meanV = meanV; stats.varV = varV; stats.skewV = skewV;
    stats.moms = mom; stats.cenmom = cenmom;
end