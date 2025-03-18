function[momz_,momz2_] = genmom(theta,n)

% DESCRIPTION OF THEORY CAN BE FOUND IN  
% Subthreshold moment analysis of neuronal populations driven by synchronous synaptic inputs
% Becker et al. 2025

%STATCORR_BETA: Theoertical Mean and Variance of V
%[meanV, varV] = statcorr_beta(theta) Computes the mean and variance of an
%                AONCB neuron. If input correlations are present, assumes a
%                Beta distribution for the deFinetti Measure. Using
%                recurrsive itterations to caluclate nth order moments
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
%
% Output:
%   momz_ - Noncentral moments
%


%Load in parameters
Ke = theta.K.ee; %# Excitatory Synapses
Ki = theta.K.ei; %# Inhibitory Synapses
lam_e = theta.r.ee; % Excitatory Rate
lam_i = theta.r.ei; % Inhibitory Rate
W = theta.Wf; We = W(1,1); Wi = W(1,2); %Synaptic Weights, excitatory and inhibitory respectively
Ve = theta.ereversal; % Excitatory Reversal
Vi = theta.ireversal; % Inhibitory Reversal
tauleak = theta.tauleak; % Leak Time Constant

% Defaults
if isfield(theta,'I'); I = theta.I; else; I = 0; end %External Current
if isfield(theta.corrinfo,'corridx'); arho = theta.corrinfo.corridx; else; arho = 0; end
if ~isfield(theta,'ir'); irflag = false; else; irflag = theta.ir; end %Default no intrinisc resistance

%If the rate was inhomogenous, take the average
if length(lam_e) > 1; lam_e = mean(lam_e(:)); end
if length(lam_i) > 1; lam_i = mean(lam_i(:)); end

% Condition set for input type (Allows computation to be quicker for E and I only cases)
if strcmp(theta.corrinfo.corr,'uncorr'); alpha_e = 0; alpha_i = 0; beta_e = 0; beta_i = 0;
else
    alpha_e = theta.corrinfo.ae1; beta_e = theta.corrinfo.ae2;  %Excitatory Beta() Paramters
    alpha_i = theta.corrinfo.ai1; beta_i = theta.corrinfo.ai2;  %Inhibitory Beta() Paramters
    if Ki == 0 || beta_i == 0; Wi = 0;
    elseif Ke == 0 || beta_e == 0; We = 0;
    end
end

%Adjusted Rates (b) from Raw: b = r*beta*(psi(beta+K)-psi(beta))
if alpha_e > 0  % Use true beta bin distribution (NEVER REALLY USED)
    ftemp = @(k) k.*exp((gammaln(Ke + 1)-gammaln(k + 1)-gammaln(Ke - k + 1)) +...
                           log(beta(k+alpha_e,(beta_e + Ke - k))./beta(alpha_e,beta_e)));
    eftemp = sum(ftemp(1:Ke),'omitnan'); 
    be = Ke*lam_e*1e-3/(eftemp);
else
    be = lam_e*beta_e*(psi(beta_e+Ke)-psi(beta_e)).*1e-3; %
end

if alpha_i > 0 % Use true beta bin distribution (NEVER REALLY USED)
    ftemp = @(k) k.*exp((gammaln(Ki + 1)-gammaln(k + 1)-gammaln(Ki - k + 1)) +...
                           log(beta(k+alpha_i,(beta_i + Ki - k))./beta(alpha_i,beta_i)));
    eftemp = sum(ftemp(1:Ki),'omitnan'); 
    bi = Ki*lam_i*1e-3/(eftemp);
else
    bi = lam_i*beta_i*(psi(beta_i+Ki)-psi(beta_i)).*1e-3;
end

%Total Rate seperated based on across correlations: b = a*max(be,bi)+(1-a)*(be+bi)
%[NOTE THIS DOES NOT REALLY WORK FOR arho ~= 0 or 1]
if arho
    b = beta_e*lam_e*(psi(beta_e+Ke+Ki)-psi(beta_e)).*1e-3;
else
    b = sum([be,bi],'omitnan');
end

%Set conditions based on a -> 0 or not. [Typically assume a -> 0]
if alpha_e == 0; den_e = psi(beta_e+Ke)-psi(beta_e); else; den_e = beta(alpha_e,beta_e); end
if alpha_i == 0; den_i = psi(beta_i+Ki)-psi(beta_i); else; den_i = beta(alpha_i,beta_i); end
if alpha_i == 0 && alpha_e == 0; den_ei = psi(beta_e+Ke+Ki)-psi(beta_e); else; den_ei = beta(alpha_e,beta_e); end

%If we want intrinsic voltage dependent conductance. Adjust tau and W
if irflag
    if ~isfield(theta,'meanv'); error("Need Desired Mean V"); end
    mv = theta.meanv;
    gv = theta.irf;
    gv = 1+gv(mv);
else
    gv = 1;
end

if Ke == 0; lam_e = 0; We = 0; end
if Ki == 0; lam_i = 0; Wi = 0; end

%Calculate Mean and Var Seperated by Case and Size of input
if strcmp(theta.corrinfo.corr,"uncorr")
    
    momz2_ = zeros(1,3);
    
    % Calc Mean
    ae1 = lam_e*1e-3*Ke*tauleak*(1-exp(-We/tauleak));
    ai1 = lam_i*1e-3*Ki*tauleak*(1-exp(-Wi/tauleak));
    meanV = (Ve*ae1 + Vi*ai1 + I)/(gv + ae1 + ai1);
    
    alphae = lam_e*1e-3*Ke*tauleak*(1-exp(-We/tauleak)).^2/2;
    alphai = lam_i*1e-3*Ki*tauleak*(1-exp(-Wi/tauleak)).^2/2;

    % Calc Var/Second moment
    varV = (alphae*(Ve-meanV).^2 + alphai*(Vi-meanV).^2)./(gv + ae1 + ai1);
    momz2_(2) = varV;
    
    %Calc Skew/Third moment
    betae = (1/3)*lam_e*1e-3*Ke*tauleak*((Ve-meanV).^3*(1-exp(-We/tauleak)).^3 + 3*(Ve-meanV)*(exp(-2*We/tauleak)-1)*(1-exp(-We/tauleak))*varV);
    betai = (1/3)*lam_i*1e-3*Ki*tauleak*((Vi-meanV).^3*(1-exp(-Wi/tauleak)).^3 + 3*(Vi-meanV)*(exp(-2*Wi/tauleak)-1)*(1-exp(-Wi/tauleak))*varV);
    momz2_(3) = (betae + betai)/(1 + ae1 + ai1);
    
    %}
    % TEMP
    re = lam_e *1e-3; ri = lam_i*1e-3;
    bnorme = re/(re+ri); bnormi = ri/(re+ri); %Coefficients for marginals
    b = (Ke*re+Ki*ri)/2;
    
    % PDFS (e only, i only, ei joint)
    epdf = zeros(Ke+1,1); epdf(2) = 1;
    ipdf = zeros(Ki+1,1); ipdf(2) = 1;

    [KE, KI] = meshgrid(0:Ke,0:Ki); %Kdist   

    fE_val = zeros(size(KE)); fI_val = zeros(size(KI)); 
    fE_val(1,:) = epdf; fI_val(:,1) = ipdf;

    R = @(k1,w1,k2,w2) (k1.*w1.*Ve + k2.*w2.*Vi)./(k1.*w1 + k2.*w2);
    x = @(k1,w1,k2,w2) exp(-(k1.*w1 + k2.*w2)/tauleak);
    
    rhs = fE_val + fI_val; % If no cross corr, all density is on the marginals
   
    % Caluclate n moments recursively
    momz = zeros(n+1,1);
    momz_ = zeros(n,1);
    for i = 1:n 
        if i == 1 % If its the first moment, calc the first NON centered moment`
            fval_1 = sum(R(KE,We,KI,Wi).^i.*(1-x(KE,We,KI,Wi)).^i.*rhs,'all','omitnan');
            momz_(i) = 0;
        else
            fval_1 = sum((R(KE,We,KI,Wi)-momz(2)).^i.*(1-x(KE,We,KI,Wi)).^i.*rhs,'all','omitnan');
        end
        fval_2 = 0;
        for j = 1:(i-2)
            fval_temp = nchoosek(i,j).*sum(x(KE,We,KI,Wi).^(i-j).*(R(KE,We,KI,Wi)-momz(2)).^j...
                                           .*(1-x(KE,We,KI,Wi)).^j.*rhs,'all','omitnan')*momz_(i-j);
            fval_2 = fval_temp + fval_2;
        end
        if i > 1
            temp = i*(momz(2)/(b*tauleak)).*momz_(i-1);
            fval_2 = fval_2 - temp;
        end
        
        fval_3 = sum((1-x(KE,We,KI,Wi).^i).*rhs,'all','omitnan');
        momz(i+1) = (fval_1 + fval_2)./(i./(b*tauleak) + fval_3);
        if i == 1
            momz_(i) = 0;
        else
            momz_(i) = momz(i+1);
        end
    end  
    momz_(1) = momz(2);

else % When there is corr
    bnorme = be/(be+bi); bnormi = bi/(be+bi); %Coefficients for marginals
    
    [KE, KI] = meshgrid(0:Ke,0:Ki); %Kdist   

    % PDFS (e only, i only, ei joint)
    epdf = @(w) exp(gammaln(Ke+1) + gammaln(beta_e+Ke-w) - gammaln(beta_e+Ke) - gammaln(Ke-w+1)-log(w))./den_e; %p_{e,k}
    ipdf = @(w) exp(gammaln(Ki+1) + gammaln(beta_i+Ki-w) - gammaln(beta_i+Ki) - gammaln(Ki-w+1)-log(w))./den_i; %p_{i,l}
    eipdf = @(we,wi) exp(gammaln(Ke+1) + gammaln(Ki+1) + gammaln(alpha_e + we + wi) + gammaln(beta_e+Ke+Ki-we-wi)...
                                           -gammaln(we+1)-gammaln(Ke-we+1)-gammaln(wi+1)-gammaln(Ki-wi+1)-gammaln(alpha_e + beta_e + Ke + Ki))./den_ei; %p_{ei,kl}
    
    fE_val = zeros(size(KE)); fI_val = zeros(size(KI)); 
    fE_val(1,:) = bnorme.*epdf(0:Ke);
    fI_val(:,1) = bnormi.*ipdf(0:Ki);

    R = @(k1,w1,k2,w2) (k1.*w1.*Ve + k2.*w2.*Vi)./(k1.*w1 + k2.*w2);
    x = @(k1,w1,k2,w2) exp(-(w1.*k1 + w2.*k2)/tauleak);
    
    if arho == 0
        rhs = fE_val + fI_val; % If no cross corr, all density is on the marginals
    elseif arho == 1 % Perfectly correlated p_{ei,kl}
        assert(lam_e == lam_i,'re ~= ri | Rates need to be the same between pops | Working on expanding this') 

        fEI_val = eipdf(KE,KI); 
        fEI_val(1,1) = 0; fEI_val(isnan(fEI_val)) = 0; % Remove any funny buisness
        
        rhs = fEI_val; 
    end
    
    % Caluclate n moments recursively
    momz = zeros(n+1,1);
    momz_ = zeros(n,1);
    for i = 1:n 
        if i == 1 % If its the first moment, calc the first NON centered moment
            fval_1 = sum((R(KE,We,KI,Wi)-I).*(1-x(KE,We,KI,Wi)).*rhs,'all','omitnan');
            momz_(i) = 0;
        else
            fval_1 = sum((R(KE,We,KI,Wi)-momz(2)).^i.*(1-x(KE,We,KI,Wi)).^i.*rhs,'all','omitnan');
        end
        fval_2 = 0;
        for j = 1:(i-2)
            fval_temp = nchoosek(i,j).*sum(x(KE,We,KI,Wi).^(i-j).*(R(KE,We,KI,Wi)-momz(2)).^j...
                                            .*(1-x(KE,We,KI,Wi)).^j.*rhs,'all','omitnan')*momz_(i-j);
            fval_2 = fval_temp + fval_2;
        end
        if i > 1
            fval_2 = fval_2 - i*(momz(2)/(b*tauleak)).*momz_(i-1);
        end
        
        fval_3 = sum((1-x(KE,We,KI,Wi).^i).*rhs,'all','omitnan');
        momz(i+1) = (fval_1 + fval_2)./(i./(b*tauleak) + fval_3);
        if i == 1
            momz_(i) = 0;
        else
            momz_(i) = momz(i+1);
        end
    end
    momz_(1) = momz(2);
    momz2_ = [];
end


