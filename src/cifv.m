function [V, y, oi] = cifv(t, Ne, Ni, input, theta)

% cifv(t,Ne,Ni,input,theta,dispflag): Simulates voltage from the
% standard conductance based integrate and fire neuron given by:
%
% tau*dV/dt = (Vl-V) + ge(Ve-V) + gi(Vi-V) + I              (1)
%
% Spikes can be generated either through a Rate based description given by
% a nonlinear conditional intensity of lam = k * V(V>0)^n or through a
% threshold method with threshold Vth. By default, the network does not
% spike and is subthreshold. 
%
% Input: 
%   t -- time vector with constant sampling frequency [1xn] (ms)
%   Ne -- scalar for number of excitatory neurons in network
%   Ni -- scalar for number of inhibitory neuron in network
%   input -- strcture giving both excitatory feedfoward inputs
%            (feedfowarde) and inhibitory feedfoward inputs (feedfowardi)
%   theta -- strcture with model paratmers:
%            tauleak -- leak time constant (ms)
%            tausyn -- synaptic time constant (ms)
%            trim -- time index for removel to be in asymtotic state [Default 500 ms]
%            ereversal -- excitatory reversal potential [Default = 60 mV]
%            ireversal -- inhibitory reversal potential [Default = -10 mV]
%            lreversal -- leak reversal potential       [Default = 0 mV]
%            ve0 -- initial excitatory voltage          [Default = 0 mV]
%            vi0 -- initial inhibitory voltage          [Default = 0 mV]
%            ref -- refractory period (ms)              [Default = 0 ms]
%            modeltype -- 'rate','thresh',sub'          [Default = 'sub']
%     
%            FOR THRESHOLD MODEL ONLY 
%            vth -- spiking threshold (for threshold model)
%            vreset -- reset voltage for threshold model
%            resettype -- 'set' (constant value) or 'average' (running average) [Default = 'set']
%            FOR RATE MODEL ONLY
%            coef -- scaling paramter for conditional intensity (for rate)
%            power -- nonlinear parameter for conditional intensity (for rate)
%
% Output:
%   V -- voltage over time given by equation 1 [mxn]
%   y -- spike events over time [mxn]
%   oi -- structure with other output information
%            rate = time varying rate
%            netinput = Ie + Ii
%            a = reccurent filteredspike info
%            ff.e = excitatory feedfoward info
%            ff.i = inhibitory feedfoward info
%            r.e = excitatory recurrent info
%            r.i = inhibitory recurrent info
%            g.e = total excitatory conductance (ff.e + r.e)
%            g.i = total inhibitory conductance (ff.i + r.i)
%            I.e = total excitatory current (g.e.*drivingforce)
%            I.i = total inhibitory current (g.i.*drivingforce)
 
% Set parameters
dt = min(diff(t)); %Time step
N = Ne + Ni; %Total Number of Neurons

tausyn = theta.tausyn; %Synaptic Time Constant
tauleak = theta.tauleak; %Leak Time Constant (vectorize for different leaks)
if length(tauleak) > 1; tauleak = [tauleak(1)*ones(1,Ne), tauleak(2)*ones(1,Ni)];
else; tauleak = [tauleak*ones(1,Ne), tauleak*ones(1,Ni)];
end

% SET DEFAULTS
if ~isfield(theta,'trim'); trimidx = 500; else; trimidx = theta.trim; end %Time trim idx
if ~isfield(theta,'model'); modeltype = 'sub'; else; modeltype = theta.model; end % Declare model type
if ~isfield(theta,'ereversal'); Ve_r = 60; else; Ve_r = theta.ereversal; end %Excitatory Reversal Potential
if ~isfield(theta,'ireversal'); Vi_r = -10; else; Vi_r = theta.ireversal; end %Excitatory Reversal Potential
if ~isfield(theta,'lreversal'); Vl = 0; else; Vl = theta.lreversal; end %Leak Voltage (default = 0)
if ~isfield(theta,'WeightRecurrent'); Jr = zeros(N); else; Jr = theta.WeightRecurrent; end %[N,N] Recurrent weight matrix (defult = none)
if ~isfield(theta,'ve0'); ve0 = zeros(1,Ne); else; ve0 = theta.ve0; end %Default inital exciatatory v = 0
if ~isfield(theta,'vi0'); vi0 = zeros(1,Ni); else; vi0 = theta.vi0; end %Default inital inhibitory v = 0
if ~isfield(theta,'ir'); irflag = false; else; irflag = theta.ir; end %Default no intrinisc resistance
if ~isfield(input,'I'); h = 0; else; h = input.I; end %External Input
if ~isfield(theta,'odemethod'); odemethod = "euler"; else; odemethod = theta.odemethod; end
if ~isfield(theta,'syndel');  syndelv = zeros(1,N); else; syndelv = theta.syndel; end
if ~isfield(theta,'memtype'); memtype = 'memsave'; else; memtype = theta.memtype; end
if ~isfield(theta,'nsave'); nsave = N; else; nsave = theta.nsave; end
if ~isfield(theta,'resettype'); resettype = 'Set'; else; resettype = theta.resettype; end

%Define the Spiking Mechanism if need
switch lower(modeltype)
    case 'thresh'
        vth = theta.vth; %Threshold for spiking mechanism
        
        %Set reset type and value
        if ~isfield(theta,'vreset'); Vr = 0; else; Vr = theta.vreset; end
        idx_avg = 1; count = 0;

        %Refractory Information
        if ~isfield(theta,'ref'); ref = 0; else; ref = theta.ref; end
        ref_idx = int32(ref/dt); %set ref time to bin idx
        %When using average, start after some period of data given by index raidx
        if ~isfield(theta,'buffersize'); buffersize = 0; else; buffersize = theta.buffersize; end
        buffer = zeros(buffersize,N);

    case 'rate'
        %Set intensity function
        if ~isfield(theta,'coef'); k = 0.3; else; k = theta.coef; end %Nonlinearity coef1
        if ~isfield(theta,'power'); n = 2; else; n = theta.power; end %Nonlinearity coef2
end

%Feedfoward Inputs
ff_e = input.feedfowarde; %Exciatory FF Input (1:Ne = EE, Ne+1:end = IE)
ff_i = input.feedfowardi; %Inhibitory FF Input (1:Ne = EI, Ne+1:end = II)

%Allocation
y = false(length(t),N); %Spike train matrix
y_ = false(1,N); yflag = 0;
ref_flag = zeros(1,N);

% How to save the memory. Save memory by only storing the important stuff
if strcmpi(memtype,'memfull')
    Nsave = N;
    netinput = zeros(length(t),N); %Network maxtrix
    leak = zeros(length(t),N); %leak matrix

    %excitatory inputs
    r_e = zeros(length(t),N); %reccurent
    g_e = zeros(length(t),N); %conductance
    I_e = zeros(length(t),N); %current

    %inhibitory inputs
    r_i = zeros(length(t),N); %recurrent
    g_i = zeros(length(t),N); %conductance
    I_i = zeros(length(t),N); %current
else
    Nsave = nsave;
end

% Preallocat
V = zeros(length(t),Nsave);
rec_e = zeros(1,N); rec_i = zeros(1,N); 

%Set initial values
V_ = zeros(1,N);
V_(1:Ne) = ve0; V_(Ne+1:end) = vi0;

% Rec Mat Check
Jridx = all(Jr == 0);

%Simulate V
for i = 2:length(t)
    if ~Jridx
        if yflag
            idxspk = find(y_>0);    
            espk = idxspk(idxspk<=Ne); ispk = idxspk(idxspk>Ne);
            
            if ~isempty(espk)
                atempe = y_(espk)*Jr(:,espk)';
            else
                atempe = 0;
            end
            
            if ~isempty(ispk)
                atempi = y_(ispk)*Jr(:,ispk)';
            else
                atempi = 0;
            end
        else
            atempe = 0; atempi = 0;
        end
        rec_e = rec_e - rec_e/tausyn*dt + atempe;
        rec_i = rec_i - rec_i/tausyn*dt + atempi;
    end

    if strcmpi(resettype,'average')
        buffer(idx_avg,ref_flag<=0) = V_(ref_flag<=0);
        idx_avg = mod(idx_avg,buffersize)+1;
        count = min(count+1,buffersize);
        runavg = sum(buffer,1)/count;
    end
    
    %total conductance
    ge_ = (ff_e(i,:) + rec_e);
    gi_ = (ff_i(i,:) + rec_i);
    
    %total current
    Ie_ = ge_.*(Ve_r-V_);
    Ii_ = gi_.*(Vi_r-V_);
    
    %total input
    netinput_ = Ie_ + Ii_;
    
    leak_ = (Vl-V_);
    if strcmpi(odemethod,'rk4')
        fV = @(v) (1./tauleak).*((Vl-v) + g_e(i,:).*(Ve_r-v) + g_i(i,:).*(Vi_r-v) + h(i));
        Vtemp = RK4step(fV,V(i-1,:),dt);
        V(i,ref_flag(i,:)) = Vtemp;
    elseif strcmpi(odemethod,'euler')
        V_(ref_flag<=0) = V_(ref_flag<=0) + (1./tauleak(ref_flag<=0)).*(leak_(ref_flag<=0) + netinput_(ref_flag<=0) + h)*dt;
    end

    switch modeltype
        case 'Rate'
            Vtemp = (V(i,:)-Vl).*(V(i,:) > Vl); %rectified
            rate(i,:) = k*Vtemp.^n; %rate
            y(i,:) = rate(i,:)*dt*1e-3 > rand(1,N); %spike train
        case 'Thresh'
            if ~isempty(find(V_ >= vth,1))
                y_(V_>=vth) = 1; y_(V_<vth) = 0;
                yflag = 1;
                ref_flag(V_>vth) = ref_idx;
                
                %MARKER FOR WHAT METHOD TO USE FOR RESET
                switch lower(resettype)
                    case lower('Set')
                        V_(y_==1) = Vr;
                
                    case lower('Average')
                        if i > buffersize
                            V_(y_==1) = runavg(y_==1);
                        else
                            V_(y_==1) = 0;
                        end
                    case lower('slope')
                        m = 0.5;
                        V(i:i+ref_idx,y(i,:)==1) = m*Vth_;
                end
            else
                y_ = zeros(1,N);
                yflag = 0;
            end
            ref_flag = ref_flag - 1;
    end
    
    y(i,:) = y_;
    V(i,:) = V_(1:nsave);
    if strcmpi(memtype,'memfull')
        netinput(i,:) = netinput_; %Network maxtrix

        %excitatory inputs
        r_e(i,:) = rec_e; %reccurent
        g_e(i,:) = g_e; %conductance
        I_e(i,:) = I_e; %current

        %inhibitory inputs
        r_i(i,:) = rec_i; %recurrent
        g_i(i,:) = g_i; %conductance
        I_i(i,:) = I_i; %current
    end
end

%If extended (do to refractory). Remove extra
if size(V,1) > length(t)
    V(length(t)+1:end,:) = [];
end

%Remove trim info for asymtotic stats
V(1:trimidx,:) = []; 
y(1:trimidx,:) = [];
if strcmpi(memtype,'memfull')
    netinput(1:trimidx,:) = []; oi.netinput = netinput;
    leak(1:trimidx,:) = []; oi.leak = leak;
    
    ff_e(1:trimidx,:) = []; oi.ff.e = ff_e;
    ff_i(1:trimidx,:) = []; oi.ff.i = ff_i;

    r_e(1:trimidx,:) = []; oi.r.e = r_e;
    r_i(1:trimidx,:) = []; oi.r.i = r_i;

    g_e(1:trimidx,:) = []; oi.g.e = g_e;
    g_i(1:trimidx,:) = []; oi.g.i = g_i;

    I_e(1:trimidx,:) = []; oi.I.e = I_e;
    I_i(1:trimidx,:) = []; oi.I.i = I_i; 
end

if strcmpi(modeltype,'sub')
    y = [];
else
    y = sparse(y);
end







