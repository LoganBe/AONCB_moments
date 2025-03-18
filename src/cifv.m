function [V, oi] = cifv(t, Ne, Ni, input, theta)

% cifv(t,Ne,Ni,input,theta,dispflag): Simulates subthreshold voltage from the
% standard conductance based integrate and fire neuron given by:
%
% tau*dV/dt = (Vl-V) + ge(Ve-V) + gi(Vi-V) + I              (1)
%
% The network does not spike and is purely subthreshold. 
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

tauleak = theta.tauleak; %Leak Time Constant (vectorize for different leaks)
if length(tauleak) > 1; tauleak = [tauleak(1)*ones(1,Ne), tauleak(2)*ones(1,Ni)];
else; tauleak = [tauleak*ones(1,Ne), tauleak*ones(1,Ni)];
end

% SET DEFAULTS
if ~isfield(theta,'trim'); trimidx = 500; else; trimidx = theta.trim; end %Time trim idx
if ~isfield(theta,'ereversal'); Ve_r = 60; else; Ve_r = theta.ereversal; end %Excitatory Reversal Potential
if ~isfield(theta,'ireversal'); Vi_r = -10; else; Vi_r = theta.ireversal; end %Excitatory Reversal Potential
if ~isfield(theta,'lreversal'); Vl = 0; else; Vl = theta.lreversal; end %Leak Voltage (default = 0)
if ~isfield(theta,'ve0'); ve0 = zeros(1,Ne); else; ve0 = theta.ve0; end %Default inital exciatatory v = 0
if ~isfield(theta,'vi0'); vi0 = zeros(1,Ni); else; vi0 = theta.vi0; end %Default inital inhibitory v = 0
if ~isfield(input,'I'); h = 0; else; h = input.I; end %External Input
if ~isfield(theta,'odemethod'); odemethod = "euler"; else; odemethod = theta.odemethod; end
if ~isfield(theta,'memtype'); memtype = 'memsave'; else; memtype = theta.memtype; end
if ~isfield(theta,'nsave'); nsave = N; else; nsave = theta.nsave; end

%Feedfoward Inputs
ff_e = input.feedfowarde; %Exciatory FF Input (1:Ne = EE, Ne+1:end = IE)
ff_i = input.feedfowardi; %Inhibitory FF Input (1:Ne = EI, Ne+1:end = II)

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
%Simulate V
for i = 2:length(t)
   
    %total conductance
    ge_ = ff_e(i,:); gi_ = ff_i(i,:);
    
    %total current
    Ie_ = ge_.*(Ve_r-V_); Ii_ = gi_.*(Vi_r-V_);
    
    %total input
    netinput_ = Ie_ + Ii_;
    
    leak_ = (Vl-V_);
    if strcmpi(odemethod,'rk4')
        fV = @(v) (1./tauleak).*((Vl-v) + g_e(i,:).*(Ve_r-v) + g_i(i,:).*(Vi_r-v) + h(i));
        Vtemp = RK4step(fV,V(i-1,:),dt);
        V(i,:) = Vtemp;
    elseif strcmpi(odemethod,'euler')
        V_ = V_ + (1./tauleak).*(leak_ + netinput_ + h)*dt;
    end

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

%Remove trim info for asymtotic stats
V(1:trimidx,:) = []; 
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







