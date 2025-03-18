function [yc,ysp,phi] = makecopinput(theta,t,dt,Ne,Ni,corrmat,nneuron,sameneuron)

% MAKECOPINPUT
%   [yc,ysp,phi] = makecopinput(theta,t,dt,npair,nneuron,corrmat)
%                  Generates correlated spiking activity given a matrix of
%                  correlation. Can be done for correlations between both within (E->I)
%                  correlations and between pairs of neurons (E1 -> E2) of neurons and
%                  repeated over nneurons.
%
%   Input:
%       theta = structure with main parameters (Ke, Ki, re, ri, alpha_e/i, beta_e/i
%       t = time vector [1xn]
%       dt = time step [1]
%       npair = number of pairs of neurons [1]
%       nneuron = number of repeats for each neuron [1]
%       corrmat = correlation matrix [2*npair,2*npair]
%
%   Output:
%       yc = spike count. Cell of size [nneuron,2*npair]. Each row in the
%            cell is for another repeat. Each Column is the pop type
%       ysp = raw spike trains. Same structure as yc
%       phi = example spike rates (only for 1 rep) [2*npair,1]
%
% *************************************************************

if nargin < 7; nneuron = 1; sameneuron = 0; end % Defualt number of repeats is 0
if nargin < 8; sameneuron = 0; end

%Mean Activity
mu_e = theta.r.ee.*dt.*1e-3; %Excitatory Rate
mu_i = theta.r.ei.*dt.*1e-3; %Inhibitory Rate

%Max Corr Paramters
b_e = theta.corrinfo.ae2; if isinf(b_e); b_e = 1e8; end
a_e = mu_e*b_e/(1-mu_e); %excitatory paramters

b_i = theta.corrinfo.ai2; if isinf(b_i); b_i = 1e8; end
a_i = mu_i*b_i/(1-mu_i); %inhibitory paramters

d = Ne+Ni;

alpha = [repmat(a_e,Ne,1);repmat(a_i,Ni,1)];
beta = [repmat(b_e,Ne,1);repmat(b_i,Ni,1)];


% Create Cov Mat from Sig Mat
D = eye(d); 
Sigma = D*corrmat*D;

% Check mat is PSD (Needed for transformation)
assert(all(eigs(corrmat) >= -1e-14) & all(isreal(eigs(corrmat))),'Corr Mat Not PSD');

T = length(t);
ysp = cell(nneuron,d); yc = cell(1,d);
if sameneuron
    for i = 1:d
        Z = mvnrnd(zeros(d,1),Sigma,T); %Take mvn samples with cov Sigma
        temp1 = normcdf(Z,zeros(size(Z)),ones(size(Z))); %Get copula given by normal cdf
        phi = betainv(temp1',alpha.*ones(size(temp1')),beta.*ones(size(temp1')))'; %Get rate through inverse beta distribution 

        if length(theta.K) > 1
            Ks = [theta.K(1).ee;theta.K(2).ee; theta.K(1).ei.*ones(npair/2,1)];
        else
            Ks = [theta.K.ee.*ones(npair/2,1); theta.K.ei.*ones(npair/2,1)];
        end
        for nn = 1:nneuron
           ysp{nn,i} = phi(:,i) > rand(length(t),Ks(i));
           yc{i} = [yc{i},sum(ysp{nn,i},2)];
        end
    end
else
    for nn = 1:nneuron
        Z = mvnrnd(zeros(d,1),Sigma,T); %Take mvn samples with cov Sigma
        temp1 = normcdf(Z,zeros(size(Z)),ones(size(Z))); %Get copula given by normal cdf
        phi = betainv(temp1',alpha.*ones(size(temp1')),beta.*ones(size(temp1')))'; %Get rate through inverse beta distribution 

        %Get spikes through Poisson process and rate phi
        if length(theta.K) > 1
            Ks = [theta.K(1).ee;theta.K(2).ee; theta.K(1).ei.*ones(npair/2,1)];
        else
            Ks = [theta.K.ee.*ones(Ne,1); theta.K.ei.*ones(Ni,1)];
        end
        for i = 1:d
           ysp{nn,i} = phi(:,i) > rand(length(t),Ks(i));
           yc{i} = [yc{i},sum(ysp{nn,i},2)];
        end
    end
end

