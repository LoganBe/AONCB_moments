function [corr_, s,covy_full] = corrlag(y,dt,win)

% [corr_, s, covy_full] = corrlag(y,d,twin) Computes the spiking
%                         correlation at time scale of win 
%
% Input:
%       y - binary spiking matrix [T,N]
%       dt - original bin size (ms)
%       win - Correlation time scale (ms)
%
% Output:
%   corr_ - average output spiking correlations
%   s - filtered spike train
%   covy_full - full correlation matrix between inputs

bin_per_newbin = floor(win/dt); % Number bins in win
num_nb = int32(floor(size(y,1) / bin_per_newbin)); 
newsize = num_nb*bin_per_newbin; % New size of spiking matrix

% Reshape spiking to window
ync_reshape = reshape(y(1:newsize,:),bin_per_newbin,num_nb,size(y,2));
s = squeeze(sum(ync_reshape,1));

% Caluclate correlation through E[s^2] - E[s]^2
A = mean(s,2);
covy_full = mean(A.^2) - mean(A).^2;
var_ = var(s,[],1);
corr_ = (mean(covy_full)-mean(var_)/size(y,2))/mean(var_); % Remove diagnols too
