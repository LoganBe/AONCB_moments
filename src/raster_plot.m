function [p] = raster_plot(spike_info,dt,color_type,offset,marker_size)
% raster_plot      Creates a raster plot of the given bined data
%
%                  [] = raster_plot(bined_data, nT ,t) creates a plot
%                  showing action potenial times across trials (raster
%                  plot).
% Input:
%   spike_info - [mxn] matrix constisting of either spike times or bined spike
%   info. Time component should match with time vector t
%   nT - Number of trials - scaler equal to m
%   t - time vector of length n
%   data_type - declare the type of data used (bined or spiketimes).
%   Default is bined
%
% Output:
%   Raster plot
%

%Default spike_info is bined
if nargin == 5
    ms = marker_size; if isempty(ms); ms = 6; end
elseif nargin == 4
    ms = 6;
elseif nargin  == 3
    ms = 6; offset = 0;
elseif nargin < 3
    color_type = "k"; ms = 6; offset = 0;
end

color_type = color_type;
nT = size(spike_info,2);
assert(size(spike_info,2) == nT,'nT is incorrect, does not match size of data')

[st_tr, st_idx] = find(spike_info');
p = plot((st_idx-1)*dt,st_tr+offset,'color',color_type,'Marker','.','Linestyle','none','MarkerSize',ms);
