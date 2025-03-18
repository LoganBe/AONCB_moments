function idx = binsearch(data,desired_value)
%
% binsearch     Find bin of data vector that is closest to a given value
%
%               idx = binsearch(data,desired_value), Find bin of data vector that
%               is closest to a given value. Value can be a vector and will find each one.
%               Code is slow if presented with a large vector of values
%       
%               Uses the minimum in the absoluted difference of the given
%               value to the data then find the first data point that
%               this min dist takes place
%
%Input:
%   data - Vector of given data set [1xn]
%   desired_value - single point or vector of points to locate [1xm] m >= 1
%Output:
%   Indexs of data that correspond to value
%*********************************
idx = zeros(1,length(desired_value));

for i = 1:length(desired_value)

    dist = abs(data-desired_value(i)); %
    minDist = min(dist);
    
    idx(i) = find(dist == minDist,1,'first');
end