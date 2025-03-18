function [corr_, s,covy_full] = corrlag(y,dt,win)
    bin_per_newbin = floor(win/dt);
    num_nb = int32(floor(size(y,1) / bin_per_newbin));
    newsize = num_nb*bin_per_newbin;
    
    ync_reshape = reshape(y(1:newsize,:),bin_per_newbin,num_nb,size(y,2));
    s = squeeze(sum(ync_reshape,1));
    A = mean(s,2);
    covy_full = mean(A.^2) - mean(A).^2;
    var_ = var(s,[],1);
    corr_ = (mean(covy_full)-mean(var_)/size(y,2))/mean(var_);
end