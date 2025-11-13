function [P_SDP_FD, P_MRT_FD, P_SDP_RS, P_MRT_RS] = Do_Cluster_BF(cluster, start_point, user_loc, boresight_gain,zmax)


X = zeros(length(cluster), 3);
for i = 1:length(cluster)
    X(i, :) =  user_loc(cluster(i), :);
end
M = size(X, 1);
if M ~= 0
    N = 100;
    row_num = sqrt(N);
    col_num = sqrt(N);
    lambda = 3e8/10e9;
    P_tx_max = 1;
    inter_dist = lambda/2;
    max_L = 1.5;
    P_tr = ones(length(X), 1).*(1/M);
    Y = Create_FD(row_num, col_num, start_point, inter_dist);
        
    channel_vec = Do_Channels(Y, X, boresight_gain, lambda);
    
    [P_SDP_FD, P_MRT_FD]  = Do_Beamforming(channel_vec, P_tx_max, P_tr);
    
    
    Y = Create_Rect_RS(start_point, inter_dist, max_L, zmax);
        
    channel_vec = Do_Channels(Y, X, boresight_gain, lambda);
    
    [P_SDP_RS, P_MRT_RS]  = Do_Beamforming(channel_vec, P_tx_max, P_tr);
else
    P_SDP_FD = 1000;
    P_MRT_FD = 1000;
    P_SDP_RS = 1000;
    P_MRT_RS = 1000;
end

end