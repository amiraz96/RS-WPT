clear
clc

boresight_gain = 2;

load('final_locs_defined_largescale_100.mat');
load('LargeScale_Combine_57_FD_SCA_Poly_CL8.mat');

Combine_Result_SCA_Vec_SDP_CL8 = zeros(ET_num, 1);
Combine_Result_SCA_Vec_MRT_CL8 = zeros(ET_num, 1);

Combine_Result_poly_Vec_SDP_CL8 = zeros(ET_num, 1);
Combine_Result_poly_Vec_MRT_CL8 = zeros(ET_num, 1);

Combine_Result_FD_Vec_SDP_CL8 = zeros(ET_num, 1);
Combine_Result_FD_Vec_MRT_CL8 = zeros(ET_num, 1);
%% Over L

P_tx_max = 1;
totit = size(M_locs, 1);
M = size(M_locs, 2);
f_c = f_fix;
lambda = 3e8/f_c;
max_L = antenna_L_fix;
inter_dist = lambda/2;
for tot = 1:totit
    xx = M_locs(tot, :, :);
    x = reshape(xx, [M, 3]);  
    for cl = 1:ET_num
        cluster_users = clustervec{cl};
        M_cluster = length(cluster_users);
        X = zeros(M_cluster, 3);
        for ii = 1:M_cluster
            X(ii, :) = x(cluster_users(ii), :);
        end
        SCA_element_loc = Y_SCAA{cl};
        FD_element_loc = Y_FDD{cl};
        Poly_element_loc = Y_polyy{cl};
    
        FD_channel_vec = Do_Channels(FD_element_loc, X, boresight_gain, lambda);
        SCA_channel_vec = Do_Channels(FD_element_loc, X, boresight_gain, lambda);
        Poly_channel_vec = Do_Channels(Poly_element_loc, X, boresight_gain, lambda);
        
        P_tr = P_SCAA{cl};
        [P_SDP, P_MRT]  = Do_Beamforming(SCA_channel_vec, P_tx_max, P_tr);
    
        Combine_Result_SCA_Vec_MRT_CL8(cl) = Combine_Result_SCA_Vec_MRT_CL8(cl) + P_MRT/totit;
        Combine_Result_SCA_Vec_SDP_CL8(cl) = Combine_Result_SCA_Vec_SDP_CL8(cl) + P_SDP/totit;
        
        P_tr = P_polyy{cl};
        [P_SDP, P_MRT]  = Do_Beamforming(Poly_channel_vec, P_tx_max, P_tr);
    
        Combine_Result_poly_Vec_MRT_CL8(cl) = Combine_Result_poly_Vec_MRT_CL8(cl) + P_MRT/totit;
        Combine_Result_poly_Vec_SDP_CL8(cl) = Combine_Result_poly_Vec_SDP_CL8(cl) + P_SDP/totit;

        P_tr = (1/M_cluster).*ones(M_cluster, 1);
        [P_SDP, P_MRT]  = Do_Beamforming(FD_channel_vec, P_tx_max, P_tr);
    
        Combine_Result_FD_Vec_MRT_CL8(cl) = Combine_Result_FD_Vec_MRT_CL8(cl) + P_MRT/totit;
        Combine_Result_FD_Vec_SDP_CL8(cl) = Combine_Result_FD_Vec_SDP_CL8(cl) + P_SDP/totit;
    end
end


save('Cluster_MonteCarlo_8.mat')
