clear
clc

boresight_gain = 2;
L_vec = [8];

% For Fixed Frequency
f_fix = 10e9;
antenna_L_vec = [0.5 1 1.5 2 2.5 3];
N_vec_overL = [36 64 100 144 169 196];
% For Fixed Length
antenna_L_fix = 1.5;
f_vec = [1e9 5e9 10e9 15e9 20e9 25e9];
N_vec_overF = [9 49 100 144 196 256];

M_vec = [7];

load('final_locs_defined_montecarlo_100.mat');
load('SGP_SCA_Results_OverF_OverL.mat');

Result_SGP_Vec_SDP_overL = zeros(length(antenna_L_vec), 1);
Result_SGP_Vec_MRT_overL = zeros(length(antenna_L_vec), 1);
Result_SCA_Vec_SDP_overL = zeros(length(antenna_L_vec), 1);
Result_SCA_Vec_MRT_overL = zeros(length(antenna_L_vec), 1);

Result_SGP_Vec_SDP_overF = zeros(length(f_vec), 1);
Result_SGP_Vec_MRT_overF = zeros(length(f_vec), 1);
Result_SCA_Vec_SDP_overF = zeros(length(f_vec), 1);
Result_SCA_Vec_MRT_overF= zeros(length(f_vec), 1);

%% Over L

P_tx_max = 1;
totit = length(M_locs);
for lll = 1:length(N_vec_overL)
    f_c = f_fix;
    lambda = 3e8/f_c;
    max_L = antenna_L_vec(lll);
    inter_dist = lambda/2;
    for tot = 1:totit
        for mmm = 1:1
            xMin = 0; xMax = L_vec(mmm);
            yMin = 0; yMax = L_vec(mmm);
            zMin = 0; zMax = 4;
            xDelta = xMax - xMin;
            yDelta = yMax - yMin; 
            zDelta = zMax/2 - zMin;
            M = M_vec(mmm);
            xx = M_locs(tot, :, :);
            x = reshape(xx, [M, 3]);
            N1 = N_vec_overL(lll);
            row_num = sqrt(N1);
            col_num = sqrt(N1);
            x_min = xMin; x_max = xMax;
            y_min = yMin; y_max = yMax;
            z_min = zMin; z_max = zMax;
            x_center = reshape(mean(x, 1), [1 3]);
            start_point = [xMax/2 yMax/2 zMax];
            SGP_element_loc = Do_Mapping(Y_Result_SGP_L{lll}, inter_dist);
            SCA_element_loc = Do_Mapping(Y_Result_SCA_L{lll}, inter_dist);

%             figure
%             scatter(SGP_element_loc(:, 1), SGP_element_loc(:, 2), 'k*')
%             hold on
%             scatter(SCA_element_loc(:, 1), SCA_element_loc(:, 2), 'b*')
%             hold on
%             scatter(x(:, 1), x(:, 2), 'r*')
%             xlim([0 8])
%             ylim([0 8])
            
            SGP_channel_vec = Do_Channels(SGP_element_loc, x, boresight_gain, lambda);
            SCA_channel_vec = Do_Channels(SCA_element_loc, x, boresight_gain, lambda);
           
            
            P_tr = P_Result_SGP_L{lll};
            [P_SDP, P_MRT]  = Do_Beamforming(SGP_channel_vec, P_tx_max, P_tr);

            Result_SGP_Vec_MRT_overL(lll) = Result_SGP_Vec_MRT_overL(lll) + P_MRT/totit;
            Result_SGP_Vec_SDP_overL(lll) = Result_SGP_Vec_SDP_overL(lll) + P_SDP/totit;
            
            P_tr = P_Result_SCA_L{lll};
            [P_SDP, P_MRT]  = Do_Beamforming(SCA_channel_vec, P_tx_max, P_tr);

            Result_SCA_Vec_MRT_overL(lll) = Result_SCA_Vec_MRT_overL(lll) + P_MRT/totit;
            Result_SCA_Vec_SDP_overL(lll) = Result_SCA_Vec_SDP_overL(lll) + P_SDP/totit;

        end
        disp('OverL')
    end
    disp('OverL Large Iter')
    disp(lll)
end

%% Over F

P_tx_max = 1;
totit = length(M_locs);
for lll = 1:length(N_vec_overF)
    f_c = f_vec(lll);
    max_L = antenna_L_fix;
    lambda = 3e8/f_c;
    inter_dist = lambda/2;
    for tot = 1:totit
        for mmm = 1:1
            xMin = 0; xMax = L_vec(mmm);
            yMin = 0; yMax = L_vec(mmm);
            zMin = 0; zMax = 4;
            xDelta = xMax - xMin;
            yDelta = yMax - yMin; 
            zDelta = zMax/2 - zMin;
            M = M_vec(mmm);
            xx = M_locs(tot, :, :);
            x = reshape(xx, [M, 3]);
            N1 = N_vec_overF(lll);
            row_num = sqrt(N1);
            col_num = sqrt(N1);
            x_min = xMin; x_max = xMax;
            y_min = yMin; y_max = yMax;
            z_min = zMin; z_max = zMax;
            x_center = reshape(mean(x, 1), [1 3]);
            start_point = [xMax/2 yMax/2 zMax];
            SGP_element_loc = Do_Mapping(Y_Result_SGP_F{lll}, inter_dist);
            SCA_element_loc = Do_Mapping(Y_Result_SCA_F{lll}, inter_dist);
            
%             figure
%             scatter(SGP_element_loc(:, 1), SGP_element_loc(:, 2), 'k*')
%             hold on
%             scatter(SCA_element_loc(:, 1), SCA_element_loc(:, 2), 'b*')
%             hold on
%             scatter(x(:, 1), x(:, 2), 'r*')
%             xlim([0 8])
%             ylim([0 8])

            SGP_channel_vec = Do_Channels(SGP_element_loc, x, boresight_gain, lambda);
            SCA_channel_vec = Do_Channels(SCA_element_loc, x, boresight_gain, lambda);
            
            P_tr = P_Result_SGP_F{lll};
            [P_SDP, P_MRT]  = Do_Beamforming(SGP_channel_vec, P_tx_max, P_tr);

            Result_SGP_Vec_MRT_overF(lll) = Result_SGP_Vec_MRT_overF(lll) + P_MRT/totit;
            Result_SGP_Vec_SDP_overF(lll) = Result_SGP_Vec_SDP_overF(lll) + P_SDP/totit;
            
            P_tr = P_Result_SCA_F{lll};
            [P_SDP, P_MRT]  = Do_Beamforming(SCA_channel_vec, P_tx_max, P_tr);

            Result_SCA_Vec_MRT_overF(lll) = Result_SCA_Vec_MRT_overF(lll) + P_MRT/totit;
            Result_SCA_Vec_SDP_overF(lll) = Result_SCA_Vec_SDP_overF(lll) + P_SDP/totit;

        end
        disp('OverF')
    end
    disp('OverF Large Iter')
    disp(lll)
end


save('SGP_SCA_MonteCarlo_100')
