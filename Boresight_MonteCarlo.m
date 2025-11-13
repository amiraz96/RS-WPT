clear
clc

boresight_gain = 2;
L_vec = [8];
boresight_vec = [2 4 8 16];

% For Fixed Frequency
f_fix = 10e9;
% For Fixed Length
antenna_L_fix = 1.5;

M_vec = [7];

load('final_locs_defined_montecarlo_100.mat');
load('Boresight_Results.mat');

Boresight_Poly_Vec_SDP = zeros(length(boresight_vec), 1);
Boresight_Poly_Vec_MRT = zeros(length(boresight_vec), 1);
Boresight_Line_Vec_SDP = zeros(length(boresight_vec), 1);
Boresight_Line_Vec_MRT = zeros(length(boresight_vec), 1);
Boresight_SGP_Vec_SDP = zeros(length(boresight_vec), 1);
Boresight_SGP_Vec_MRT = zeros(length(boresight_vec), 1);
Boresight_SCA_Vec_SDP = zeros(length(boresight_vec), 1);
Boresight_SCA_Vec_MRT = zeros(length(boresight_vec), 1);
Boresight_FD_Vec_SDP = zeros(length(boresight_vec), 1);
Boresight_FD_Vec_MRT = zeros(length(boresight_vec), 1);
Boresight_RS_Vec_SDP = zeros(length(boresight_vec), 1);
Boresight_RS_Vec_MRT = zeros(length(boresight_vec), 1);


%% Over B

P_tx_max = 1;
totit = length(M_locs);
for bbb = 1:length(boresight_vec)
    boresight_gain = boresight_vec(bbb);
    f_c = f_fix;
    lambda = 3e8/f_c;
    max_L = antenna_L_fix;
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
            x_min = xMin; x_max = xMax;
            y_min = yMin; y_max = yMax;
            z_min = zMin; z_max = zMax;
            x_center = reshape(mean(x, 1), [1 3]);
            start_point = [xMax/2 yMax/2 zMax];
            Poly_element_loc = Y_Result_Poly_B{bbb};
            Line_element_loc = Y_Result_Line_B{bbb};
            SGP_element_loc = Do_Mapping(Y_Result_SGP_B{bbb}, inter_dist);
            SCA_element_loc = Do_Mapping(Y_Result_SCA_B{bbb}, inter_dist);
            N1 = size(Poly_element_loc, 1);
            row_num = sqrt(N1);
            col_num = sqrt(N1);
            FD_element_loc = Create_FD(row_num, col_num, start_point, inter_dist);
            RS_element_loc = Create_Rect_RS(start_point, inter_dist, max_L, zMax); 

            Poly_channel_vec = Do_Channels(Poly_element_loc, x, boresight_gain, lambda);
            Line_channel_vec = Do_Channels(Line_element_loc, x, boresight_gain, lambda);
            SGP_channel_vec = Do_Channels(SGP_element_loc, x, boresight_gain, lambda);
            SCA_channel_vec = Do_Channels(SCA_element_loc, x, boresight_gain, lambda);
            FD_channel_vec = Do_Channels(FD_element_loc, x, boresight_gain, lambda);
            RS_channel_vec = Do_Channels(RS_element_loc, x, boresight_gain, lambda);

            P_tr = P_Result_Poly_B{bbb};
            [P_SDP, P_MRT]  = Do_Beamforming(Poly_channel_vec, P_tx_max, P_tr);

            Boresight_Poly_Vec_MRT(bbb) = Boresight_Poly_Vec_MRT(bbb) + P_MRT/totit;
            Boresight_Poly_Vec_SDP(bbb) = Boresight_Poly_Vec_SDP(bbb) + P_SDP/totit;
            
            P_tr = P_Result_Line_B{bbb};
            [P_SDP, P_MRT]  = Do_Beamforming(Line_channel_vec, P_tx_max, P_tr);

            Boresight_Line_Vec_MRT(bbb) = Boresight_Line_Vec_MRT(bbb) + P_MRT/totit;
            Boresight_Line_Vec_SDP(bbb) = Boresight_Line_Vec_SDP(bbb) + P_SDP/totit;

            P_tr = P_Result_SGP_B{bbb};
            [P_SDP, P_MRT]  = Do_Beamforming(SGP_channel_vec, P_tx_max, P_tr);

            Boresight_SGP_Vec_MRT(bbb) = Boresight_SGP_Vec_MRT(bbb) + P_MRT/totit;
            Boresight_SGP_Vec_SDP(bbb) = Boresight_SGP_Vec_SDP(bbb) + P_SDP/totit;
            
            if bbb < 4 
                P_tr = P_Result_SCA_B{bbb};
                [P_SDP, P_MRT]  = Do_Beamforming(SCA_channel_vec, P_tx_max, P_tr);
    
                Boresight_SCA_Vec_MRT(bbb) = Boresight_SCA_Vec_MRT(bbb) + P_MRT/totit;
                Boresight_SCA_Vec_SDP(bbb) = Boresight_SCA_Vec_SDP(bbb) + P_SDP/totit;
            end

            P_tr = (1/M).*ones(M, 1);
            [P_SDP, P_MRT]  = Do_Beamforming(FD_channel_vec, P_tx_max, P_tr);

            Boresight_FD_Vec_MRT(bbb) = Boresight_FD_Vec_MRT(bbb) + P_MRT/totit;
            Boresight_FD_Vec_SDP(bbb) = Boresight_FD_Vec_SDP(bbb) + P_SDP/totit;

            P_tr = (1/M).*ones(M, 1);
            [P_SDP, P_MRT]  = Do_Beamforming(RS_channel_vec, P_tx_max, P_tr);

            Boresight_RS_Vec_MRT(bbb) = Boresight_RS_Vec_MRT(bbb) + P_MRT/totit;
            Boresight_RS_Vec_SDP(bbb) = Boresight_RS_Vec_SDP(bbb) + P_SDP/totit;

        end
        disp(tot)
    end
end



save('Boresight_MonteCarlo_100.mat')
