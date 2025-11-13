clear
clc
boresight_gain = 2;
L_vec = [8];
f_vec = [5e9 10e9 15e9 20e9 25e9];
antenna_L_vec = [0.4 0.6 0.8 1 1.2];
M_vec = [7];
load('Final_Shape_Line_Y_P.mat')
load('final_locs_defined_montecarlo_100.mat');
Result_OPT_Shape_Vec_avg = zeros(length(antenna_L_vec), length(f_vec), length(M_vec));
Result_OPT_Shape_Vec_MRT = zeros(length(antenna_L_vec), length(f_vec), length(M_vec));
Result_OPT_Line_Vec_avg_new = zeros(length(antenna_L_vec), length(f_vec), length(M_vec));
Result_OPT_Line_Vec_MRT_new = zeros(length(antenna_L_vec), length(f_vec), length(M_vec));
P_tx_max = 1;

%% BF
totit = length(M_locs);
for lll = 1:length(antenna_L_vec)
    for fff = 1:length(f_vec) 
        f_c = f_vec(fff);
        lambda = 3e8/f_vec(fff);
        for mmm = 1:1
            xMin = 0; xMax = L_vec(mmm);
            yMin = 0; yMax = L_vec(mmm);
            zMin = 0; zMax = 4;
            xDelta = xMax - xMin;
            yDelta = yMax - yMin; 
            zDelta = zMax/2 - zMin;                
            for tot = 1:totit
                M = M_vec(mmm);
                xx = M_locs(tot, :, :);
                x = reshape(xx, [M, 3]);
                Y_opt_f = Y_Result_Shape_Vec_avg{lll, fff, mmm};
                P_tr = P_tr_vec_Shape{lll, fff, mmm};
                
                channel_vec_opt = Do_Channels(Y_opt_f, x, boresight_gain, lambda);
                [P_SDP, P_MRT]  = Do_Beamforming(channel_vec_opt, P_tx_max, P_tr);
           
                Result_OPT_Shape_Vec_avg(lll, fff, mmm) = Result_OPT_Shape_Vec_avg(lll, fff, mmm) + P_SDP/totit;
                Result_OPT_Shape_Vec_MRT(lll, fff, mmm) = Result_OPT_Shape_Vec_MRT(lll, fff, mmm) + P_MRT/totit;

                Y_opt_f = Y_Result_Line_Vec_avg_new{lll, fff, mmm};
                P_tr = P_tr_vec_Line_new{lll, fff, mmm};

                channel_vec_opt = Do_Channels(Y_opt_f, x, boresight_gain, lambda);
                [P_SDP, P_MRT]  = Do_Beamforming(channel_vec_opt, P_tx_max, P_tr);
                
                Result_OPT_Line_Vec_avg_new(lll, fff, mmm) = Result_OPT_Line_Vec_avg_new(lll, fff, mmm) + P_SDP/totit;
                Result_OPT_Line_Vec_MRT_new(lll, fff, mmm) = Result_OPT_Line_Vec_MRT_new(lll, fff, mmm) + P_MRT/totit;

            end

        end
        disp(fff)
    end
    disp(lll)
end


save('New_Shape_Line_100iters')