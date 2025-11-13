clear all
clc
cvx_solver mosek

boresight_gain = 2;
ET_vec = 1:8;

f_fix = 10e9;
antenna_L_fix = 1.5;
N_FD = 64;
max_L = antenna_L_fix;
c = 3e8;
f_c = f_fix;
lambda = c/f_c;
inter_dist = 0.5*lambda;

Y_FD = cell(length(ET_vec), 1);
Y_SCA = cell(length(ET_vec), 1);
Y_poly = cell(length(ET_vec), 1);
P_SCA = cell(length(ET_vec), 1);
P_poly = cell(length(ET_vec), 1);


load('final_locs_defined_largescale_100.mat');

IT_loc = size(M_locs, 1);
Cheby_obj = zeros(length(ET_vec), IT_loc);
opt_obj = zeros(length(ET_vec), IT_loc);
Cheby_locs = cell(length(ET_vec), IT_loc);
Cheby_P = cell(length(ET_vec), IT_loc);
Opt_locs = cell(length(ET_vec), IT_loc);
Opt_P = cell(length(ET_vec), IT_loc);


for iii = 1:3
    
    ET_num = ET_vec(iii);
    
    Y_FDD = cell(length(ET_num), 1);
    Y_SCAA = cell(length(ET_num), 1);
    Y_polyy = cell(length(ET_num), 1);
    P_SCAA = cell(length(ET_num), 1);
    P_polyy = cell(length(ET_num), 1);


    ET_loc = zeros(ET_num, 3);
    for i = 1:ET_num
        ET_loc(i, 1) = (i)*xmax/(ET_num + 1);
        ET_loc(i, 2) = ymax/2;
    end
    ET_loc(:, 3) = zmax;

    for ttt = 57:57
        x = reshape(M_locs(ttt, :), [M_vec 3]);
        user_num = size(x, 1);
        %% Chebychev
    
        [ET_loc_opt, r_vec, clustervec] = K_Cheby_Cluster(L_vec, L_vec, zmax, x, ET_num, ET_loc);
        ET_loc_cheby = ET_loc_opt;
        clustervec_cheby = clustervec;

        %% Fair Part
    
        [ET_loc_opt, bvecopt, obj_init, obj_final, bvec_init] = K_Cheby_Cluster_Fair_min_newapproach(xmax, ymax, zmax, x, ET_num,ET_loc, 2);
        clustervec = cell(ET_num, 1);
        for i = 1:ET_num
            for j = 1:user_num
                if bvecopt(i, j) == 1
                    clustervec{i}(end + 1) = j;
                end
            end
        end

        for i = 1:ET_num
            disp('start for cluster number')
            disp(i)
            [Y_FDD{i}, Y_SCAA{i}, P_SCAA{i}, Y_polyy{i}, P_polyy{i}] = ...
                Do_Deploy(clustervec{i}, ET_loc_opt(i, :), x, boresight_gain, zmax, max_L, inter_dist, N_FD);
        end
        Y_FD{iii} = Y_FDD;
        Y_SCA{iii} = Y_SCAA;
        Y_poly{iii} = Y_polyy;
        P_SCA{iii} = P_SCAA;
        P_poly{iii} = P_polyy;

        


    end
    disp('clusters set')
    disp(iii)
    filename = strcat('Final_Cluster_Deploy_57_OPT_ETnum_', string(ET_num), '.mat');
    save(filename);
end

