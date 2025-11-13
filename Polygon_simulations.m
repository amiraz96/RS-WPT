clear 
clc

boresight_gain = 2;
load('final_locs_defined.mat');
L_vec = [8];


% For Fixed Length
antenna_L_fix = 1.5;
f_vec = [1e9 5e9 10e9 15e9 20e9 25e9];

M_vec = [7];
Y_Result_Ploy_F = cell(length(f_vec), 1);
P_Result_Poly_F = cell(length(f_vec), 1);
c = 3e8;
sgp_iters = 30;
sgp_thr = 1.1;
P_tx_max = 1;
opt_thr = 1e-7;
deltaaa = 1.01;


for lll = 1:4
    max_L = antenna_L_fix;
    f_c = f_vec(lll);
    lambda = c/f_c;
    inter_dist = 0.5*lambda;
    for mmm = 1:1
        xMin = 0; xMax = L_vec(mmm);
        yMin = 0; yMax = L_vec(mmm);
        zMin = 0; zMax = 4;
        xDelta = xMax - xMin;
        yDelta = yMax - yMin; 
        zDelta = zMax/2 - zMin;
        M = M_vec(mmm);
        xx = x_total{mmm};
        x = reshape(xx(1, :), [M, 1, 3]); % x coordinates of Poisson points
        x = reshape(x, [M 3]);
        L_row = max_L/4;
        x_center = reshape(mean(x, 1), [1 3]);
        start_point = [x_center(1)-L_row/2 x_center(2)-L_row/2 zMax];
        RS_loc_heu = Create_Rect_RS(start_point, inter_dist, max_L, zMax);
        N = length(RS_loc_heu);
        start_point = [x_center(1) x_center(2)];
        max_L = (N - 1)*inter_dist;
        antenna_length = max_L;
        Y00 = start_point;
        r0 = (inter_dist/2)/(sin(pi/N));
        teta = 2*pi/N;
        Y_opt = zeros(N, 3);
        for j = 1:N
            Y_opt(j, 1) = Y00(1, 1) + r0*cos((j-1)*teta);
            Y_opt(j, 2) = Y00(1, 2) + r0*sin((j-1)*teta);
            Y_opt(j, 3) = zMax;
        end
        d00 = zeros(M, N);
        for i = 1:M
            for j = 1:N
                Xi = reshape(x(i, :), [1 3]);
                d00(i,j) = norm(Y_opt(j, :) - Xi);
            end
        end


        [Y_opt, P_tr] = Do_SGP_Polygon(x, Y00, d00, sgp_iters, ...
            sgp_thr, opt_thr, boresight_gain, inter_dist, xMax, yMax, zMax);
        
        
        Y_Result_Ploy_F{lll} = Y_opt;
        P_Result_Poly_F{lll} = P_tr;
        


        disp(strcat('Done Ploygon Simulation For L ---> :', string(max_L*10),...
            '-- For F Iter ---> :', string(f_c/1e9), ...
            '-- For N Iter ---> :', string(N), ...
            '-- For M Iter ---> :', string(M)))

            
        filename = strcat('Ploygon_Results_overF', '_L', ...
                string(max_L*10), '_F', string(f_c/1e9), ...
                '_M', string(M),'.mat');
        save(filename);
    end
end
        
    




