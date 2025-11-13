clear 
clc

boresight_gain = 2;
load('final_locs_defined.mat');
L_vec = [8];


% For Fixed Length
antenna_L_fix = 1.5;
f_vec = [1e9 5e9 10e9 15e9 20e9 25e9];

M_vec = [7];
Y_Result_AngleLine_F = cell(length(f_vec), 1);
P_Result_AngleLine_F = cell(length(f_vec), 1);
c = 3e8;
sgp_iters = 10;
sgp_thr = 1.1;
P_tx_max = 1;
opt_thr = 1e-7;
phi_res = 8;
phivec = pi/phi_res:pi/phi_res:pi;

for lll = 1:length(f_vec)
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
        start_point = [x_center(1) x_center(2) zMax];
        max_L = (N - 1)*inter_dist;
        antenna_length = max_L;
        Y00 = start_point(1:2);
        
        obj_opt = 0;
        for pp = 1:length(phivec)
            d00 = zeros(M, N);
            for i = 1:M
                Xi = reshape(x(i, :), [1 3]);
                for j = 1:N
                    ytmp = [Y00(1) + - (floor(N/2) - j)*inter_dist*cos(phivec(pp)) ...
                        Y00(2) + - (floor(N/2) - j)*inter_dist*sin(phivec(pp)) zMax];
                    d00(i, j) = norm(ytmp - Xi);
                end
            end
            [Y_opt_temp, P_tr_temp, obj_temp] = Do_SGP_Line_Angle(x, Y00, d00, phivec(pp), sgp_iters, ...
                sgp_thr, opt_thr, boresight_gain, inter_dist, xMax, yMax, zMax, N);
            if obj_temp > obj_opt
                Y_opt = Y_opt_temp;
                P_tr = P_tr_temp;
                obj_opt = obj_temp;
            end
        end
        
        Y_Result_AngleLine_F{lll} = Y_opt;
        P_Result_AngleLine_F{lll} = P_tr;
        


        disp(strcat('Done Angle Line Simulation For L ---> :', string(max_L*10),...
            '-- For F Iter ---> :', string(f_c/1e9), ...
            '-- For N Iter ---> :', string(N), ...
            '-- For M Iter ---> :', string(M)))

            
        filename = strcat('Angle_Line_Results_overF', '_L', ...
                string(max_L*10), '_F', string(f_c/1e9), ...
                '_M', string(M),'.mat');
        save(filename);
    end
end
        
    




