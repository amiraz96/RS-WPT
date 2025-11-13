function [Y_FD, Y_SCA, P_tr_SCA, Y_poly, P_tr_poly] = Do_Deploy(cluster, x_center, user_loc, boresight_gain, ...
    zmax, max_L, inter_dist, N_FD)

sgp_iters = 30;
sgp_thr = 1.1;
P_tx_max = 1;
opt_thr = 1e-5;
deltaaa = 1.01;
rho = 0.1;

x = zeros(length(cluster), 3);
for i = 1:length(cluster)
    x(i, :) =  user_loc(cluster(i), :);
end
M = size(x, 1);

%% SCA

L_row = max_L/4;
start_point = [x_center(1)-L_row/2 x_center(2)-L_row/2 zmax];
RS_loc_heu = Create_Rect_RS(start_point, inter_dist, max_L, zmax);
N = length(RS_loc_heu);
max_L = (N - 1)*inter_dist;
Y00 = reshape(RS_loc_heu(:, 1:2), [N, 2]);
d00 = zeros(M, N);
for i = 1:M
    for j = 1:N
        Xi = x(i, :);
        d00(i,j) = norm(RS_loc_heu(j, :) - Xi);
    end
end
gamma_var_00 = zeros(N, N);
for j = 1:N
    for n = 1:N
        gamma_var_00(j,n) = norm(Y00(j, :) - Y00(n, :));
    end
end


% [Y_opt, P_tr_SCA] = Do_SCA(x, Y00, d00, sgp_iters, ...
%     opt_thr, rho, boresight_gain, inter_dist, max_L, zmax, zmax, zmax);

% Y_SCA = Do_Mapping(Y_opt, inter_dist);
Y_SCA = [];
P_tr_SCA =[];
disp('SCA done')


start_point = [x_center(1) x_center(2) zmax];
Y_FD = Create_FD(sqrt(N_FD), sqrt(N_FD), start_point, inter_dist);

disp('FD done')

start_point = [x_center(1) x_center(2)];
Y00 = start_point;
r0 = (inter_dist/2)/(sin(pi/N));
teta = 2*pi/N;
Y_opt = zeros(N, 3);
for j = 1:N
    Y_opt(j, 1) = Y00(1, 1) + r0*cos((j-1)*teta);
    Y_opt(j, 2) = Y00(1, 2) + r0*sin((j-1)*teta);
    Y_opt(j, 3) = zmax;
end
d00 = zeros(M, N);
for i = 1:M
    for j = 1:N
        Xi = reshape(x(i, :), [1 3]);
        d00(i,j) = norm(Y_opt(j, :) - Xi);
    end
end


[Y_poly, P_tr_poly] = Do_SGP_Polygon(x, Y00, d00, sgp_iters, ...
    sgp_thr, opt_thr, boresight_gain, inter_dist, zmax, zmax, zmax);
disp('ploygon done')

end