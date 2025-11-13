function [Y_opt, P_tr] = Do_SGP_Line(x, Y00, d00, sgp_iters, ...
    sgp_thr, opt_thr, boresight_gain, inter_dist, xMax, yMax, zMax, N)

f_opt = 0;
iter = 1;
M = size(x, 1); 

while iter < sgp_iters 
    dz = zeros(M, 1);
    for i = 1:M
        Xi = reshape(x(i, :), [1 3]);
        dz(i) = norm(zMax - Xi(3));
    end
    % Solve the convex approximation problem
    cvx_begin gp quiet
        variables t Y(1, 3) P_tr(M, 1) d(M, 1)
        maximize(t)
        subject to
            sum(P_tr) <= 1
            for m = 1:M
                f0 = sum(d00(m, :).^(-(boresight_gain + 2)));
                temp_cons = 1;
                for j = 1:1
                    d0 = d00(m, j);
                    alpha = (d0/f0)*(-(boresight_gain + 2))*(d0^(-(boresight_gain + 3)));
                    temp_cons = temp_cons*((d(m, j)./d0)^(alpha));
                end
                (dz(m).^(-boresight_gain))*t*(1/P_tr(m)) <= f0*(temp_cons)
            end
            for i = 1:M
                for j = 1:1
                    d0 = d00(i,j);
                    f0 = 1 + (d00(i,j)^(-2))*(2*Y00(j, 1)*x(i, 1) + 2*Y00(j, 2)*x(i, 2) + 2*Y00(j, 3)*x(i, 3));
                    alphad = -2*(d0/f0)*(d0^(-3))*(2*Y00(j, 1)*x(i, 1) + 2*Y00(j, 2)*x(i, 2) + 2*Y00(j, 3)*x(i, 3));
                    alpha1 = 2*(Y00(j,1)/f0)*(d0^(-2))*x(i,1);
                    alpha2 = 2*(Y00(j,2)/f0)*(d0^(-2))*x(i,2);
                    alpha3 = 2*(Y00(j,3)/f0)*(d0^(-2))*x(i,3);
                    f0*((d(i,j)/d0)^alphad)*((Y(j,1)/Y00(j,1))^alpha1)*((Y(j,2)/Y00(j,2))^alpha2)*((Y(j,3)/Y00(j,3))^alpha3) ...
                        >= (d(i,j)^(-2))*(Y(j, 1)^2 + Y(j, 2)^2 + Y(j, 3)^2 + ...
                            x(i, 1)^2 + x(i, 2)^2 + x(i, 3)^2)
                end
            end
            (1/sgp_thr).*d00 <= d <= sgp_thr.*d00
            (1/sgp_thr).*Y00 <= Y <= sgp_thr.*Y00
            (1/yMax).*Y(:, 2) <= ones(2, 1)
            (1/xMax).*Y(:, 1) <= ones(2, 1)
            Y(:, 3) == zMax*ones(2, 1)
        cvx_end
        Y00 = Y; d00 = d;
        disp(iter)
        disp(t)
        iter = iter + 1;
        opt_eps = abs(f_opt - t);
        if iter > 4 && opt_eps < opt_thr
            break
        end
        f_opt = t;
end

tetvec = 0:pi/50:2*pi;
Y_opt = zeros(N, 3);
Y_opt_temp = zeros(N, 3);
temp_obj = 0;
for tt = 1:length(tetvec)
    teta = tetvec(tt);
    for j = 1:floor(N/2)
        Y_opt_temp(j, 1) = Y(1, 1) - (j-1)*inter_dist*cos(teta);
        Y_opt_temp(j, 2) = Y(1, 2) - (j-1)*inter_dist*sin(teta);
        Y_opt_temp(j, 3) = zMax;
    end
    for j = floor(N/2) + 1:N
        Y_opt_temp(j, 1) = Y(1, 1) + (j - floor(N/2))*inter_dist*cos(teta);
        Y_opt_temp(j, 2) = Y(1, 2) + (j - floor(N/2))*inter_dist*sin(teta);
        Y_opt_temp(j, 3) = zMax;
    end
    obj_vec = zeros(M, 1);
    for ii = 1:M
        tempval = 0;
        for jj = 1:N
            tempval = tempval + (norm(Y_opt_temp(jj, :) - x(ii, :)))^(-(boresight_gain + 2));
        end
        obj_vec(ii) = ((zMax - x(ii, 3))^(boresight_gain))*P_tr(ii)*tempval;
    end
    objval = min(obj_vec);
    if objval > temp_obj
        temp_obj = objval;
        Y_opt = Y_opt_temp;
    end
end


end