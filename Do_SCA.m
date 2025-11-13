function [Y_opt, P_tr] = Do_SCA(x, Y00, d00, sgp_iters, ...
    opt_thr, rho, boresight_gain, inter_dist, max_L, xMax, yMax, zMax)

f_opt = 0;
iter2 = 1;
M = size(x, 1); N = size(Y00, 1);
P_tr = 1/M.*ones(M, 1);
while iter2 < sgp_iters
    iter = 1;
    while iter < sgp_iters 
        dz = zeros(M, 1);
        for i = 1:M
            Xi = reshape(x(i, :), [1 3]);
            dz(i) = norm(zMax - Xi(3));
        end
        % Solve the convex approximation problem
        cvx_begin quiet
            variables t Y(N, 2) d(M, N)
            maximize(t)
            subject to
                for m = 1:M
                    temp_cons = 0;
                    for j = 1:N
                        norm(Y(j, :) - Y00(j, :)) <= rho
                        d0 = d00(m, j);
                        d(m,j) >= norm([Y(j, :) zMax] - x(m, :))
                        alpha = d0^(-(boresight_gain + 2));
                        temp_cons = temp_cons + alpha - (boresight_gain + 2)*(d0^(-(boresight_gain + 3))).*(d(m, j) - d0);
                    end
                    t*(1/P_tr(m))*((dz(m))^(-boresight_gain)) <= temp_cons
                end
                use_L = 0;
                for j = 1:N-1
                    use_L = use_L + norm(Y(j, :) - Y(j+1, :));
                end
                use_L <= max_L;
                for j = 1:N-1
                    for n = j+1:N
                        norm(Y00(j, :) - Y00(n, :))^2 + 2*(Y00(j, :) - Y00(n, :))*...
                            transpose(Y(j, :) - Y00(j, :)) - 2*(Y00(j, :) - Y00(n, :))*...
                            transpose(Y(n, :) - Y00(n, :)) >= inter_dist^2
                    end
                end
%                 (1/yMax).*Y(:, 2) <= ones(N, 1)
%                 (1/xMax).*Y(:, 1) <= ones(N, 1)
                Y(:, 1) >= 0; Y(:, 2) >= 0
            cvx_end
            Y00 = Y; d00 = d;
            disp(iter)
            disp(t)
            opt_eps = norm(f_opt - t);
            if iter > 4 && opt_eps < opt_thr
                break
            end
            iter = iter + 1;
            f_opt = t;
    end
    Y_temp = zeros(N, 3);
    Y_temp(:, 1:2) = Y;
    Y_temp(:, 3) = zMax;
    cvx_begin quiet
        variables t P_tr(M, 1)
        maximize(t)
        subject to
            for m = 1:M
                temp_cons = 0;
                for j = 1:N
                    temp_cons = temp_cons + norm(Y_temp(j, :) - x(m, :))^(-(boresight_gain + 2));
                end
                t*(1/(dz(m ,1)^(boresight_gain))) <= temp_cons*P_tr(m);
            end
            sum(P_tr) <= 1;
    cvx_end
%     disp('Now Alter')
%     disp(iter2)
%     disp(t)
    opt_eps = norm(f_opt - t);
    if iter2 > 2 && opt_eps < opt_thr
        break
    end
    if t < f_opt
        break
    end
    iter2 = iter2 + 1;
    f_opt = t;
end

Y_opt = zeros(N, 3);
Y_opt(:, 1:2) = Y;
Y_opt(:, 3) = zMax;

end