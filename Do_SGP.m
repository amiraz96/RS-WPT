function [Y_opt, P_tr] = Do_SGP(x, Y00, d00, gamma_var_00, sgp_iters, ...
    sgp_thr, opt_thr, deltaaa, boresight_gain, inter_dist, max_L, xMax, yMax, zMax)

f_opt = 0;
iter = 1;
M = size(x, 1); N = size(Y00, 1);
while iter < sgp_iters 
    dz = zeros(M, 1);
    for i = 1:M
        Xi = reshape(x(i, :), [1 3]);
        dz(i) = norm(zMax - Xi(3));
    end
    % Solve the convex approximation problem
    cvx_begin gp quiet
        variables t Y(N, 2) P_tr(M, 1) gamma_var(N, N) d(M, N)
        maximize(t)
        subject to
            sum(P_tr) <= 1
            for m = 1:M
                f0 = sum(d00(m, :).^(-(boresight_gain + 2)));
                temp_cons = 1;
                for j = 1:N
                    d0 = d00(m, j);
                    alpha = (d0/f0)*(-(boresight_gain + 2)*(d0^(-(boresight_gain + 3))));
                    temp_cons = temp_cons*((d(m, j)./d0)^(alpha));
                end
                (dz(m).^(-boresight_gain))*t*(1/P_tr(m)) <= f0*(temp_cons)
            end
            use_L = 0;
            for j = 1:N-1
                use_L = use_L + gamma_var(j, j+1);
            end
            (1/max_L)*use_L <= 1;
            for j = 1:N-1
                for n = j+1:N
                    inter_dist*(gamma_var(j, n)^(-1)) <= 1

                    f0 = (gamma_var_00(j, n)^(-2))*(Y00(j, 1)^2 + Y00(j, 2)^2 + Y00(n, 1)^2 + Y00(n, 2)^2);
                    alphag1 = (gamma_var_00(j,n)/f0)*(-2*(gamma_var_00(j, n)^(-3)))*...
                        (Y00(j, 1)^2 + Y00(j, 2)^2 + Y00(n, 1)^2 + Y00(n, 2)^2);
                    alpha1 = (gamma_var_00(j, n)^(-2))*(Y00(j, 1)/f0)*(2)*(Y00(j, 1));
                    alpha2 = (gamma_var_00(j, n)^(-2))*(Y00(j, 2)/f0)*(2)*(Y00(j, 2));
                    alpha3 = (gamma_var_00(j, n)^(-2))*(Y00(n, 1)/f0)*(2)*(Y00(n, 1));
                    alpha4 = (gamma_var_00(j, n)^(-2))*(Y00(n, 2)/f0)*(2)*(Y00(n, 2));
                    f02 = 1 +  (gamma_var_00(j, n)^(-2))*2*(Y00(j, 1)*Y00(n, 1) + Y00(j, 2)*Y00(n, 2));
                    alphag2 = (gamma_var_00(j,n)/f02)*(-2*gamma_var_00(j, n)^(-3))...
                        *2*((Y00(j, 1)*Y00(n, 1) + Y00(j, 2)*Y00(n, 2)));
                    alpha12 = (Y00(j, 1)/f02)*(2)*(gamma_var_00(j, n)^(-2))*(Y00(n, 1));
                    alpha22 = (Y00(j, 2)/f02)*(2)*(gamma_var_00(j, n)^(-2))*(Y00(n, 2));
                    alpha32 = (Y00(n, 1)/f02)*(2)*(gamma_var_00(j, n)^(-2))*(Y00(j, 1));
                    alpha42 = (Y00(n, 2)/f02)*(2)*(gamma_var_00(j, n)^(-2))*(Y00(j, 2));
                    
                    f0*((Y(j, 1)/Y00(j, 1))^(alpha1))*...
                          ((Y(j, 2)/Y00(j, 2))^(alpha2))*...
                          ((Y(n, 1)/Y00(n, 1))^(alpha3))*...
                          ((Y(n, 1)/Y00(n, 1))^(alpha4))* ...
                          ((gamma_var(j, n)/gamma_var_00(j,n))^alphag1) == ...
                          (f02*((gamma_var(j, n)/gamma_var_00(j, n)))^(alphag2)*...
                            (Y(j, 1)/Y00(j, 1))^(alpha12)*...
                          (Y(j, 2)/Y00(j, 2))^(alpha22)*...
                          (Y(n, 1)/Y00(n, 1))^(alpha32)*...
                          (Y(n, 1)/Y00(n, 1))^(alpha42))

                end
            end
            for i = 1:M
                for j = 1:N
                    d0 = d00(i,j);
                    f0 = 1 + (d00(i,j)^(-2))*(2*Y00(j, 1)*x(i, 1) + 2*Y00(j, 2)*x(i, 2) + 2*zMax*x(i, 3));
                    alphad = -2*(d0/f0)*(d0^(-3))*(2*Y00(j, 1)*x(i, 1) + 2*Y00(j, 2)*x(i, 2) + 2*zMax*x(i, 3));
                    alpha1 = 2*(Y00(j,1)/f0)*(d0^(-2))*x(i,1);
                    alpha2 = 2*(Y00(j,2)/f0)*(d0^(-2))*x(i,2);
                    f0*((d(i,j)/d0)^alphad)*((Y(j,1)/Y00(j,1))^alpha1)*((Y(j,2)/Y00(j,2))^alpha2) ...
                        >= (d(i,j)^(-2))*(Y(j, 1)^2 + Y(j, 2)^2 + zMax^2 + ...
                            x(i, 1)^2 + x(i, 2)^2 + x(i, 3)^2)
                end
            end
            (1/sgp_thr).*d00 <= d <= sgp_thr.*d00
            (1/sgp_thr).*Y00 <= Y <= sgp_thr.*Y00
            (1/sgp_thr).*gamma_var_00 <= gamma_var <= sgp_thr.*gamma_var_00
            (1/yMax).*Y(:, 2) <= ones(N, 1)
            (1/xMax).*Y(:, 1) <= ones(N, 1)
        cvx_end
        if t > 100
            Y = Y00;
            break
        end
        disp(iter)
        disp(t)
        opt_eps = norm(f_opt - t);
        if iter > 2 && opt_eps < opt_thr
            break
        end
        if t < f_opt
            break
        end
        iter = iter + 1;
        f_opt = t;
end

Y_opt = zeros(N, 3);
Y_opt(:, 1:2) = Y;
Y_opt(:, 3) = zMax;


end