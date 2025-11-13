function [Y_opt, P_tr] = Do_SGP_Polygon(x, Y00, d00, sgp_iters, ...
    sgp_thr, opt_thr, boresight_gain, inter_dist, xMax, yMax, zMax)

f_opt = 0;
iter = 1;
M = size(x, 1); N = size(d00, 2);
teta = 2*pi/N;
r0 = (inter_dist/2)/(sin(pi/N));

while iter < sgp_iters 
    dz = zeros(M, 1);
    for i = 1:M
        Xi = reshape(x(i, :), [1 3]);
        dz(i) = norm(zMax - Xi(3));
    end
    cvx_begin gp quiet
        variables t Y(1, 2) P_tr(M, 1) d(M, N)
        maximize(t)
        subject to
            sum(P_tr) <= 1
            for m = 1:M
                f0 = sum(d00(m, :).^(-(boresight_gain + 2)));
                temp_cons = 1;
                for j = 1:N
                    d0 = d00(m, j);
                    alpha = (d0/f0)*(-(boresight_gain + 2))*(d0^(-(boresight_gain + 3)));
                    temp_cons = temp_cons*((d(m, j)./d0)^(alpha));
                end
                 (dz(m).^(-boresight_gain))*t*(1/P_tr(m)) <= f0*(temp_cons)
            end
            for i = 1:M
                for j = 1:N
                    d0 = d00(i,j);
                    f0 = 1 + (d00(i,j)^(-2))*(2*Y00(1, 1)*(x(i, 1) - r0*cos((j-1)*teta)) + ...
                        2*Y00(1, 2)*(x(i, 2) - r0*sin((j-1)*teta)) + 2*zMax*x(i, 3));
                    alphad = -2*(d0/f0)*(d0^(-3))*(2*Y00(1, 1)*(x(i, 1) - r0*cos((j-1)*teta)) + ...
                        2*Y00(1, 2)*(x(i, 2) - r0*sin((j-1)*teta)) + 2*zMax*x(i, 3));
                    alpha1 = 2*(Y00(1,1)/f0)*(d0^(-2))*(x(i, 1) - r0*cos((j-1)*teta));
                    alpha2 = 2*(Y00(1,2)/f0)*(d0^(-2))*(x(i, 2) - r0*sin((j-1)*teta));
                    f0*((d(i,j)/d0)^alphad)*((Y(1,1)/Y00(1,1))^alpha1)*((Y(1,2)/Y00(1,2))^alpha2) ...
                        >= (d(i,j)^(-2))*(Y(1, 1)^2 + Y(1, 2)^2 + zMax^2 + ...
                            (x(i, 1) - r0*cos((j-1)*teta))^2 + (x(i, 2) - r0*sin((j-1)*teta))^2 + x(i, 3)^2)
                end
            end
            (1/sgp_thr).*d00 <= d <= sgp_thr.*d00
            (1/sgp_thr).*Y00 <= Y <= sgp_thr.*Y00
%             (1/yMax).*Y(1, 2) <= ones(1, 1)
%             (1/xMax).*Y(1, 1) <= ones(1, 1)
        cvx_end
        Y00 = Y; d00 = d;
        disp(iter)
        disp(t)
        opt_eps = norm(1 - f_opt/t);
        if iter > 2 && opt_eps < opt_thr
            break
        end
        iter = iter + 1;
        f_opt = t;
end

Y_opt = zeros(N, 3);
r0 = (inter_dist/2)/(sin(pi/N));
teta = 2*pi/N;
for j = 1:N
    Y_opt(j, 1) = Y(1, 1) + r0*cos((j-1)*teta);
    Y_opt(j, 2) = Y(1, 2) + r0*sin((j-1)*teta);
    Y_opt(j, 3) = zMax;
end


end