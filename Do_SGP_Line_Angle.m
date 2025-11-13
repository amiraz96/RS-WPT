function [Y_opt, P_tr, objval] = Do_SGP_Line_Angle(x, Y00, d00, phi, sgp_iters, ...
    sgp_thr, opt_thr, boresight_gain, inter_dist, xMax, yMax, zMax, N)

cos1 = cos(phi);
sin1 = sin(phi);
f_opt = 0;
iter = 1;
M = size(x, 1); 
N = size(d00, 2);
while iter < sgp_iters 
    dz = zeros(M, 1);
    for i = 1:M
        Xi = reshape(x(i, :), [1 3]);
        dz(i) = norm(zMax - Xi(3));
    end
    % Solve the convex approximation problem
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
%             for i = 1:M
%                 for j = 1:N
%                     d0 = d00(i,j);
%                     a2 = (zMax* - x(i, 3))^2;
%                     f0 = 1 + 2*(d00(i,j)^(-2))*(Y00(1)*x(i, 1) + Y00(2)*x(i, 2) ...
%                         + (j - 1)*x(i, 1)*inter_dist*omega00(1) + (j - 1)*x(i, 2)*inter_dist*omega00(2));
%                     alphad = -4*(d0/f0)*(d0^(-3))*(Y00(1)*x(i, 1) + Y00(2)*x(i, 2) ...
%                         + (j - 1)*x(i, 1)*inter_dist*omega00(1) + (j - 1)*x(i, 2)*inter_dist*omega00(2));
%                     alpha1 = 2*(Y00(1)/f0)*(d0^(-2))*x(i,1);
%                     alpha2 = 2*(Y00(2)/f0)*(d0^(-2))*x(i,2);
%                     alpha3 = 2*(omega00(1)/f0)*(d0^(-2))*x(i,1)*(j - 1)*inter_dist;
%                     alpha4 = 2*(omega00(2)/f0)*(d0^(-2))*x(i,2)*(j - 1)*inter_dist;
%                     f0*((d(i,j)/d0)^alphad)*((Y(1)/Y00(1))^alpha1)*((Y(2)/Y00(2))^alpha2)*....
%                         ((omega(1)/omega00(1))^alpha3)*((omega(2)/omega00(2))^alpha4) ...
%                         >= (d(i,j)^(-2))*(Y(1)^2 + Y(2)^2 + ...
%                             x(1)^2 + x(2)^2 + a2 + 2*(j-1)*inter_dist*Y(1)*omega(1) ...
%                             + 2*(j-1)*inter_dist*Y(2)*omega(2) + ...
%                             ((j - 1)^2)*(inter_dist^2)*(omega(1)^2) + ...
%                             ((j - 1)^2)*(inter_dist^2)*(omega(2)^2))
%                 end
%             end
            for i = 1:M
                for j = 1:N
                    d0 = d00(i,j);
                    f0 = 1 + (d00(i,j)^(-2))*(2*Y00(1, 1)*(x(i, 1) + (floor(N/2) - j)*inter_dist*cos1) + ...
                        2*Y00(1, 2)*(x(i, 2) + (floor(N/2) - j)*inter_dist*sin1) + 2*zMax*x(i, 3));
                    alphad = -2*(d0/f0)*(d0^(-3))*(2*Y00(1, 1)*(x(i, 1) + (floor(N/2) - j)*inter_dist*cos1) + ...
                        2*Y00(1, 2)*(x(i, 2) + (floor(N/2) - j)*inter_dist*sin1) + 2*zMax*x(i, 3));
                    alpha1 = 2*(Y00(1,1)/f0)*(d0^(-2))*(x(i, 1) + (floor(N/2) - j)*inter_dist*cos1);
                    alpha2 = 2*(Y00(1,2)/f0)*(d0^(-2))*(x(i, 2) + (floor(N/2) - j)*inter_dist*sin1);
                    f0*((d(i,j)/d0)^alphad)*((Y(1,1)/Y00(1,1))^alpha1)*((Y(1,2)/Y00(1,2))^alpha2) ...
                        >= (d(i,j)^(-2))*(Y(1, 1)^2 + Y(1, 2)^2 + zMax^2 + ...
                            (x(i, 1) + (floor(N/2) - j)*inter_dist*cos1)^2 + (x(i, 2) + (floor(N/2) - j)*inter_dist*sin1)^2 + x(i, 3)^2)
                end
            end
%             f0 = omega00(1)^2 + omega00(2)^2;
%             alpha1 = 2*(omega00(1)^2)/f0;
%             alpha2 = 2*(omega00(2)^2)/f0;
%             f0*((omega(1)/omega00(1))^alpha1)*((omega(2)/omega00(2))^alpha2) <= 1

            (1/sgp_thr).*d00 <= d <= sgp_thr.*d00
            (1/sgp_thr).*Y00 <= Y <= sgp_thr.*Y00
%             (1/sgp_thr).*omega00 <= omega <= sgp_thr.*omega00
            (1/yMax).*Y(2) <= 1
            (1/xMax).*Y(1) <= 1
%             omega <= ones(1, 2)
%             omega >= ones(1, 2).*0.001
        cvx_end
        Y00 = Y; d00 = d;
        disp(iter)
        disp(t)
        iter = iter + 1;
        opt_eps = abs(f_opt - t);
        if iter > 3 && opt_eps < opt_thr
            break
        end
        f_opt = t;
end

Y_opt = zeros(N, 3);
for j = 1:N
    Y_opt(j, 3) = zMax;
    Y_opt(j, 1) = Y(1) - (floor(N/2) - j)*inter_dist*cos1;
    Y_opt(j, 2) = Y(2) - (floor(N/2) - j)*inter_dist*sin1;
end

obj_vec = zeros(M, 1);
for ii = 1:M
    tempval = 0;
    for jj = 1:N
        tempval = tempval + (norm(Y_opt(jj, :) - x(ii, :)))^(-(boresight_gain + 2));
    end
    obj_vec(ii) = ((zMax - x(ii, 3))^(boresight_gain))*P_tr(ii)*tempval;
end
objval = min(obj_vec);

end