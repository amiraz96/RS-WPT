function [P_SDP, P_MRT]  = Do_Beamforming(channel_vec, P_tx_max, P_tr)

a = channel_vec;
N = size(channel_vec, 1);
M = size(channel_vec, 2);
Bk = zeros(M, N, N);
for k = 1:M
    bk = a(:, k);
    Bk(k, :, :) = transpose(bk*bk');
end
cvx_begin sdp quiet
    variable W_opt(N,N) hermitian semidefinite
    variable p_r_opt_shape
    maximize(p_r_opt_shape)
    subject to
        for kk = 1:M
            p_r_opt_shape <= real(trace(transpose(W_opt)*reshape(Bk(kk, :, :), [N N])));
        end
        real(trace(W_opt)) <= P_tx_max;
cvx_end
p_r_opt_shape_MRT = zeros(M,1);

w = zeros(N, M);
for m = 1:M
    w(:, m) = sqrt(P_tr(m)*P_tx_max).*a(:, m)./norm(a(:, m));
end
for kk = 1:M
    for m = 1:M
        p_r_opt_shape_MRT(kk) = p_r_opt_shape_MRT(kk) + norm(a(:, kk)'*w(:, m))^2;
    end
end

P_MRT = min(p_r_opt_shape_MRT);
P_SDP = p_r_opt_shape;
end