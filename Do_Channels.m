function channel_vec = Do_Channels(Y, X, boresight_gain, lambda)

N = size(Y, 1); M = size(X, 1);

channel_vec = zeros(N, M);

for i = 1:N
    for m = 1:M
        d3D  = sqrt(sum((reshape(Y(i, :), [1, 3]) - X(m, :)) .^ 2));
        dz = norm(Y(i, 3) - X(m, 3));
        channel_rad = 2*(boresight_gain+1)*...
            ((dz/d3D)^boresight_gain);
        A_channel = sqrt(channel_rad) * (lambda/(4*pi*d3D));
        phase_channel = d3D/lambda;
        channel_vec(i, m) = A_channel*exp(-1i*2*pi*phase_channel);
    end
end

end