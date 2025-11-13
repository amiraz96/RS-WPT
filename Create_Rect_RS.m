function Y = Create_Rect_RS(start_point, inter_dist, max_L, zMax)

L_row = max_L/4;
num_N = ceil(L_row/(inter_dist));
start_point = start_point - [L_row/2 L_row/2 0];
Y = [];
Y(1, :) = start_point;
i = 1;
for j = 2:num_N - 1
    Y(j, :) = start_point + [(i)*inter_dist 0 0];
    i = i + 1;
end
Y(num_N, :) = start_point + [(num_N-1)*inter_dist 0 0];
Y(num_N + 1, :) = start_point + [(num_N)*inter_dist 0 0];
i = 1;
for j = num_N + 2:2*num_N - 1
    Y(j, :) = start_point + [(num_N)*inter_dist (i)*inter_dist 0];
    i = i + 1;
end
Y(num_N*2, :) = start_point + [(num_N)*inter_dist (num_N - 1)*inter_dist 0];
Y(num_N*2 + 1, :) = start_point + [(num_N)*inter_dist (num_N)*inter_dist 0];
Y(num_N*2 + 2, :) = start_point + [(num_N - 1)*inter_dist (num_N)*inter_dist 0];
i = 1;
for j = 2*num_N + 3:3*num_N
    Y(j, :) = start_point + [(num_N - 1)*inter_dist - (i)*inter_dist (num_N)*inter_dist 0];
    i = i + 1;
end
Y(3*num_N + 1, :) = start_point + [0 (num_N)*inter_dist 0];
i = 1;
for j = 3*num_N + 2:4*num_N
    Y(j, :) = start_point + [0 (num_N)*inter_dist - (i)*inter_dist 0];
    i = i + 1;
end
Y(:, 3) = zMax;

end