function Y = Create_FD(row_num, col_num, start_point, inter_dist)

DMA_loc = [start_point(1) - (col_num/2)*inter_dist ...
                    start_point(2) - (row_num/2)*inter_dist start_point(3)];
Y = zeros(row_num*col_num, 3);

for i = 1:row_num
    for j = 1:col_num
        Y((i-1)*col_num + j, :) = [DMA_loc(1) + (j - 1)*inter_dist ...
            DMA_loc(2) + (i - 1)*inter_dist DMA_loc(3)];
    end
end


end