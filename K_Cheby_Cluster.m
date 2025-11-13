function [ET_loc_opt, r_vec, clustervec] = K_Cheby_Cluster(xmax, ymax, zmax, user_loc, ET_num, ET_loc)

%% Area and System Parameters
bor = 2;
user_num = size(user_loc, 1); % Set The Number of Users U Have, The Code Will Take Care of The Rest!
xmin = 0; ymin = 0;


%% K-Means with Chebyshev



r_vec = zeros(ET_num, 1); 
ET_loc_opt = ET_loc;
ET_loc_opt_2 = zeros(ET_num, 3);
inc_end = 20;
inc = 1;
while sum(ET_loc_opt ~= ET_loc_opt_2, 'all') ~= 0 && inc <= inc_end% Continiue Until Locations Do Not Change !!

    ET_loc_opt_2 = ET_loc_opt; % Keep The Previous Result in Mind !

    % Assign The Data Points To The Nearest Cluster Head & Generate A User-
    % Set For Each Cluster !

    d_vec = zeros(ET_num, user_num);
    clustervec = cell(ET_num, 1);
    

    for j = 1:user_num
        for i = 1:ET_num
            d_vec(i, j) = norm(ET_loc_opt(i, :) - user_loc(j, :));
        end
        [a, b] = min(d_vec(:, j));
        clustervec{b}(end + 1) = j;
    end


    % Recompute The Location Of ETs And The Radius For Each Cluster !

    for i = 1:ET_num % Optimization for each ET Using CVX Convex Programming !
        cvx_begin quiet
            variable r 
            variable loc_opt(1, 3)
            minimize r
            subject to
                loc_opt(3) == zmax
                xmin <= loc_opt(1) <= xmax
                ymin <= loc_opt(2) <= ymax
                for j = clustervec{i} % Check Which Users Are In The ET_i Set !
                    norm(loc_opt - user_loc(j, :)) - r <= 0 % Now Take Care of The Constraints !
                end
        cvx_end
        r_vec(i) = r; % Store Best Radius for ET i
        ET_loc_opt(i, :) = loc_opt; % Store Best Location for ET i
    end
%     disp(inc)
    inc = inc + 1;
end



end