function [ET_loc_opt_final, bvec_final, obj_init, obj_obj_star, bvec_init] = K_Cheby_Cluster_Fair_min_newapproach(xmax, ymax, zmax, user_loc, ET_num,ET_loc, bor)

%% Area and System Parameters
user_num = size(user_loc, 1); % Set The Number of Users U Have, The Code Will Take Care of The Rest!
xmin = 0; ymin = 0;

%% K-Means with Chebyshev

incend = 50;
inc = 1;

[ET_loc_opt, r_vec, clustervec] = K_Cheby_Cluster(xmax, ymax, zmax, user_loc, ET_num, ET_loc);
r_vec1 = r_vec;

bvec_init = zeros(ET_num, user_num);
loss_vec = zeros(ET_num, user_num);

for i = 1:ET_num
    usidc = clustervec{i};
    for j = 1:length(usidc)
        idx = usidc(j);
        bvec_init(i, idx) = 1;
    end
end
max_it = 10;
for i = 1:ET_num
    for j = 1:user_num
        loss_vec(i, j) = ((norm(ET_loc_opt(i, :) - user_loc(j, :)))^(bor + 2))/((norm(zmax - user_loc(j, 3)))^bor);
    end
end
% objvec = [];
% bbar = [];
% for i = 1:ET_num
%     bbar(i) = sum(bvec(i, :));
%     objvec(i) = bbar(i).*max(loss_vec(i, :).*bvec(i, :));
% end
loss_vec = loss_vec./max(max(loss_vec));

obj_init = max(sum(loss_vec.*bvec_init, 2));
obj_obj_star = obj_init;
ET_loc_opt_final = ET_loc_opt;
bvec_final = bvec_init;


if ET_num > 1 
    tstar = 10000000;
    inc_it = 0;
    while true % Continiue Until Locations Do Not Change !!
            
        cvx_begin quiet
            variable t
            variable bvec(ET_num, user_num) binary
            minimize(t)
            subject to
                for i = 1:ET_num
                    t>=sum(bvec(i,:).*loss_vec(i,:));
                end
                for j = 1:user_num
                    sum(bvec(:, j)) >= 1;
                end
        cvx_end

        bvec = full(bvec);

        for i = 1:ET_num % Optimization for each ET Using CVX Convex Programming !
            cvx_begin quiet
                variable loc_opt(1, 3)
                consval = 0;
                for j = 1:user_num % Check Which Users Are In The ET_i Set !
                    if bvec(i,j) == 1
                        consval = consval + (norm(loc_opt - user_loc(j, :))...
                            /((norm(zmax - user_loc(j, 3)))^(bor/(bor + 2)))); 
                    end
                end
                minimize(consval)
                subject to
                    loc_opt(3) == zmax
                    xmin <= loc_opt(1) <= xmax
                    ymin <= loc_opt(2) <= ymax
            cvx_end
            ET_loc_opt(i, :) = loc_opt; % Store Best Location for ET i
        end


        for i = 1:ET_num
            for j = 1:user_num
                loss_vec(i, j) = ((norm(ET_loc_opt(i, :) - user_loc(j, :)))^(bor + 2))/((norm(zmax - user_loc(j, 3)))^bor);
            end
        end
        loss_vec = loss_vec./max(max(loss_vec));
        


        obj_obj = max(sum(loss_vec.*bvec, 2));


        if norm(1 - obj_obj/tstar) <= 1e-4 || inc_it == max_it
            break
        end

        if obj_obj < obj_obj_star
            obj_obj_star = obj_obj;
            ET_loc_opt_final = ET_loc_opt;
            bvec_final = bvec;
        end
        tstar = obj_obj;

        inc_it = inc_it + 1;
        



    end
end
end