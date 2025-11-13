function Y_opt_f = Do_Mapping(Y_opt_f, inter_dist)
ccount = 0;
N = size(Y_opt_f, 1);
while ccount < 100
%     disp(ccount)
    for j = 1:N-1
        for n = j+1:N
            xtet = Y_opt_f(n, 1) - Y_opt_f(j, 1);
            ytet = Y_opt_f(n, 2) - Y_opt_f(j, 2);
            rtet = norm(Y_opt_f(j, :) - Y_opt_f(n, :));
            a11 = [xtet ytet 0];
            a22 = [xtet 0 0];
            tet = acos((a11*transpose(a22))./(norm(a11)*norm(a22)));
            if rtet < inter_dist 
                if xtet < 0
                    Y_opt_f(n, 1) = Y_opt_f(j, 1) - inter_dist*(cos(tet));
                else
                    Y_opt_f(n, 1) = Y_opt_f(j, 1) + inter_dist*(cos(tet));
                end
                if ytet < 0
                    Y_opt_f(n, 2) = Y_opt_f(j, 2) - inter_dist*(sin(tet));
                else
                    Y_opt_f(n, 2) = Y_opt_f(j, 2) + inter_dist*(sin(tet));
                end
            end
        end
        xtet = Y_opt_f(j + 1, 1) - Y_opt_f(j, 1);
        ytet = Y_opt_f(j + 1, 2) - Y_opt_f(j, 2);
        rtet = norm(Y_opt_f(j, :) - Y_opt_f(j+1, :));
        sintet = ytet/rtet;
        a11 = [xtet ytet 0];
        a22 = [xtet 0 0];
        tet = acos((a11*transpose(a22))./(norm(a11)*norm(a22)));
        if mod(ccount, 10) == 0
            addtet = (pi/25);
            tet = acos((a11*transpose(a22))./(norm(a11)*norm(a22))) + addtet;
        end
        if rtet ~= inter_dist 
            if xtet < 0
                Y_opt_f(j+1, 1) = Y_opt_f(j, 1) - inter_dist*(cos(tet));
            else
                Y_opt_f(j+1, 1) = Y_opt_f(j, 1) + inter_dist*(cos(tet));
            end
            if ytet < 0
                Y_opt_f(j+1, 2) = Y_opt_f(j, 2) - inter_dist*(sin(tet));
            else
                Y_opt_f(j+1, 2) = Y_opt_f(j, 2) + inter_dist*(sin(tet));
            end
        end
    end
    ccount = ccount + 1;
    %% Nothing
    
    count = 0;
    for j = 1:N -1
        if (norm(Y_opt_f(j, :) - Y_opt_f(j+1, :))) > 1.01*inter_dist || ...
                (norm(Y_opt_f(j, :) - Y_opt_f(j+1, :))) < 0.99*inter_dist
            count = count + 1;
        end
    end
    for j = 1:N -1
        for l = j + 1:N
            if (norm(Y_opt_f(j, :) - Y_opt_f(l, :))) < 0.99*inter_dist
                count = count + 1;
            end
        end
    end
    if count == 0
        break
    end
end

end