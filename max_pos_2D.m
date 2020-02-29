function [R,V] = max_pos_2D(g_matrix)

    if size(g_matrix,1)==1 || size(g_matrix,2)==1
        max_value = max(abs(g_matrix));
        RV = find(abs(g_matrix) == max_value);
        R = RV;
        V = RV;
    else
        max_value = max(max(abs(g_matrix)));
        [V,R] = find(abs(g_matrix) == max_value);
        V=V(1);
        R=R(1);

        if R==1
            R = -inf;
        end
        if V==1
            V = -inf;
        end
    end
    
end