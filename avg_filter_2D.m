% Author: Lin Junyang
% Date: 2019-10-27
% 2D average filter

function mat_ret = avg_filter_2D(matrix_2d, marg)
    
    [M,N] = size(matrix_2d);
    matrix_fil = zeros([M+2*marg,N+2*marg]);
    matrix_fil(1+marg:M+marg,1+marg:N+marg) = matrix_2d;
    mat_ret = matrix_fil;
    
    for m=1+marg:M+marg
        for n=1+marg:N+marg
            mat_ret(m,n) = mean(matrix_fil(m,n-marg:n+marg));
        end
    end
    
    matrix_fil = mat_ret;
    
    for m=1+marg:M+marg
        for n=1+marg:N+marg
            mat_ret(m,n) = mean(matrix_fil(m-marg:m+marg,n));
        end
    end

    mat_ret = mat_ret(1+marg:M+marg,1+marg:N+marg);

end