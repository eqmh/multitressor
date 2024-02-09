function [trend, z, p_value] = mann_kendall(pix_data)
    n = length(pix_data);
    S = 0;
    
    for i = 1:n-1
        for j = i+1:n
            S = S + sign(pix_data(j) - pix_data(i));
        end
    end
    
    var_S = (n * (n - 1) * (2 * n + 5)) / 18;
    
    if S > 0
        z = (S - 1) / sqrt(var_S);
    elseif S < 0
        z = (S + 1) / sqrt(var_S);
    else
        z = 0;
    end
    
    p_value = 2 * (1 - normcdf(abs(z)));
    
    if p_value <= 0.05
        trend = 'Significant';
    else
        trend = 'Not Significant';
    end
end
