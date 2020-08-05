function [sum] = square_norm_sum(X, dx, dy)

    column_num = size(X, 2);
    
    sum = 0;
    
    for k = 1 : column_num
%         sum = sum + norm(X(:, k)) ^ 2;
        sum = sum + dx * dy * (X(:, k).' * conj(X(:, k)));
    end
end