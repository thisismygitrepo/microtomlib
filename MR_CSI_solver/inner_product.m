function [inner_production] = inner_product(X, Y)

    if size(X, 1) ~= size(Y, 1) || size(X, 2) ~= size(Y, 2)
        error('The size of X and Y have to be the same\n')
    end
   
    column_num = size(X, 2);
    
    inner_production = zeros(column_num, 1);
       
    for k = 1 : column_num
        inner_production(k) = (X(:, k).' * conj(Y(:, k)));
    end
end
          