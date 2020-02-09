% Benson Pan
% this function creates random points
% with correct covariance and mean using matlab's
% randn function and chol (for covariance decomposition)
function points = generate_clusters(n, m, cov)
    points = zeros(2, n);
    for i = 1:n
        points(:, i) = (chol(cov) * randn(2, 1)) + m;
    end
end