% Benson Pan
% this function finds the euclidean distance for
% k nearest neighbours for a set of points
function dist = k_nearest(points, x, k)
    e_dist = zeros(1, size(points, 2));
    
    % calc the e distance of all points
    for i = 1:size(points, 2)
        e_dist(i) = sqrt(sum((x' - points(:,i)).^2));
    end
    
    [~, idx] = sort(e_dist);
    k_avg = mean(points(:, idx(1:k)), 2);
    dist = sqrt(sum((x' - k_avg).^2));
end