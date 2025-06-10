function [nd_points_x, nd_points_y, ind] = find_nondominated(X, Y)
    % Returns near non-dominated points with tolerance epsilon in both X and Y.
    % Assumes minimization problem in 2D.
    
     epsilon = 0;
    % if nargin < 3
    %     epsilon = 0;
    % end

    if size(X, 2) > size(X, 1)
        X = X';
        Y = Y';
    end

    n = length(X);
    is_dominated = false(n, 1);

    % Compare each point to all others
    for i = 1:n
        xi = X(i);
        yi = Y(i);

        % Mark as dominated if any point is better by more than epsilon
        dominates = (X <= xi - epsilon) & (Y <= yi - epsilon);
        dominates(i) = false;

        if any(dominates)
            is_dominated(i) = true;
        end
    end

    % Keep near non-dominated points
    keep = ~is_dominated;
    nd_X = X(keep);
    nd_Y = Y(keep);
    nd_ind = find(keep);

    % Sort the kept points by X, then Y for smooth plotting
    nd_points = [nd_X, nd_Y];
    [sorted_points, sort_idx] = sortrows(nd_points, [1, 2]);

    nd_points_x = sorted_points(:,1);
    nd_points_y = sorted_points(:,2);
    ind = nd_ind(sort_idx);
end
