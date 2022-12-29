function hough_img = generateHoughAccumulator(img, theta_num_bins, rho_num_bins)
    function rho = getRho(x0, y0, theta)
        %
        % N: number of points
        % T: number of possible theta values = theta_num_bins
        % 
        % x0    : N x 1
        % y0    : N x 1
        % theta : 1 x T
        %
        sinusoid = [sin(theta); cos(theta)];    % 2 x T
        xy_vecs = [-x0, y0];                    % N x 2
        rho = xy_vecs * sinusoid;               % N x T after matrix multiplication
    end
    
    function i_rho = getRhoIndex(rho, rho_min, rho_max)
        %
        % rho               : N x T (output of getRho above)
        % rho_min, rho_max  : scalars
        %
        rho_spacing = (rho_max - rho_min) / (rho_num_bins - 1);
        % index should be between 1 and rho_num_bins, inclusive
        i_rho = 1 + round((rho - rho_min) / rho_spacing);       % still N x T
    end
    
    % set up coordinate system
    % here, origin = middle of image, not top-left corner
    % you can use the top-left corner too, but your bins may be different
    [H, W] = size(img);
    centre_x = floor(W/2);
    centre_y = floor(H/2);
    [x, y] = meshgrid(1:W, 1:H);
    x = x - centre_x;
    y = y - centre_y;
    
    % img is an edge image
    x_edge = x(img > 0);
    y_edge = y(img > 0);

    thetas = linspace(-pi/2, pi/2, theta_num_bins);

    rho_min = -sqrt(H^2 + W^2) / 2;
    rho_max = sqrt(H^2 + W^2) / 2;
    rhos_edges = getRho(x_edge, y_edge, thetas);
    rho_bin_votes = getRhoIndex(rhos_edges, rho_min, rho_max);

    hough_img = zeros(rho_num_bins, theta_num_bins);
    
    % take votes from each edge pixel
    N_edge_px = size(rho_bin_votes, 1);
    theta_bin_idxs = 1:theta_num_bins;
    for i = 1:N_edge_px
        inds = sub2ind(size(hough_img), rho_bin_votes(i,:), theta_bin_idxs);
        hough_img(inds) = hough_img(inds) + 1;
    end
end