function cropped_line_img = lineSegmentFinder(orig_img, hough_img, hough_threshold)
    % repeated a few steps from lineFinder, because we can't change the
    % signature.
    % it's okay...
    
    % remove weak votes upfront
    hough_img(hough_img < hough_threshold) = 0;

    [H, W] = size(orig_img);
    [N_rho, N_theta] = size(hough_img);
    % non-max suppression window size: chosen experimentally
    % starting point was ~1px and ~1deg
    half_window_size = [int64((3 * N_rho) / max(H, W)), int64((9 * N_theta) / 360)];
    peaky_hough_img = nonMaximalSuppression(hough_img, half_window_size);

    % need to convert bin index to pixels to be able to draw
    centre_x = floor(W / 2);
    centre_y = floor(H / 2);
    rho_min = -sqrt(H^2 + W^2) / 2;
    rho_max = sqrt(H^2 + W^2) / 2;
    rho_spacing = (rho_max - rho_min) / (N_rho - 1);
    theta_min = -pi/2;
    theta_max = pi/2;
    theta_spacing = (theta_max - theta_min) / (N_theta - 1);

    % Need this for line segment detection
    % would have been nice to be able to pass the threshold in as a parameter...
    % well, next time.
    edges = edge(orig_img, 'canny', 0.1);
    % edges = edge(orig_img, 'canny', 0.25);
    
    % bounds on centred coordinates
    xc_min = 1 - centre_x;
    xc_max = W - centre_x;
    yc_min = 1 - centre_y;
    yc_max = H - centre_y;

    % Now draw.
    fig = figure();
    imshow(orig_img);
    for i = 1:N_rho
        for j = 1:N_theta
            if peaky_hough_img(i, j) > 0
                % convert bin index to line parameters
                rho = rho_min + (i-1) * rho_spacing;
                theta = theta_min + (j-1) * theta_spacing;

                % Create 2 points on this line.
                % These 2 should align with the top and bottom horizontal
                % borders of the image.
                if theta ~= 0           % non-horizontal line
                    % bottom edge
                    y0_c = yc_max;
                    x0_c = (y0_c * cos(theta) - rho) / sin(theta);
                    
                    % top edge
                    y1_c = yc_min;
                    x1_c = (y1_c * cos(theta) - rho) / sin(theta);
                else
                    % There's no top or bottom.
                    % Left and right then.
                    %
                    % right edge
                    y0_c = rho;
                    x0_c = xc_max;
                    
                    % left edge
                    y1_c = rho;
                    x1_c = xc_min;
                end
                
                % Keep their pixel counterparts as well
                x0_p = x0_c + centre_x;
                y0_p = y0_c + centre_y;
                x1_p = x1_c + centre_x;
                y1_p = y1_c + centre_y;
                
                % These points might not be in the image. In that case, we
                % find the topmost point which is.
                if theta ~= 0
                    if x0_c < xc_min || x0_c > xc_max
                        % Find a new y0_c such that x0_c becomes 
                        % xc_min or xc_max
                        if x0_c < xc_min
                            x0_c = xc_min;
                        else
                            x0_c = xc_max;
                        end
                        
                        % if the line were vertical and visible, its x
                        % coordinates cannot overshoot
                        assert(cos(theta) ~= 0, ...
                            "Vertical line cannot overshoot in x!");
                        y0_c = (x0_c * sin(theta) + rho) / cos(theta);

                        % update pixel coordinates of end-point
                        x0_p = x0_c + centre_x;
                        y0_p = y0_c + centre_y;
                    end
                    if x1_c < xc_min || x1_c > xc_max
                        % Find a new y0 such that x1 becomes 1 or W
                        if x1_c < xc_min
                            x1_c = xc_min;
                        else
                            x1_c = xc_max;
                        end
                        
                        assert(cos(theta) ~= 0, ...
                            "vertical line cannot overshoot in x!");
                        y1_c = (x1_c * sin(theta) + rho) / cos(theta);
                        
                        % update pixel coordinates of end point
                        y1_p = y1_c + centre_y;
                        x1_p = x1_c + centre_x;
                    end
                else
                    % We already are on the left and right edges.
                    % Only the y-coordinates can overshoot.
                    % but then how is the line visible at all?
                    % 
                    % there's something wrong with this line.
                    fprintf(2, "Overshooting line with theta = 0.\n");
                    continue;   % ignore this line
                end
                
                % use (x0_p, y0_p) and (x1_p, y1_p) to draw the line
                %
                % let's see the edge map along the line joining these
                % points
                %
                d_01 = sqrt((y1_p - y0_p) ^ 2 + (x1_p - x0_p) ^ 2);
                nl = ceil(d_01);
                t = linspace(0, 1, nl);
                pts = zeros(nl, 2);
                % generates 2D locations along the line segment between the
                % 2 end-points
                for i_t = 1:nl
                    pts(i_t,:) = round((1 - t(i_t)) * [x0_p, y0_p] + t(i_t) * [x1_p, y1_p]);
                end
                % pull out edge map values from the generated locations
                line_inds = sub2ind(size(edges), pts(:,2), pts(:,1));
                line_edgemap = edges(line_inds);
                
                % median filter replaces a value with the median of its 
                % neighborhood.
                % with size = 3, it effectively silences isolated values
                % e.g.,
                %       1 1 0   -->     1 1 0       (remains the same)
                %  but  0 1 0   -->     0 0 0
                %
                % We can use this to remove pixels coming from a line at a
                % different angle than what we are looking at now. They
                % will only intersect with it at a single, isolated, point,
                % which the median filter can remove.
                % (at least in theory)
                % line_edgemap = medfilt1(single(line_edgemap), 3);
                %
                % this value seems to work better in practice, at least on
                % hough_3
                line_edgemap = medfilt1(single(line_edgemap), 11);
                
                % now, to draw the line segment
                % it's a subset of the locations we generated above, but
                % which subset?
                % i use a single loop: 
                %   go inwards from the end-points of the line
                %   stop if the density of edge points within the
                %   end-points is above a threshold
                % 
                % the idea is to try finding the longest kind-of-unbroken
                % line with this rho and theta.
                % we capture "unbroken-ness" through the segment's density
                D = zeros(nl, nl);
                have_line = false;
                for k = 1:nl
                    if line_edgemap(k)
                        for l = k+1:nl
                            if line_edgemap(l)
                                segment = line_edgemap(k:l);
                                density = sum(segment) / (l - k);
                                % can try very high densities (0.95) for
                                % reliable segments
                                % but need a way to have multiple segments
                                % from the same line then
                                if density > 0.6
                                    have_line = true;   % there's at least one line
                                    D(k, l) = l - k;    % proportional to length of line
                                end
                            end
                        end
                    end
                end
                [~, ind_max] = max(D(:));
                [kmax, lmax] = ind2sub(size(D), ind_max);
                x0_p = pts(kmax, 1);
                x1_p = pts(lmax, 1);
                y0_p = pts(kmax, 2);
                y1_p = pts(lmax, 2);
                if have_line
                    hold on;
                    line([x0_p; x1_p], [y0_p; y1_p], 'LineWidth', 5);
                end
            end
        end
    end
    % provided in demoMATLABTricksFun.m
    cropped_line_img = saveAnnotatedImg(fig);
    close(fig);
end
