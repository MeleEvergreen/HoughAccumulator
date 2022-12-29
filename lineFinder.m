function line_detected_img = lineFinder(orig_img, hough_img, hough_threshold)
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
    
    % bounds on centred coordinates
    xc_min = 1 - centre_x;
    xc_max = W - centre_x;
    yc_min = 1 - centre_y;
    yc_max = H - centre_y;

    fig = figure();
    imshow(orig_img);
    for i = 1:N_rho
        for j = 1:N_theta
            if peaky_hough_img(i, j) > 0
                % convert bin index to line parameters
                rho = rho_min + i * rho_spacing;
                theta = theta_min + j * theta_spacing;
                
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
                
                hold on;
                line([x0_p; x1_p], [y0_p; y1_p], 'LineWidth', 5);                
            end
        end
    end
    % provided in demoMATLABTricksFun.m
    line_detected_img = saveAnnotatedImg(fig);
    close(fig);
end
