function [peaky_img] = nonMaximalSuppression(img, half_window_size)
    [H, W] = size(img);
    hw = half_window_size(1); vw = half_window_size(2);
    peaks = zeros(size(img));
    peaky_img = img;
    for i = 1:H
        for j = 1:W
            if img(i, j) > 0
                up = int64(max(i - vw, 1));
                down = int64(min(i + vw, H));
                left = int64(max(j - hw, 1));
                right = int64(min(j + hw, W));
                window = img(up:down, left:right);
                peaks(i, j) = isPeak(window, i - up + 1, j - left + 1);
                % 
                % Idea: absorb neighboring values into the detected peak.
                % Didn't work out too well for me.
                %
                % if peaks(i, j)
                %     peaky_img(i, j) = sum(window(:));
                % end
            end
        end
    end
    % non-max suppression
    peaky_img(~peaks) = 0;
end

