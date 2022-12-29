function peaky = isPeak(window, centre_row, centre_col)
    centre_val = window(centre_row, centre_col);
    % make sure we don't mess up with all-zero windows
    if centre_val == 0
        % defintely can't call something a peak if the value is zero
        peaky = false;
        return;
    end
    % TODO: figure out what to do with ties
    peaky = (numel(find(window > centre_val)) == 0); % nothing is greater
end
