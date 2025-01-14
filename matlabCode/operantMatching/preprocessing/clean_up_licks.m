function [licksL_cleaned, licksR_cleaned] = clean_up_licks(licksL, licksR, varargin)
    % Function to remove elements from licksL if preceded by smaller licksR within a threshold,
    % and vice versa. Returns cleaned vectors licksL_cleaned and licksR_cleaned.
    %
    % Inputs:
    %   licksL - vector of lick times for the left side (in ms)
    %   licksR - vector of lick times for the right side (in ms)
    %   crosstalk_thresh - time threshold (in ms) for detecting crosstalk
    %
    % Outputs:
    %   licksL_cleaned - cleaned vector of lick times for the left side
    %   licksR_cleaned - cleaned vector of lick times for the right side
    p = inputParser;
    % default parameters if none given
    p.addParameter('crosstalk_thresh', 100);
    p.addParameter('rebound_thresh', 20);
    p.addParameter('plot', 0);
    p.parse(varargin{:});


    % Sort inputs to ensure time order
    licksL = sort(licksL);
    licksR = sort(licksR);

    % Initialize cleaned vectors
    licksL_cleaned = licksL;
    licksR_cleaned = licksR;
    % crosstalk
    % Remove elements in licksL if preceded by smaller elements in licksR within the threshold
    licksL_cleaned = licksL_cleaned(~arrayfun(@(x) any((licksR < x) & (x - licksR <= p.Results.crosstalk_thresh)), licksL_cleaned));

    % Remove elements in licksR if preceded by smaller elements in licksL within the threshold
    licksR_cleaned = licksR_cleaned(~arrayfun(@(x) any((licksL < x) & (x - licksL <= p.Results.crosstalk_thresh)), licksR_cleaned));
    % rebound
    mask_1 = [true, diff(licksL_cleaned)>p.Results.rebound_thresh];
    licksL_cleaned = licksL_cleaned(mask_1);
    mask_1 = [true, diff(licksR_cleaned)>p.Results.rebound_thresh];
    licksR_cleaned = licksR_cleaned(mask_1);
    

    if p.Results.plot
        figure2Wide;
        bins_same = linspace(0, 300, 30);
        bins_diff = linspace(0, 300, 30);
        % before clean up
        all_licks = [licksL, licksR];
        all_licks_id = [zeros(size(licksL)), ones(size(licksR))];
        [all_licks, sorted_ind] = sort(all_licks);
        all_licks_diff = diff(all_licks);
        all_licks_id = all_licks_id(sorted_ind);
        all_licks_id_pre = all_licks_id(1:end-1);
        all_licks_id_post = all_licks_id(2:end);
        

        subplot(2,4,1)
        histogram(all_licks_diff(all_licks_id_post==0 & all_licks_id_pre==0), bins_same, EdgeColor="none");
        title('L_ILI', 'Interpreter','none')
        ylabel('Before clean-up')
        subplot(2,4,2)
        histogram(all_licks_diff(all_licks_id_post==1 & all_licks_id_pre==1), bins_same, EdgeColor="none");
        title('R_ILI', 'Interpreter','none')
        subplot(2,4,3)
        histogram(all_licks_diff(all_licks_id_post==0 & all_licks_id_pre==1), bins_diff, EdgeColor="none");
        title('L-R_ILI', 'Interpreter','none')
        subplot(2,4,4)
        histogram(all_licks_diff(all_licks_id_post==1 & all_licks_id_pre==0), bins_diff, EdgeColor="none");
        title('R-L_ILI', 'Interpreter','none')
        % after clean up
        all_licks = [licksL_cleaned, licksR_cleaned];
        all_licks_id = [zeros(size(licksL_cleaned)), ones(size(licksR_cleaned))];
        [all_licks, sorted_ind] = sort(all_licks);
        all_licks_diff = diff(all_licks);
        all_licks_id = all_licks_id(sorted_ind);
        all_licks_id_pre = all_licks_id(1:end-1);
        all_licks_id_post = all_licks_id(2:end);
        

        subplot(2,4,5)
        histogram(all_licks_diff(all_licks_id_post==0 & all_licks_id_pre==0), bins_same, EdgeColor="none");
        title('L_ILI', 'Interpreter','none')
        ylabel('After clean-up')
        subplot(2,4,6)
        histogram(all_licks_diff(all_licks_id_post==1 & all_licks_id_pre==1), bins_same, EdgeColor="none");
        title('R_ILI', 'Interpreter','none')
        subplot(2,4,7)
        histogram(all_licks_diff(all_licks_id_post==0 & all_licks_id_pre==1), bins_diff, EdgeColor="none");
        title('L-R_ILI', 'Interpreter','none')
        subplot(2,4,8)
        histogram(all_licks_diff(all_licks_id_post==1 & all_licks_id_pre==0), bins_diff, EdgeColor="none");
        title('R-L_ILI', 'Interpreter','none')
    end



end

