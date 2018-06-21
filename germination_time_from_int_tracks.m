function [A_fits, germ_times, germ_bool] = germination_time_from_int_tracks(intensity_tracks, position, h, user_profile, framenumbers, ts)

% This function takes intensity tracks and returns the points of germination 
% based on the inflection point of a sigmoidal curve fitted to the tracks.
% It also decides whether a track describes a germinating or
% non-germinating spore based on the drop_threshold.
%
% input:
%           intensity_tracks:                           a n*m matrix containing n intensity tracks of length m 
%
%           position, h, user_profile, framenumbers:    objects containing indices and algorithmic parameters
%
%           ts:                                         vector of timepoints after induction of germination (in minutes)
%
% output:
%           fits:                                       vector of size n of function handles to sigmoidal functions fitting the data
%
%           A_fits:                                     vector of size n of parameter sets of these functions
%
%           germ_times:                                 vector of size n of the extracted times of germination (in min)
%
%           germ_bool:                                  boolean vector of size n, indicating whether spore is classified as germinating or non-germinating


drop_threshold = user_profile.cfg.DROP_THRESH;
graph_output = user_profile.cfg.GRAPH_OUTPUT;
movie = user_profile.cfg.position_numbers(position);
no_of_tracks = size(intensity_tracks,1);


%% Reduce the data to the frames chosen for analysis
if string(user_profile.cfg.frames_to_use) == "all"
    intensity_tracks = intensity_tracks;
    members = ones(1,length(framenumbers));
else
    frames_to_use = str2num(user_profile.cfg.frames_to_use);
    members = ismember(framenumbers, frames_to_use);
    intensity_tracks(:, ~members) = nan;
    ts = ts(members);
end


%% Normalize tracks to range [min(track)/max(track) : 1]
i = 1;
for track = intensity_tracks'
    highest = max(track);
    intensity_tracks(i,:) = track/highest;
    i = i+1;
end

%% Smooth all tracks with weighted mean (1:2:1)
smooth_tracks = zeros(size(intensity_tracks));
for i = 1:size(intensity_tracks,2)
    if i == 1
        before = intensity_tracks(:,i);
    else
        before = intensity_tracks(:,i-1);
    end
    if i == size(intensity_tracks,2)
        after = intensity_tracks(:,i);
    else
        after = intensity_tracks(:,i+1);
    end
    mid = intensity_tracks(:,i);
    if all(isnan([before mid after]))
        smooth_tracks(:,i) = nan;
    else
        smooth_tracks(:,i) = mean([before, mid, mid, after],2, 'omitnan');
    end
end




%% Fit a sigmoidal curve to each track, determine point of half max and define as germination point


% adjust options for fitting algorithm
options = optimoptions('fmincon');
options.Algorithm = 'sqp';
options.MaxIterations = 4000;
options.MaxFunctionEvaluations = 4000;
options.ConstraintTolerance = 1e-13;
options.OptimalityTolerance = 1e-13;
options.StepTolerance = 1e-13;
options.Display = 'off';

% initialize container for germination times, 
% differences between fit plateaus (for classification germ vs nongerm),
% fit functions and
% the fitting function parameters
germ_times = nan(no_of_tracks, 1);
diffs = zeros(no_of_tracks, 1);
A_fits = cell(size(smooth_tracks, 1),1);



for track_no = 1:size(smooth_tracks, 1)
    track = smooth_tracks(track_no,:);
    track = track(find(members));
    % estimate starting parameters for fit
    early_asymp_estimate = max(track);
    late_asymp_estimate = min(track);
    switch_estimate = max(ts)/2;
    growth_rate_estimate = 1;
    A0 = [early_asymp_estimate, late_asymp_estimate, growth_rate_estimate, switch_estimate];
    
    % set lower boundaries for sigmoidal parameters
    A_lb = [0 0 0.1 0];
    % set upper boundaries for sigmoidal parameters
    A_ub = [2 2 5 max(ts)*1.05];
    
    % define function to minimize (sum of squared errors fit - data)
    f = @(A)optifunc(A, ts, track);
    
    % fit the curve
    [A_fit, fval, exitflag] = fmincon(f, A0, [], [], [], [], A_lb, A_ub, [], options);
    
    % store fit function parameters, fit functions and plateau differences
    A_fits{track_no} = A_fit;
    %fits{track_no} = @(x) (A_fit(1) + ((A_fit(2) - A_fit(1)) ./ (1 + exp( -A_fit(3) .* (x - A_fit(4)))).^(1 ./ 0.5)));
        
    diffs(track_no) = A_fit(1) - A_fit(2);
    
    % define half-max value
    half_intensity = (A_fit(1) + A_fit(2))/2;
    
    % extract time of half-max intensity
    
        % calculate with inverse sigmoidal function (only necessary, if a fit parameter V is introduced such that the
        
        % generalized logistic function is defined as
        % y = A(1) + ((A(2) - A(1)) ./ (1 + exp( -A(3) .* (x - A(4)))).^1/A(5));
        % with V = A(5);
        %germ_times(track_no) = inv_sigfunc(A_fit, half_intensity);
          
        %extract from A_fit if generalized logistic function is simplified
        %(V == 1)
        germ_times(track_no) = A_fit(3);
        
    % for visualisation purposes
%     plot(ts, track, 'xb', 'LineWidth', 2)
%     legend('show');
%     hold on
%     plot(ts, sigfunc(A_fit, ts), '-r')
%     plot(x_half_intensity, half_intensity, 'k.', 'MarkerSize', 20, 'LineWidth', 2)
%     legend('track' ,'fit', 'germination point')
%     title(['Track #' num2str(track_no) ' with fit, germination at ' num2str(x_half_intensity) ' minutes'])
%     xlabel('time [min]')
%     ylabel('Intensity')
%     hold off
%     saveas(gcf,['movie_' num2str(movie) '_fit_track_no_' num2str(track_no) '_' original_value '_constrained.svg']);

end

% classification germinating vs non-germinating (should be developed
% further, since drop alone seems not to be sufficient for precise
% classification in noisy data

germ_bool = diffs > drop_threshold;

% optional: graphical output of tracks
if graph_output
    out_folder = h.session_graph_output_directory;
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    if any(~germ_bool)
        plot(ts, intensity_tracks(~germ_bool,find(members))', 'ro-', 'MarkerSize', 4);
    end
    xlim([min(ts)-1 max(ts)+1]);ylim([0 1.05]);
    hold on
    if any(germ_bool)
        plot(ts, intensity_tracks(germ_bool,find(members))', 'bo-', 'MarkerSize', 4);
    end
    h = zeros(2,1);
    h(1) = plot(NaN, NaN, 'r-o');
    h(2) = plot(NaN, NaN, 'b-o');
    legend(h, 'nongerminating spores', 'germinating spores', 'Location', 'SouthWest');
    ylabel('normalized intensity');xlabel('minutes after germination induction');
    title(['Intensity tracks position ' num2str(movie)]);
    hold off
    saveas(gcf, [out_folder, filesep, 'tracks_plot_position_' num2str(movie)], 'png');
    close(gcf)
end

    %% Optimization function describing the sum of squared errors between sigmoidal function sigfunc and data exp_data
    function out = optifunc(A, x, exp_data)
       
        fit_data = sigfunc(A,x);
        errors = exp_data - fit_data;
        out = errors * errors';
    end

    %% Generalized sigmoidal function
    function y = sigfunc(A,x)
            
        y = A(1) + ((A(2) - A(1)) ./ (1 + exp( -A(3) .* (x - A(4)))));
        
        % variable names (A, K, B...) as seen in https://en.wikipedia.org/wiki/Generalised_logistic_function
        % A(1): lower asymptote A
        % A(2): upper asymptote K
        % A(3): growth rate B
        % A(4): relates to a left-right shifting M
    end

    %% The inverse of the generalized sigmoidal function, used to determine point of germination
    function x = inv_sigfunc(A,y)
        
        x = (log( ((A(2) - A(1)) ./ (y - A(1))) - 1)) ./ -A(3) + A(4);
    end

end


