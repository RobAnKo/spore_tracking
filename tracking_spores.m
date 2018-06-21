function tracking_spores(h)
% This function is called by GUI_spore_tracking and extracts intensity values of spores and 
% subsequently their germination times from .dv-files (defined in the GUI)
%
%input:     h                   - complete GUI handle structure containing the input/output/algorithm parameter values
%
%output:    report              - contains:
%                                           intensity tracks,
%                                           fits to those tracks, 
%                                           germination times,
%                                           fluorescence intensity values
%                                           and additional info
%
%           framovie(optional)  - contains:
%                                           bf and fl image raw data(framovie.bf and framovie.fl), 
%                                           the segmentation masks (framovie.spore), 
%                                           and the relabeled tracked spores (framovie.cluster)
%
%           plots(optional)     - shows intensity tracks of tracked spores

%% Transfer algorithmic parameters to user_profile structure to use as input for subfunctions
user_profile = h.user_profile;
user_profile.cfg.positions_to_use = h.Germination_Analysis.uicontrols.PositionsField.String;
user_profile.cfg.dvfiles_to_use = h.Germination_Analysis.uicontrols.DVFilesField.String;
user_profile.cfg.frames_to_use = h.Germination_Analysis.uicontrols.FramesField.String;
user_profile.cfg.MIN_SIZE = str2double(h.Germination_Analysis.uicontrols.SporeMinSize.String);
user_profile.cfg.MAX_SIZE = str2double(h.Germination_Analysis.uicontrols.SporeMaxSize.String);
user_profile.cfg.INT_THRESH = str2double(h.Germination_Analysis.uicontrols.SporeInt.String);
user_profile.cfg.PREPROCESS_DARK_CELLS = str2double(h.Germination_Analysis.uicontrols.DarkCells.String);
user_profile.cfg.PREPROCESS_DOME_HEIGHT = str2double(h.Germination_Analysis.uicontrols.DomeHeight.String);
user_profile.cfg.PREPROCESS_CONTRAST = str2double(h.Germination_Analysis.uicontrols.Contrast.String);
user_profile.cfg.BF_CHANNEL = str2num(h.Germination_Analysis.uicontrols.BFField.String);
user_profile.cfg.FL_CHANNELS = str2num(h.Germination_Analysis.uicontrols.FlField.String);
user_profile.cfg.EXP_VAR = str2double(h.Germination_Analysis.uicontrols.ExpVarInt.String);
user_profile.cfg.DROP_THRESH = str2double(h.Germination_Analysis.uicontrols.DropThresh.String);
user_profile.cfg.GRAPH_OUTPUT = h.Germination_Analysis.uicontrols.GraphOutput.Value;
user_profile.cfg.SAVE_FRAMOVIE = h.Germination_Analysis.uicontrols.FramovieOutput.Value;
user_profile.cfg.STORE_BF_TRACKS = h.Germination_Analysis.uicontrols.BFTracksOutput.Value;
user_profile.cfg.STORE_FL_TRACKS = h.Germination_Analysis.uicontrols.FlTracksOutput.Value;
user_profile.cfg.STORE_FITS = h.Germination_Analysis.uicontrols.FitFunctionOutput.Value;
user_profile.cfg.start_time = datenum(h.Germination_Analysis.uicontrols.TimeField.String);

% necessary for ReturnPlanes() to function properly
user_profile.FrameRange = 0;
% necessary for assign_fquant() to function properly
user_profile.Overwrite = 0;
% necessary for bf_quant() to function properly
user_profile.cfg.quant_type =   struct('flat_field', 0,...
                                       'camera_bg_subtract', 0,...
                                       'bgsubtract', 1,...
                                       'absoluteVals', 0);


% create a cell array of all dv-files per position analyzed
files_per_position = all_dv_to_position_lists(h.session_input_directory, user_profile.cfg.positions_to_use, user_profile.cfg.dvfiles_to_use);
% remove positions with no files
files_per_position = files_per_position(cellfun(@(x) ~isempty(x), files_per_position));
savefiles_parts = cellfun(@(x) split(x{1}, {'_' '.'}), files_per_position, 'UniformOutput', false);
savefiles = cellfun(@(x) join(x([1 3 4]), '_'), savefiles_parts, 'UniformOutput', false);
user_profile.cfg.position_numbers = cellfun(@(x) str2num(x{3}), savefiles_parts);

%% Produce a report and optionally track plots/movies for each position by calling 'movie_report' on every position

for position = 1:size(files_per_position, 2)
    savefile_report = fullfile(h.session_output_directory,[savefiles{position}{1} '_report.mat']);
    [report, framovie] = movie_report(files_per_position, position, h, user_profile);
    save(savefile_report,'report','-v7.3');
    disp(['report saved: ', savefile_report]); 
    if user_profile.cfg.SAVE_FRAMOVIE
            savefile_framovie = fullfile(h.session_framovie_output_directory, [savefiles{position}{1} '_framovie.mat']);
            save(savefile_framovie, 'framovie', '-v7.3');
            disp(['framovie saved: ', savefile_framovie]);
    end
end




    function [report, framovie] = movie_report(dvfiles, position, h, user_profile)
    % This function produces a report for one position 'position' from data stored in 'dvfiles', 
    % using parameters defined in the handle structure 'h' and the configuration structure 'user_profile'.
    
        store_bf_tracks = user_profile.cfg.STORE_BF_TRACKS;
        store_fl_tracks = user_profile.cfg.STORE_FL_TRACKS;
        store_fits = user_profile.cfg.STORE_FITS;
        dvfiles = dvfiles{position};
        dv_base = h.session_input_directory;
        fqs = {};
        framovie = {};
         
        
        %% Test whether there are multiple dv-files?
        
        if iscell(dvfiles)
            multiple_dvs = true;
            len = length(dvfiles);
        else
            multiple_dvs = false;
            len = 1;
        end
        
        %% Iterate over all dv-files to store data in fra-files and quantify their content 
        
        for i = 1:len
            if multiple_dvs
                dvpath = fullfile(dv_base, dvfiles{i});
            else
                dvpath = fullfile(dv_base, dvfiles);
            end
            % make framovie and quantify content of dv-files (stored in fqs)
            [fq, framovie] = dv2fq(framovie, dvpath, user_profile);
            fqs = [fqs fq];
        end
        
        split_names = cellfun(@(x) split(x(1).general.filename, {':', '_', 'lapse'}), fqs, 'UniformOutput', false);
        framenumbers = cellfun(@(x) str2num(x{end-3}) + str2num(x{end})-1, split_names);
        
        %% Extract time values from fqs        
        
        split_dates = cellfun(@(x) split(x.general.logdata.Date, ' '), fqs, 'UniformOutput', false);
        for i = 1:length(split_dates)
            if length(split_dates{i}{4}) < 2
                split_dates{i}{4} = ['0' split_dates{i}{4}];
            end
        end
        dates = cellfun(@(x) join({x{end}, x{2}, x{4}}, '-'), split_dates, 'UniformOutput', false);
        times = cellfun(@(x) x{5}, split_dates, 'UniformOutput', false);
        datimes = cellfun(@(x,y) join({x{1}, y}, ' '), dates, times, 'UniformOutput', false);
        ts = cellfun(@(x) datenum(x), datimes);
        ts = ts - user_profile.cfg.start_time;
        ts = minutes(days(ts));
        
       
        %% Construct correspondence matrices between two subsequent frames for subsequent mapping (= creating tracks)
        
        disp(['Extraction and quantification done, now starting mapping of spores between frames']);
        no_of_frames = length(fqs);
        Ms = cell(1,no_of_frames-1);
        for frame = 1:no_of_frames-1
            if frame == 1
                [optimizer, metric] = imregconfig('monomodal');
                img_1 = dip_array(framovie{frame}.spore);
                %img_1 = dip_array(framovie{frame}.bf_original);
                img_2 = dip_array(framovie{frame+1}.spore);
                tform = imregtform(img_2, img_1, 'translation', optimizer, metric);
                t_initial = tform.T(3,1:2);
                %img_2 = dip_array(framovie{frame+1}.bf_original);
            else
                t_initial = [0,0];
            end
            fq1 = fqs{frame};
            fq2 = fqs{frame+1};
            [Mbin,~] = directmatch(fq1, fq2,user_profile.cfg, t_initial);
            Ms{frame} = Mbin;
        end
        
        
        %% Relabel the spores in the images and store it in framovie{}.cluster
        disp(['Mapping done, relabeling spores']);
        framovie = relabel_spores(framovie, Ms);
        
               
        %% construct tracks of bf and fluorescence intensity
        
        num_objects = fqs{1}.spore.bf.numobjects;
        mapping = nan(num_objects, length(fqs));
        mapping(:, 1) = 1:num_objects;
        
        % iterate over all tracks
        for track = 1:num_objects
            index = track;
            % iterate over all frames
            for frame = 2:length(fqs)
                % use respective correspondence matrix
                corr_matrix = Ms{frame-1};
                % check whether object has a correspondence in next frame
                if find(corr_matrix(index,:))
                    mapping(track,frame) = find(corr_matrix(index,:));
                    index = mapping(track,frame);
                else
                    break
                end
            end
        end
        
        % remove incomplete mappings
        mapping = mapping(~any(isnan(mapping),2),:);
        num_complete_tracks = size(mapping, 1);
        
        %% Extract the tracks
        
        disp('Extracting brightfield tracks');
        % bf
        intensity_tracks = zeros(size(mapping));
        for track = 1:size(intensity_tracks,1)
            for frame = 1:size(intensity_tracks,2)
                intensity_tracks(track,frame) = fqs{frame}.spore.bf.mean(mapping(track,frame));
            end
        end
        
        % fl
        channels = fieldnames(fqs{1}.spore.fl);
        fluorescence_tracks = struct();
        for channel = channels'
            curr_fluorescence_tracks = zeros(size(mapping));
            for track = 1:size(mapping,1)
                for frame = 1:size(mapping,2)
                    curr_fluorescence_tracks(track,frame) = fqs{frame}.spore.fl.(channel{1}).int(mapping(track,frame));
                end
            end
            fluorescence_tracks.(channel{1}) = curr_fluorescence_tracks;
        end
        
        %% Extract the germination points from inclination points of sigmoidal fits to intensity tracks
        
        disp('Fitting and calculating germination points');
        [A_fits, germ_points, germ_bool] = germination_time_from_int_tracks(intensity_tracks, position, h, user_profile, framenumbers, ts);
        num_germ_spores = sum(germ_bool); 
        
        
        %% Extract fluorescence data
        
        if ~store_fl_tracks
            % take the mean intensity as fluorescence value
            %fluorescence_values = cellfun(@(x) mean(x,2), struct2cell(fluorescence_tracks), 'UniformOutput', false);

            % alternatively: take the first frames intensity as fluorescence value
            fluorescence_values = cellfun(@(x) x(:,1), struct2cell(fluorescence_tracks), 'UniformOutput', false)

            fluorescence_values = cell2struct(fluorescence_values, channels);
        else
            fluorescence_values = fluorescence_tracks;
        end

        
        %% Combine all the parts into a report structure
        
        disp('Create report');
        report = struct('times_after_revival_signal', [],...
                        'bf_track', [],...                        
                        'fit', [],...
                        'fit_func_parameters',[],...
                        'germination_time', [],...
                        'germination_classification', [],...
                        'fluorescence_values', struct,...
                        'number_total_spores', [],...
                        'number_tracked_spores', [],...
                        'number_germinated_spores', []);
        for channel = channels'
            report(1).fluorescence_values.(channel{1}) = [];
        end
        report.times_after_revival_signal = ts;
        report.number_total_spores = num_objects;
        report.number_tracked_spores = num_complete_tracks;
        report.number_germinated_spores = num_germ_spores;            
        if store_fits
            report.fit = @(A,x) A(1) + ((A(2) - A(1)) ./ (1 + exp( -A(3) .* (x - A(4)))).^(1 ./ 0.5));            
        end
        for i = 1:size(germ_points,1)
            if store_bf_tracks
                report(i).bf_track = intensity_tracks(i,:)/max(intensity_tracks(i,:));
            end
            if store_fits
                report(i).fit_func_parameters = A_fits{i};
            end
            report(i).germination_time = germ_points(i);
            report(i).germination_classification = germ_bool(i);
            for channel = channels'
                if ~store_fl_tracks
                    report(i).fluorescence_values.(channel{1}) = fluorescence_values.(channel{1})(i);
                else
                    report(i).fluorescence_values.(channel{1}) = fluorescence_values.(channel{1})(i,:);
                end
            end
        end
    end
end