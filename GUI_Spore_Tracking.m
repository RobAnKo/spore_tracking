function GUI_Spore_Tracking
%This GUI starts the germination time analysis of B.subtilis spores.
%
%How is this done? -> type 'doc tracking_spores'

%% Initialization of the tabbed windows

% get user screen size
SC = get(0, 'ScreenSize');
MaxMonitorX = SC(3);
MaxMonitorY = SC(4);

% set the figure window size values
MainFigScale = .8;          % Change this value to adjust the figure size
MaxWindowX = round(MaxMonitorX*MainFigScale);
MaxWindowY = round(MaxMonitorY*MainFigScale);
XBorder = (MaxMonitorX-MaxWindowX)/2;
YBorder = (MaxMonitorY-MaxWindowY)/2;
TabOffset = 0;              % This value offsets the tabs inside the figure.
ButtonHeight = 40;
PanelWidth = MaxWindowX-2*TabOffset+4;
PanelHeight = MaxWindowY-ButtonHeight-2*TabOffset;

% set the color variables
White = [1  1  1];   % White - Selected tab color
Grey = .9*White;     % Light Grey - Background color
BGColor = .8*White;

% initiate structure to hold all handles to all windows/tabs
h = struct( 'Germination_Analysis', struct());%,...
    %'Algorithm',struct(),...
    %'Config',struct());
fnames = fieldnames(h);

% import User Profile into h
h.user_profile=evalin('base','user_profile');
h.session_input_directory = h.user_profile.base_folder;
h.session_output_directory = h.user_profile.output_folder;
h.session_graph_output_directory =[h.session_output_directory filesep 'plots'];
h.session_framovie_output_directory =[h.session_output_directory filesep 'movies'];
% add fields for storage of parameters
h.parameters = struct;
h.parameters.BFChannel = '1';
h.parameters.FlChannel = '2,3';





% specific size/content settings for windows
NumberOfTabs = [2];               % Number of tabs to be generated for each of the three windows
TabLabels = {{'Input_Output'; 'Data_Analysis'}};
for entry = 1:length(NumberOfTabs)
    if size(TabLabels{entry},1) ~= NumberOfTabs(entry)
        errordlg('Number of tabs and tab labels must be the same','Setup Error');
        return
    end
end

% initiate all the figures
for fig_num = 1:length(fnames)
    
    ButtonWidth = round((PanelWidth-NumberOfTabs(fig_num))/NumberOfTabs(fig_num));
    h.(fnames{fig_num}).fig = figure('Units', 'pixels',...
        'Toolbar', 'none',...
        'Position',[ XBorder, YBorder, MaxWindowX, MaxWindowY ],...
        'NumberTitle', 'off',...
        'Name', fnames{fig_num},...
        'MenuBar', 'none',...
        'Resize', 'off',...
        'DockControls', 'off',...
        'Visible', 'off',... % Only the first window should be visible from the start
        'Color', White);
    
    % initiate all the tabs
    for tab_number = 1:NumberOfTabs(fig_num)
        
        % create a UIPanel
        h.(fnames{fig_num}).panels.([TabLabels{fig_num}{tab_number}])...
            = uipanel('Parent',h.(fnames{fig_num}).fig,...
            'Units', 'pixels', ...
            'Visible', 'off', ...
            'Backgroundcolor', Grey, ...
            'BorderWidth',1, ...
            'Position', [TabOffset TabOffset ...
            PanelWidth PanelHeight]);
        
        % create a selection pushbutton
        h.(fnames{fig_num}).tabbuttons.([TabLabels{fig_num}{tab_number}])...
            = uicontrol('Parent',h.(fnames{fig_num}).fig,...
            'Style', 'pushbutton',...
            'Units', 'pixels', ...
            'BackgroundColor', BGColor, ...
            'Position', [TabOffset+(tab_number-1)*ButtonWidth PanelHeight+TabOffset...
            ButtonWidth ButtonHeight], ...
            'String', TabLabels{fig_num}{tab_number},...
            'HorizontalAlignment', 'center',...
            'FontName', 'arial',...
            'FontWeight', 'bold',...
            'FontSize', 10);
    end
end

% set Callbacks for the tab selection buttons retroactively for completeness of h
for fig_num = 1:length(fnames)
    for tab_number = 1:NumberOfTabs(fig_num)
        set(h.(fnames{fig_num}).tabbuttons.([TabLabels{fig_num}{tab_number}]),...
            'Callback', {@TabSelectCallback, fig_num, tab_number})
    end
end


%% Fill the tabs with buttons/checkboxes etc

P = h.Germination_Analysis.panels.Input_Output;   %define parent panel,
%then add children uicontrols

h.Germination_Analysis.uicontrols.InputText =          uicontrol('Parent',P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Position', [.15 .9 .3 .07],...
    'String', 'Input',...
    'Fontsize', 16);

h.Germination_Analysis.uicontrols.OutputText =          uicontrol('Parent',P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Position', [.65 .9 .3 .07],...
    'String', 'Output',...
    'Fontsize', 16);

h.Germination_Analysis.uicontrols.InputButton =   uicontrol('Parent', P,...
    'Style', 'pushbutton',...
    'Units', 'normalized',...
    'BackgroundColor', BGColor,...
    'Position', [.05 .8 .1 .07],...
    'String', 'Input Folder');

h.Germination_Analysis.uicontrols.InputField =    uicontrol( 'Parent', P,...
    'Style', 'edit',...
    'Units', 'normalized',...
    'BackgroundColor', White,...
    'Position', [.18 .8 .3 .07],...
    'String', h.session_input_directory);

h.Germination_Analysis.uicontrols.OutputButton =  uicontrol('Parent', P,...
    'Style', 'pushbutton',...
    'Units', 'normalized',...
    'BackgroundColor', BGColor,...
    'Position', [.55 .8 .1 .07],...
    'String', 'Output Folder');

h.Germination_Analysis.uicontrols.OutputField =   uicontrol('Parent', P,...
    'Style', 'edit',...
    'Units', 'normalized',...
    'BackgroundColor', White,...
    'Position', [.68 .8 .3 .07],...
    'String', h.session_output_directory);

h.Germination_Analysis.uicontrols.PositionsText =          uicontrol('Parent',P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Position', [.05 .58 .1 .07],...
    'String', 'Which positions to use?',...
    'Fontsize', 10);

h.Germination_Analysis.uicontrols.PositionsField =   uicontrol('Parent', P,...
    'Style', 'edit',...
    'Units', 'normalized',...
    'BackgroundColor', White,...
    'Position', [.18 .6 .3 .07],...
    'String', 'all');

h.Germination_Analysis.uicontrols.DVFilesText =          uicontrol('Parent',P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Position', [.05 .48 .1 .07],...
    'String', 'Which DV-lapse files to use?',...
    'Fontsize', 10);

h.Germination_Analysis.uicontrols.DVFilesField =   uicontrol('Parent', P,...
    'Style', 'edit',...
    'Units', 'normalized',...
    'BackgroundColor', White,...
    'Position', [.18 .5 .3 .07],...
    'String', 'all');

h.Germination_Analysis.uicontrols.FramesText =          uicontrol('Parent',P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Position', [.05 .38 .1 .07],...
    'String', 'Which frames to use?',...
    'Fontsize', 10);

h.Germination_Analysis.uicontrols.FramesField =   uicontrol('Parent', P,...
    'Style', 'edit',...
    'Units', 'normalized',...
    'BackgroundColor', White,...
    'Position', [.18 .4 .3 .07],...
    'String', 'all');

h.Germination_Analysis.uicontrols.GraphOutput =          uicontrol('Parent',P,...
    'Style', 'checkbox',...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Position', [.75 .68 .15 .07],...
    'String', 'Plots of intensity tracks',...
    'Fontsize', 10,...
    'Callback', @GraphOutputCallback);

h.Germination_Analysis.uicontrols.GraphOutputField =   uicontrol('Parent', P,...
    'Style', 'edit',...
    'Units', 'normalized',...
    'BackgroundColor', White,...
    'Visible', 'off',...
    'Position', [.68 .6 .3 .07],...
    'String', [h.Germination_Analysis.uicontrols.OutputField.String filesep 'plots']); 

h.Germination_Analysis.uicontrols.GraphOutputButton = uicontrol('Parent', P,...
    'Style', 'pushbutton',...
    'Units', 'normalized',...
    'BackgroundColor', BGColor,...
    'Visible', 'off',...
    'Position', [.55 .6 .1 .07],...
    'String', 'Plots Output Folder');

h.Germination_Analysis.uicontrols.FramovieOutput =          uicontrol('Parent',P,...
    'Style', 'checkbox',...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Position', [.75 .48 .15 .07],...
    'String', 'Store framovies',...
    'Fontsize', 10,...
    'Callback', @FramovieOutputCallback);

h.Germination_Analysis.uicontrols.FramovieOutputField =   uicontrol('Parent', P,...
    'Style', 'edit',...
    'Units', 'normalized',...
    'BackgroundColor', White,...
    'Visible', 'off',...
    'Position', [.68 .4 .3 .07],...
    'String', [h.Germination_Analysis.uicontrols.OutputField.String filesep 'movies']); 

h.Germination_Analysis.uicontrols.FramovieOutputButton = uicontrol('Parent', P,...
    'Style', 'pushbutton',...
    'Units', 'normalized',...
    'BackgroundColor', BGColor,...
    'Visible', 'off',...
    'Position', [.55 0.4 .1 .07],...
    'String', 'Movies Output Folder');

h.Germination_Analysis.uicontrols.IncludeText =          uicontrol('Parent',P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Position', [.63 .21 .15 .07],...
    'String', 'Include in report:',...
    'Fontsize', 10);

h.Germination_Analysis.uicontrols.FlTracksOutput =          uicontrol('Parent',P,...
    'Style', 'checkbox',...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Position', [.75 .23 .20 .07],...
    'String', 'Fl intensity tracks (instead of one value)',...
    'Fontsize', 10);

h.Germination_Analysis.uicontrols.BFTracksOutput =          uicontrol('Parent',P,...
    'Style', 'checkbox',...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Position', [.75 .18 .15 .07],...
    'String', 'BF intensity tracks',...
    'Fontsize', 10);

h.Germination_Analysis.uicontrols.FitFunctionOutput =          uicontrol('Parent',P,...
    'Style', 'checkbox',...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Position', [.75 .13 .15 .07],...
    'String', 'BF fitting function',...
    'Fontsize', 10);

h.Germination_Analysis.uicontrols.InfoText =          uicontrol('Parent',P,...
    'Style', 'text',...
    'HorizontalAlignment', 'left',...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Position', [.18 .0 .3 .35],...
    'String', ['For these fields, use the standard MATLAB indexing to choose the position/dv-file/frame IDs:' newline newline 'only 7 --> 7,' newline newline '2, 3, 4 and 5 --> [2:5],' newline newline '9, 13, 24 --> [9, 13, 24],' newline '2, 7, 8, 9, 12 --> [2, 7:9, 12],' newline...
    'all even: --> [2:2:end],' newline newline 'etc...' newline newline 'If you want all positions/dvs/frames, use "all"'],...
    'Fontsize', 10);






P = h.Germination_Analysis.panels.Data_Analysis; %define parent panel,
%then add children uicontrols


h.Germination_Analysis.uicontrols.DataText =          uicontrol('Parent',P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Position', [.05 .9 .1 .07],...
    'String', 'Data',...
    'HorizontalAlignment', 'left',...
    'Fontsize', 16);

h.Germination_Analysis.uicontrols.BFText =   uicontrol('Parent', P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Position', [.05 .78 .1 .07],...
    'String', 'BF Channel',...
    'FontSize', 10);

h.Germination_Analysis.uicontrols.BFField =    uicontrol( 'Parent', P,...
    'Style', 'edit',...
    'Units', 'normalized',...
    'BackgroundColor', White,...
    'Position', [.18 .8 .3 .07],...
    'String', h.parameters.BFChannel);

h.Germination_Analysis.uicontrols.FlText =  uicontrol('Parent', P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Position', [.05 .68 .1 .07],...
    'FontSize', 10,...
    'String', 'Fluorescence Channel(s)',...
    'Callback', @BFchannelFieldCallback);

h.Germination_Analysis.uicontrols.FlField =   uicontrol('Parent', P,...
    'Style', 'edit',...
    'Units', 'normalized',...
    'BackgroundColor', White,...
    'Position', [.18 .7 .3 .07],...
    'String', h.parameters.FlChannel,...
    'Callback', @FlchannelFieldCallback);

h.Germination_Analysis.uicontrols.ParameterText =          uicontrol('Parent',P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Position', [.05 .5 .1 .07],...
    'String', 'Parameters',...
    'HorizontalAlignment', 'left',...
    'Fontsize', 16);

h.Germination_Analysis.uicontrols.TimeText =   uicontrol('Parent', P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Position', [.53 .76 .22 .1],...
    'String', ['At which time was germination induced?',newline, ' Format like "yyyy-mmm-dd HH:MM:SS"' newline ' Example: "2017-Sep-29 14:33:00"'],...
    'FontSize', 10);

h.Germination_Analysis.uicontrols.TimeField =    uicontrol( 'Parent', P,...
    'Style', 'edit',...
    'Units', 'normalized',...
    'BackgroundColor', White,...
    'Position', [.75 .8 .2 .07],...
    'String', '2017-Oct-05 13:51:00');

h.Germination_Analysis.uicontrols.TrackingButton =     uicontrol('Parent',P,...
    'Style', 'pushbutton',...
    'Units', 'normalized',...
    'BackgroundColor', 'green',...
    'Position', [.8 .1 .12 .07],...
    'String', 'Start analysis ->',...
    'callback',@TrackingButtonCallback);


%% Create subpanels for parameter sets

h.Germination_Analysis.uicontrols.subpanels.Preprocess =               uipanel('Parent', P,...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Visible','On',...
    'Position', [.05 .15 .2 .3]);

h.Germination_Analysis.uicontrols.subpanels.Segmentation =               uipanel('Parent', P,...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Visible','On',...
    'Position', [.3 .15 .2 .3]);

h.Germination_Analysis.uicontrols.subpanels.Tracking =               uipanel('Parent', P,...
    'Units', 'normalized',...
    'BackgroundColor', Grey,...
    'Visible','On',...
    'Position', [.55 .15 .2 .3]);

P = h.Germination_Analysis.uicontrols.subpanels.Preprocess;

h.Germination_Analysis.uicontrols.PreprocessText =         uicontrol('Parent', P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'string', 'Preprocessing:',...
    'FontWeight', 'Bold',...
    'HorizontalAlignment', 'left',...
    'Backgroundcolor',Grey,...
    'Visible', 'On',...
    'Position', [.05 .8 1 .1]);

h.Germination_Analysis.uicontrols.DarkCellsText =         uicontrol('Parent', P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'string', 'DARK_CELLS:',...
    'Backgroundcolor',Grey,...
    'Visible', 'On',...
    'Position', [0 .7 .5 .05]);

h.Germination_Analysis.uicontrols.DarkCells =         uicontrol('Parent', P,...
    'Style', 'edit',...
    'Units', 'normalized',...
    'string', 0,...
    'Visible', 'On',...
    'Position', [.5 .65 .3 .166]);

h.Germination_Analysis.uicontrols.DomeHeightText =         uicontrol('Parent', P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'string', 'DOME_HEIGHT:',...
    'Backgroundcolor',Grey,...
    'Visible', 'On',...
    'Position', [0 .54 .5 .05]);

h.Germination_Analysis.uicontrols.DomeHeight =         uicontrol('Parent', P,...
    'Style', 'edit',...
    'Units', 'normalized',...
    'string', 0.05,...
    'Visible', 'On',...
    'Position', [.5 .48 .3 .166]);


h.Germination_Analysis.uicontrols.ContrastText =         uicontrol('Parent', P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'string', 'CONTRAST:',...
    'Backgroundcolor',Grey,...
    'Visible', 'On',...
    'Position', [0 .37 .5 .05]);

h.Germination_Analysis.uicontrols.Contrast =         uicontrol('Parent', P,...
    'Style', 'edit',...
    'Units', 'normalized',...
    'string', 8,...
    'Visible', 'On',...
    'Position', [.5 .31 .3 .166]);

P = h.Germination_Analysis.uicontrols.subpanels.Segmentation;

h.Germination_Analysis.uicontrols.SegmentationText =       uicontrol('Parent', P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'string', 'Segmentation:',...
    'FontWeight', 'Bold',...
    'HorizontalAlignment', 'left',...
    'Backgroundcolor',Grey,...
    'Visible', 'On',...
    'Position', [.05 .8 1 .1]);

h.Germination_Analysis.uicontrols.SporeIntText =       uicontrol('Parent', P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'string', 'SPORE_INT:',...
    'Backgroundcolor',Grey,...
    'Visible', 'On',...
    'Position', [0 .7 .5 .05]);

h.Germination_Analysis.uicontrols.SporeInt =           uicontrol('Parent', P,...
    'Style', 'edit',...
    'Units', 'normalized',...
    'string', 2.5,...
    'Visible', 'On',...
    'Position', [.5 .65 .3 .166]);

h.Germination_Analysis.uicontrols.SporeMinSizeText =       uicontrol('Parent', P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'string', 'SPORE_MIN_SIZE:',...
    'Backgroundcolor',Grey,...
    'Visible', 'On',...
    'Position', [0 .54 .5 .05]);

h.Germination_Analysis.uicontrols.SporeMinSize =           uicontrol('Parent', P,...
    'Style', 'edit',...
    'Units', 'normalized',...
    'string', 10,...
    'Visible', 'On',...
    'Position', [.5 .48 .3 .166]);
h.Germination_Analysis.uicontrols.SporeMaxSizeText =       uicontrol('Parent', P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'string', 'SPORE_MAX_SIZE:',...
    'Backgroundcolor',Grey,...
    'Visible', 'On',...
    'Position', [0 .37 .5 .05]);

h.Germination_Analysis.uicontrols.SporeMaxSize =           uicontrol('Parent', P,...
    'Style', 'edit',...
    'Units', 'normalized',...
    'string', 100,...
    'Visible', 'On',...
    'Position', [.5 .31 .3 .166]);

P = h.Germination_Analysis.uicontrols.subpanels.Tracking;

h.Germination_Analysis.uicontrols.TrackingText =       uicontrol('Parent', P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'string', 'Tracking:',...
    'FontWeight', 'Bold',...
    'HorizontalAlignment', 'left',...
    'Backgroundcolor',Grey,...
    'Visible', 'On',...
    'Position', [.05 .8 1 .1]);

h.Germination_Analysis.uicontrols.ExpVarText =       uicontrol('Parent', P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'string', 'EXP_VARIANCE:',...
    'Backgroundcolor',Grey,...
    'Visible', 'On',...
    'Position', [0 .7 .5 .05]);

h.Germination_Analysis.uicontrols.ExpVarInt =           uicontrol('Parent', P,...
    'Style', 'edit',...
    'Units', 'normalized',...
    'string', 2,...
    'Visible', 'On',...
    'Position', [.5 .65 .3 .166]);

h.Germination_Analysis.uicontrols.DropThreshText =       uicontrol('Parent', P,...
    'Style', 'text',...
    'Units', 'normalized',...
    'string', 'DROP_THRESH:',...
    'Backgroundcolor',Grey,...
    'Visible', 'On',...
    'Position', [0 .54 .5 .05]);

h.Germination_Analysis.uicontrols.DropThresh =           uicontrol('Parent', P,...
    'Style', 'edit',...
    'Units', 'normalized',...
    'string', 0.2,...
    'Visible', 'On',...
    'Position', [.5 .48 .3 .166]);


% set callbacks only now to assure that h is already completely filled
set(h.Germination_Analysis.uicontrols.InputButton, 'Callback', @InputButtonCallback);
set(h.Germination_Analysis.uicontrols.InputField, 'Callback', @InputFieldCallback);
set(h.Germination_Analysis.uicontrols.OutputButton, 'Callback', @OutputButtonCallback);
set(h.Germination_Analysis.uicontrols.OutputField, 'Callback', @OutputFieldCallback);
set(h.Germination_Analysis.uicontrols.GraphOutput, 'Callback', @GraphOutputCallback);
set(h.Germination_Analysis.uicontrols.GraphOutputButton, 'Callback', @GraphOutputButtonCallback);
set(h.Germination_Analysis.uicontrols.GraphOutputField, 'Callback', @GraphOutputFieldCallback);
set(h.Germination_Analysis.uicontrols.FramovieOutput, 'Callback', @FramovieOutputCallback);
set(h.Germination_Analysis.uicontrols.FramovieOutputButton, 'Callback', @FramovieOutputButtonCallback);
set(h.Germination_Analysis.uicontrols.FramovieOutputField, 'Callback', @FramovieOutputFieldCallback);
% make first figure and first panel(tab) visible
set(h.Germination_Analysis.fig, 'Visible', 'On')
TabSelectCallback(0,0,1,1)



%% Callback functions

% general callback for selection of tabs
    function TabSelectCallback(~,~,active_figure, selected_tab) 
        % set active figure visible
        set(h.(fnames{active_figure}).fig, 'Visible', 'On');
        % turn all tabs off
        pnames = fieldnames(h.(fnames{active_figure}).panels);
        for p = 1:numel(pnames)
            set(h.(fnames{active_figure}).panels.(pnames{p}), 'Visible', 'off');
            set(h.(fnames{active_figure}).tabbuttons.(pnames{p}), 'BackgroundColor', BGColor);
        end
        
        %   Enable the selected tab
        set(h.(fnames{active_figure}).panels.(pnames{selected_tab}), 'Visible', 'on');
        set(h.(fnames{active_figure}).tabbuttons.(pnames{selected_tab}), 'BackgroundColor', White);
    end

% window: "Image", panel: "Input_Output"
    function InputButtonCallback(~,~)
        directory=uigetdir(h.session_input_directory);
        if ~directory
            return
        else
            h.session_input_directory = directory;
            set(h.Germination_Analysis.uicontrols.InputField, 'string', directory)
        end
    end

    function InputFieldCallback(~,~)
        directory = get(gcbo,'String');
        h.session_input_directory = directory;
    end

    function OutputButtonCallback(~,~)
        if ~exist(h.session_output_directory,'dir')
            mkdir(h.session_output_directory);
        end
        directory=uigetdir(h.session_output_directory);
        if ~directory
            return
        else
            h.session_output_directory = directory;
            set(h.Germination_Analysis.uicontrols.OutputField, 'string', directory)
        end
        h.Germination_Analysis.uicontrols.GraphOutputField.String = [directory filesep 'plots'];
        h.session_graph_output_directory = [directory filesep 'plots'];
        h.Germination_Analysis.uicontrols.FramovieOutputField.String = [directory filesep 'movies'];
        h.session_framovie_output_directory = [directory filesep 'movies'];
    end

    

    function OutputFieldCallback(~,~)
        directory = get(gcbo,'String');
        if ~exist(directory,'dir')
            mkdir(directory);
        end
        h.session_output_directory = directory;
        h.session_framovie_output_directory =[h.session_output_directory filesep 'movies'];
        h.session_graph_output_directory =[h.session_output_directory filesep 'plots'];

    end


    function GraphOutputCallback(~,~)
        active = get(gcbo, 'Value');
        if active
            set(h.Germination_Analysis.uicontrols.GraphOutputField, 'visible', 'on');
            set(h.Germination_Analysis.uicontrols.GraphOutputField, 'string', [h.Germination_Analysis.uicontrols.OutputField.String filesep 'plots']);
            set(h.Germination_Analysis.uicontrols.GraphOutputButton, 'visible', 'on');
            h.session_graph_output_directory = h.Germination_Analysis.uicontrols.GraphOutputField.String;
        else
            set(h.Germination_Analysis.uicontrols.GraphOutputField, 'visible', 'off');
            set(h.Germination_Analysis.uicontrols.GraphOutputButton, 'visible', 'off');
            h.session_graph_output_directory = h.Germination_Analysis.uicontrols.GraphOutputField.String;
        end
        directory = h.session_graph_output_directory;
        if ~exist(directory, 'dir')
            mkdir(directory);
        end
        if exist(directory, 'dir') && ~active
            files = dir(fullfile(directory, '*.png'));
            for k = 1:length(files)
                basefilename = files(k).name;
                delete(fullfile(directory, basefilename));
            end
            rmdir(directory, 's');
        end
    end

    function GraphOutputButtonCallback(~,~)
        directory = get(h.Germination_Analysis.uicontrols.GraphOutputField, 'string');
        if ~exist(directory,'dir')
            mkdir(directory);
        end
        directory=uigetdir(directory);
        if ~directory 
            return
        end
        set(h.Germination_Analysis.uicontrols.GraphOutputField, 'string', directory);
        h.session_graph_output_directory = directory;
    end

    function GraphOutputFieldCallback(~,~)
        directory = get(gcbo,'String');
        if ~exist(directory,'dir')
            mkdir(directory);
        end
        h.session_graph_output_directory = directory;
    end

function FramovieOutputCallback(~,~)
        active = get(gcbo, 'Value');
        if active
            set(h.Germination_Analysis.uicontrols.FramovieOutputField, 'visible', 'on');
            set(h.Germination_Analysis.uicontrols.FramovieOutputField, 'string', [h.Germination_Analysis.uicontrols.OutputField.String filesep 'movies']);
            set(h.Germination_Analysis.uicontrols.FramovieOutputButton, 'visible', 'on');
            h.session_framovie_output_directory = h.Germination_Analysis.uicontrols.FramovieOutputField.String;
        else
            set(h.Germination_Analysis.uicontrols.FramovieOutputField, 'visible', 'off');
            set(h.Germination_Analysis.uicontrols.FramovieOutputButton, 'visible', 'off');
            h.session_framovie_output_directory = h.Germination_Analysis.uicontrols.FramovieOutputField.String;
        end
        directory = h.session_framovie_output_directory;
        if ~exist(directory, 'dir')
            mkdir(directory);
        end
        if exist(directory, 'dir') && ~active
            files = dir(fullfile(directory, '*.mat'));
            for k = 1:length(files)
                basefilename = files(k).name;
                delete(fullfile(directory, basefilename));
            end
            rmdir(directory,'s');
        end
    end

    function FramovieOutputButtonCallback(~,~)
        directory = get(h.Germination_Analysis.uicontrols.FramovieOutputField, 'string');
        if ~exist(directory,'dir')
            mkdir(directory);
        end
        directory=uigetdir(directory);
        if ~directory 
            return
        end
        set(h.Germination_Analysis.uicontrols.FramovieOutputField, 'string', directory);
        h.session_framovie_output_directory = directory;
    end

    function FramovieOutputFieldCallback(~,~)
        directory = get(gcbo,'String');
        if ~exist(directory,'dir')
            mkdir(directory);
        end
        h.session_framovie_output_directory = directory;
    end


% window 'Germination_Analysis', panel 'Data_Analysis'
    
    
    function BFchannelFieldCallback(~,~)
        str = get(gcbo, 'string');
        h.parameters.BFChannel = str;
    end

    function FlchannelFieldCallback(~,~)
        str = get(gcbo, 'string');
        h.parameters.FlChannel = str;
    end


  

    function TrackingButtonCallback(~,~)
        overall_time = tic;
        tracking_spores(h);
        overall_time = toc(overall_time);
        disp(['processing took ' num2str(overall_time/3600) ' hours']);
    end

   
end
