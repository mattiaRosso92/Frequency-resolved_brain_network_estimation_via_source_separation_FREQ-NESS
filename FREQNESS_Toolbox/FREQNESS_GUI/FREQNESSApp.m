classdef FREQNESSApp < handle
    %FREQNESSAPP Programmatic graphical interface for the FREQ-NESS Toolbox.

    properties (SetAccess = private)
        UIFigure
    end

    properties (Access = private)
        Dataset
        NetworkSet
        Config
        IsRunning = false
        PageGroup
        NetworkTab
        SecondaryTab
        MainGrid

        DatasetDropLabel
        DatasetPathField
        DatasetBrowseButton
        DatasetSummaryLabel
        DatasetTable

        MNI
        MNIDropLabel
        MNIPathField
        MNIBrowseButton
        MNISummaryLabel
        MNIAxes

        FrequencyRangeSlider
        FrequencyMinField
        FrequencyMaxField
        FrequencyStepField
        FrequencySummaryLabel
        FWHMAxes
        SamplingRateField
        DurationField
        FilterWidthField
        FilterTypeDropDown
        RegularisationField
        ComponentsField
        BadSegmentsField
        RescaleCheckBox
        RecomputeCheckBox
        NetworkOutputField
        RunButton
        SecondaryPageButton
        RunStatusLabel
        ProgressTextArea

        NetworkDropLabel
        NetworkPathField
        NetworkBrowseButton
        AnalysisOutputField
        SecondarySummaryLabel
        SecondaryParticipantButton
        SecondarySelectedParticipants = {}
        SecondaryFrequencyLabel
        SecondaryComponentsLabel
        SecondaryMNILabel
        SecondarySourceLabel
        SecondaryTree
        SecondaryConfigPanel
        SecondaryConfigGrid
        SecondaryOptionalPanel
        SecondaryAdvancedButton
        SecondaryFields = struct()
        SecondarySchema
        SecondarySelectedModuleId = ''
        SecondaryReadinessTextArea
        SecondaryModuleOutputField
        SecondaryValidationLabel
        SecondaryValidateButton
        SecondaryBackgroundCheckBox
        SecondaryRunButton
        SecondaryProgressTextArea
        SecondaryProcess
        SecondaryProcessTimer
        SecondaryProcessControl
        SecondaryLastProgressSequence = 0
        SecondaryCancellationFile = ''
        SecondaryCancelRequested = false
        SecondaryCancelTic
        SecondaryRunInBackground = false
        SecondaryRunConfiguration

        FooterStatusLabel
    end

    methods
        function app = FREQNESSApp(varargin)
            app.Config = freqnessgui.defaultConfig();
            app.createInterface();
            setappdata(app.UIFigure,'FREQNESSApp',app);
            app.loadDefaultMNI();
            app.registerFileDropTargets();
            drawnow
            try
                scroll(app.MainGrid,'top');
            catch
                % Older releases retain their default top scroll position.
            end

            if nargin >= 1 && ~isempty(varargin{1})
                inputPath = char(varargin{1});
                [~,inputName] = fileparts(inputPath);
                networkPrefix = 'FREQ_Networks';
                if isfolder(inputPath) && numel(inputName) >= numel(networkPrefix) && ...
                        strcmpi(inputName(1:numel(networkPrefix)),networkPrefix)
                    app.importNetworkFolder(inputPath);
                    app.PageGroup.SelectedTab = app.SecondaryTab;
                else
                    app.importDataset(inputPath);
                end
            end
        end

        function delete(app)
            app.shutdownSecondaryExecution();
            if ~isempty(app.UIFigure) && isvalid(app.UIFigure)
                app.UIFigure.CloseRequestFcn = '';
                delete(app.UIFigure);
            end
        end
    end

    methods (Access = private)
        function createInterface(app)
            colors.background = [0.955 0.965 0.975];
            colors.navy = [0.075 0.145 0.235];
            colors.blue = [0.055 0.415 0.690];
            colors.softBlue = [0.925 0.955 0.985];
            colors.green = [0.125 0.545 0.365];
            colors.border = [0.78 0.82 0.87];
            colors.muted = [0.36 0.41 0.47];

            app.UIFigure = uifigure( ...
                'Name','FREQ-NESS Toolbox', ...
                'Tag','FREQNESS_GUI_MainFigure', ...
                'Color',colors.background, ...
                'Position',[60 20 1360 980], ...
                'CloseRequestFcn',@(~,~)delete(app));

            rootGrid = uigridlayout(app.UIFigure,[1 1]);
            rootGrid.Padding = [0 0 0 0];
            app.PageGroup = uitabgroup(rootGrid, ...
                'SelectionChangedFcn',@(~,~)app.updateSecondaryReadiness());
            app.NetworkTab = uitab(app.PageGroup, ...
                'Title','1  Network estimation');
            app.SecondaryTab = uitab(app.PageGroup, ...
                'Title','2  Secondary analyses');

            app.MainGrid = uigridlayout(app.NetworkTab,[4 1]);
            app.MainGrid.RowHeight = {72,250,455,26};
            app.MainGrid.Padding = [18 14 18 10];
            app.MainGrid.RowSpacing = 10;
            app.MainGrid.Scrollable = 'on';
            mainGrid = app.MainGrid;

            headerGrid = uigridlayout(mainGrid,[2 1]);
            headerGrid.Layout.Row = 1;
            headerGrid.RowHeight = {38,24};
            headerGrid.Padding = [12 3 12 2];
            headerGrid.RowSpacing = 0;
            headerGrid.BackgroundColor = colors.navy;

            titleLabel = uilabel(headerGrid, ...
                'Text','FREQ-NESS Toolbox', ...
                'FontSize',25, ...
                'FontWeight','bold', ...
                'FontColor',[1 1 1]);
            titleLabel.Layout.Row = 1;
            subtitleLabel = uilabel(headerGrid, ...
                'Text','Frequency-resolved brain network estimation via source separation', ...
                'FontSize',12, ...
                'FontColor',[0.80 0.87 0.94]);
            subtitleLabel.Layout.Row = 2;

            app.createDataLayer(mainGrid,colors);
            app.createCoreLayer(mainGrid,colors);

            app.FooterStatusLabel = uilabel(mainGrid, ...
                'Text','Ready', ...
                'FontSize',11, ...
                'FontColor',colors.muted);
            app.FooterStatusLabel.Layout.Row = 4;

            app.createSecondaryPage(app.SecondaryTab,colors);
        end

        function createDataLayer(app,parent,colors)
            panel = uipanel(parent, ...
                'Title','1  Data and anatomical import', ...
                'FontSize',14, ...
                'FontWeight','bold', ...
                'ForegroundColor',colors.navy, ...
                'BackgroundColor',[1 1 1]);
            panel.Layout.Row = 2;

            grid = uigridlayout(panel,[1 2]);
            grid.ColumnWidth = {'1x','1x'};
            grid.Padding = [12 8 12 10];
            grid.ColumnSpacing = 12;

            datasetPanel = uipanel(grid, ...
                'Title','Participant data', ...
                'FontWeight','bold', ...
                'BackgroundColor',[1 1 1]);
            datasetPanel.Layout.Column = 1;
            datasetGrid = uigridlayout(datasetPanel,[3 3]);
            datasetGrid.RowHeight = {34,20,'1x'};
            datasetGrid.ColumnWidth = {170,'1x',86};
            datasetGrid.Padding = [8 6 8 8];
            datasetGrid.RowSpacing = 6;
            datasetGrid.ColumnSpacing = 8;

            app.DatasetDropLabel = uilabel(datasetGrid, ...
                'Text',sprintf('Drop a Dataset* folder\nor its MAT-files here'), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','center', ...
                'FontWeight','bold', ...
                'FontColor',colors.blue, ...
                'BackgroundColor',colors.softBlue);
            app.DatasetDropLabel.Layout.Row = [1 3];
            app.DatasetDropLabel.Layout.Column = 1;

            app.DatasetPathField = uieditfield(datasetGrid,'text', ...
                'Editable','off', ...
                'Placeholder','No Dataset* folder selected');
            app.DatasetPathField.Layout.Row = 1;
            app.DatasetPathField.Layout.Column = 2;

            app.DatasetBrowseButton = uibutton(datasetGrid,'push', ...
                'Text','Browse...', ...
                'ButtonPushedFcn',@(~,~)app.browseDataset());
            app.DatasetBrowseButton.Layout.Row = 1;
            app.DatasetBrowseButton.Layout.Column = 3;

            app.DatasetSummaryLabel = uilabel(datasetGrid, ...
                'Text','Waiting for a Dataset* folder', ...
                'FontColor',colors.muted);
            app.DatasetSummaryLabel.Layout.Row = 2;
            app.DatasetSummaryLabel.Layout.Column = [2 3];

            app.DatasetTable = uitable(datasetGrid, ...
                'Data',cell(0,4), ...
                'ColumnName',{'Participant','Variable','Matrix size','Validation'}, ...
                'ColumnWidth',{125,105,85,'auto'}, ...
                'RowName',{});
            app.DatasetTable.Layout.Row = 3;
            app.DatasetTable.Layout.Column = [2 3];

            mniPanel = uipanel(grid, ...
                'Title','MNI coordinates', ...
                'FontWeight','bold', ...
                'BackgroundColor',[1 1 1]);
            mniPanel.Layout.Column = 2;
            mniGrid = uigridlayout(mniPanel,[3 3]);
            mniGrid.RowHeight = {34,20,'1x'};
            mniGrid.ColumnWidth = {170,'1x',86};
            mniGrid.Padding = [8 6 8 8];
            mniGrid.RowSpacing = 6;
            mniGrid.ColumnSpacing = 8;

            app.MNIDropLabel = uilabel(mniGrid, ...
                'Text',sprintf('Drop an MNI coordinate\nMAT-file here'), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','center', ...
                'FontWeight','bold', ...
                'FontColor',colors.blue, ...
                'BackgroundColor',colors.softBlue);
            app.MNIDropLabel.Layout.Row = [1 3];
            app.MNIDropLabel.Layout.Column = 1;

            app.MNIPathField = uieditfield(mniGrid,'text', ...
                'Editable','off', ...
                'Placeholder','No MNI coordinate file selected');
            app.MNIPathField.Layout.Row = 1;
            app.MNIPathField.Layout.Column = 2;

            app.MNIBrowseButton = uibutton(mniGrid,'push', ...
                'Text','Browse...', ...
                'ButtonPushedFcn',@(~,~)app.browseMNIFile());
            app.MNIBrowseButton.Layout.Row = 1;
            app.MNIBrowseButton.Layout.Column = 3;

            app.MNISummaryLabel = uilabel(mniGrid, ...
                'Text','No MNI coordinates loaded', ...
                'FontColor',colors.muted);
            app.MNISummaryLabel.Layout.Row = 2;
            app.MNISummaryLabel.Layout.Column = [2 3];

            app.MNIAxes = uiaxes(mniGrid);
            app.MNIAxes.Layout.Row = 3;
            app.MNIAxes.Layout.Column = [2 3];
            app.MNIAxes.Interactions = rotateInteraction;
            app.showEmptyMNIPreview();
        end

        function createCoreLayer(app,parent,colors)
            panel = uipanel(parent, ...
                'Title','2  Core analysis — FREQNESS_NetworkEstimation', ...
                'FontSize',14, ...
                'FontWeight','bold', ...
                'ForegroundColor',colors.navy, ...
                'BackgroundColor',[1 1 1]);
            panel.Layout.Row = 3;

            grid = uigridlayout(panel,[1 3]);
            grid.ColumnWidth = {520,330,'1x'};
            grid.Padding = [12 10 12 12];
            grid.ColumnSpacing = 12;

            mandatoryPanel = uipanel(grid, ...
                'Title','Mandatory inputs', ...
                'FontWeight','bold', ...
                'BackgroundColor',[1 1 1]);
            mandatoryPanel.Layout.Column = 1;
            mandatoryGrid = uigridlayout(mandatoryPanel,[8 1]);
            mandatoryGrid.RowHeight = {20,34,54,0,28,'1x',14,42};
            mandatoryGrid.Padding = [10 8 10 8];
            mandatoryGrid.RowSpacing = 4;

            label = uilabel(mandatoryGrid,'Text','Frequency range and resolution (Hz)', ...
                'FontWeight','bold');
            label.Layout.Row = 1;

            frequencyInputGrid = uigridlayout(mandatoryGrid,[1 8]);
            frequencyInputGrid.Layout.Row = 2;
            frequencyInputGrid.ColumnWidth = {24,'1x',28,'1x',34,'1x',38,'1x'};
            frequencyInputGrid.Padding = [0 0 0 0];
            frequencyInputGrid.ColumnSpacing = 5;
            label = uilabel(frequencyInputGrid,'Text','Min');
            label.Layout.Column = 1;
            app.FrequencyMinField = uieditfield(frequencyInputGrid,'numeric', ...
                'Value',app.Config.network.frequencies(1), ...
                'Limits',[eps Inf], ...
                'ValueDisplayFormat','%.4g Hz', ...
                'ValueChangedFcn',@(~,~)app.commitFrequencyFields());
            app.FrequencyMinField.Layout.Column = 2;
            label = uilabel(frequencyInputGrid,'Text','Max');
            label.Layout.Column = 3;
            app.FrequencyMaxField = uieditfield(frequencyInputGrid,'numeric', ...
                'Value',app.Config.network.frequencies(end), ...
                'Limits',[eps Inf], ...
                'ValueDisplayFormat','%.4g Hz', ...
                'ValueChangedFcn',@(~,~)app.commitFrequencyFields());
            app.FrequencyMaxField.Layout.Column = 4;
            label = uilabel(frequencyInputGrid,'Text','Step');
            label.Layout.Column = 5;
            defaultStep = median(diff(app.Config.network.frequencies));
            app.FrequencyStepField = uieditfield(frequencyInputGrid,'numeric', ...
                'Value',defaultStep, ...
                'Limits',[eps Inf], ...
                'ValueDisplayFormat','%.4g Hz', ...
                'ValueChangedFcn',@(~,~)app.commitFrequencyFields());
            app.FrequencyStepField.Layout.Column = 6;

            label = uilabel(frequencyInputGrid,'Text','Srate');
            label.Layout.Column = 7;
            app.SamplingRateField = uieditfield(frequencyInputGrid,'numeric', ...
                'Value',app.Config.network.samplingRate, ...
                'Limits',[eps Inf], ...
                'ValueDisplayFormat','%.4g Hz', ...
                'ValueChangedFcn',@(~,~)app.samplingRateChanged());
            app.SamplingRateField.Layout.Column = 8;

            sliderCeiling = app.frequencySliderCeiling( ...
                app.Config.network.frequencies(end),defaultStep);
            app.FrequencyRangeSlider = uislider(mandatoryGrid,'range', ...
                'Limits',[0 sliderCeiling], ...
                'Value',[app.Config.network.frequencies(1), ...
                    app.Config.network.frequencies(end)], ...
                'Step',defaultStep, ...
                'MinorTicks',[], ...
                'ValueChangingFcn',@(~,event) ...
                    app.previewFrequencyRange(event.Value), ...
                'ValueChangedFcn',@(~,event) ...
                    app.commitFrequencyRange(event.Value));
            app.FrequencyRangeSlider.Layout.Row = 3;

            app.FrequencySummaryLabel = uilabel(mandatoryGrid, ...
                'Text','', ...
                'FontColor',colors.blue, ...
                'HorizontalAlignment','center', ...
                'Visible','off');
            app.FrequencySummaryLabel.Layout.Row = 4;

            app.FWHMAxes = uiaxes(mandatoryGrid);
            app.FWHMAxes.Layout.Row = 6;
            app.FWHMAxes.FontSize = 9;
            app.FWHMAxes.Box = 'on';
            app.FWHMAxes.XGrid = 'on';
            app.FWHMAxes.YGrid = 'on';
            xlabel(app.FWHMAxes,'Frequency (Hz)');
            ylabel(app.FWHMAxes,'Normalized gain');

            note = uilabel(mandatoryGrid, ...
                'Text',sprintf(['Participant matrices are read as voxels x time.\n' ...
                'Each .mat file becomes one persisted FREQ result.']), ...
                'FontColor',colors.muted, ...
                'VerticalAlignment','top');
            note.Layout.Row = 8;

            optionalPanel = uipanel(grid, ...
                'Title','Optional inputs', ...
                'FontWeight','bold', ...
                'BackgroundColor',[1 1 1]);
            optionalPanel.Layout.Column = 2;
            optionalGrid = uigridlayout(optionalPanel,[7 2]);
            optionalGrid.RowHeight = {30,30,30,30,30,30,30};
            optionalGrid.ColumnWidth = {155,'1x'};
            optionalGrid.Padding = [10 8 10 8];
            optionalGrid.RowSpacing = 5;

            app.DurationField = app.addTextSetting(optionalGrid,1, ...
                'Duration (seconds)','', 'Empty uses the full recording');
            app.FilterWidthField = app.addTextSetting(optionalGrid,2, ...
                'Filter width (FWHM)','', 'Empty uses automatic widths');
            app.FilterWidthField.ValueChangedFcn = @(~,~)app.updateFWHMPreview();

            label = uilabel(optionalGrid,'Text','Filter progression');
            label.Layout.Row = 3;
            label.Layout.Column = 1;
            app.FilterTypeDropDown = uidropdown(optionalGrid, ...
                'Items',{'logarithmic','linear'}, ...
                'Value',app.Config.network.filter, ...
                'ValueChangedFcn',@(~,~)app.updateFWHMPreview());
            app.FilterTypeDropDown.Layout.Row = 3;
            app.FilterTypeDropDown.Layout.Column = 2;

            label = uilabel(optionalGrid,'Text','Regularisation');
            label.Layout.Row = 4;
            label.Layout.Column = 1;
            app.RegularisationField = uieditfield(optionalGrid,'numeric', ...
                'Value',app.Config.network.regularisation, ...
                'Limits',[0 1], ...
                'LowerLimitInclusive','on', ...
                'UpperLimitInclusive','off');
            app.RegularisationField.Layout.Row = 4;
            app.RegularisationField.Layout.Column = 2;

            label = uilabel(optionalGrid,'Text','Retained components');
            label.Layout.Row = 5;
            label.Layout.Column = 1;
            app.ComponentsField = uieditfield(optionalGrid,'numeric', ...
                'Value',app.Config.network.ncomps, ...
                'Limits',[1 Inf], ...
                'RoundFractionalValues','on');
            app.ComponentsField.Layout.Row = 5;
            app.ComponentsField.Layout.Column = 2;

            app.BadSegmentsField = app.addTextSetting(optionalGrid,6, ...
                'Bad sample indices','', 'Example: 1:100 350:420');

            checkGrid = uigridlayout(optionalGrid,[1 2]);
            checkGrid.Layout.Row = 7;
            checkGrid.Layout.Column = [1 2];
            checkGrid.ColumnWidth = {'1x','1x'};
            checkGrid.Padding = [0 0 0 0];
            app.RescaleCheckBox = uicheckbox(checkGrid, ...
                'Text','Rescale low-amplitude data', ...
                'Value',app.Config.network.rescale);
            app.RescaleCheckBox.Layout.Column = 1;
            app.RecomputeCheckBox = uicheckbox(checkGrid, ...
                'Text','Recompute existing outputs', ...
                'Value',app.Config.execution.recomputeExisting);
            app.RecomputeCheckBox.Layout.Column = 2;

            outputPanel = uipanel(grid, ...
                'Title','Participant outputs', ...
                'FontWeight','bold', ...
                'BackgroundColor',[1 1 1]);
            outputPanel.Layout.Column = 3;
            outputGrid = uigridlayout(outputPanel,[8 1]);
            outputGrid.RowHeight = {20,34,48,28,48,44,32,'1x'};
            outputGrid.Padding = [10 8 10 8];
            outputGrid.RowSpacing = 5;

            label = uilabel(outputGrid,'Text','Derived output folder', ...
                'FontWeight','bold');
            label.Layout.Row = 1;
            app.NetworkOutputField = uieditfield(outputGrid,'text', ...
                'Editable','off', ...
                'Placeholder','FREQ_Networks* is derived from Dataset*');
            app.NetworkOutputField.Layout.Row = 2;

            note = uilabel(outputGrid, ...
                'Text',sprintf(['Files: <participant>_FREQ.mat\n' ...
                'Index: FREQNESS_Manifest.mat']), ...
                'FontColor',colors.muted, ...
                'VerticalAlignment','top');
            note.Layout.Row = 3;

            resumeLabel = uilabel(outputGrid, ...
                'Text','Existing results are retained by default, enabling resumable runs.', ...
                'FontColor',colors.green);
            resumeLabel.Layout.Row = 4;

            app.RunStatusLabel = uilabel(outputGrid, ...
                'Text','Import a valid Dataset* folder to begin.', ...
                'FontColor',colors.muted, ...
                'VerticalAlignment','center');
            app.RunStatusLabel.Layout.Row = 5;

            app.RunButton = uibutton(outputGrid,'push', ...
                'Text','Run Network Estimation', ...
                'Enable','off', ...
                'FontSize',14, ...
                'FontWeight','bold', ...
                'FontColor',[1 1 1], ...
                'BackgroundColor',colors.blue, ...
                'ButtonPushedFcn',@(~,~)app.runCoreAnalysis());
            app.RunButton.Layout.Row = 6;

            app.SecondaryPageButton = uibutton(outputGrid,'push', ...
                'Text','Continue to Secondary Analyses  →', ...
                'FontWeight','bold', ...
                'FontColor',colors.blue, ...
                'ButtonPushedFcn',@(~,~)app.openSecondaryPage());
            app.SecondaryPageButton.Layout.Row = 7;

            app.ProgressTextArea = uitextarea(outputGrid, ...
                'Editable','off', ...
                'Value',{'Progress messages will appear here.'}, ...
                'FontName','Courier New', ...
                'FontSize',10);
            app.ProgressTextArea.Layout.Row = 8;

            app.applyFrequencyRange([app.Config.network.frequencies(1), ...
                app.Config.network.frequencies(end)]);
        end

        function createSecondaryPage(app,parent,colors)
            pageGrid = uigridlayout(parent,[4 1]);
            pageGrid.RowHeight = {72,142,'1x',26};
            pageGrid.Padding = [18 14 18 10];
            pageGrid.RowSpacing = 10;
            pageGrid.BackgroundColor = colors.background;

            headerGrid = uigridlayout(pageGrid,[2 1]);
            headerGrid.Layout.Row = 1;
            headerGrid.RowHeight = {38,24};
            headerGrid.Padding = [12 3 12 2];
            headerGrid.RowSpacing = 0;
            headerGrid.BackgroundColor = colors.navy;
            titleLabel = uilabel(headerGrid, ...
                'Text','FREQ-NESS Secondary Analyses', ...
                'FontSize',25, ...
                'FontWeight','bold', ...
                'FontColor',[1 1 1]);
            titleLabel.Layout.Row = 1;
            subtitleLabel = uilabel(headerGrid, ...
                'Text','Select one analysis, configure only its relevant inputs, and review readiness.', ...
                'FontSize',12, ...
                'FontColor',[0.80 0.87 0.94]);
            subtitleLabel.Layout.Row = 2;

            contextPanel = uipanel(pageGrid, ...
                'Title','Analysis context', ...
                'FontWeight','bold', ...
                'ForegroundColor',colors.navy, ...
                'BackgroundColor',[1 1 1]);
            contextPanel.Layout.Row = 2;
            contextGrid = uigridlayout(contextPanel,[2 7]);
            contextGrid.RowHeight = {36,'1x'};
            contextGrid.ColumnWidth = {210,'1x',95,160,160,160,185};
            contextGrid.Padding = [10 7 10 8];
            contextGrid.RowSpacing = 7;
            contextGrid.ColumnSpacing = 8;

            app.NetworkDropLabel = uilabel(contextGrid, ...
                'Text',sprintf('Drop a FREQ_Networks* folder here\nor use Browse'), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','center', ...
                'FontWeight','bold', ...
                'FontColor',colors.blue, ...
                'BackgroundColor',colors.softBlue);
            app.NetworkDropLabel.Layout.Row = [1 2];
            app.NetworkDropLabel.Layout.Column = 1;

            app.NetworkPathField = uieditfield(contextGrid,'text', ...
                'Editable','off', ...
                'Placeholder','No FREQ_Networks* folder selected');
            app.NetworkPathField.Layout.Row = 1;
            app.NetworkPathField.Layout.Column = 2;
            app.NetworkBrowseButton = uibutton(contextGrid,'push', ...
                'Text','Browse...', ...
                'ButtonPushedFcn',@(~,~)app.browseNetworkFolder());
            app.NetworkBrowseButton.Layout.Row = 1;
            app.NetworkBrowseButton.Layout.Column = 3;

            app.AnalysisOutputField = uieditfield(contextGrid,'text', ...
                'Editable','off', ...
                'Placeholder','Derived FREQ_Analyses* output folder');
            app.AnalysisOutputField.Layout.Row = 1;
            app.AnalysisOutputField.Layout.Column = [4 7];

            app.SecondarySummaryLabel = app.addContextBadge( ...
                contextGrid,2,2,'Participants','Waiting for FREQ results',colors);
            app.SecondaryParticipantButton = uibutton(contextGrid,'push', ...
                'Text','Select...', ...
                'Enable','off', ...
                'ButtonPushedFcn',@(~,~)app.selectSecondaryParticipants());
            app.SecondaryParticipantButton.Layout.Row = 2;
            app.SecondaryParticipantButton.Layout.Column = 3;
            app.SecondaryFrequencyLabel = app.addContextBadge( ...
                contextGrid,2,4,'Frequencies','—',colors);
            app.SecondaryComponentsLabel = app.addContextBadge( ...
                contextGrid,2,5,'Components','—',colors);
            app.SecondaryMNILabel = app.addContextBadge( ...
                contextGrid,2,6,'MNI coordinates','—',colors);
            app.SecondarySourceLabel = app.addContextBadge( ...
                contextGrid,2,7,'Source data','—',colors);

            workspaceGrid = uigridlayout(pageGrid,[1 3]);
            workspaceGrid.Layout.Row = 3;
            workspaceGrid.ColumnWidth = {260,570,'1x'};
            workspaceGrid.Padding = [0 0 0 0];
            workspaceGrid.ColumnSpacing = 12;

            browserPanel = uipanel(workspaceGrid, ...
                'Title','Analysis browser', ...
                'FontWeight','bold', ...
                'BackgroundColor',[1 1 1]);
            browserPanel.Layout.Column = 1;
            browserGrid = uigridlayout(browserPanel,[2 1]);
            browserGrid.RowHeight = {42,'1x'};
            browserGrid.Padding = [8 7 8 8];
            browserGrid.RowSpacing = 6;
            note = uilabel(browserGrid, ...
                'Text','Choose one module. Only its settings appear.', ...
                'FontColor',colors.muted, ...
                'WordWrap','on', ...
                'VerticalAlignment','top');
            note.Layout.Row = 1;
            app.SecondaryTree = uitree(browserGrid, ...
                'SelectionChangedFcn',@(~,event)app.secondaryTreeSelection(event));
            app.SecondaryTree.Layout.Row = 2;
            app.populateSecondaryTree();

            app.SecondaryConfigPanel = uipanel(workspaceGrid, ...
                'Title','Analysis configuration', ...
                'FontWeight','bold', ...
                'BackgroundColor',[1 1 1]);
            app.SecondaryConfigPanel.Layout.Column = 2;

            readinessPanel = uipanel(workspaceGrid, ...
                'Title','Readiness and execution', ...
                'FontWeight','bold', ...
                'BackgroundColor',[1 1 1]);
            readinessPanel.Layout.Column = 3;
            readinessGrid = uigridlayout(readinessPanel,[9 1]);
            readinessGrid.RowHeight = {20,112,20,34,44,24,40,34,'1x'};
            readinessGrid.Padding = [10 8 10 8];
            readinessGrid.RowSpacing = 6;

            label = uilabel(readinessGrid,'Text','Input readiness', ...
                'FontWeight','bold');
            label.Layout.Row = 1;
            app.SecondaryReadinessTextArea = uitextarea(readinessGrid, ...
                'Editable','off', ...
                'Value',{'Import a FREQ_Networks* folder to begin.'});
            app.SecondaryReadinessTextArea.Layout.Row = 2;
            label = uilabel(readinessGrid,'Text','Function output folder', ...
                'FontWeight','bold');
            label.Layout.Row = 3;
            app.SecondaryModuleOutputField = uieditfield(readinessGrid,'text', ...
                'Editable','off', ...
                'Placeholder','One subfolder per analysis function');
            app.SecondaryModuleOutputField.Layout.Row = 4;
            app.SecondaryValidationLabel = uilabel(readinessGrid, ...
                'Text','Select a module and import network results.', ...
                'FontColor',colors.muted, ...
                'WordWrap','on');
            app.SecondaryValidationLabel.Layout.Row = 5;
            backgroundAvailable = freqnessgui.secondaryBackgroundAvailable();
            app.SecondaryBackgroundCheckBox = uicheckbox(readinessGrid, ...
                'Text','Run in background (keeps the GUI responsive)', ...
                'Value',backgroundAvailable, ...
                'Enable','on', ...
                'Tooltip',['Recommended for long analyses. Clear this box ' ...
                'to run synchronously in the current MATLAB session.']);
            if ~backgroundAvailable
                app.SecondaryBackgroundCheckBox.Text = ...
                    'Background execution unavailable in this MATLAB installation';
                app.SecondaryBackgroundCheckBox.Enable = 'off';
            end
            app.SecondaryBackgroundCheckBox.Layout.Row = 6;
            app.SecondaryValidateButton = uibutton(readinessGrid,'push', ...
                'Text','Validate Configuration', ...
                'Enable','off', ...
                'FontWeight','bold', ...
                'ButtonPushedFcn',@(~,~)app.validateSecondaryConfiguration());
            app.SecondaryValidateButton.Layout.Row = 7;
            app.SecondaryRunButton = uibutton(readinessGrid,'push', ...
                'Text','Run Analysis', ...
                'Enable','off', ...
                'FontWeight','bold', ...
                'FontColor',[1 1 1], ...
                'BackgroundColor',colors.blue, ...
                'Tooltip','Run the selected analysis and persist its outputs.', ...
                'ButtonPushedFcn',@(~,~)app.handleSecondaryRunButton());
            app.SecondaryRunButton.Layout.Row = 8;
            app.SecondaryProgressTextArea = uitextarea(readinessGrid, ...
                'Editable','off', ...
                'Value',{'Validate a configuration, then run the analysis.', ...
                    'Participant and group outputs will be saved automatically.'}, ...
                'FontName','Courier New', ...
                'FontSize',10);
            app.SecondaryProgressTextArea.Layout.Row = 9;

            footer = uilabel(pageGrid, ...
                'Text','Secondary analyses use one derived output subfolder per backend function.', ...
                'FontSize',11, ...
                'FontColor',colors.muted);
            footer.Layout.Row = 4;

            app.selectSecondaryModule('entropy');
        end

        function badge = addContextBadge(~,parent,row,column,titleText,valueText,colors)
            badge = uilabel(parent, ...
                'Text',sprintf('%s\n%s',titleText,valueText), ...
                'FontSize',11, ...
                'FontColor',colors.muted, ...
                'BackgroundColor',colors.softBlue, ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','center');
            badge.Layout.Row = row;
            badge.Layout.Column = column;
        end

        function populateSecondaryTree(app)
            modules = freqnessgui.moduleRegistry();
            categories = unique({modules.category},'stable');
            firstLeaf = [];
            for categoryi = 1:numel(categories)
                categoryNode = uitreenode(app.SecondaryTree, ...
                    'Text',categories{categoryi}, ...
                    'NodeData','');
                categoryModules = modules(strcmp({modules.category}, ...
                    categories{categoryi}));
                for modulei = 1:numel(categoryModules)
                    leaf = uitreenode(categoryNode, ...
                        'Text',categoryModules(modulei).name, ...
                        'NodeData',categoryModules(modulei).id);
                    if isempty(firstLeaf)
                        firstLeaf = leaf;
                    end
                end
            end
            try
                expand(app.SecondaryTree,'all');
            catch
                % Tree categories can still be expanded manually.
            end
            if ~isempty(firstLeaf)
                app.SecondaryTree.SelectedNodes = firstLeaf;
            end
        end

        function secondaryTreeSelection(app,event)
            if isempty(event.SelectedNodes)
                return
            end
            moduleId = event.SelectedNodes(1).NodeData;
            if ischar(moduleId) && ~isempty(moduleId)
                app.selectSecondaryModule(moduleId);
            elseif isstring(moduleId) && isscalar(moduleId) && strlength(moduleId) > 0
                app.selectSecondaryModule(char(moduleId));
            end
        end

        function module = secondaryModule(~,moduleId)
            modules = freqnessgui.moduleRegistry();
            moduleIndex = find(strcmp({modules.id},moduleId),1);
            if isempty(moduleIndex)
                error('FREQNESS:GUI:UnknownSecondaryModule', ...
                    'Unknown secondary-analysis module: %s.',moduleId);
            end
            module = modules(moduleIndex);
        end

        function context = secondaryContext(app)
            context = struct();
            if isempty(app.NetworkSet)
                context.frequencies = app.Config.network.frequencies;
                context.nComponents = app.Config.network.ncomps;
                context.samplingRate = app.Config.network.samplingRate;
            else
                context.frequencies = app.NetworkSet.frequencies;
                context.nComponents = app.NetworkSet.nComponents;
                context.samplingRate = app.NetworkSet.samplingRate;
                if isempty(context.frequencies)
                    context.frequencies = app.Config.network.frequencies;
                end
                if isempty(context.nComponents)
                    context.nComponents = app.Config.network.ncomps;
                end
            end
        end

        function selectSecondaryModule(app,moduleId)
            module = app.secondaryModule(moduleId);
            app.SecondarySelectedModuleId = module.id;
            app.SecondarySchema = freqnessgui.secondaryModuleSchema( ...
                module.id,app.secondaryContext());
            app.buildSecondaryConfiguration(module);

            if ~isempty(app.NetworkSet)
                moduleFolder = fullfile(app.NetworkSet.analysisFolder, ...
                    module.folderName);
                app.SecondaryModuleOutputField.Value = moduleFolder;
                app.SecondaryModuleOutputField.Tooltip = moduleFolder;
            else
                app.SecondaryModuleOutputField.Value = '';
                app.SecondaryModuleOutputField.Placeholder = ...
                    ['FREQ_Analyses*/' module.folderName];
            end
            app.updateSecondaryReadiness();
        end

        function buildSecondaryConfiguration(app,module)
            delete(app.SecondaryConfigPanel.Children);
            app.SecondaryFields = struct();
            app.SecondaryConfigPanel.Title = ['Configuration — ' module.name];

            app.SecondaryConfigGrid = uigridlayout( ...
                app.SecondaryConfigPanel,[6 1]);
            app.SecondaryConfigGrid.RowHeight = {62,24, ...
                app.secondaryPanelHeight(app.SecondarySchema.mandatory), ...
                34,0,'1x'};
            app.SecondaryConfigGrid.Padding = [10 8 10 8];
            app.SecondaryConfigGrid.RowSpacing = 7;
            app.SecondaryConfigGrid.Scrollable = 'on';

            description = uilabel(app.SecondaryConfigGrid, ...
                'Text',module.description, ...
                'FontSize',13, ...
                'FontColor',[0.20 0.24 0.29], ...
                'WordWrap','on', ...
                'VerticalAlignment','top');
            description.Layout.Row = 1;

            requirements = {'FREQ results'};
            if module.requiresMNI
                requirements{end+1} = 'MNI coordinates';
            end
            if module.requiresSourceData
                requirements{end+1} = 'original source data';
            end
            if module.requiresEvents
                requirements{end+1} = 'events file';
            end
            requirementLabel = uilabel(app.SecondaryConfigGrid, ...
                'Text',['Requires: ' strjoin(requirements,'  •  ')], ...
                'FontWeight','bold', ...
                'FontColor',[0.055 0.415 0.690]);
            requirementLabel.Layout.Row = 2;

            app.createSecondaryFieldPanel(app.SecondaryConfigGrid,3, ...
                'Mandatory inputs',app.SecondarySchema.mandatory,true);

            app.SecondaryAdvancedButton = uibutton( ...
                app.SecondaryConfigGrid,'state', ...
                'Text','Show optional settings  ▾', ...
                'Value',false, ...
                'ValueChangedFcn',@(~,~)app.toggleSecondaryOptions());
            app.SecondaryAdvancedButton.Layout.Row = 4;

            app.SecondaryOptionalPanel = app.createSecondaryFieldPanel( ...
                app.SecondaryConfigGrid,5,'Optional settings', ...
                app.SecondarySchema.optional,false);
            app.SecondaryOptionalPanel.Visible = 'off';
            if isempty(app.SecondarySchema.optional)
                app.SecondaryAdvancedButton.Text = 'No optional settings';
                app.SecondaryAdvancedButton.Enable = 'off';
            end

            hint = uilabel(app.SecondaryConfigGrid, ...
                'Text',['Values are constrained by the imported frequencies ' ...
                    'and retained components.'], ...
                'FontColor',[0.36 0.41 0.47], ...
                'VerticalAlignment','top');
            hint.Layout.Row = 6;
        end

        function height = secondaryPanelHeight(app,fields)
            if isempty(fields)
                height = 68;
                return
            end
            fieldHeights = arrayfun(@(field) ...
                app.secondaryFieldHeight(field),fields);
            height = max(68,34+sum(fieldHeights)+6*(numel(fields)-1));
        end

        function height = secondaryFieldHeight(~,field)
            if ismember(field.type,{'frequency','frequencyRange'})
                height = 72;
            else
                height = 31;
            end
        end

        function panel = createSecondaryFieldPanel(app,parent,row,titleText,fields,isRequired)
            panel = uipanel(parent, ...
                'Title',titleText, ...
                'FontWeight','bold', ...
                'BackgroundColor',[1 1 1]);
            panel.Layout.Row = row;

            if isempty(fields)
                panelGrid = uigridlayout(panel,[1 1]);
                panelGrid.Padding = [10 7 10 7];
                label = uilabel(panelGrid, ...
                    'Text','No additional inputs are required.', ...
                    'FontColor',[0.36 0.41 0.47]);
                label.Layout.Row = 1;
                return
            end

            panelGrid = uigridlayout(panel,[numel(fields) 2]);
            panelGrid.ColumnWidth = {190,'1x'};
            panelGrid.RowHeight = arrayfun(@(field) ...
                app.secondaryFieldHeight(field),fields, ...
                'UniformOutput',false);
            panelGrid.Padding = [10 7 10 7];
            panelGrid.RowSpacing = 6;
            panelGrid.ColumnSpacing = 8;
            for fieldi = 1:numel(fields)
                app.addSecondaryControl(panelGrid,fieldi,fields(fieldi),isRequired);
            end
        end

        function addSecondaryControl(app,parent,row,field,isRequired)
            label = uilabel(parent,'Text',field.label, ...
                'Tooltip',field.help);
            label.Layout.Row = row;
            label.Layout.Column = 1;

            control = struct('type',field.type,'handles',{{}}, ...
                'required',isRequired,'label',field.label, ...
                'frequencyValues',[]);
            switch field.type
                case {'frequency','frequencyRange'}
                    context = app.secondaryContext();
                    frequencies = context.frequencies(:)';
                    [defaultIndices,~] = freqnessgui.mapFrequencySelection( ...
                        frequencies,field.default);
                    [majorTicks,majorLabels] = ...
                        app.secondaryFrequencyTicks(frequencies);
                    sliderGrid = uigridlayout(parent,[2 1]);
                    sliderGrid.Layout.Row = row;
                    sliderGrid.Layout.Column = 2;
                    sliderGrid.RowHeight = {45,20};
                    sliderGrid.Padding = [0 0 0 0];
                    sliderGrid.RowSpacing = 0;
                    exactLabel = uilabel(sliderGrid, ...
                        'FontColor',[0.055 0.415 0.690], ...
                        'HorizontalAlignment','center');
                    exactLabel.Layout.Row = 2;
                    isRange = strcmp(field.type,'frequencyRange');
                    if isscalar(frequencies)
                        sliderLimits = [1 2];
                        sliderEnabled = 'off';
                    else
                        sliderLimits = [1 numel(frequencies)];
                        sliderEnabled = 'on';
                    end
                    if isRange
                        frequencySlider = uislider(sliderGrid,'range', ...
                            'Limits',sliderLimits, ...
                            'Value',defaultIndices, ...
                            'Step',1, ...
                            'MajorTicks',majorTicks, ...
                            'MajorTickLabels',majorLabels, ...
                            'MinorTicks',[], ...
                            'Enable',sliderEnabled, ...
                            'Tooltip',field.help, ...
                            'ValueChangingFcn',@(~,event) ...
                                app.previewSecondaryFrequencySlider( ...
                                event.Value,exactLabel,frequencies,true), ...
                            'ValueChangedFcn',@(source,event) ...
                                app.commitSecondaryFrequencySlider( ...
                                source,event.Value,exactLabel,frequencies,true));
                    else
                        frequencySlider = uislider(sliderGrid, ...
                            'Limits',sliderLimits, ...
                            'Value',defaultIndices, ...
                            'MajorTicks',majorTicks, ...
                            'MajorTickLabels',majorLabels, ...
                            'MinorTicks',[], ...
                            'Enable',sliderEnabled, ...
                            'Tooltip',field.help, ...
                            'ValueChangingFcn',@(~,event) ...
                                app.previewSecondaryFrequencySlider( ...
                                event.Value,exactLabel,frequencies,false), ...
                            'ValueChangedFcn',@(source,event) ...
                                app.commitSecondaryFrequencySlider( ...
                                source,event.Value,exactLabel,frequencies,false));
                    end
                    frequencySlider.Layout.Row = 1;
                    app.previewSecondaryFrequencySlider(defaultIndices, ...
                        exactLabel,frequencies,isRange);
                    label.VerticalAlignment = 'top';
                    control.handles = {frequencySlider,exactLabel};
                    control.frequencyValues = frequencies;

                case {'component','number','integer'}
                    limits = [-Inf Inf];
                    if strcmp(field.type,'component')
                        context = app.secondaryContext();
                        limits = [1 max(1,context.nComponents)];
                    elseif ismember(field.type,{'number','integer'})
                        limits = [0 Inf];
                    end
                    numericField = uieditfield(parent,'numeric', ...
                        'Value',field.default, ...
                        'Limits',limits, ...
                        'Tooltip',field.help, ...
                        'ValueChangedFcn',@(~,~)app.updateSecondaryReadiness());
                    if ismember(field.type,{'component','integer'})
                        numericField.RoundFractionalValues = 'on';
                    end
                    numericField.Layout.Row = row;
                    numericField.Layout.Column = 2;
                    control.handles = {numericField};

                case 'timeRange'
                    rangeGrid = uigridlayout(parent,[1 3]);
                    rangeGrid.Layout.Row = row;
                    rangeGrid.Layout.Column = 2;
                    rangeGrid.ColumnWidth = {'1x',24,'1x'};
                    rangeGrid.Padding = [0 0 0 0];
                    rangeGrid.ColumnSpacing = 5;
                    lowerField = uieditfield(rangeGrid,'numeric', ...
                        'Value',field.default(1), ...
                        'ValueChangedFcn',@(~,~)app.updateSecondaryReadiness());
                    lowerField.Layout.Column = 1;
                    connector = uilabel(rangeGrid,'Text','to', ...
                        'HorizontalAlignment','center');
                    connector.Layout.Column = 2;
                    upperField = uieditfield(rangeGrid,'numeric', ...
                        'Value',field.default(2), ...
                        'ValueChangedFcn',@(~,~)app.updateSecondaryReadiness());
                    upperField.Layout.Column = 3;
                    lowerField.Tooltip = field.help;
                    upperField.Tooltip = field.help;
                    control.handles = {lowerField,upperField};

                case {'componentVector','numberOrEmpty'}
                    textField = uieditfield(parent,'text', ...
                        'Value',char(string(field.default)), ...
                        'Tooltip',field.help, ...
                        'ValueChangedFcn',@(~,~)app.updateSecondaryReadiness());
                    textField.Layout.Row = row;
                    textField.Layout.Column = 2;
                    control.handles = {textField};

                case 'logical'
                    checkBox = uicheckbox(parent, ...
                        'Text','Enabled', ...
                        'Value',field.default, ...
                        'Tooltip',field.help, ...
                        'ValueChangedFcn',@(~,~)app.updateSecondaryReadiness());
                    checkBox.Layout.Row = row;
                    checkBox.Layout.Column = 2;
                    control.handles = {checkBox};

                case 'file'
                    fileGrid = uigridlayout(parent,[1 2]);
                    fileGrid.Layout.Row = row;
                    fileGrid.Layout.Column = 2;
                    fileGrid.ColumnWidth = {'1x',82};
                    fileGrid.Padding = [0 0 0 0];
                    fileGrid.ColumnSpacing = 6;
                    fileField = uieditfield(fileGrid,'text', ...
                        'Value',field.default, ...
                        'Placeholder','Select a MAT-file', ...
                        'ValueChangedFcn',@(~,~)app.updateSecondaryReadiness());
                    fileField.Layout.Column = 1;
                    browseButton = uibutton(fileGrid,'push', ...
                        'Text','Browse...', ...
                        'ButtonPushedFcn',@(~,~)app.browseSecondaryFile(field.key));
                    browseButton.Layout.Column = 2;
                    control.handles = {fileField,browseButton};

                otherwise
                    error('FREQNESS:GUI:UnknownSecondaryFieldType', ...
                        'Unsupported secondary field type: %s.',field.type);
            end
            app.SecondaryFields.(field.key) = control;
        end

        function [ticks,labels] = secondaryFrequencyTicks(~,frequencies)
            tickCount = min(7,numel(frequencies));
            ticks = unique(round(linspace(1,numel(frequencies),tickCount)));
            labels = arrayfun(@(index) ...
                sprintf('%.4g',frequencies(index)),ticks, ...
                'UniformOutput',false);
        end

        function previewSecondaryFrequencySlider(~,rawValue,label, ...
                frequencies,isRange)
            indices = round(rawValue);
            indices = max(1,min(numel(frequencies),indices));
            if isRange
                indices = sort(indices);
                if indices(1) == indices(2)
                    label.Text = sprintf('%.4g Hz  |  1 frequency', ...
                        frequencies(indices(1)));
                else
                    label.Text = sprintf('%.4g–%.4g Hz  |  %d frequencies', ...
                        frequencies(indices(1)),frequencies(indices(2)), ...
                        indices(2)-indices(1)+1);
                end
            else
                label.Text = sprintf('Selected: %.4g Hz',frequencies(indices));
            end
        end

        function commitSecondaryFrequencySlider(app,slider,rawValue,label, ...
                frequencies,isRange)
            indices = round(rawValue);
            indices = max(1,min(numel(frequencies),indices));
            if isRange
                indices = sort(indices);
            end
            slider.Value = indices;
            app.previewSecondaryFrequencySlider(indices,label,frequencies,isRange);
            app.updateSecondaryReadiness();
        end

        function toggleSecondaryOptions(app)
            if app.SecondaryAdvancedButton.Value
                app.SecondaryAdvancedButton.Text = 'Hide optional settings  ▴';
                app.SecondaryOptionalPanel.Visible = 'on';
                rowHeights = app.SecondaryConfigGrid.RowHeight;
                rowHeights{5} = app.secondaryPanelHeight( ...
                    app.SecondarySchema.optional);
                app.SecondaryConfigGrid.RowHeight = rowHeights;
            else
                app.SecondaryAdvancedButton.Text = 'Show optional settings  ▾';
                app.SecondaryOptionalPanel.Visible = 'off';
                rowHeights = app.SecondaryConfigGrid.RowHeight;
                rowHeights{5} = 0;
                app.SecondaryConfigGrid.RowHeight = rowHeights;
            end
        end

        function browseSecondaryFile(app,fieldKey)
            [selectedFile,selectedFolder] = uigetfile( ...
                {'*.mat','MAT-files (*.mat)'}, ...
                'Select a secondary-analysis input file');
            if isequal(selectedFile,0)
                return
            end
            control = app.SecondaryFields.(fieldKey);
            control.handles{1}.Value = fullfile(selectedFolder,selectedFile);
            app.updateSecondaryReadiness();
        end

        function openSecondaryPage(app)
            if isempty(app.NetworkSet) && ~isempty(app.NetworkOutputField.Value) && ...
                    isfolder(app.NetworkOutputField.Value)
                try
                    app.importNetworkFolder(app.NetworkOutputField.Value);
                catch
                    % Page 2 remains available for importing an existing folder.
                end
            end
            app.PageGroup.SelectedTab = app.SecondaryTab;
            app.updateSecondaryReadiness();
        end

        function selectSecondaryParticipants(app)
            if isempty(app.NetworkSet)
                return
            end
            dialog = uifigure( ...
                'Name','Select participants', ...
                'Position',[420 260 420 480], ...
                'WindowStyle','modal');
            dialogGrid = uigridlayout(dialog,[3 1]);
            dialogGrid.RowHeight = {48,'1x',40};
            dialogGrid.Padding = [14 12 14 12];
            prompt = uilabel(dialogGrid, ...
                'Text','Choose the participant results included in this analysis.', ...
                'WordWrap','on');
            prompt.Layout.Row = 1;
            participantList = uilistbox(dialogGrid, ...
                'Items',app.NetworkSet.participantIds, ...
                'Multiselect','on', ...
                'Value',app.SecondarySelectedParticipants);
            participantList.Layout.Row = 2;
            buttonGrid = uigridlayout(dialogGrid,[1 2]);
            buttonGrid.Layout.Row = 3;
            buttonGrid.ColumnWidth = {'1x','1x'};
            buttonGrid.Padding = [0 0 0 0];
            cancelButton = uibutton(buttonGrid,'push', ...
                'Text','Cancel', ...
                'ButtonPushedFcn',@(~,~)delete(dialog));
            cancelButton.Layout.Column = 1;
            applyButton = uibutton(buttonGrid,'push', ...
                'Text','Apply selection', ...
                'FontWeight','bold', ...
                'ButtonPushedFcn',@(~,~) ...
                    app.applySecondaryParticipantSelection( ...
                    dialog,participantList));
            applyButton.Layout.Column = 2;
        end

        function applySecondaryParticipantSelection(app,dialog,participantList)
            selected = participantList.Value;
            if ischar(selected) || isstring(selected)
                selected = cellstr(selected);
            end
            if isempty(selected)
                uialert(dialog,'Select at least one participant.', ...
                    'Empty participant selection');
                return
            end
            app.SecondarySelectedParticipants = selected(:)';
            delete(dialog);
            app.updateSecondaryContext();
        end

        function updateSecondaryContext(app)
            if isempty(app.SecondarySummaryLabel) || ...
                    ~isvalid(app.SecondarySummaryLabel)
                return
            end
            if isempty(app.NetworkSet)
                app.SecondarySummaryLabel.Text = sprintf( ...
                    'Participants\nWaiting for FREQ results');
                app.SecondaryFrequencyLabel.Text = sprintf('Frequencies\n—');
                app.SecondaryComponentsLabel.Text = sprintf('Components\n—');
                app.SecondarySourceLabel.Text = sprintf('Source data\n—');
                app.SecondaryParticipantButton.Enable = 'off';
            else
                app.SecondarySummaryLabel.Text = sprintf( ...
                    'Participants\n%d/%d selected', ...
                    numel(app.SecondarySelectedParticipants), ...
                    app.NetworkSet.nParticipants);
                app.SecondaryParticipantButton.Enable = 'on';
                if isempty(app.NetworkSet.frequencies)
                    frequencyText = 'unknown';
                else
                    frequencyText = sprintf('%.4g–%.4g Hz', ...
                        app.NetworkSet.frequencies(1), ...
                        app.NetworkSet.frequencies(end));
                end
                app.SecondaryFrequencyLabel.Text = sprintf( ...
                    'Frequencies\n%s',frequencyText);
                if isempty(app.NetworkSet.nComponents)
                    componentText = 'unknown';
                else
                    componentText = sprintf('%d retained', ...
                        app.NetworkSet.nComponents);
                end
                app.SecondaryComponentsLabel.Text = sprintf( ...
                    'Components\n%s',componentText);
                if app.NetworkSet.sourceDataAvailable
                    sourceText = 'available';
                else
                    sourceText = 'not linked';
                end
                app.SecondarySourceLabel.Text = sprintf( ...
                    'Source data\n%s',sourceText);
            end
            if isempty(app.MNI)
                mniText = 'not loaded';
            else
                mniText = sprintf('%d points',app.MNI.nPoints);
            end
            app.SecondaryMNILabel.Text = sprintf('MNI coordinates\n%s',mniText);
            app.updateSecondaryReadiness();
        end

        function [configuration,message] = collectSecondaryConfiguration(app)
            if isempty(app.NetworkSet)
                error('FREQNESS:GUI:MissingNetworkResults', ...
                    'Import a FREQ_Networks* folder.');
            end
            module = app.secondaryModule(app.SecondarySelectedModuleId);
            if module.requiresMNI && isempty(app.MNI)
                error('FREQNESS:GUI:MissingMNI', ...
                    'This analysis requires MNI coordinates.');
            end
            if module.requiresMNI && ~isempty(app.NetworkSet.nVoxels) && ...
                    app.MNI.nPoints ~= app.NetworkSet.nVoxels
                error('FREQNESS:GUI:IncompatibleMNI', ...
                    ['MNI coordinates contain %d points, but the imported ' ...
                    'FREQ results contain %d voxels.'], ...
                    app.MNI.nPoints,app.NetworkSet.nVoxels);
            end
            if module.requiresSourceData && ~app.NetworkSet.sourceDataAvailable
                error('FREQNESS:GUI:MissingSourceData', ...
                    'This analysis requires the original Dataset* source data.');
            end

            values = struct();
            selectedFrequencies = struct();
            fieldNames = fieldnames(app.SecondaryFields);
            context = app.secondaryContext();
            for fieldi = 1:numel(fieldNames)
                key = fieldNames{fieldi};
                control = app.SecondaryFields.(key);
                handles = control.handles;
                switch control.type
                    case 'frequency'
                        index = round(handles{1}.Value);
                        index = max(1,min(numel(control.frequencyValues),index));
                        value = control.frequencyValues(index);
                        selectedFrequencies.(key) = value;
                    case 'frequencyRange'
                        indices = round(handles{1}.Value);
                        indices = sort(max(1,min( ...
                            numel(control.frequencyValues),indices)));
                        value = control.frequencyValues(indices);
                        selectedFrequencies.(key) = ...
                            control.frequencyValues(indices(1):indices(2));
                    case {'component','number','integer'}
                        value = handles{1}.Value;
                        if ~isfinite(value)
                            error('FREQNESS:GUI:InvalidSecondaryConfig', ...
                                '%s must be finite.',control.label);
                        end
                        if strcmp(control.type,'component') && ...
                                (value < 1 || value > context.nComponents)
                            error('FREQNESS:GUI:InvalidSecondaryConfig', ...
                                '%s exceeds the retained component count.',control.label);
                        end
                    case 'timeRange'
                        value = [handles{1}.Value handles{2}.Value];
                        if any(~isfinite(value)) || value(1) > value(2)
                            error('FREQNESS:GUI:InvalidSecondaryConfig', ...
                                '%s must contain ascending finite endpoints.',control.label);
                        end
                    case 'componentVector'
                        value = freqnessgui.parseNumericVector( ...
                            handles{1}.Value,false,control.label);
                        if any(value < 1) || any(value ~= round(value)) || ...
                                any(value > context.nComponents)
                            error('FREQNESS:GUI:InvalidSecondaryConfig', ...
                                ['%s must contain positive integer indices no ' ...
                                'greater than %d.'],control.label,context.nComponents);
                        end
                    case 'numberOrEmpty'
                        value = freqnessgui.parseNumericVector( ...
                            handles{1}.Value,true,control.label);
                        if ~isempty(value) && (~isscalar(value) || value <= 0)
                            error('FREQNESS:GUI:InvalidSecondaryConfig', ...
                                '%s must be empty or one positive value.',control.label);
                        end
                    case 'logical'
                        value = handles{1}.Value;
                    case 'file'
                        value = handles{1}.Value;
                        if control.required && ~isfile(value)
                            error('FREQNESS:GUI:InvalidSecondaryConfig', ...
                                '%s must be an existing MAT-file.',control.label);
                        end
                    otherwise
                        error('FREQNESS:GUI:UnknownSecondaryFieldType', ...
                            'Unsupported secondary field type: %s.',control.type);
                end
                values.(key) = value;
            end

            if strcmp(module.id,'induced_responses')
                if values.baselineWindow(1) < values.epochWindow(1) || ...
                        values.baselineWindow(2) > values.epochWindow(2)
                    error('FREQNESS:GUI:InvalidSecondaryConfig', ...
                        'The baseline window must lie inside the epoch window.');
                end
            elseif strcmp(module.id,'cross_coupling') && ...
                    values.lfoFrequency >= values.carrierRange(2)
                error('FREQNESS:GUI:InvalidSecondaryConfig', ...
                    'The LFO frequency must be below the carrier-range maximum.');
            elseif strcmp(module.id,'comp_gradients')
                expectedComponents = values.components(1):values.components(end);
                if ~isequal(values.components,expectedComponents)
                    error('FREQNESS:GUI:InvalidSecondaryConfig', ...
                        'Components to model must define one ascending contiguous range.');
                end
            end

            configuration = struct();
            configuration.schemaVersion = 1;
            configuration.moduleId = module.id;
            configuration.functionName = module.functionName;
            configuration.networkFolder = app.NetworkSet.folder;
            configuration.outputFolder = fullfile( ...
                app.NetworkSet.analysisFolder,module.folderName);
            configuration.values = values;
            configuration.selectedFrequencies = selectedFrequencies;
            if isempty(app.SecondarySelectedParticipants)
                error('FREQNESS:GUI:MissingParticipants', ...
                    'Select at least one participant result.');
            end
            configuration.participants = app.SecondarySelectedParticipants;
            message = sprintf('%s configuration is valid.',module.name);
        end

        function updateSecondaryReadiness(app)
            if isempty(app.SecondaryReadinessTextArea) || ...
                    ~isvalid(app.SecondaryReadinessTextArea) || ...
                    isempty(app.SecondarySelectedModuleId)
                return
            end
            module = app.secondaryModule(app.SecondarySelectedModuleId);
            lines = cell(0,1);
            if isempty(app.NetworkSet)
                lines{end+1,1} = '✗ FREQ results — not imported';
            else
                lines{end+1,1} = sprintf('✓ FREQ results — %d/%d participants', ...
                    numel(app.SecondarySelectedParticipants), ...
                    app.NetworkSet.nParticipants);
            end
            if module.requiresMNI
                if isempty(app.MNI)
                    lines{end+1,1} = '✗ MNI coordinates — required';
                elseif ~isempty(app.NetworkSet) && ...
                        ~isempty(app.NetworkSet.nVoxels) && ...
                        app.MNI.nPoints ~= app.NetworkSet.nVoxels
                    lines{end+1,1} = sprintf( ...
                        '✗ MNI coordinates — %d points / %d voxels', ...
                        app.MNI.nPoints,app.NetworkSet.nVoxels);
                else
                    lines{end+1,1} = sprintf('✓ MNI coordinates — %d points', ...
                        app.MNI.nPoints);
                end
            else
                lines{end+1,1} = '— MNI coordinates — not required';
            end
            if module.requiresSourceData
                if ~isempty(app.NetworkSet) && app.NetworkSet.sourceDataAvailable
                    lines{end+1,1} = '✓ Original source data — linked';
                else
                    lines{end+1,1} = '✗ Original source data — required';
                end
            else
                lines{end+1,1} = '— Original source data — not required';
            end
            if module.requiresEvents
                eventsReady = isfield(app.SecondaryFields,'eventsFile') && ...
                    isfile(app.SecondaryFields.eventsFile.handles{1}.Value);
                if eventsReady
                    lines{end+1,1} = '✓ Events file — selected';
                else
                    lines{end+1,1} = '✗ Events file — required';
                end
            else
                lines{end+1,1} = '— Events file — not required';
            end

            try
                [~,message] = app.collectSecondaryConfiguration();
                isReady = true;
                lines{end+1,1} = '✓ Configuration — valid';
                app.SecondaryValidationLabel.Text = message;
                app.SecondaryValidationLabel.FontColor = [0.125 0.545 0.365];
            catch exception
                isReady = false;
                lines{end+1,1} = '✗ Configuration — incomplete';
                app.SecondaryValidationLabel.Text = exception.message;
                app.SecondaryValidationLabel.FontColor = [0.78 0.24 0.16];
            end
            app.SecondaryReadinessTextArea.Value = lines;
            if isReady && ~app.IsRunning
                app.SecondaryValidateButton.Enable = 'on';
                app.SecondaryRunButton.Enable = 'on';
            else
                app.SecondaryValidateButton.Enable = 'off';
                app.SecondaryRunButton.Enable = 'off';
            end
        end

        function validateSecondaryConfiguration(app)
            try
                [configuration,message] = app.collectSecondaryConfiguration();
                app.SecondaryValidationLabel.Text = message;
                app.SecondaryValidationLabel.FontColor = [0.125 0.545 0.365];
                app.SecondaryProgressTextArea.Value = { ...
                    sprintf('%s is ready.',configuration.functionName), ...
                    sprintf('Output: %s',configuration.outputFolder), ...
                    'Press Run Analysis to execute and persist the results.'};
            catch exception
                app.SecondaryValidationLabel.Text = exception.message;
                app.SecondaryValidationLabel.FontColor = [0.78 0.24 0.16];
                uialert(app.UIFigure,exception.message, ...
                    'Invalid secondary-analysis configuration');
            end
        end

        function handleSecondaryRunButton(app)
            if app.IsRunning
                if app.SecondaryRunInBackground
                    app.requestSecondaryCancellation();
                end
                return
            end
            app.runSecondaryAnalysis();
        end

        function runSecondaryAnalysis(app)
            if app.IsRunning
                return
            end
            try
                [configuration,~] = app.collectSecondaryConfiguration();
            catch exception
                uialert(app.UIFigure,exception.message, ...
                    'Invalid secondary-analysis configuration');
                return
            end

            mniCoordinates = [];
            if ~isempty(app.MNI)
                mniCoordinates = app.MNI.coordinates;
            end
            useBackground = app.SecondaryBackgroundCheckBox.Value && ...
                freqnessgui.secondaryBackgroundAvailable();
            if useBackground
                app.startSecondaryBackgroundAnalysis( ...
                    configuration,mniCoordinates);
            else
                app.runSecondarySynchronously(configuration,mniCoordinates);
            end
        end

        function runSecondarySynchronously(app,configuration,mniCoordinates)
            app.setSecondaryRunningState(true,false);
            runningCleanup = onCleanup( ...
                @()app.setSecondaryRunningState(false));
            try
                app.updateSecondaryProgress(0, ...
                    ['Running synchronously; the GUI will resume when ' ...
                    'the backend function returns.']);
                report = freqnessgui.runSecondaryAnalysis( ...
                    app.NetworkSet,configuration,mniCoordinates, ...
                    @app.updateSecondaryProgress);
                clear runningCleanup
                app.showSecondaryReport(report);
            catch exception
                clear runningCleanup
                app.showSecondaryError(exception);
            end
        end

        function startSecondaryBackgroundAnalysis( ...
                app,configuration,mniCoordinates)
            app.SecondaryCancelRequested = false;
            app.SecondaryRunConfiguration = configuration;
            try
                app.SecondaryCancellationFile = ...
                    freqnessgui.prepareSecondaryCancellation( ...
                    configuration.outputFolder);
                app.setSecondaryRunningState(true,true);
                [process,control] = ...
                    freqnessgui.launchSecondaryAnalysisProcess( ...
                    app.NetworkSet,configuration,mniCoordinates, ...
                    app.SecondaryCancellationFile);
                app.SecondaryProcess = process;
                app.SecondaryProcessControl = control;
                pollTimer = timer( ...
                    'Name','FREQNESS secondary-analysis monitor', ...
                    'ExecutionMode','fixedSpacing', ...
                    'Period',0.25, ...
                    'BusyMode','drop', ...
                    'TimerFcn',@(~,~)app.pollSecondaryProcess());
                app.SecondaryProcessTimer = pollTimer;
                app.SecondaryLastProgressSequence = 0;
                app.updateSecondaryProgress(0, ...
                    ['Background MATLAB process started. The GUI remains responsive; ' ...
                    'press Cancel Analysis to stop safely.']);
                start(pollTimer);
            catch exception
                if app.isSecondaryProcessAlive()
                    try
                        app.SecondaryProcess.destroyForcibly();
                        app.SecondaryProcess.waitFor();
                    catch
                        % Continue cleanup and report the startup error.
                    end
                end
                app.cleanupSecondaryProcessState();
                app.setSecondaryRunningState(false);
                app.showSecondaryError(exception);
            end
        end

        function requestSecondaryCancellation(app)
            if ~app.IsRunning || ~app.SecondaryRunInBackground || ...
                    app.SecondaryCancelRequested
                return
            end
            app.SecondaryCancelRequested = true;
            app.SecondaryRunButton.Text = 'Cancelling...';
            app.SecondaryRunButton.Enable = 'off';
            app.updateSecondaryProgress(1, ...
                ['Cancellation requested. Completed atomic outputs will be ' ...
                'retained; pending work will be marked cancelled.']);
            try
                freqnessgui.requestSecondaryCancellation( ...
                    app.SecondaryCancellationFile);
            catch exception
                app.updateSecondaryProgress(1,sprintf( ...
                    'Could not create the cancellation marker: %s', ...
                    exception.message));
            end
            app.SecondaryCancelTic = tic;
        end

        function pollSecondaryProcess(app)
            if isempty(app.SecondaryProcess)
                return
            end
            app.relaySecondaryProcessProgress();
            processIsAlive = app.isSecondaryProcessAlive();
            if app.SecondaryCancelRequested && processIsAlive && ...
                    ~isempty(app.SecondaryCancelTic)
                cancellationWait = toc(app.SecondaryCancelTic);
                if cancellationWait > 2
                    try
                        app.SecondaryProcess.destroyForcibly();
                    catch
                        % Poll again; the process may have just stopped.
                    end
                elseif cancellationWait > 0.75
                    try
                        app.SecondaryProcess.destroy();
                    catch
                        % Poll again; the process may have just stopped.
                    end
                end
                processIsAlive = app.isSecondaryProcessAlive();
            end
            if processIsAlive
                return
            end
            cancellationRequested = app.SecondaryCancelRequested;
            configuration = app.SecondaryRunConfiguration;
            control = app.SecondaryProcessControl;
            workerResult = struct();
            if isstruct(control) && isfield(control,'resultFile') && ...
                    isfile(control.resultFile)
                try
                    loaded = load(control.resultFile,'workerResult');
                    workerResult = loaded.workerResult;
                catch
                    workerResult = struct();
                end
            end
            if cancellationRequested && ...
                    (~isfield(workerResult,'status') || ...
                    ~strcmp(workerResult.status,'completed'))
                try
                    manifest = freqnessgui.markSecondaryAnalysisCancelled( ...
                        configuration, ...
                        'Cancelled by the user from the FREQ-NESS GUI.');
                catch manifestException
                    app.cleanupSecondaryProcessState();
                    app.setSecondaryRunningState(false);
                    app.showSecondaryError(manifestException);
                    return
                end
                app.cleanupSecondaryProcessState();
                app.setSecondaryRunningState(false);
                app.showSecondaryCancellation(manifest);
                return
            end
            if ~isfield(workerResult,'status')
                message = app.secondaryProcessFailureMessage(control);
                app.cleanupSecondaryProcessState();
                app.setSecondaryRunningState(false);
                app.showSecondaryError(MException( ...
                    'FREQNESS:GUI:BackgroundProcessFailed','%s',message));
                return
            end
            switch workerResult.status
                case {'completed','cancelled'}
                    report = workerResult.report;
                    app.cleanupSecondaryProcessState();
                    app.setSecondaryRunningState(false);
                    app.showSecondaryReport(report);
                otherwise
                    errorInfo = workerResult.error;
                    app.cleanupSecondaryProcessState();
                    app.setSecondaryRunningState(false);
                    app.showSecondaryError(MException( ...
                        errorInfo.identifier,'%s',errorInfo.message));
            end
        end

        function relaySecondaryProcessProgress(app)
            control = app.SecondaryProcessControl;
            if ~isstruct(control) || ~isfield(control,'progressFile') || ...
                    ~isfile(control.progressFile)
                return
            end
            try
                loaded = load(control.progressFile,'progress');
                progress = loaded.progress;
            catch
                return
            end
            if ~isstruct(progress) || ~isfield(progress,'sequence') || ...
                    progress.sequence <= app.SecondaryLastProgressSequence
                return
            end
            app.SecondaryLastProgressSequence = progress.sequence;
            app.updateSecondaryProgress(progress.fraction,progress.message);
        end

        function isAlive = isSecondaryProcessAlive(app)
            isAlive = false;
            if isempty(app.SecondaryProcess)
                return
            end
            try
                isAlive = logical(app.SecondaryProcess.isAlive());
            catch
                % A process handle that can no longer be queried is finished.
            end
        end

        function message = secondaryProcessFailureMessage(~,control)
            message = ['The background MATLAB process stopped before it ' ...
                'reported a result.'];
            if ~isstruct(control) || ~isfield(control,'logFile') || ...
                    ~isfile(control.logFile)
                return
            end
            try
                logText = strtrim(fileread(control.logFile));
            catch
                return
            end
            if isempty(logText)
                return
            end
            maxCharacters = 1200;
            if numel(logText) > maxCharacters
                logText = logText(end-maxCharacters+1:end);
            end
            message = sprintf('%s\n\nWorker log:\n%s',message,logText);
        end

        function showSecondaryReport(app,report)
            if isfield(report,'status') && strcmp(report.status,'cancelled')
                loaded = load(report.manifestFile,'manifest');
                app.showSecondaryCancellation(loaded.manifest);
                return
            end
            summary = sprintf( ...
                'Finished: %d participant output(s), %d figure(s).', ...
                report.nCompleted,numel(report.figureFiles));
            app.SecondaryValidationLabel.Text = summary;
            app.SecondaryValidationLabel.FontColor = [0.125 0.545 0.365];
            app.FooterStatusLabel.Text = summary;
        end

        function showSecondaryCancellation(app,manifest)
            statuses = {manifest.participants.status};
            nCompleted = sum(strcmp(statuses,'completed'));
            nCancelled = sum(strcmp(statuses,'cancelled'));
            summary = sprintf( ...
                'Cancelled: %d completed output(s) retained; %d cancelled.', ...
                nCompleted,nCancelled);
            app.SecondaryValidationLabel.Text = summary;
            app.SecondaryValidationLabel.FontColor = [0.72 0.42 0.05];
            app.FooterStatusLabel.Text = summary;
            app.updateSecondaryProgress(1,summary);
            app.SecondaryValidationLabel.Text = summary;
            app.SecondaryValidationLabel.FontColor = [0.72 0.42 0.05];
        end

        function showSecondaryError(app,exception)
            app.SecondaryValidationLabel.Text = exception.message;
            app.SecondaryValidationLabel.FontColor = [0.78 0.24 0.16];
            app.FooterStatusLabel.Text = exception.message;
            uialert(app.UIFigure,exception.message, ...
                'Secondary-analysis error');
        end

        function updateSecondaryProgress(app,fraction,message)
            percentage = max(0,min(100,round(100*fraction)));
            app.SecondaryValidationLabel.Text = sprintf( ...
                '%s  (%d%%)',message,percentage);
            app.SecondaryValidationLabel.FontColor = [0.055 0.415 0.690];
            app.FooterStatusLabel.Text = message;
            currentLines = app.SecondaryProgressTextArea.Value;
            if ischar(currentLines) || isstring(currentLines)
                currentLines = cellstr(currentLines);
            end
            timeText = char(datetime('now','Format','HH:mm:ss'));
            currentLines{end+1,1} = sprintf('[%s] %s',timeText,message);
            maxVisibleLines = 12;
            if numel(currentLines) > maxVisibleLines
                currentLines = currentLines(end-maxVisibleLines+1:end);
            end
            app.SecondaryProgressTextArea.Value = currentLines;
            drawnow limitrate
        end

        function setSecondaryRunningState(app,isRunning,isBackground)
            if nargin < 3
                isBackground = false;
            end
            app.IsRunning = isRunning;
            app.SecondaryRunInBackground = isRunning && isBackground;
            if isRunning
                app.SecondaryBackgroundCheckBox.Enable = 'off';
                if isBackground
                    app.SecondaryRunButton.Enable = 'on';
                    app.SecondaryRunButton.Text = 'Cancel Analysis';
                    app.SecondaryRunButton.BackgroundColor = [0.78 0.24 0.16];
                    app.SecondaryRunButton.Tooltip = ...
                        'Cancel the background analysis and retain completed outputs.';
                else
                    app.SecondaryRunButton.Enable = 'off';
                    app.SecondaryRunButton.Text = 'Running...';
                end
                app.SecondaryValidateButton.Enable = 'off';
                app.SecondaryTree.Enable = 'off';
                app.SecondaryParticipantButton.Enable = 'off';
                app.NetworkBrowseButton.Enable = 'off';
                app.DatasetBrowseButton.Enable = 'off';
                app.MNIBrowseButton.Enable = 'off';
                app.RunButton.Enable = 'off';
                app.SecondaryProgressTextArea.Value = { ...
                    'Starting secondary analysis...'};
            else
                app.SecondaryRunButton.Text = 'Run Analysis';
                app.SecondaryRunButton.BackgroundColor = [0.055 0.415 0.690];
                app.SecondaryRunButton.Tooltip = ...
                    'Run the selected analysis and persist its outputs.';
                if freqnessgui.secondaryBackgroundAvailable()
                    app.SecondaryBackgroundCheckBox.Enable = 'on';
                end
                app.SecondaryTree.Enable = 'on';
                app.NetworkBrowseButton.Enable = 'on';
                app.DatasetBrowseButton.Enable = 'on';
                app.MNIBrowseButton.Enable = 'on';
                if ~isempty(app.NetworkSet)
                    app.SecondaryParticipantButton.Enable = 'on';
                end
                if ~isempty(app.Dataset) && app.Dataset.isValid
                    app.RunButton.Enable = 'on';
                end
                app.updateSecondaryReadiness();
            end
        end

        function cleanupSecondaryProcessState(app)
            if ~isempty(app.SecondaryProcessTimer) && ...
                    isvalid(app.SecondaryProcessTimer)
                try
                    stop(app.SecondaryProcessTimer);
                catch
                    % The timer may already have stopped itself.
                end
                delete(app.SecondaryProcessTimer);
            end
            app.SecondaryProcessTimer = [];
            app.SecondaryProcess = [];
            if ~isempty(app.SecondaryCancellationFile) && ...
                    isfile(app.SecondaryCancellationFile)
                delete(app.SecondaryCancellationFile);
            end
            if isstruct(app.SecondaryProcessControl) && ...
                    isfield(app.SecondaryProcessControl,'controlFolder') && ...
                    isfolder(app.SecondaryProcessControl.controlFolder)
                try
                    rmdir(app.SecondaryProcessControl.controlFolder,'s');
                catch
                    % Temporary control files can be removed on the next run.
                end
            end
            app.SecondaryProcessControl = [];
            app.SecondaryLastProgressSequence = 0;
            app.SecondaryCancellationFile = '';
            app.SecondaryRunConfiguration = [];
            app.SecondaryCancelRequested = false;
            app.SecondaryCancelTic = [];
        end

        function shutdownSecondaryExecution(app)
            configuration = app.SecondaryRunConfiguration;
            if app.isSecondaryProcessAlive()
                try
                    freqnessgui.requestSecondaryCancellation( ...
                        app.SecondaryCancellationFile);
                catch
                    % Continue with direct process cancellation.
                end
                try
                    app.SecondaryProcess.destroyForcibly();
                    app.SecondaryProcess.waitFor();
                catch
                    % Application shutdown must continue.
                end
                if ~isempty(configuration)
                    try
                        freqnessgui.markSecondaryAnalysisCancelled( ...
                            configuration, ...
                            'Cancelled because the FREQ-NESS GUI was closed.');
                    catch
                        % Do not block application shutdown on manifest I/O.
                    end
                end
            end
            app.cleanupSecondaryProcessState();
        end

        function field = addTextSetting(~,parent,row,labelText,defaultValue,tooltipText)
            label = uilabel(parent,'Text',labelText);
            label.Layout.Row = row;
            label.Layout.Column = 1;
            field = uieditfield(parent,'text', ...
                'Value',defaultValue, ...
                'Tooltip',tooltipText);
            field.Layout.Row = row;
            field.Layout.Column = 2;
        end

        function registerFileDropTargets(app)
            try
                uiFileDnD(app.DatasetDropLabel, ...
                    @(~,dropData)app.handleDatasetDrop(dropData));
                uiFileDnD(app.MNIDropLabel, ...
                    @(~,dropData)app.handleMNIDrop(dropData));
                uiFileDnD(app.NetworkDropLabel, ...
                    @(~,dropData)app.handleNetworkDrop(dropData));
            catch exception
                app.FooterStatusLabel.Text = sprintf( ...
                    'Drag-and-drop unavailable (%s). Browse buttons remain available.', ...
                    exception.message);
            end
        end

        function handleDatasetDrop(app,dropData)
            names = dropData.names;
            if ischar(names) || isstring(names)
                names = cellstr(names);
            end
            folder = '';
            if isscalar(names) && isfolder(names{1})
                folder = names{1};
            elseif ~isempty(names) && all(cellfun(@isfile,names))
                parentFolders = cell(size(names));
                validExtensions = false(size(names));
                for namei = 1:numel(names)
                    [parentFolders{namei},~,extension] = fileparts(names{namei});
                    validExtensions(namei) = strcmpi(extension,'.mat');
                end
                if all(validExtensions) && ...
                        isscalar(unique(parentFolders))
                    folder = parentFolders{1};
                end
            end
            if isempty(folder)
                uialert(app.UIFigure, ...
                    ['Drop one Dataset* folder, or MAT-files that all ' ...
                    'belong to the same Dataset* folder.'], ...
                    'Invalid dataset drop');
                return
            end
            if ~isempty(folder)
                app.importDataset(folder);
            end
        end

        function handleNetworkDrop(app,dropData)
            folder = app.singleDroppedFolder(dropData,'FREQ_Networks*');
            if ~isempty(folder)
                app.importNetworkFolder(folder);
            end
        end

        function handleMNIDrop(app,dropData)
            names = dropData.names;
            if ischar(names) || isstring(names)
                names = cellstr(names);
            end
            if isscalar(names)
                [~,~,extension] = fileparts(names{1});
            else
                extension = '';
            end
            if numel(names) ~= 1 || ~isfile(names{1}) || ...
                    ~strcmpi(extension,'.mat')
                uialert(app.UIFigure, ...
                    'Drop exactly one MNI coordinate MAT-file.', ...
                    'Invalid MNI drop');
                return
            end
            app.importMNIFile(names{1});
        end

        function folder = singleDroppedFolder(app,dropData,expectedName)
            folder = '';
            names = dropData.names;
            if ischar(names) || isstring(names)
                names = cellstr(names);
            end
            if numel(names) ~= 1 || ~isfolder(names{1})
                uialert(app.UIFigure, ...
                    sprintf('Drop exactly one %s folder.',expectedName), ...
                    'Invalid drop');
                return
            end
            folder = names{1};
        end

        function browseDataset(app)
            selectedFolder = uigetdir(pwd,'Select a Dataset* folder');
            if isequal(selectedFolder,0)
                return
            end
            app.importDataset(selectedFolder);
        end

        function browseMNIFile(app)
            [selectedFile,selectedFolder] = uigetfile( ...
                {'*.mat','MAT-files (*.mat)'}, ...
                'Select an MNI coordinate file');
            if isequal(selectedFile,0)
                return
            end
            app.importMNIFile(fullfile(selectedFolder,selectedFile));
        end

        function loadDefaultMNI(app)
            appFile = which('FREQNESSApp');
            toolboxFolder = fileparts(fileparts(appFile));
            defaultFile = freqnessgui.findDefaultMNI(toolboxFolder);
            if isempty(defaultFile)
                app.Config.mniFile = '';
                app.showEmptyMNIPreview();
                return
            end
            app.importMNIFile(defaultFile,false);
        end

        function importMNIFile(app,mniFile,showAlert)
            if nargin < 3
                showAlert = true;
            end
            try
                mni = freqnessgui.loadMNICoordinates(mniFile);
                app.MNI = mni;
                app.Config.mniFile = mni.path;
                app.MNIPathField.Value = mni.path;
                app.MNIPathField.Tooltip = mni.path;
                app.plotMNIPreview();
                app.updateMNICompatibility();
                app.FooterStatusLabel.Text = sprintf( ...
                    'Loaded %s with %d MNI coordinates.',mni.name,mni.nPoints);
            catch exception
                if showAlert
                    app.FooterStatusLabel.Text = 'MNI coordinate import failed.';
                    uialert(app.UIFigure,exception.message,'MNI import error');
                else
                    app.MNI = [];
                    app.Config.mniFile = '';
                    app.MNIPathField.Value = '';
                    app.showEmptyMNIPreview();
                end
            end
        end

        function showEmptyMNIPreview(app)
            cla(app.MNIAxes);
            axis(app.MNIAxes,'off');
            text(app.MNIAxes,0.5,0.5,'No MNI coordinates loaded', ...
                'Units','normalized', ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'Color',[0.36 0.41 0.47]);
            app.MNISummaryLabel.Text = 'No MNI coordinates loaded';
            app.MNISummaryLabel.FontColor = [0.36 0.41 0.47];
            app.updateSecondaryContext();
        end

        function plotMNIPreview(app)
            coordinates = app.MNI.coordinates;
            cla(app.MNIAxes);
            scatter3(app.MNIAxes,coordinates(:,1),coordinates(:,2), ...
                coordinates(:,3),3,'k','filled', ...
                'MarkerFaceAlpha',0.28, ...
                'MarkerEdgeColor','none');
            axis(app.MNIAxes,'equal');
            axis(app.MNIAxes,'vis3d');
            grid(app.MNIAxes,'on');
            view(app.MNIAxes,3);
            xlabel(app.MNIAxes,'X');
            ylabel(app.MNIAxes,'Y');
            zlabel(app.MNIAxes,'Z');
            app.MNIAxes.FontSize = 8;
            app.MNIAxes.Interactions = rotateInteraction;
        end

        function updateMNICompatibility(app)
            if isempty(app.MNI)
                app.updateSecondaryContext();
                return
            end

            if isempty(app.Dataset)
                app.MNISummaryLabel.Text = sprintf( ...
                    '%d coordinates — drag the plot to rotate', ...
                    app.MNI.nPoints);
                app.MNISummaryLabel.FontColor = [0.055 0.415 0.690];
                app.updateSecondaryContext();
                return
            end

            voxelCounts = arrayfun( ...
                @(participant)participant.dataSize(1), ...
                app.Dataset.participants);
            if all(voxelCounts == app.MNI.nPoints)
                app.MNISummaryLabel.Text = sprintf( ...
                    '%d coordinates — matches participant voxels', ...
                    app.MNI.nPoints);
                app.MNISummaryLabel.FontColor = [0.125 0.545 0.365];
            else
                app.MNISummaryLabel.Text = sprintf( ...
                    '%d coordinates — participant voxel count differs', ...
                    app.MNI.nPoints);
                app.MNISummaryLabel.FontColor = [0.78 0.38 0.05];
            end
            app.updateSecondaryContext();
        end

        function importDataset(app,datasetFolder)
            try
                dataset = freqnessgui.inspectDataset(datasetFolder);
                app.Dataset = dataset;
                app.Config.datasetFolder = dataset.folder;
                app.Config.outputFolder = dataset.networkFolder;

                app.DatasetPathField.Value = dataset.folder;
                app.DatasetPathField.Tooltip = dataset.folder;
                app.NetworkOutputField.Value = dataset.networkFolder;
                app.NetworkOutputField.Tooltip = dataset.networkFolder;

                tableData = cell(dataset.nParticipants,4);
                for participanti = 1:dataset.nParticipants
                    participant = dataset.participants(participanti);
                    if numel(participant.dataSize) >= 2
                        sizeText = sprintf('%d x %d', ...
                            participant.dataSize(1),participant.dataSize(2));
                    else
                        sizeText = '-';
                    end
                    tableData(participanti,:) = {participant.id, ...
                        participant.variableName,sizeText,participant.message};
                end
                app.DatasetTable.Data = tableData;
                app.DatasetSummaryLabel.Text = sprintf('%d/%d participants ready', ...
                    dataset.nValid,dataset.nParticipants);
                app.updateMNICompatibility();

                if dataset.isValid
                    app.RunButton.Enable = 'on';
                    app.RunStatusLabel.Text = ...
                        'Dataset validated. Review settings, then run.';
                    app.FooterStatusLabel.Text = sprintf( ...
                        'Imported %s with %d participants.', ...
                        dataset.name,dataset.nParticipants);
                else
                    app.RunButton.Enable = 'off';
                    app.RunStatusLabel.Text = ...
                        'Resolve invalid participant files before running.';
                    app.FooterStatusLabel.Text = 'Dataset validation failed.';
                end
            catch exception
                app.RunButton.Enable = 'off';
                app.FooterStatusLabel.Text = 'Dataset import failed.';
                uialert(app.UIFigure,exception.message,'Dataset import error');
            end
        end

        function browseNetworkFolder(app)
            selectedFolder = uigetdir(pwd,'Select a FREQ_Networks* folder');
            if isequal(selectedFolder,0)
                return
            end
            app.importNetworkFolder(selectedFolder);
        end

        function importNetworkFolder(app,networkFolder)
            try
                networkSet = freqnessgui.inspectNetworkFolder(networkFolder);
                app.NetworkSet = networkSet;
                app.SecondarySelectedParticipants = networkSet.participantIds(:)';
                app.NetworkPathField.Value = networkSet.folder;
                app.NetworkPathField.Tooltip = networkSet.folder;
                app.AnalysisOutputField.Value = networkSet.analysisFolder;
                app.AnalysisOutputField.Tooltip = networkSet.analysisFolder;
                if ~isempty(networkSet.mniFile) && isfile(networkSet.mniFile) && ...
                        (isempty(app.MNI) || ...
                        ~strcmp(app.MNI.path,networkSet.mniFile))
                    app.importMNIFile(networkSet.mniFile,false);
                end
                app.updateSecondaryContext();
                if ~isempty(app.SecondarySelectedModuleId)
                    app.selectSecondaryModule(app.SecondarySelectedModuleId);
                end
                app.FooterStatusLabel.Text = sprintf( ...
                    'Imported %s for secondary analyses.',networkSet.name);
            catch exception
                app.FooterStatusLabel.Text = 'Network-result import failed.';
                uialert(app.UIFigure,exception.message,'Network import error');
            end
        end

        function config = collectConfig(app)
            config = app.Config;
            config.datasetFolder = app.Dataset.folder;
            config.outputFolder = app.Dataset.networkFolder;
            config.mniFile = app.Config.mniFile;
            config.network.frequencies = app.currentFrequencies();
            config.network.samplingRate = app.SamplingRateField.Value;
            config.network.duration = freqnessgui.parseNumericVector( ...
                app.DurationField.Value,true,'Duration');
            config.network.fwidth = freqnessgui.parseNumericVector( ...
                app.FilterWidthField.Value,true,'Filter width');
            config.network.filter = app.FilterTypeDropDown.Value;
            config.network.regularisation = app.RegularisationField.Value;
            config.network.ncomps = app.ComponentsField.Value;
            config.network.badSegments = freqnessgui.parseNumericVector( ...
                app.BadSegmentsField.Value,true,'Bad sample indices');
            config.network.rescale = app.RescaleCheckBox.Value;
            config.execution.recomputeExisting = app.RecomputeCheckBox.Value;
        end

        function ceiling = frequencySliderCeiling(app,requestedMaximum,step)
            if isempty(app.SamplingRateField)
                samplingRate = app.Config.network.samplingRate;
            else
                samplingRate = app.SamplingRateField.Value;
            end
            nyquist = samplingRate/2;
            tolerance = max(16*eps(nyquist),nyquist*1e-12);
            maximumGridValue = floor((nyquist-tolerance)/step)*step;
            if maximumGridValue < step
                error('FREQNESS:GUI:InvalidFrequencyStep', ...
                    ['Frequency resolution must be smaller than the ' ...
                    'Nyquist frequency (%.4g Hz).'],nyquist);
            end

            preferredCeiling = max(requestedMaximum,min(40,maximumGridValue));
            ceiling = ceil(preferredCeiling/step)*step;
            ceiling = min(maximumGridValue,max(step,ceiling));
        end

        function commitFrequencyFields(app)
            previousStep = app.FrequencyRangeSlider.Step;
            previousRange = app.FrequencyRangeSlider.Value;
            try
                app.applyFrequencyRange([app.FrequencyMinField.Value, ...
                    app.FrequencyMaxField.Value]);
            catch exception
                app.FrequencyStepField.Value = previousStep;
                app.FrequencyMinField.Value = previousRange(1);
                app.FrequencyMaxField.Value = previousRange(2);
                uialert(app.UIFigure,exception.message, ...
                    'Invalid frequency settings');
            end
        end

        function commitFrequencyRange(app,requestedRange)
            try
                app.applyFrequencyRange(requestedRange);
            catch exception
                uialert(app.UIFigure,exception.message, ...
                    'Invalid frequency settings');
            end
        end

        function applyFrequencyRange(app,requestedRange)
            step = app.FrequencyStepField.Value;
            [~,snappedRange] = freqnessgui.buildFrequencyVector( ...
                requestedRange,step);
            ceiling = app.frequencySliderCeiling(snappedRange(2),step);
            snappedRange = min(snappedRange,ceiling);
            snappedRange(1) = min(snappedRange(1),snappedRange(2));
            [frequencies,snappedRange] = freqnessgui.buildFrequencyVector( ...
                snappedRange,step);

            app.FrequencyRangeSlider.Limits = [0 ceiling];
            app.FrequencyRangeSlider.Step = step;
            app.FrequencyRangeSlider.Value = snappedRange;
            app.FrequencyMinField.Value = snappedRange(1);
            app.FrequencyMaxField.Value = snappedRange(2);
            app.Config.network.frequencies = frequencies;
            app.renderFrequencyPreview(frequencies,snappedRange,step);
        end

        function previewFrequencyRange(app,requestedRange)
            try
                [frequencies,snappedRange] = ...
                    freqnessgui.buildFrequencyVector(requestedRange, ...
                    app.FrequencyStepField.Value);
                app.renderFrequencyPreview(frequencies,snappedRange, ...
                    app.FrequencyStepField.Value);
            catch
                % The committed slider value remains the source of truth.
            end
        end

        function frequencies = currentFrequencies(app)
            [frequencies,~] = freqnessgui.buildFrequencyVector( ...
                [app.FrequencyMinField.Value,app.FrequencyMaxField.Value], ...
                app.FrequencyStepField.Value);
            if any(frequencies >= app.SamplingRateField.Value/2)
                error('FREQNESS:GUI:InvalidFrequencyRange', ...
                    'Every frequency must be below the Nyquist frequency.');
            end
        end

        function renderFrequencyPreview(app,frequencies,frequencyRange,step)
            app.FrequencySummaryLabel.Text = sprintf( ...
                'Filter bank  |  %.4g–%.4g Hz  |  Δf %.4g Hz  |  %d Gaussian filters', ...
                frequencyRange(1),frequencyRange(2),step,numel(frequencies));
            app.updateFWHMPreview(frequencies);
        end

        function updateFWHMPreview(app,frequencies)
            if nargin < 2
                try
                    frequencies = app.currentFrequencies();
                catch exception
                    app.showInvalidFWHMPreview(exception.message);
                    return
                end
            end

            try
                fwidth = freqnessgui.parseNumericVector( ...
                    app.FilterWidthField.Value,true,'Filter width');
                fwhm = freqnessgui.computeFWHM(frequencies,fwidth, ...
                    app.FilterTypeDropDown.Value);
                plotLimits = app.FrequencyRangeSlider.Limits;
                [frequencyAxis,responses] = ...
                    freqnessgui.computeFilterResponses( ...
                    frequencies,fwhm,plotLimits);
                lowColor = [0.48 0.80 0.91];
                highColor = [0.035 0.18 0.42];
                if isscalar(frequencies)
                    filterColors = [0.055 0.415 0.690];
                else
                    colorPosition = linspace(0,1,numel(frequencies))';
                    filterColors = (1-colorPosition).*lowColor+ ...
                        colorPosition.*highColor;
                end

                cla(app.FWHMAxes);
                hold(app.FWHMAxes,'on');
                patch(app.FWHMAxes, ...
                    [frequencies(1) frequencies(end) ...
                    frequencies(end) frequencies(1)], ...
                    [0 0 1.05 1.05],[0.055 0.415 0.690], ...
                    'FaceAlpha',0.055,'EdgeColor','none');
                for frequencyi = 1:numel(frequencies)
                    plot(app.FWHMAxes,frequencyAxis,responses(frequencyi,:), ...
                        'Color',filterColors(frequencyi,:), ...
                        'LineWidth',1.15);
                end
                scatter(app.FWHMAxes,frequencies,ones(size(frequencies)), ...
                    13,filterColors,'filled');
                yline(app.FWHMAxes,0.5,':', ...
                    'Color',[0.48 0.52 0.57], ...
                    'LineWidth',0.8);
                text(app.FWHMAxes,0.985,0.93, ...
                    'light → dark = low → high', ...
                    'Units','normalized', ...
                    'HorizontalAlignment','right', ...
                    'VerticalAlignment','top', ...
                    'FontSize',8, ...
                    'Color',[0.36 0.41 0.47], ...
                    'BackgroundColor',[1 1 1], ...
                    'Margin',1);
                hold(app.FWHMAxes,'off');
                axis(app.FWHMAxes,'on');
                app.FWHMAxes.Box = 'on';
                app.FWHMAxes.XGrid = 'on';
                app.FWHMAxes.YGrid = 'on';
                xlabel(app.FWHMAxes,'Frequency (Hz)');
                ylabel(app.FWHMAxes,'Normalized gain');
                xlim(app.FWHMAxes,plotLimits);
                ylim(app.FWHMAxes,[0 1.05]);
                app.FWHMAxes.YTick = [0 0.5 1];
            catch exception
                app.showInvalidFWHMPreview(exception.message);
            end
        end

        function showInvalidFWHMPreview(app,message)
            cla(app.FWHMAxes);
            axis(app.FWHMAxes,'off');
            text(app.FWHMAxes,0.5,0.5,message, ...
                'Units','normalized', ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'Color',[0.78 0.24 0.16], ...
                'Interpreter','none');
        end

        function samplingRateChanged(app)
            previousSamplingRate = app.Config.network.samplingRate;
            app.Config.network.samplingRate = app.SamplingRateField.Value;
            try
                app.applyFrequencyRange([app.FrequencyMinField.Value, ...
                    app.FrequencyMaxField.Value]);
            catch exception
                app.Config.network.samplingRate = previousSamplingRate;
                app.SamplingRateField.Value = previousSamplingRate;
                app.applyFrequencyRange([app.FrequencyMinField.Value, ...
                    app.FrequencyMaxField.Value]);
                uialert(app.UIFigure,exception.message, ...
                    'Invalid sampling rate');
            end
        end

        function runCoreAnalysis(app)
            if isempty(app.Dataset) || app.IsRunning
                return
            end

            try
                config = app.collectConfig();
                freqnessgui.validateNetworkConfig(config,app.Dataset);
            catch exception
                uialert(app.UIFigure,exception.message,'Invalid analysis settings');
                return
            end

            app.Config = config;
            app.setRunningState(true);
            runningCleanup = onCleanup(@()app.setRunningState(false));

            try
                report = freqnessgui.runNetworkEstimation( ...
                    app.Dataset,config,@app.updateProgress);
                app.importNetworkFolder(report.outputFolder);

                if report.nFailed == 0
                    summary = sprintf( ...
                        'Finished: %d completed, %d retained.', ...
                        report.nCompleted,report.nSkipped);
                else
                    summary = sprintf( ...
                        'Finished with failures: %d completed, %d retained, %d failed.', ...
                        report.nCompleted,report.nSkipped,report.nFailed);
                end
                app.RunStatusLabel.Text = summary;
                app.FooterStatusLabel.Text = summary;
            catch exception
                app.RunStatusLabel.Text = 'Network estimation stopped.';
                app.FooterStatusLabel.Text = exception.message;
                uialert(app.UIFigure,exception.message,'Network estimation error');
            end

            clear runningCleanup
        end

        function updateProgress(app,fraction,message)
            percentage = max(0,min(100,round(100*fraction)));
            app.RunStatusLabel.Text = sprintf('%s  (%d%%)',message,percentage);
            app.FooterStatusLabel.Text = message;
            currentLines = app.ProgressTextArea.Value;
            if ischar(currentLines) || isstring(currentLines)
                currentLines = cellstr(currentLines);
            end
            timeText = char(datetime('now','Format','HH:mm:ss'));
            currentLines{end+1,1} = sprintf('[%s] %s',timeText,message);
            maxVisibleLines = 12;
            if numel(currentLines) > maxVisibleLines
                currentLines = currentLines(end-maxVisibleLines+1:end);
            end
            app.ProgressTextArea.Value = currentLines;
            drawnow limitrate
        end

        function setRunningState(app,isRunning)
            app.IsRunning = isRunning;
            if isRunning
                app.RunButton.Enable = 'off';
                app.RunButton.Text = 'Running...';
                app.DatasetBrowseButton.Enable = 'off';
                app.MNIBrowseButton.Enable = 'off';
                app.NetworkBrowseButton.Enable = 'off';
                app.ProgressTextArea.Value = {'Starting network estimation...'};
            else
                app.RunButton.Text = 'Run Network Estimation';
                app.DatasetBrowseButton.Enable = 'on';
                app.MNIBrowseButton.Enable = 'on';
                app.NetworkBrowseButton.Enable = 'on';
                if ~isempty(app.Dataset) && app.Dataset.isValid
                    app.RunButton.Enable = 'on';
                else
                    app.RunButton.Enable = 'off';
                end
            end
            drawnow
        end
    end
end
