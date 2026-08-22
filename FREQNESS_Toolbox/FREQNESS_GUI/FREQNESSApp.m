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
        MainGrid

        DatasetDropLabel
        DatasetPathField
        DatasetBrowseButton
        DatasetSummaryLabel
        DatasetTable

        FrequenciesField
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
        RunStatusLabel
        ProgressTextArea

        NetworkDropLabel
        NetworkPathField
        NetworkBrowseButton
        AnalysisOutputField
        SecondarySummaryLabel
        ModuleTable

        FooterStatusLabel
    end

    methods
        function app = FREQNESSApp(varargin)
            app.Config = freqnessgui.defaultConfig();
            app.createInterface();
            setappdata(app.UIFigure,'FREQNESSApp',app);
            app.registerFileDropTargets();
            drawnow
            try
                scroll(app.MainGrid,'top');
            catch
                % Older releases retain their default top scroll position.
            end

            if nargin >= 1 && ~isempty(varargin{1})
                app.importDataset(varargin{1});
            end
        end

        function delete(app)
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

            app.MainGrid = uigridlayout(app.UIFigure,[5 1]);
            app.MainGrid.RowHeight = {72,168,398,190,26};
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
            app.createSecondaryLayer(mainGrid,colors);

            app.FooterStatusLabel = uilabel(mainGrid, ...
                'Text','Ready', ...
                'FontSize',11, ...
                'FontColor',colors.muted);
            app.FooterStatusLabel.Layout.Row = 5;
        end

        function createDataLayer(app,parent,colors)
            panel = uipanel(parent, ...
                'Title','1  Data import', ...
                'FontSize',14, ...
                'FontWeight','bold', ...
                'ForegroundColor',colors.navy, ...
                'BackgroundColor',[1 1 1]);
            panel.Layout.Row = 2;

            grid = uigridlayout(panel,[2 4]);
            grid.RowHeight = {58,'1x'};
            grid.ColumnWidth = {220,'1x',110,260};
            grid.Padding = [12 8 12 10];
            grid.RowSpacing = 8;
            grid.ColumnSpacing = 10;

            app.DatasetDropLabel = uilabel(grid, ...
                'Text',sprintf('Drop a Dataset* folder here\nor use Browse'), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','center', ...
                'FontWeight','bold', ...
                'FontColor',colors.blue, ...
                'BackgroundColor',colors.softBlue);
            app.DatasetDropLabel.Layout.Row = [1 2];
            app.DatasetDropLabel.Layout.Column = 1;

            app.DatasetPathField = uieditfield(grid,'text', ...
                'Editable','off', ...
                'Placeholder','No Dataset* folder selected');
            app.DatasetPathField.Layout.Row = 1;
            app.DatasetPathField.Layout.Column = 2;

            app.DatasetBrowseButton = uibutton(grid,'push', ...
                'Text','Browse...', ...
                'ButtonPushedFcn',@(~,~)app.browseDataset());
            app.DatasetBrowseButton.Layout.Row = 1;
            app.DatasetBrowseButton.Layout.Column = 3;

            app.DatasetSummaryLabel = uilabel(grid, ...
                'Text','Waiting for a Dataset* folder', ...
                'FontColor',colors.muted, ...
                'HorizontalAlignment','right');
            app.DatasetSummaryLabel.Layout.Row = 1;
            app.DatasetSummaryLabel.Layout.Column = 4;

            app.DatasetTable = uitable(grid, ...
                'Data',cell(0,4), ...
                'ColumnName',{'Participant','Variable','Matrix size','Validation'}, ...
                'ColumnWidth',{170,150,110,'auto'}, ...
                'RowName',{});
            app.DatasetTable.Layout.Row = 2;
            app.DatasetTable.Layout.Column = [2 4];
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
            grid.ColumnWidth = {310,390,'1x'};
            grid.Padding = [12 10 12 12];
            grid.ColumnSpacing = 12;

            mandatoryPanel = uipanel(grid, ...
                'Title','Mandatory inputs', ...
                'FontWeight','bold', ...
                'BackgroundColor',[1 1 1]);
            mandatoryPanel.Layout.Column = 1;
            mandatoryGrid = uigridlayout(mandatoryPanel,[6 1]);
            mandatoryGrid.RowHeight = {20,32,20,32,50,'1x'};
            mandatoryGrid.Padding = [10 8 10 8];
            mandatoryGrid.RowSpacing = 4;

            label = uilabel(mandatoryGrid,'Text','Frequencies (Hz)', ...
                'FontWeight','bold');
            label.Layout.Row = 1;
            app.FrequenciesField = uieditfield(mandatoryGrid,'text', ...
                'Value','1.2:1.2:24', ...
                'Tooltip','Examples: 1.2:1.2:24 or 2 4 8 12');
            app.FrequenciesField.Layout.Row = 2;

            label = uilabel(mandatoryGrid,'Text','Sampling rate (Hz)', ...
                'FontWeight','bold');
            label.Layout.Row = 3;
            app.SamplingRateField = uieditfield(mandatoryGrid,'numeric', ...
                'Value',app.Config.network.samplingRate, ...
                'Limits',[eps Inf]);
            app.SamplingRateField.Layout.Row = 4;

            note = uilabel(mandatoryGrid, ...
                'Text',sprintf(['Participant matrices are read as voxels x time.\n' ...
                'Each .mat file becomes one persisted FREQ result.']), ...
                'FontColor',colors.muted, ...
                'VerticalAlignment','top');
            note.Layout.Row = 5;

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

            label = uilabel(optionalGrid,'Text','Filter progression');
            label.Layout.Row = 3;
            label.Layout.Column = 1;
            app.FilterTypeDropDown = uidropdown(optionalGrid, ...
                'Items',{'logarithmic','linear'}, ...
                'Value',app.Config.network.filter);
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
            outputGrid = uigridlayout(outputPanel,[7 1]);
            outputGrid.RowHeight = {20,34,48,28,48,44,'1x'};
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

            app.ProgressTextArea = uitextarea(outputGrid, ...
                'Editable','off', ...
                'Value',{'Progress messages will appear here.'}, ...
                'FontName','Courier New', ...
                'FontSize',10);
            app.ProgressTextArea.Layout.Row = 7;
        end

        function createSecondaryLayer(app,parent,colors)
            panel = uipanel(parent, ...
                'Title','3  Secondary analyses', ...
                'FontSize',14, ...
                'FontWeight','bold', ...
                'ForegroundColor',colors.navy, ...
                'BackgroundColor',[1 1 1]);
            panel.Layout.Row = 4;

            grid = uigridlayout(panel,[2 5]);
            grid.RowHeight = {56,'1x'};
            grid.ColumnWidth = {220,'1x',110,'1x',250};
            grid.Padding = [12 8 12 10];
            grid.RowSpacing = 8;
            grid.ColumnSpacing = 10;

            app.NetworkDropLabel = uilabel(grid, ...
                'Text',sprintf('Drop a FREQ_Networks* folder here\nor use Browse'), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','center', ...
                'FontWeight','bold', ...
                'FontColor',colors.blue, ...
                'BackgroundColor',colors.softBlue);
            app.NetworkDropLabel.Layout.Row = [1 2];
            app.NetworkDropLabel.Layout.Column = 1;

            app.NetworkPathField = uieditfield(grid,'text', ...
                'Editable','off', ...
                'Placeholder','No FREQ_Networks* folder selected');
            app.NetworkPathField.Layout.Row = 1;
            app.NetworkPathField.Layout.Column = 2;

            app.NetworkBrowseButton = uibutton(grid,'push', ...
                'Text','Browse...', ...
                'ButtonPushedFcn',@(~,~)app.browseNetworkFolder());
            app.NetworkBrowseButton.Layout.Row = 1;
            app.NetworkBrowseButton.Layout.Column = 3;

            app.AnalysisOutputField = uieditfield(grid,'text', ...
                'Editable','off', ...
                'Placeholder','FREQ_Analyses* is derived from FREQ_Networks*');
            app.AnalysisOutputField.Layout.Row = 1;
            app.AnalysisOutputField.Layout.Column = 4;

            app.SecondarySummaryLabel = uilabel(grid, ...
                'Text','One output subfolder per function', ...
                'FontColor',colors.muted, ...
                'HorizontalAlignment','right');
            app.SecondarySummaryLabel.Layout.Row = 1;
            app.SecondarySummaryLabel.Layout.Column = 5;

            modules = freqnessgui.moduleRegistry();
            moduleData = cell(numel(modules),4);
            for modulei = 1:numel(modules)
                if modules(modulei).requiresMNI && ...
                        modules(modulei).requiresSourceData
                    requirements = 'FREQ + MNI + source data';
                elseif modules(modulei).requiresMNI
                    requirements = 'FREQ + MNI';
                elseif modules(modulei).requiresSourceData
                    requirements = 'FREQ + source data';
                else
                    requirements = 'FREQ';
                end
                moduleData(modulei,:) = {modules(modulei).name, ...
                    modules(modulei).functionName,modules(modulei).folderName,requirements};
            end

            app.ModuleTable = uitable(grid, ...
                'Data',moduleData, ...
                'ColumnName',{'Analysis','Backend function','Output subfolder','Requires'}, ...
                'ColumnWidth',{190,200,150,'auto'}, ...
                'RowName',{});
            app.ModuleTable.Layout.Row = 2;
            app.ModuleTable.Layout.Column = [2 5];
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
                uiFileDnD(app.NetworkDropLabel, ...
                    @(~,dropData)app.handleNetworkDrop(dropData));
            catch exception
                app.FooterStatusLabel.Text = sprintf( ...
                    'Drag-and-drop unavailable (%s). Browse buttons remain available.', ...
                    exception.message);
            end
        end

        function handleDatasetDrop(app,dropData)
            folder = app.singleDroppedFolder(dropData,'Dataset*');
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
                app.NetworkPathField.Value = networkSet.folder;
                app.NetworkPathField.Tooltip = networkSet.folder;
                app.AnalysisOutputField.Value = networkSet.analysisFolder;
                app.AnalysisOutputField.Tooltip = networkSet.analysisFolder;
                app.SecondarySummaryLabel.Text = sprintf( ...
                    '%d participant results ready',networkSet.nParticipants);
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
            config.network.frequencies = freqnessgui.parseNumericVector( ...
                app.FrequenciesField.Value,false,'Frequencies');
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
                app.NetworkBrowseButton.Enable = 'off';
                app.ProgressTextArea.Value = {'Starting network estimation...'};
            else
                app.RunButton.Text = 'Run Network Estimation';
                app.DatasetBrowseButton.Enable = 'on';
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
