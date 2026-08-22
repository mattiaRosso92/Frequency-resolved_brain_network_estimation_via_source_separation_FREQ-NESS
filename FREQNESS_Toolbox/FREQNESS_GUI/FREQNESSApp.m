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
            app.loadDefaultMNI();
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
            app.MainGrid.RowHeight = {72,250,455,190,26};
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
            mandatoryGrid.RowHeight = {20,34,54,20,10,'1x',14,42};
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
                'HorizontalAlignment','center');
            app.FrequencySummaryLabel.Layout.Row = 4;

            app.FWHMAxes = uiaxes(mandatoryGrid);
            app.FWHMAxes.Layout.Row = 6;
            app.FWHMAxes.FontSize = 9;
            app.FWHMAxes.Box = 'on';
            app.FWHMAxes.XGrid = 'on';
            app.FWHMAxes.YGrid = 'on';
            ylabel(app.FWHMAxes,'FWHM (Hz)');

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

            app.applyFrequencyRange([app.Config.network.frequencies(1), ...
                app.Config.network.frequencies(end)]);
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
                return
            end

            if isempty(app.Dataset)
                app.MNISummaryLabel.Text = sprintf( ...
                    '%d coordinates — drag the plot to rotate', ...
                    app.MNI.nPoints);
                app.MNISummaryLabel.FontColor = [0.055 0.415 0.690];
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
                'FWHM preview  |  %.4g–%.4g Hz  |  Δf %.4g Hz  |  %d frequencies', ...
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
                cla(app.FWHMAxes);
                plot(app.FWHMAxes,frequencies,fwhm,'-o', ...
                    'Color',[0.055 0.415 0.690], ...
                    'MarkerFaceColor',[0.055 0.415 0.690], ...
                    'MarkerSize',3, ...
                    'LineWidth',1.2);
                axis(app.FWHMAxes,'on');
                app.FWHMAxes.Box = 'on';
                app.FWHMAxes.XGrid = 'on';
                app.FWHMAxes.YGrid = 'on';
                ylabel(app.FWHMAxes,'FWHM (Hz)');
                xlim(app.FWHMAxes,app.FrequencyRangeSlider.Limits);
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
