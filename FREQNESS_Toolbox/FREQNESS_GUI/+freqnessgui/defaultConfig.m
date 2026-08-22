function config = defaultConfig()
%DEFAULTCONFIG Return the default configuration for network estimation.

config = struct();
config.schemaVersion = 1;
config.datasetFolder = '';
config.outputFolder = '';
config.mniFile = '';

config.network = struct();
config.network.frequencies = 1.2:1.2:24;
config.network.samplingRate = 250;
config.network.duration = [];
config.network.fwidth = [];
config.network.filter = 'logarithmic';
config.network.regularisation = 0.01;
config.network.ncomps = 30;
config.network.badSegments = [];
config.network.rescale = false;

config.execution = struct();
config.execution.recomputeExisting = false;

end
