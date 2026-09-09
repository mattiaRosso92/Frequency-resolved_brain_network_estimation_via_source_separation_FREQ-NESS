function modules = moduleRegistry()
%MODULEREGISTRY Describe secondary-analysis modules and output folders.

rows = {
    'entropy',          'Entropy landscape',       'Characterize network structure',      'FREQNESS_EntropyLandscape', 'EntropyLandscape', false, false, false, 'Quantify frequency-resolved network diversity and effective dimensionality.'
    'exponential_dk',   'Exponential decay',        'Characterize network structure',      'FREQNESS_ExponentialDK',    'ExponentialDK',    false, false, false, 'Model decay across the ordered component eigenspectrum.'
    'freq_gradients',   'Frequency gradients',      'Spatial organization',                'FREQNESS_FreqGradients',    'FreqGradients',    true,  false, false, 'Model how spatial activation patterns vary across frequencies.'
    'comp_gradients',   'Component gradients',      'Spatial organization',                'FREQNESS_CompGradients',    'CompGradients',    true,  false, false, 'Model spatial organization across ordered network components.'
    'cross_coupling',   'Cross-frequency coupling', 'Network dynamics',                    'FREQNESS_CrossCoupling',    'CrossCoupling',    false, false, false, 'Estimate phase-amplitude coupling between low and carrier frequencies.'
    'induced_responses','Induced responses',         'Network dynamics',                    'FREQNESS_InducedResponses', 'InducedResponses', false, false, true,  'Estimate event-related induced power in FREQ-NESS component time series.'
    'backprojection',   'Backprojection',            'Reconstruction and manipulation',    'FREQNESS_BackProjection',   'BackProjection',   false, false, false, 'Reconstruct selected components in voxel space.'
    'network_removal',  'Network removal',           'Reconstruction and manipulation',    'FREQNESS_NetworkRemoval',   'NetworkRemoval',   false, true,  false, 'Remove selected network contributions from the original broadband data.'
    'visualizer',       'Visualization & maps',      'Visualize and export',                'FREQNESS_Visualizer',       'Visualizer',       true,  false, false, 'Create network-landscape figures, spatial maps, and optional NIFTI exports.'
    };

template = struct('id','','name','','functionName','','folderName','', ...
    'category','','description','','requiresMNI',false, ...
    'requiresSourceData',false,'requiresEvents',false);
modules = repmat(template,size(rows,1),1);
for rowi = 1:size(rows,1)
    modules(rowi).id = rows{rowi,1};
    modules(rowi).name = rows{rowi,2};
    modules(rowi).category = rows{rowi,3};
    modules(rowi).functionName = rows{rowi,4};
    modules(rowi).folderName = rows{rowi,5};
    modules(rowi).requiresMNI = rows{rowi,6};
    modules(rowi).requiresSourceData = rows{rowi,7};
    modules(rowi).requiresEvents = rows{rowi,8};
    modules(rowi).description = rows{rowi,9};
end

end
