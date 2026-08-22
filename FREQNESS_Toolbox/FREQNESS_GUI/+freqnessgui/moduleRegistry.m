function modules = moduleRegistry()
%MODULEREGISTRY Describe secondary-analysis modules and output folders.

rows = {
    'visualizer',       'Visualization & spatial maps', 'FREQNESS_Visualizer',       'Visualizer',       true,  false
    'entropy',          'Entropy landscape',            'FREQNESS_EntropyLandscape', 'EntropyLandscape', false, false
    'exponential_dk',   'Exponential decay',             'FREQNESS_ExponentialDK',    'ExponentialDK',    false, false
    'freq_gradients',   'Frequency gradients',           'FREQNESS_FreqGradients',    'FreqGradients',    true,  false
    'comp_gradients',   'Component gradients',           'FREQNESS_CompGradients',    'CompGradients',    true,  false
    'cross_coupling',   'Cross-frequency coupling',      'FREQNESS_CrossCoupling',    'CrossCoupling',    false, false
    'backprojection',   'Backprojection',                'FREQNESS_BackProjection',   'BackProjection',   false, false
    'network_removal',  'Network removal',               'FREQNESS_NetworkRemoval',   'NetworkRemoval',   false, true
    'induced_responses','Induced responses',              'FREQNESS_InducedResponses', 'InducedResponses', false, false
    };

template = struct('id','','name','','functionName','','folderName','', ...
    'requiresMNI',false,'requiresSourceData',false);
modules = repmat(template,size(rows,1),1);
for rowi = 1:size(rows,1)
    modules(rowi).id = rows{rowi,1};
    modules(rowi).name = rows{rowi,2};
    modules(rowi).functionName = rows{rowi,3};
    modules(rowi).folderName = rows{rowi,4};
    modules(rowi).requiresMNI = rows{rowi,5};
    modules(rowi).requiresSourceData = rows{rowi,6};
end

end
