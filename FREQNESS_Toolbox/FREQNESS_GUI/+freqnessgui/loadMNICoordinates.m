function mni = loadMNICoordinates(mniFile)
%LOADMNICOORDINATES Validate and load one N-by-3 coordinate matrix.

if ~(ischar(mniFile) || (isstring(mniFile) && isscalar(mniFile))) || ...
        ~isfile(mniFile)
    error('FREQNESS:GUI:InvalidMNIFile', ...
        'Select an existing MNI coordinate MAT-file.');
end

mniFile = char(mniFile);
variables = whos('-file',mniFile);
isCandidate = false(size(variables));
for variablei = 1:numel(variables)
    variableSize = variables(variablei).size;
    isCandidate(variablei) = ismember(variables(variablei).class, ...
        {'double','single','int8','uint8','int16','uint16', ...
        'int32','uint32','int64','uint64'}) && ...
        numel(variableSize) == 2 && any(variableSize == 3) && ...
        all(variableSize > 0) && ~variables(variablei).complex;
end

candidateIndices = find(isCandidate);
if numel(candidateIndices) ~= 1
    error('FREQNESS:GUI:InvalidMNIFile', ...
        ['The MAT-file must contain exactly one real numeric coordinate ' ...
        'matrix sized N-by-3 or 3-by-N.']);
end

variableName = variables(candidateIndices).name;
loaded = load(mniFile,variableName);
coordinates = loaded.(variableName);
if size(coordinates,2) ~= 3
    coordinates = coordinates.';
end
coordinates = double(coordinates);
if size(coordinates,2) ~= 3 || any(~isfinite(coordinates(:)))
    error('FREQNESS:GUI:InvalidMNIFile', ...
        'MNI coordinates must be a finite N-by-3 numeric matrix.');
end

[~,fileName,extension] = fileparts(mniFile);
mni = struct();
mni.schemaVersion = 1;
mni.path = mniFile;
mni.name = [fileName extension];
mni.variableName = variableName;
mni.coordinates = coordinates;
mni.nPoints = size(coordinates,1);

end
