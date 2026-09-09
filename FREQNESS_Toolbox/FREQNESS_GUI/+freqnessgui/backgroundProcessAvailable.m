function available = backgroundProcessAvailable()
%BACKGROUNDPROCESSAVAILABLE Whether a separate MATLAB can be launched.

available = false;
try
    matlabExecutable = fullfile(matlabroot,'bin','matlab');
    if ispc && ~isfile(matlabExecutable)
        matlabExecutable = [matlabExecutable '.exe'];
    end
    available = ~isdeployed && usejava('jvm') && ...
        isfile(matlabExecutable);
catch
    % Synchronous execution remains available on unsupported releases.
end

end
