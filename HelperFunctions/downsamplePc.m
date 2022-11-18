function pcDownsampled = downsamplePc(pc,downsampleRate)
% Downsample point cloud with ratio specified in downsampleRate
    stepSize = round(1/downsampleRate);
    indices = 1:stepSize:pc.Count;
    pcDownsampled = select(pc,indices);
end

