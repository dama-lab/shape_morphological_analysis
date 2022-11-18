function CorNorm = normalizeCoordinate(Cor)
    %% Ref: https://www.mathworks.com/matlabcentral/fileexchange/47490-surface-reconstruction-using-implicit-b-splines-fast
    MinX = min(Cor(:,1))-0.1*abs(min(Cor(:,1)));
    MaxX = max(Cor(:,1))+0.1*abs(min(Cor(:,1)));
    MinY = min(Cor(:,2))-0.1*abs(min(Cor(:,2)));
    MaxY = max(Cor(:,2))+0.1*abs(min(Cor(:,2)));
    MinZ = min(Cor(:,3))-0.1*abs(min(Cor(:,3)));
    MaxZ = max(Cor(:,3))+0.1*abs(min(Cor(:,3)));

    CorNorm = Cor;
    CorNorm(:,1) = (CorNorm(:,1)-MinX)/((MaxX-MinX));
    CorNorm(:,2) = (CorNorm(:,2)-MinY)/((MaxY-MinY));
    CorNorm(:,3) = (CorNorm(:,3)-MinZ)/((MaxZ-MinZ));
end