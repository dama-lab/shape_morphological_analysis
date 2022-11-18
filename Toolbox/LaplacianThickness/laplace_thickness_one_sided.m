function [L, epsilon_array] = laplace_thickness_one_sided(p, direction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculate layer thickness map from the laplace gradient field
%  Version: 0.4 (Jul 05 2013)
%  Author: Da Ma ( d.ma.11@ucl.ac.uk ) CMIC/CABI UCL
%  Reference:
%  [1] Yezzi, A. J., & Prince, J. L. (2002). A PDE approach for thickness,
%       correspondence, and gridding of annular tissues.
%  [2] Rocha, K. R., Yezzi, A. J., & Prince, J. L. (2007).
%       A Hybrid Eulerian�CLagrangian Approach for Thickness, Correspondence,
%       and Gridding of Annular Tissues.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% direction: 1=upwind; -1=downwind

%% initialize
p = double(p);
layer = double((p>0).*(p<1));
[x,y,z]=size(layer);

%% Step 1: Create normalize (unit tangental) vector field for (gx,gy,gz)
[gx,gy,gz]=gradient(1-p);
d=(gx.^2+gy.^2+gz.^2).^0.5; % calculate denominator
vx=layer .* gx./d; % normalize x
vy=layer .* gy./d; % normalize y
vz=layer .* gz./d; % normalize z

%% set non-layer region to zero
vx(isnan(vx))=0;
vy(isnan(vy))=0;
vz(isnan(vz))=0;

vx(isinf(vx))=0;
vy(isinf(vy))=0;
vz(isinf(vz))=0;

%%
vx = vx * direction;
vy = vy * direction;
vz = vz * direction;

%% Step 2: calculate layer thickness on every voxel with PDF (Yezzi et
% al. 2003)
niter=500;            % number of maximum iterations
epsilon_min=10^-5;     % minimum energy (convergence ratio) threshold = 1e-4
epsilon_array=[];

%%======== update L ========
% fprintf('update L ...');
L=zeros(x,y,z); % compute L map towards upper layer (upwind)
epsilon=1;            % initialize convergence ratio

for it=1:niter
    % determine convergence point
    if ( epsilon < epsilon_min )
        break;
    end
    % or no further improvement over epsilon
    if (length(epsilon_array) > 2) && ((epsilon_array(end-1)-epsilon_array(end)) < 10^-8)
        break;
    end
    
    L_old=L;
    L=length_map(vy,vx,vz,L_old,layer);     
    % impose boundar condition
    % (may not neccessary as already constrained in fuction length_map)
    L=L.*layer;

    % calculate the total energy for L;
    epsilon=sum(abs(L(:)-L_old(:)))/sum(L_old(:));    
    epsilon_array = [epsilon_array; epsilon];
    % display for every 20 iterations
    if (rem(it,20) == 0) 
        fprintf('[%d] %.2d ',it,epsilon);
        if (rem(it,100) == 0) 
            fprintf('\n');
        end
    end
end

end
