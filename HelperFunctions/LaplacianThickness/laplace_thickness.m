function [T, L0, L1, epsilon_array] = laplace_thickness(p)
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
% for retinal thickness, only L0 is enough, no need to calculate both sides

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

%% Step 2: calculate layer thickness on every voxel with PDF (Yezzi et
% al. 2003)
niter=500;            % number of maximum iterations
epsilon_min=10^-5;     % minimum energy (convergence ratio) threshold = 1e-4
epsilon_array=[];

%======== update L0 ========
fprintf('update L0 ...');
L0=zeros(x,y,z); % compute L0 map towards upper layer (upwind)
epsilon=1;            % initialize convergence ratio

for it=1:niter
    % determine convergence point
    if ( epsilon < epsilon_min )
        break;
    end
    
    L0_old=L0;
    L0=length_map(vy,vx,vz,L0_old,layer);     
    % impose boundar condition
    % (may not neccessary as already constrained in fuction length_map)
    L0=L0.*layer;
    
    % calculate the total energy for L0;
    epsilon=sum(abs(L0(:)-L0_old(:)))/sum(L0_old(:));    
    epsilon_array = [epsilon_array; epsilon];
    % display for every 20 iterations
    if (rem(it,20) == 0) 
        fprintf('[%d] %.2d ',it,epsilon);
        if (rem(it,100) == 0) 
            fprintf('\n');
        end
    end
end

%======== update L1 ========
fprintf(' update L1 ...')
L1=zeros(x,y,z); % compute L1 map towards WM (downwind)
epsilon=1;            % initialize convergence ratio
for it=1:niter
    % determine convergence point
    if ( epsilon < epsilon_min )
        break;
    end
    
    L1_old=L1;
    L1=length_map(-vy,-vx,-vz,L1_old,layer);    
    % impose boundar condition
    % (may not neccessary as already constrained in fuction length_map)
    L1=L1.*layer;
    
    % calculate the total energy for L0;
    epsilon=sum(abs(L1(:)-L1_old(:)))/sum(L1_old(:));    
    
    % display for every 20 iterations
    if (rem(it,20) == 0) 
        fprintf('[%d] %.2d ',it,epsilon);
        if (rem(it,100) == 0) 
            fprintf('\n');
        end
    end
end

%% Thickness map
T = L0 + L1;

function L=length_map(vx,vy,vz,L_old,layer)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate length map based on normalized tangent vector field (gx,gy,gz) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L =  thickness map along one direction
[x,y,z]=size(vx);

%% new fast implementation
% create 9 3D matrixs (or 6 for 2D) to calculate each pixel in L from L_old at once
% x
sx_plus=(vx>0);
sx_minus=(vx<0);
sx_zero=(vx==0);
x_plus=cat(1,L_old(2:x,:,:),zeros(1,y,z)).*sx_plus;
x_minus=cat(1,zeros(1,y,z),L_old(1:x-1,:,:)).*sx_minus;
x_zero=L_old.*sx_zero;

% y
sy_plus=(vy>0);
sy_minus=(vy<0);
sy_zero=(vy==0);
y_plus=cat(2,L_old(:,2:y,:),zeros(x,1,z)).*sy_plus;
y_minus=cat(2,zeros(x,1,z),L_old(:,1:y-1,:)).*sy_minus;
y_zero=L_old.*sy_zero;

% z
sz_plus=(vz>0);
sz_minus=(vz<0);
sz_zero=(vz==0);
z_plus=cat(3,L_old(:,:,2:z),zeros(x,y,1)).*sz_plus;
z_minus=cat(3,zeros(x,y,1),L_old(:,:,1:z-1)).*sz_minus;
z_zero=L_old.*sz_zero;


% PDE [equation (8) in Yezzi, 2002; equation (6)(7) in Rocha, 2007] 
ax=abs(vx).*(x_plus+x_minus+x_zero); % +x_zero
ay=abs(vy).*(y_plus+y_minus+y_zero); % +y_zero
az=abs(vz).*(z_plus+z_minus+z_zero); % +z_zero
L=layer.*(ones(x,y,z)+ax+ay+az)./(abs(vx)+abs(vy)+abs(vz));

% ax=(vx).*(x_plus+x_minus); % +x_zero
% ay=(vy).*(y_plus+y_minus); % +y_zero
% az=(vz).*(z_plus+z_minus); % +z_zero
% L=layer.*(ones(x,y,z)+ax+ay+az)./((vx)+(vy)+(vz));

L(isnan(L))=0;
L(isinf(L))=0;



%%