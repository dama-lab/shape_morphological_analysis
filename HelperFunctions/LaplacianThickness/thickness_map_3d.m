%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculate cortical thickness map from the laplace gradient field
%  Version: 0.4 (Jul 05 2013)
%  Author: Da Ma ( d.ma.11@ucl.ac.uk ) CMIC/CABI UCL
%  Reference:
%  [1] Yezzi, A. J., & Prince, J. L. (2002). A PDE approach for thickness,
%       correspondence, and gridding of annular tissues.
%  [2] Rocha, K. R., Yezzi, A. J., & Prince, J. L. (2007).
%       A Hybrid Eulerian¨CLagrangian Approach for Thickness, Correspondence,
%       and Gridding of Annular Tissues.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function t=thickness_map_3d_fast(p,cortical,three_layer_nii)
% calculate thickness map from:
% cortical:  cortical region (label=1 for cortical)
% p:         potential field derived from Laplace equation

function [L0,L1,t]=thickness_map_3d(nifti_name,result_folder)
% nifti_name: file name of the cortical gradient file (Laplace field)
% result_folder: folder to put result images

%start counting execution time
time=cputime;

if ~exist(result_folder,'dir') % 7=folder
    mkdir(result_folder);
end

% extract file name from path string
[~, input_name_nii, input_ext] = fileparts(nifti_name);
[~, input_name, input_ext_nii] = fileparts(input_name_nii);
% print laplace map file name for debug
fprintf('create cortical thickness map for: %s\n',input_name);

laplace_nii=load_untouch_nii(nifti_name); %
result=laplace_nii;
p=double(laplace_nii.img); % p: potential field derived from Laplace equation
cortical=double((p>0).*(p<1)); % cortical_label=1
% Step 1: Create normalize (unit tangental) vector field for (gx,gy,gz)
[gx,gy,gz]=gradient(1-p);
% gx,gy,gz:gradient along x, y, z axis, can improve for fast implementation;
d=(gx.^2+gy.^2+gz.^2).^0.5; % calculate denominator
vx=gx./d; % normalize x
vy=gy./d; % normalize y
vz=gz./d; % normalize z

vx=vx.*cortical;
vy=vy.*cortical;
vz=vz.*cortical;

% set non-cortical region to zero
vx(isnan(vx))=0;
vy(isnan(vy))=0;
vz(isnan(vz))=0;

vx(isinf(vx))=0;
vy(isinf(vy))=0;
vz(isinf(vz))=0;

% for debug purpose only
% vx_nii=make_nii(vx);
% view_nii(vx_nii);
% vy_nii=make_nii(vy);
% view_nii(vy_nii);


% % ==========================================================
% % old: create unit vector field for L1
% [gx_1,gy_1,gz_1]=gradient(1-p);
% d=(gx_1.^2+gy_1.^2+gz_1.^2).^0.5;
% vx_1=gx_1./d; % normalize x_1
% vy_1=gy_1./d; % normalize y_1
% vz_1=gz_1./d; % normalize z_1
% % =========================================================


% Step 2: calculate cortical thickness on every voxel with PDF (Yezzi et
% al. 2003)
niter=500;            % number of maximum iterations
epsilon_min=10^-5;     % minimum energy (convergence ratio) threshold = 1e-4

% initialize L0 and L1
[x,y,z]=size(cortical);

%======== update L1 ========
L1=zeros(x,y,z); % compute L1 map towards WM (downwind)
epsilon=1;            % initialize convergence ratio
for it=1:niter
    % determine convergence point
    if ( epsilon < epsilon_min )
        break;
    end
    
    L1_old=L1;
    L1=length_map(-vy,-vx,-vz,L1_old,cortical);    
    % impose boundar condition
    % (may not neccessary as already constrained in fuction length_map)
    L1=L1.*cortical;
    
    % calculate the total energy for L0;
    epsilon=sum(abs(L1(:)-L1_old(:)))/sum(L1_old(:));    
    
    if (it/10==floor(it/10)) 
        fprintf('L1: %d, %d\n',it,epsilon);
    end
end

% derive White matter distance (WMD)
if ~exist([result_folder,'/WMD/'],'dir') % 7=folder
    mkdir([result_folder,'/WMD/']);
end 
result.img=L1;
nifti_WMD=strcat(result_folder,'/WMD/',input_name,input_ext_nii,input_ext);
save_untouch_nii(result,nifti_WMD);

% for debug purpose only
%     L1_nii=make_nii(L1);
%    close(gcf);
%    view_nii(L1_nii);

%======== update L0 ========
L0=zeros(x,y,z); % compute L0 map towards CSF (upwind)
epsilon=1;            % initialize convergence ratio
for it=1:niter
    % determine convergence point
    if ( epsilon < epsilon_min )
        break;
    end
    
    L0_old=L0;
    L0=length_map(vy,vx,vz,L0_old,cortical);    
    % impose boundar condition
    % (may not neccessary as already constrained in fuction length_map)
    L0=L0.*cortical;
    
    % calculate the total energy for L0;
    epsilon=sum(abs(L0(:)-L0_old(:)))/sum(L0_old(:));    

    % display for every 10 iterations
    if (it/10==floor(it/10)) 
        fprintf('L0: %d, %d\n',it,epsilon);
    end
end

% derive CSF distance (CSFD)
if ~exist([result_folder,'/CSFD/'],'dir') % 7=folder
    mkdir([result_folder,'/CSFD/']);
end
result.img=L0;
nifti_CSFD=strcat(result_folder,'/CSFD/',input_name,input_ext_nii,input_ext);
save_untouch_nii(result,nifti_CSFD);

% for debug purpose only
%     L0_nii=make_nii(L0);
%      close(gcf);
%    view_nii(L0_nii);

% Preparing thickness folders
if ~exist([result_folder,'/CS/'],'dir') % 7=folder
    mkdir([result_folder,'/CS/']);
end
if ~exist([result_folder,'/thickness_map/'],'dir') % 7=folder
    mkdir([result_folder,'/thickness_map/']);
end
if ~exist([result_folder,'/CS_thickness/'],'dir') % 7=folder
    mkdir([result_folder,'/CS_thickness/']);
end


% ========= central surface ===========
result.img=(abs(L0-L1)<=1).*(laplace_nii.img>0).*(laplace_nii.img<1);
nifti_CS=strcat(result_folder,'/CS/',input_name,input_ext_nii,input_ext);
save_untouch_nii(result,nifti_CS);

% ========= thickness map ===========
result.img=L0+L1;
nifti_thickness=strcat(result_folder,'/thickness_map/',input_name,input_ext_nii,input_ext);
save_untouch_nii(result,nifti_thickness);

% ========= thickness map only on central surface ===========
result.img=(L0+L1).*(abs(L0-L1)<=1).*(laplace_nii.img>0).*(laplace_nii.img<1);
nifti_thickness=strcat(result_folder,'/CS_thickness/',input_name,input_ext_nii,input_ext);
save_untouch_nii(result,nifti_thickness);

% view_nii(three_layer_nii);
time=cputime-time; % calculate execute time
minutes=time/60;
fprintf('Total execute time: %d minutes\n',minutes);


function L=length_map(vx,vy,vz,L_old,cortical)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate length map based on normalized tangent vector field (gx,gy,gz) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L =  thickness map along one direction
[x,y,z]=size(vx);

%% new fast implementation
% create 6 (9) matrix to calculate each pixel in L from L_old at once
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
L=cortical.*(ones(x,y,z)+ax+ay+az)./(abs(vx)+abs(vy)+abs(vz));

% ax=(vx).*(x_plus+x_minus); % +x_zero
% ay=(vy).*(y_plus+y_minus); % +y_zero
% az=(vz).*(z_plus+z_minus); % +z_zero
% L=cortical.*(ones(x,y,z)+ax+ay+az)./((vx)+(vy)+(vz));

L(isnan(L))=0;
L(isinf(L))=0;

% following function is for debug, never called in the main function
function L=length_map_L1(vx,vy,vz,L_old,cortical) % for debug
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate length map based on normalized tangent vector field (gx,gy,gz) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L =  thickness map along one direction
[x,y,z]=size(vx);

%% new fast implementation
% create 6 (9) matrix to calculate each pixel in L from L_old at once
% x
sx_plus=(vx<0);
sx_minus=(vx>0);
% sx_zero=(vx==0);
x_plus=cat(1,L_old(2:x,:,:),zeros(1,y,z)).*sx_plus;
x_minus=cat(1,zeros(1,y,z),L_old(1:x-1,:,:)).*sx_minus;
% x_zero=L_old.*sx_zero;

% y
sy_plus=(vy<0);
sy_minus=(vy>0);
% sy_zero=(vy==0);
y_plus=cat(2,L_old(:,2:y,:),zeros(x,1,z)).*sy_plus;
y_minus=cat(2,zeros(x,1,z),L_old(:,1:y-1,:)).*sy_minus;
% y_zero=L_old.*sy_zero;

% z
sz_plus=(vz<0);
sz_minus=(vz>0);
% sz_zero=(vz==0);
z_plus=cat(3,L_old(:,:,2:z),zeros(x,y,1)).*sz_plus;
z_minus=cat(3,zeros(x,y,1),L_old(:,:,1:z-1)).*sz_minus;
% z_zero=L_old.*sz_zero;

% PDE (equation 8) in Yezzi
ax=abs(vx).*(x_plus+x_minus); % +x_zero
ay=abs(vy).*(y_plus+y_minus); % +y_zero
az=abs(vz).*(z_plus+z_minus); % +z_zero
L=cortical.*(ones(x,y,z)+ax+ay+az)./(abs(vx)+abs(vy)+abs(vz));
L(isnan(L))=0;
L(isinf(L))=0;

%% original slow implementation
% create sign matrix for Tangental field vx,vy,vz
% sx=sign(vx);
% sy=sign(vy);
% sz=sign(vz);
% ax=zeros(x,y,z); % first term in the numerator to update L0/L1
% ay=zeros(x,y,z); % second term in the numerator to update L0/L1
% az=zeros(x,y,z); % third term in the numerator to update L0/L1
% % PDE (equation 8/9) in Yezzi
% for i=2:x-1
%     for j=2:y-1
%         for k=2:z-1
%             ax(i,j,k)=abs(vx(i,j,k))*L_old(i+sx(i,j,k),j,k);
%             ay(i,j,k)=abs(vy(i,j,k))*L_old(i,j+sy(i,j,k),k);
%             az(i,j,k)=abs(vy(i,j,k))*L_old(i,j,k+sz(i,j,k));
%             if ( cortical(i,j,k)==1 ) % only update the cortical region
%                 denominator=abs(vx(i,j,k))+abs(vy(i,j,k))+abs(vz(i,j,k));
%                 if denominator==0
%                     L(i,j,k)=(L_old(i+1,j,k)+L_old(i-1,j,k)+L_old(i,j-1,k)+L_old(i,j,k+1)+L_old(i,j,k-1))/6; % is this right?
%                 else
%                     L(i,j,k)=(1+ax(i,j,k)+ay(i,j,k)+az(i,j,k))/(abs(vx(i,j,k))+abs(vy(i,j,k))+abs(vz(i,j,k)));
%                 end
%             end
%             % for debug purpose only
%             % fprintf('%d, %d, %d\n',i,j,k);
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Change Log:
% 2013.07.05 Version 0.4 use original nii file head instead of make_nii
%    to do: include laplace_equation_3d in the function