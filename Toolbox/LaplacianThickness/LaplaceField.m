function [LaplacianField,LaplacianRegion] = LaplaceField(vol,l_middle,l_lower,verbose)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculate laplace gradient 
%  Version: 0.5 (Dec 1 2019)
%  Author: Da Ma ( dma73@sfu.ca, d.ma.11@ucl.ac.uk )
%  Reference:
%  Jones, S. E., Buchbinder, B. R., & Aharon, I. (2000).
%  Three-dimensional mapping of cortical thickness using Laplace's Equation.
%
%  Input:
%     l_upper:    (0) label number for upper boundary (outer csf+backgound)
%     l_middle:   (2) label number for middle slab (grey matter - cerebellum)
%     l_lower:    (1) label number for lower boundary (white matter + cerebellum)
%     verbose:     0): no verbose; 1) show progress
%  Output:
%     LaplacianField
%     LaplacianRegion: Region of Laplacian Field (0< LaplacianRegion <1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Python approach: just a filter in scipy:
%  import scipy.ndimage.filters as filters
%  lap = filters.laplace(image)
% Ref: https://stackoverflow.com/questions/1764859/how-to-compute-laplacian-of-a-field/1770788

% By defaul, don't display progress
if ~exist('verbose','var'); verbose=0;end 

% % start counting execution time
% time=cputime;

% %% Initialize layer information
% vol = vol_layer;
% l_middle = laplace_layer;
% l_lower =  laplace_layer+1;

%% Initialize boundary condition
[x,y,z]=size(vol);
% convert to double type
cortical = double(vol == l_middle);
boundary = double(vol == l_lower);
p = zeros(x,y,z); % initialize potential field

%% Initialize solver parameters
i=2:x-1;
j=2:y-1;
k=2:z-1;
epsilon_ratio_min=10^-5;	% energy convergence value = 1e-5
epsilon_ratio=1;            % Initialize energy convergence term
epsilon=1000;               % Initialize energy
niter=500;                % Number of iterations 

%% Iteratively solveing Laplace equation
for it=1:niter
    if ( epsilon_ratio < epsilon_ratio_min )
        break;
    end
    %%
    epsilon_old=epsilon;
    % Jacobei's method to update potential field p (Jones, 2000)
    p(i,j,k)=(p(i+1,j,k)+p(i-1,j,k)+p(i,j+1,k)+p(i,j-1,k)+p(i,j,k+1)+p(i,j,k-1))/6;
    
    % update Laplacian field
    % p=p.*cortical+boundary; % Boundary conditions (Dirichlet Condition?)
    p=p.*cortical;
    
    % Boundary conditions (Dirichlet Condition?)
    p(boundary==1) = 1; 
    
    %% Calculate the total energy term
    [gx,gy,gz]=gradient(p);
    [gxx,gxy,gxz]=gradient(gx);
    [gyx,gyy,gyz]=gradient(gy);
    [gzx,gzy,gzz]=gradient(gz);

    % diff is 3x faster than gradient, but will loss boundary voxels
%     gx  = diff(p,1,1);
%     gy  = diff(p,1,2);
%     gz  = diff(p,1,3);
%     gxx = diff(gx,1,1);
%     gxy = diff(gx,1,2);
%     gxz = diff(gx,1,3);
%     gyx = diff(gy,1,1);
%     gyy = diff(gy,1,2);
%     gyz = diff(gy,1,3);
%     gzx = diff(gz,1,1);
%     gzy = diff(gz,1,2);
%     gzz = diff(gz,1,3);

    % Potential solution: https://www.mathworks.com/matlabcentral/answers/142346-how-can-i-differentiate-without-decreasing-the-length-of-a-vector
    % 1. Fit your vectors with polyfit, use polyder to calculate the derivatives of the polynomial function, and then use polyval with the results of polyder to calculate the actual values at the x-values you choose;
    % 2. Take the derivative using diff, then use interp1 (with the 'extrap' option if necessary) to interpolate (and extrapolate) the derivative.
    % These approaches both assume your data are smooth and noise-free. Taking the derivative of a noisy signal is generally not recommended because the noise will predominate in the derivative.

    %% potential improvement of coding
%     for dim = 1:3
%         d
%     g
%     energy = sumsqr(gxx)
    
    energy=gxx.^2+gxy.^2+gxz.^2+gyx.^2+gyy.^2+gyz.^2+gzx.^2+gzy.^2+gzz.^2;
    epsilon=sum(energy(:));
    epsilon_ratio=abs(epsilon_old-epsilon)/epsilon_old;

    % If verbose, display for every 10 iterations
    if verbose == 1
        if (rem(it,10) == 0) 
            fprintf('[%d]%.2d ',it, epsilon_ratio);
            if (rem(it,100) == 0) 
                fprintf('\n');
            end
        end
    end
end
%%
LaplacianField = p;
LaplacianRegion = single(p>0.*p<1);