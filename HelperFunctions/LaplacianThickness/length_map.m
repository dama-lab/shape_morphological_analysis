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