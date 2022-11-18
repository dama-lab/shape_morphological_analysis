function pt2surf(vertex)
    





%% %%% Some additional comments to useIBSL3_3D_T()
% Ref: https://www.mathworks.com/matlabcentral/fileexchange/47490-surface-reconstruction-using-implicit-b-splines-fast#overview_tab
% as far as I have figured it out the indices of faces have to be zero-based. These indices will lateron be elevated in the function:
% [in out,Normal]=NewOffset(r,pnt, tri) --- line25: tri=tri+1;
% I solved this issue like this ;):
% Tri = Tri-1;
% before calling IBSL3_3D_T().
% 
% Another thing that I had to fix: If the coordinates have their lower bounds to zero and their upper bounds at one it might result an error when the offsets "s&t" are created. There are values possible that are lower zeros. This results in an error lateron in function:
% IND=activeINDEX(P,SET)
% 
% So you wanna make sure that your 1x1x1 cube has a small offset to your point-coordinates. There might be easier ways for normalizing the point-coordinates but this one works (Cor are the coordinates [NPx3]):
% 
% MinX = min(Cor(:,1))-0.1*abs(min(Cor(:,1)));
% MaxX = max(Cor(:,1))+0.1*abs(min(Cor(:,1)));
% MinY = min(Cor(:,2))-0.1*abs(min(Cor(:,2)));
% MaxY = max(Cor(:,2))+0.1*abs(min(Cor(:,2)));
% MinZ = min(Cor(:,3))-0.1*abs(min(Cor(:,3)));
% MaxZ = max(Cor(:,3))+0.1*abs(min(Cor(:,3)));
% 
% CorNorm = Cor;
% CorNorm(:,1) = (CorNorm(:,1)-MinX)/((MaxX-MinX));
% CorNorm(:,2) = (CorNorm(:,2)-MinY)/((MaxY-MinY));
% CorNorm(:,3) = (CorNorm(:,3)-MinZ)/((MaxZ-MinZ));
% 
% BTW: Great work Mohammad! Really like your surface tool! =)

%% For IBSL3_3DTRI: Ref: https://www.mathworks.com/matlabcentral/fileexchange/44654-surface-reconstruction-using-implicit-b-splines
% 
% I solved the usage for your own data set. Here is the code I am using:
% 
% %(x,y,z) is my point-cloud
% N = 1000; %for large set of data -> use every N-th point
% tri = delaunay(x(1:N:end),y(1:N:end));
% tri = tri-1; %needed for this code
% 
% %normalization
% MinX = min(x)-0.1*abs(min(x));
% MaxX = max(x)+0.1*abs(min(x));
% MinY = min(y)-0.1*abs(min(y));
% MaxY = max(y)+0.1*abs(min(y));
% MinZ = min(z)-0.1*abs(min(z));
% MaxZ = max(z)+0.1*abs(min(z));
% xn = x; yn = y; zn = z;
% xn = (xn-MinX)/((MaxX-MinX));
% yn = (yn-MinY)/((MaxY-MinY));
% zn = (zn-MinZ)/((MaxZ-MinZ));
% 
% L=10;
% P=IBSL3_3DTRI(.01,20,L,[xn(1:N:end),yn(1:N:end),zn(1:N:end)], tri);
% IBSLevelSurf(P,[.5 .6 .8],0.03);

