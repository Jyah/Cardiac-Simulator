clc;clear;
%% Model Shape
lambda_i = 0.35;% inner radius
lambda_o = 0.55;% outer radius
delta = 4; % focal radius, cm
Nphi = 200;
phi = linspace(0,360,Nphi)/180*pi;
Nl = 100;
lambda = linspace(lambda_i,lambda_o,Nl);
Ni = 100;
ita = linspace(0,120,Ni)/180*pi;
[PHI,LAMBDA,ITA] = ndgrid(phi,lambda,ita);
x = delta*sinh(LAMBDA).*sin(ITA).*cos(PHI);
y = delta*sinh(LAMBDA).*sin(ITA).*sin(PHI);
z = delta*cosh(LAMBDA).*cos(ITA);
figure;scatter3(x(:),y(:),z(:),'.k');axis image;
% K = boundary(x(:),y(:),z(:),0.9);
% trisurf(K,x(:),y(:),z(:),'Facecolor','red','Facealpha',0.1,'Edgecolor','none');axis image;
% xlabel('x (cm)');ylabel('y (cm)');zlabel('z (cm)');

%% create a 2D slice in original iamge and add mask
dr = 0.1;
[rx,ry,rz] = ndgrid(-2.5:dr:2.5,-2.5:dr:2.5,-2.5:dr:4.5);
r1_INV = sqrt(rx.^2+ry.^2+(rz+delta).^2);
r2_INV = sqrt(rx.^2+ry.^2+(rz-delta).^2);
lambda_INV = acosh((r1_INV+r2_INV)/(2*delta));
ita_INV = acos((r1_INV-r2_INV)/(2*delta));
phi_INV = atan(ry./rx);

%% Add MRI Tags
D0 = 300; % AU*
TR = 10;% sec
TE = 0.03;% sec
T1 = 0.6;% sec
T2 = 0.1; % sec
kx = 8; % rad/cm
ky = 8; % rad/cm
theta = 45; % degrees
Td = 0.017;% Td = ti-to; sec
mu = zeros(size(rx));
for layer = 1:size(rx,3)
  rxl = rx(:,:,layer);
  ryl = ry(:,:,layer);
  EPS_r0 = funEPS(rxl,ryl,theta,kx,ky);
  mu(:,:,layer) = D0*exp(-TE/T2)*(1+((1-exp(-(TR-Td)/T1))*EPS_r0-1)*exp(-Td/T2));
end
mask = zeros(size(rx));
mask(lambda_INV>=lambda_i & lambda_INV<=lambda_o & ita_INV>=0 & ita_INV<=(120/180*pi))=1;
% figure;scatter3(rx(mask>0),ry(mask>0),rz(mask>0),1,mu(mask>0));axis image;
% xlabel('x (cm)');ylabel('y (cm)');zlabel('z (cm)');

%% Motion Model
a = mean(cosh(lambda)./sinh(lambda)); % correcitonal parameter
% the wall volume of the model LV
Vw = pi*delta^3/4*(3*(cosh(lambda_o)-cosh(lambda_i))+4*((cosh(lambda_o))^3-(cosh(lambda_i))^3));
fpath = 'C:\Users\Jyahway\Git-code\Data\';
kpath = [fpath,'ks_csv'];
Ks = convert13ks(kpath);
mxc = rx(mask>0);
px = mxc(:);
myc = ry(mask>0);
py = myc(:);
mzc = rz(mask>0);
pz = mzc(:);
new_x = cell(60,1);
new_y = cell(60,1);
new_z = cell(60,1);
new_mu = mu(mask>0);
parfor k = 1:60
    [new_x{k},new_y{k},new_z{k}] = transK(a,Ks(:,k),px,py,pz,Vw);
end
xmin = min(cell2mat(new_x));
xmax = max(cell2mat(new_x));
ymin = min(cell2mat(new_y));
ymax = max(cell2mat(new_y));
zmin = min(cell2mat(new_z));
zmax = max(cell2mat(new_z));

%% Write VTK to Paraview
[xq,yq,zq] = ndgrid(xmin:dr:xmax,ymin:dr:ymax,zmin:dr:zmax);
slice = cell(60,1);
parfor k = 1:60
    F = scatteredInterpolant(new_x{k},new_y{k},new_z{k},new_mu,'natural','none');
    slice{k} = F(xq,yq,zq);
    fname = sprintf('%s%sslice_3d%02.0f.vtk',fpath,filesep,k);
    Mat2VTK(fname,slice{k},'binary');
end
%% function p = invTran(r,k,a)
lyr = 1;
r_x = xq(:);
r_y = yq(:);
r_z = zq(:);
kk = Ks(:,lyr);
F0 = diag([a^(1/3),a^(1/3),a^(-2/3),1]);
Fs = FaF6F5F4F3(kk,a);

%% Iteratively find the reference map
N = length(r_x);
p_0 = zeros(N,3);
parfor ir = 1:N
    x0 = [r_x(ir),r_y(ir),r_z(ir)];
    objfun = @(p)newPF(p(1),p(2),p(3),x0,Vw,F0,Fs,a,kk);
%     options = optimset('Display','iter','PlotFcns',@optimplotfval);
    p_0(ir,:) = fminsearch(objfun, x0);
end

%% Transform reference map to all frames
% figure;scatter3(p_0(:,1),p_0(:,2),p_0(:,3),1,slice{1}(:));axis image;
px = p_0(:,1);
py = p_0(:,2);
pz = p_0(:,3);
new_lx = cell(60,1);
new_ly = cell(60,1);
new_lz = cell(60,1);

parfor k = 1:60
    [new_lx{k},new_ly{k},new_lz{k}] = transK(a,Ks(:,k),px,py,pz,Vw);
end
X = cell2mat(new_lx');
Y = cell2mat(new_ly');
Z = cell2mat(new_lz');

%% Displacement Field
DX = X-X(:,1);
DY = Y-Y(:,1);
DZ = Z-Z(:,1);
n = 40;
sa_mask = ~isnan(slice{n});
vx = xq(sa_mask>0);
vy = yq(sa_mask>0);
vz = zq(sa_mask>0);
dx = reshape(DX(:,n),size(sa_mask));
dy = reshape(DY(:,n),size(sa_mask));
dz = reshape(DZ(:,n),size(sa_mask));
quiver3(vx(:),vy(:),vz(:),dx(sa_mask>0),dy(sa_mask>0),dz(sa_mask>0),'Color','black');




