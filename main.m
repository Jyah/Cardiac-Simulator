clc;clear;close all;
%% Model Shape
lambda_i = 0.35;% inner radius
lambda_o = 0.55;% outer radius
delta = 4; % focal radius, cm
Nphi = 100;
phi = linspace(0,360,Nphi)/180*pi;
Nl = 50;
lambda = linspace(lambda_i,lambda_o,Nl);
Ni = 50;
ita = linspace(0,120,Ni)/180*pi;
[PHI,LAMBDA,ITA] = ndgrid(phi,lambda,ita);
x = delta*sinh(LAMBDA).*sin(ITA).*cos(PHI);
y = delta*sinh(LAMBDA).*sin(ITA).*sin(PHI);
z = delta*cosh(LAMBDA).*cos(ITA);
% figure;scatter3(x(:),y(:),z(:),'.k');axis image;
% K = boundary(x(:),y(:),z(:),0.9);
% trisurf(K,x(:),y(:),z(:),'Facecolor','red','Facealpha',0.1,'Edgecolor','none');axis image;
% xlabel('x (cm)');ylabel('y (cm)');zlabel('z (cm)');


r1_INV = sqrt(x.^2+y.^2+(z+delta).^2);
r2_INV = sqrt(x.^2+y.^2+(z-delta).^2);
lambda_INV = acosh((r1_INV+r2_INV)/(2*delta));
ita_INV = acos((r1_INV-r2_INV)/(2*delta));
phi_INV = atan(y./x);
mask = zeros(size(x));
mask(lambda_INV>=lambda_i & lambda_INV<=lambda_o & ita_INV>=0 & ita_INV<=(120/180*pi))=1;
%% create a 2D slice in original iamge and add mask
xt = min(x(:)):dr:max(x(:));
yt = min(y(:)):dr:max(y(:));
zt = 2;
[xt,yt,zt] = ndgrid(xt,yt,zt);
r1_INV = sqrt(xt.^2+yt.^2+(zt+delta).^2);
r2_INV = sqrt(xt.^2+yt.^2+(zt-delta).^2);
lambda_INV = acosh((r1_INV+r2_INV)/(2*delta));
ita_INV = acos((r1_INV-r2_INV)/(2*delta));
phi_INV = atan(yt./xt);
maskt = zeros(size(xt));
maskt(lambda_INV>=lambda_i & lambda_INV<=lambda_o & ita_INV>=0 & ita_INV<=(120/180*pi))=1;

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
mu = zeros(size(x));
for layer = 1:size(x,3)
  rxl = x(:,:,layer);
  ryl = y(:,:,layer);
  EPS_r0 = EPS(rxl,ryl,theta,kx,ky);
  mu(:,:,layer) = D0*exp(-TE/T2)*(1+((1-exp(-(TR-Td)/T1))*EPS_r0-1)*exp(-Td/T2));
end

% figure;scatter3(x(mask>0),y(mask>0),z(mask>0),1,mu(mask>0));axis image;
% xlabel('x (cm)');ylabel('y (cm)');zlabel('z (cm)');

%% Motion Model
a = mean(cosh(lambda)./sinh(lambda)); % correcitonal parameter
% the wall volume of the model LV
Vw = pi*delta^3/4*(3*(cosh(lambda_o)-cosh(lambda_i))+4*((cosh(lambda_o))^3-(cosh(lambda_i))^3));
fpath = 'C:\Users\Jyahway\Git-code\Data\';
kpath = [fpath,'ks_csv'];
Ks = convert13ks(kpath);
mxc = x(mask>0);
px = mxc(:);
myc = y(mask>0);
py = myc(:);
mzc = z(mask>0);
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
dr = 0.1;
% idx = find(z>1.99&z<2.01);
% xq = x(idx);
% yq = y(idx);
% zq = z(idx);
% figure;scatter3(xq(:),yq(:),zq(:),'.k');axis image;
[xq,yq,zq] = ndgrid(xmin:dr:xmax,ymin:dr:ymax,2);

r1q = sqrt(xq.^2+yq.^2+(zq+delta).^2);
r2q = sqrt(xq.^2+yq.^2+(zq-delta).^2);
lambdaq = acosh((r1q+r2q)/(2*delta));
itaq = acos((r1q-r2q)/(2*delta));
phiq = atan(yq./xq);

parfor k = 1:60
    F = scatteredInterpolant(new_x{k},new_y{k},new_z{k},new_mu,'natural','none');
    Fxyz = F(xq,yq,zq);
    
    slice{k} = Fxyz;
%     fname = sprintf('%s%sslice%02.0f.vtk',fpath,filesep,k);
%     Mat2VTK(fname,slice{k},'binary');
end

%% Disaplacement field



%% functions
function EPS_r = EPS(rxl,ryl,theta,kx,ky)
part1 = (cos(theta))^2-(sin(theta))^2*cos(kx*rxl);
part2 = (cos(theta))^2-(sin(theta))^2*cos(ky*ryl);
EPS_r = part1.*part2;
end
