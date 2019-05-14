clc;clear;
%% Model Shape
lambda_i = 0.35;% inner radius
lambda_o = 0.55;% outer radius
delta = 4; % focal radius, cm
Nphi = 100;
phi = linspace(0,360,Nphi)/180*pi;
Nl = 50;
lambda = linspace(lambda_i,lambda_o,Nl);
Ni = 100;
ita = linspace(0,120,Ni)/180*pi;
[PHI,LAMBDA,ITA] = ndgrid(phi,lambda,ita);
x = delta*sinh(LAMBDA).*sin(ITA).*cos(PHI);
y = delta*sinh(LAMBDA).*sin(ITA).*sin(PHI);
z = delta*cosh(LAMBDA).*cos(ITA);
% figure;scatter3(x(:),y(:),z(:),'.k');axis image;
% K = boundary(x(:),y(:),z(:),0.9);
% trisurf(K,x(:),y(:),z(:),'Facecolor','red','Facealpha',0.1,'Edgecolor','none');axis image;
% xlabel('x (cm)');ylabel('y (cm)');zlabel('z (cm)');
% inner layer
[PHI_i,LAMBDA_i,ITA_i] = ndgrid(phi,lambda_i,ita);
xi = delta*sinh(LAMBDA_i).*sin(ITA_i).*cos(PHI_i);
yi = delta*sinh(LAMBDA_i).*sin(ITA_i).*sin(PHI_i);
zi = delta*cosh(LAMBDA_i).*cos(ITA_i);
% outer layer
[PHI_o,LAMBDA_o,ITA_o] = ndgrid(phi,lambda_o,ita);
xo = delta*sinh(LAMBDA_o).*sin(ITA_o).*cos(PHI_o);
yo = delta*sinh(LAMBDA_o).*sin(ITA_o).*sin(PHI_o);
zo = delta*cosh(LAMBDA_o).*cos(ITA_o);

%% check lambda_INV to get mask 
r1_INV = sqrt(x.^2+y.^2+(z+delta).^2);
r2_INV = sqrt(x.^2+y.^2+(z-delta).^2);
lambda_INV = acosh((r1_INV+r2_INV)/(2*delta));
ita_INV = acos((r1_INV-r2_INV)/(2*delta));
phi_INV = atan(y./x);
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
  EPS_r0 = funEPS(rxl,ryl,theta,kx,ky);
  mu(:,:,layer) = D0*exp(-TE/T2)*(1+((1-exp(-(TR-Td)/T1))*EPS_r0-1)*exp(-Td/T2));
end
mask = zeros(size(x));
mask(lambda_INV>=lambda_i & lambda_INV<=lambda_o & ita_INV>=0 & ita_INV<=(120/180*pi))=1;
% figure;scatter3(rx(mask>0),ry(mask>0),rz(mask>0),1,mu(mask>0));axis image;
% xlabel('x (cm)');ylabel('y (cm)');zlabel('z (cm)');
%% Motion Model
a = mean(cosh(lambda)./sinh(lambda)); % correcitonal parameter
% the wall volume of the model LV
Vw = pi*delta^3/4*(3*(cosh(lambda_o)-cosh(lambda_i))+4*((cosh(lambda_o))^3-(cosh(lambda_i))^3));
fpath = 'Data\';
kpath = [fpath,'ks_csv'];
N = 60;
Ks = convert13ks(kpath,N);
mxc = x(mask>0);
px = mxc(:);
myc = y(mask>0);
py = myc(:);
mzc = z(mask>0);
pz = mzc(:);
new_x = cell(N,1);
new_y = cell(N,1);
new_z = cell(N,1);
new_mu = mu(mask>0);
parfor k = 1:N
    [new_x{k},new_y{k},new_z{k}] = transK(a,Ks(:,k),px,py,pz,Vw);
end
xmin = min(cell2mat(new_x));
xmax = max(cell2mat(new_x));
ymin = min(cell2mat(new_y));
ymax = max(cell2mat(new_y));
zmin = min(cell2mat(new_z));
zmax = max(cell2mat(new_z));

%% Write VTK to Paraview
for k = 1:N
    scatter3(new_x{k},new_y{k},new_z{k},2,new_mu,'filled');
    axis image;
    set(gca,'Zdir','reverse');
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    zlim([zmin zmax]);
    pause(0.1);
    filename = sprintf('%s%sslice_3d%02.0f',fpath,filesep,k);
    p = [new_x{k},new_y{k},new_z{k}];
    t = delaunayn(unique(p));
    u = new_mu;
    writeVTK(filename,t,p,u);
end
