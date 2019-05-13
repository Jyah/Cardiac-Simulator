% function p = invTran(r,k,a)
lyr = 1;
rx = xt(:);
ry = yt(:);
rz = zt(:);
kk = Ks(:,lyr);
F0 = diag([a^(1/3),a^(1/3),a^(-2/3),1]);
Fs = FaF6F5F4F3(kk,a);

%% Iteratively find the reference map
N = length(rx);
p_0 = zeros(N,3);
parfor ir = 1:N
    x0 = [rx(ir),ry(ir),rz(ir)];
    objfun = @(p)fun(p(1),p(2),p(3),x0,Vw,F0,Fs,a,kk);
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
n = 30;
sa_mask = ~isnan(slice{n});
vx = xq(sa_mask>0);
vy = yq(sa_mask>0);
vz = zq(sa_mask>0);
dx = reshape(DX(:,n),size(sa_mask));
dy = reshape(DY(:,n),size(sa_mask));
dz = reshape(DZ(:,n),size(sa_mask));
quiver3(vx(:),vy(:),vz(:),dx(sa_mask>0),dy(sa_mask>0),dz(sa_mask>0),'Color','black');
figure;
quiver(vx(:),vy(:),dx(sa_mask>0),dy(sa_mask>0),'Color','black');
function diF = fun(g1,g2,g3,rc,Vw,F0,Fs,a,k)
rc = [rc,1]';
g = [g1,g2,g3,1]';
eps = (1+3*k(1)*Vw/(4*pi*(norm(F0*g)).^3)).^(1/3);
F1 = diag([eps,eps,eps,1]);
r1 = F1*F0*g;
F2 = [cos(a*k(2)*r1(3)/norm(r1)),-sin(a*k(2)*r1(3)/norm(r1)),0,0;
    sin(a*k(2)*r1(3)/norm(r1)),cos(a*k(2)*r1(3)/norm(r1)),0,0;
    0,0,1,0;
    0,0,0,1];
r2 = F2*F1*F0*g;
R2 = [cos(a*k(2)*r2(3))/norm(r2),sin(a*k(2)*r2(3))/norm(r2),0,0;
    -sin(a*k(2)*r2(3))/norm(r2),cos(a*k(2)*r2(3))/norm(r2),0,0;
    0,0,1,0;
    0,0,0,1];
r0 = F0*g;
norm_r1 = norm(r0)*(1+3*k(1)*Vw/(4*pi*(norm(r1))^3))^(1/3);
epsc = (1-3*k(1)*Vw/(4*pi*(norm_r1)^3))^(1/3);
R1 = diag([epsc,epsc,epsc,1]);
p = inv(F0)*R1*R2*inv(Fs)*rc;
diF = norm(p-g);
end

