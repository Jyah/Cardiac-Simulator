% function p = invTran(r,k,a)
lyr = 1;
r_x = new_x{lyr};
r_y = new_y{lyr};
r_z = new_z{lyr};
kk = Ks(:,lyr);
F0 = diag([a^(1/3),a^(1/3),a^(-2/3),1]);
Fs = FaF6F5F4F3(kk,a);

%% Iteratively find the reference map
N = length(r_x);
p_0 = zeros(N,3);

x0 = [r_x(ir),r_y(ir),r_z(ir)];
objfun = @(p)fun(p(1),p(2),p(3),x0,Vw,F0,Fs,a,kk);
options = optimset('Display','iter','PlotFcns',@optimplotfval);
p_min = fminsearch(objfun, x0,options);
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
