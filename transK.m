function [new_x,new_y,new_z] = transK(a,k,px,py,pz,Vw)
r = zeros(4,length(px));
F0 = diag([a^(1/3),a^(1/3),a^(-2/3),1]);
F3 = [a^(-1/3)*exp(k(4)-k(3)/2),0,0,0;
  0,a^(-1/3)*exp(-k(4)-k(3)/2),0,0;
  0,0,a^(2/3)*exp(k(3)),0;
  0,0,0,1];
F4 = [1,k(5),0,0;
  k(5),1+(k(5))^2,0,0;
  0,0,1,0;
  0,0,0,1];
F5 = [1,0,k(6),0;
  0,1,0,0;
  k(6),0,1+(k(6))^2,0;
  0,0,0,1];
F6 = [1,0,0,0;
  0,1,k(7),0;
  0,k(7),1+(k(7))^2,0;
  0,0,0,1];
A1 = [1,0,0,0;
  0,cos(k(8)),-sin(k(8)),0;
  0,sin(k(8)),cos(k(8)),0;
  0,0,0,1];
A2 = [cos(k(9)),0,sin(k(9)),0;
  0,1,0,0;
  -sin(k(9)),0,cos(k(9)),0;
  0,0,0,1];
A3 = [cos(k(10)),-sin(k(10)),0,0;
  sin(k(10)),cos(k(10)),0,0;
  0,0,1,0;
  0,0,0,1];
A4 = [1,0,0,k(11);
  0,1,0,k(12);
  0,0,1,k(13);
  0,0,0,1];
Fa = A4*A3*A2*A1;

for iP = 1:length(px)
  P = [px(iP);py(iP);pz(iP);1];
  eps = (1+3*k(1)*Vw/(4*pi*(norm(F0*P)).^3)).^(1/3);
  F1 = diag([eps,eps,eps,1]);
  r1 = F1*F0*P;% =[x1;y1;z1;1];
  z1 = r1(3,1);
  F2 = [cos(a*k(2)*z1/norm(r1)),-sin(a*k(2)*z1/norm(r1)),0,0;
    sin(a*k(2)*z1/norm(r1)),cos(a*k(2)*z1/norm(r1)),0,0;
    0,0,1,0;
    0,0,0,1];
  r(:,iP) = Fa*F6*F5*F4*F3*F2*F1*F0*P;%[rx,ry,rz,1];
end
new_x = reshape(r(1,:),size(px));
new_y = reshape(r(2,:),size(py));
new_z = reshape(r(3,:),size(pz));
end
