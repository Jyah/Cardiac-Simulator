function EPS_r = funEPS(rxl,ryl,theta,kx,ky)
part1 = (cos(theta))^2-(sin(theta))^2*cos(kx*rxl);
part2 = (cos(theta))^2-(sin(theta))^2*cos(ky*ryl);
EPS_r = part1.*part2;
end
