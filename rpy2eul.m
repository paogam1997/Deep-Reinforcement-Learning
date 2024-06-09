function euler = rpy2eul(r_deg, p_deg, y_deg)

zyx = [r_deg, p_deg, y_deg] .* (pi/180);
R = eul2rotm(zyx, "ZYX");
euler = rotm2eul(R, "ZYZ") .* (180/pi);