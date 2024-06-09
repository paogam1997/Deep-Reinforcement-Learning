function rpy_deg = eul2rpy(z_deg, y_deg, z2_deg)

eul = [z_deg, y_deg, z2_deg] .* (pi/180);
R = eul2rotm(eul, "ZYZ");
rpy_deg = rotm2eul(R, "ZYX") .* (180/pi);