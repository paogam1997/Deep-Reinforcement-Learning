function info = CheckCollisionRandOrient(roll_deg, pitch_deg, yaw_deg, CI_z, toll)

phi = roll_deg*pi/180;
theta = pitch_deg*pi/180;
psi = yaw_deg*pi/180;

C_b_i = [ cos(theta)*cos(psi), cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi), cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi);
          cos(theta)*sin(psi), sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi), sin(psi)*sin(theta)*cos(phi)-sin(phi)*cos(psi);
          -sin(theta),          cos(theta)*sin(phi),                            cos(theta)*cos(phi)];

OAb = [0.235;   0.125; 0.316];
OBb = [-0.235;  0.125; 0.316];
OCb = [-0.235; -0.125; 0.316];
ODb = [0.235;  -0.125; 0.316];

OAi = C_b_i * OAb;
OBi = C_b_i * OBb;
OCi = C_b_i * OCb;
ODi = C_b_i * ODb;

if OAi(3) + toll >= CI_z
    info = "!! COLLISIONE SU ANTERIORE DESTRO!!";
elseif OBi(3) + toll >= CI_z
    info = "!! COLLISIONE SU POSTERIORE DESTRO !!";
elseif OCi(3)+ toll >= CI_z
    info = "!! COLLISIONE SU POSTERIORE SINISTRO !!";
elseif ODi(3)+ toll >= CI_z
    info = "!! COLLISIONE SU ANTERIORE SINISTRO !!";
else
    info = " OK randomizzazione angoli Eulero non porta a collisione ";
end


end