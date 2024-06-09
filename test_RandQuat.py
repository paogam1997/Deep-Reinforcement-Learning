import torch
import math

def degrees_to_radians(degrees):
    return degrees * math.pi / 180

def quaternion_to_euler_angle(q):
    w, x, y, z = q
    ysqr = y * y

    t0 = +2.0 * (w * x + y * z)
    t1 = +1.0 - 2.0 * (x * x + ysqr)
    roll_x = math.atan2(t0, t1)

    t2 = +2.0 * (w * y - z * x)
    t2 = +1.0 if t2 > +1.0 else t2
    t2 = -1.0 if t2 < -1.0 else t2
    pitch_y = math.asin(t2)

    t3 = +2.0 * (w * z + x * y)
    t4 = +1.0 - 2.0 * (ysqr + z * z)
    yaw_z = math.atan2(t3, t4)

    return roll_x, pitch_y, yaw_z # in radians

def euler_angle_to_quaternion(roll, pitch, yaw):
    qx = torch.sin(roll/2) * torch.cos(pitch/2) * torch.cos(yaw/2) - torch.cos(roll/2) * torch.sin(pitch/2) * torch.sin(yaw/2)
    qy = torch.cos(roll/2) * torch.sin(pitch/2) * torch.cos(yaw/2) + torch.sin(roll/2) * torch.cos(pitch/2) * torch.sin(yaw/2)
    qz = torch.cos(roll/2) * torch.cos(pitch/2) * torch.sin(yaw/2) - torch.sin(roll/2) * torch.sin(pitch/2) * torch.cos(yaw/2)
    qw = torch.cos(roll/2) * torch.cos(pitch/2) * torch.cos(yaw/2) + torch.sin(roll/2) * torch.sin(pitch/2) * torch.sin(yaw/2)
    return torch.tensor([qw, qx, qy, qz])

def torch_rand_quaternion(roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, env_ids):
    random_values = torch.rand(env_ids, 3)
    random_values[:, 0] = random_values[:, 0] * (degrees_to_radians(roll_max - roll_min)) + degrees_to_radians(roll_min)
    random_values[:, 1] = random_values[:, 1] * (degrees_to_radians(pitch_max - pitch_min)) + degrees_to_radians(pitch_min)
    random_values[:, 2] = random_values[:, 2] * (degrees_to_radians(yaw_max - yaw_min)) + degrees_to_radians(yaw_min)

    random_quaternions = torch.zeros(env_ids, 4)
    for i in range(env_ids):
        random_quaternions[i] = euler_angle_to_quaternion(random_values[i, 0], random_values[i, 1], random_values[i, 2])

    return random_quaternions

# Esempio di utilizzo
env_ids = 4096  # Numero di ambienti
roll_min = -5  # Angolo minimo di roll in gradi
roll_max = 5  # Angolo massimo di roll in gradi
pitch_min = -10  # Angolo minimo di pitch in gradi
pitch_max = 10  # Angolo massimo di pitch in gradi
yaw_min = -15  # Angolo minimo di yaw in gradi
yaw_max = 15  # Angolo massimo di yaw in gradi

# Genera quaternioni casuali nel range specificato
quat_randomization = torch_rand_quaternion(roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, env_ids)

# Verifica se i quaternioni hanno norma unitaria
err_norm = 0
err_range = 0
for j in range(env_ids):
    norm = torch.norm(quat_randomization[j])
    if not math.isclose(norm.item(), 1.0, abs_tol=1e-6):
        print(f"!!!!!Norma q_rand non unitaria ({j}):", norm ,"!!!!!!!!!!")
        err_norm = +1
    # else:
    #     print(f"OK, norma = {norm}")
        # Verifica se i quaternioni generati sono nel range specificato
        for j in range(env_ids):
            roll, pitch, yaw = quaternion_to_euler_angle(quat_randomization[j])
            roll = math.degrees(roll)
            pitch = math.degrees(pitch)
            yaw = math.degrees(yaw)
            if (
                    roll < roll_min or roll > roll_max or pitch < pitch_min or pitch > pitch_max or yaw < yaw_min or yaw > yaw_max):
                print(f"!!!!!!Quaternione fuori dal range ({j}): Roll={roll}, Pitch={pitch}, Yaw={yaw}!!!!!!!!!")
                err_range = +1
            # else:
            #     print(f"OK roll={roll}, pitch={pitch}, yaw={yaw}")
print("-------------------------------------")
print(f"errori in rannge = {err_range}")
print(f"errori nella norma= {err_norm}")