import torch
import numpy as np
import matplotlib.pyplot as plt


# Definisci la funzione FadingFilter
def FadingFilter(beta, previous_x_hat, measured_x, filter_order, Ts=0):
    if filter_order == 1:
        # Ridifinisci gli input per semplicità
        x_hat_p = previous_x_hat
        x_star = measured_x
        # Definisci il parametro del filtro
        G = 1 - beta
        # Applica il Fading Filter
        res = x_star - x_hat_p
        x_hat = x_hat_p + G * res
        dx_hat = torch.tensor([0.0])
    elif filter_order == 2:
        # Ridifinisci gli input per semplicità
        x_hat_p = previous_x_hat[0]
        dx_hat_p = previous_x_hat[1]
        x_star = measured_x[0]
        dx_star = measured_x[1]
        # Definisci i parametri del filtro
        G  = 1 - beta**2
        H = (1 - beta)**2
        # Applica il Fading Filter di 2° ordine
        res = x_star - (x_hat_p + dx_hat_p * Ts)
        x_hat = (x_hat_p + dx_hat_p * Ts) + G * res
        dx_hat = (dx_hat_p) + H/Ts * res

    return x_hat, dx_hat

# Definisci i parametri
t = torch.linspace(0, 3, 500)
hip_0 = 15 * np.pi / 180
knee_0 = 10 * np.pi / 180
Hip1 = 15 * np.pi / 180
Hip2 = 5 * np.pi / 180
Hip3 = 1 * np.pi / 180
w_hip1 = 2 * np.pi / 2
w_hip2 = 2 * np.pi / 0.5
w_hip3 = 2 * np.pi / 0.2
Knee1 = 20 * np.pi / 180
Knee2 = 15 * np.pi / 180
Knee3 = 3 * np.pi / 180
w_knee1 = 2 * np.pi / 1
w_knee2 = 2 * np.pi / 0.3
w_knee3 = 2 * np.pi / 0.1

# Calcola il movimento
hip = hip_0 + Hip1 * torch.sin(w_hip1 * t) + Hip2 * torch.sin(w_hip2 * t) + Hip3 * torch.sin(w_hip3 * t)
knee = knee_0 + Knee1 * torch.sin(w_knee1 * t) + Knee2 * torch.sin(w_knee2 * t) + Knee3 * torch.sin(w_knee3 * t)
theta = torch.stack((hip, knee))
dhip = Hip1 * w_hip1 * torch.cos(w_hip1 * t) + Hip2 * w_hip2 * torch.cos(w_hip2 * t) + Hip3 * w_hip3 * torch.cos(w_hip3 * t)
dknee = Knee1 * w_knee1 * torch.cos(w_knee1 * t) + Knee2 * w_knee2 * torch.cos(w_knee2 * t) + Knee3 * w_knee3 * torch.cos(w_knee3 * t)
omega = torch.stack((dhip, dknee))

# Aggiungi rumore
devstd_hip = 0.25
devstd_knee = 0.25
devstd_dhip = 2.5
devstd_dknee = 2.5
noise_hip = devstd_hip * torch.randn_like(hip)
noise_knee = devstd_knee * torch.randn_like(knee)
noise_dhip = devstd_dhip * torch.randn_like(dhip)
noise_dknee = devstd_dknee * torch.randn_like(dknee)
policy_hip = hip + noise_hip
policy_knee = knee + noise_knee
policy_dhip = dhip + noise_dhip
policy_dknee = dknee + noise_dknee

# Applica il Fading Filter
beta_hip = 0.9
beta_hip2 = 0.95
beta_hip2inv = 0.8
beta_knee = 0.9
beta_knee2 = 0.95
beta_knee2inv = 0.8
beta_dhip = 0.9
beta_dknee = 0.9
hip_hat = torch.zeros_like(t)
hip_hat2 = torch.zeros_like(t)
hip_hat2inv = torch.zeros_like(t)
knee_hat = torch.zeros_like(t)
knee_hat2 = torch.zeros_like(t)
knee_hat2inv = torch.zeros_like(t)
dhip_hat = torch.zeros_like(t)
dhip_hat2 = torch.zeros_like(t)
dhip_hat2inv = torch.zeros_like(t)
dknee_hat = torch.zeros_like(t)
dknee_hat2 = torch.zeros_like(t)
dknee_hat2inv = torch.zeros_like(t)

for i in range(1, len(t)):
    # 1° ordine
    hip_hat[i], _ = FadingFilter(beta_hip, hip_hat[i-1], policy_hip[i], 1)
    knee_hat[i], _ = FadingFilter(beta_knee, knee_hat[i-1], policy_knee[i], 1)
    dhip_hat[i], _ = FadingFilter(beta_dhip, dhip_hat[i-1], policy_dhip[i], 1)
    dknee_hat[i], _ = FadingFilter(beta_dknee, dknee_hat[i-1], policy_dknee[i], 1)
    # 2° ordine
    hip_hat2[i], dhip_hat2[i] = FadingFilter(beta_hip2, torch.tensor([hip_hat2[i-1], dhip_hat2[i-1]]), torch.tensor([policy_hip[i], policy_dhip[i]]), 2, 0.01)
    knee_hat2[i], dknee_hat2[i] = FadingFilter(beta_knee2, torch.tensor([knee_hat2[i-1], dknee_hat2[i-1]]), torch.tensor([policy_knee[i], policy_dknee[i]]), 2, 0.01)


# Create the figure and subplots
fig, axs = plt.subplots(2, 2, figsize=(10, 10))
fig.suptitle("Hip & Knee joint Angles & Velocities")

# Plot the hip angle
axs[0, 0].plot(t, policy_hip, color='#0072BD')
axs[0, 0].plot(t, hip, color='#D95319', linestyle='--', linewidth=1.5)
axs[0, 0].plot(t, hip_hat, color='#77AC30', linewidth=1)
axs[0, 0].plot(t, hip_hat2, color='#7E2F8E', linestyle='-.', linewidth=0.75)
axs[0, 0].grid(True)
axs[0, 0].set_title("Hip Angle")
axs[0, 0].set_xlabel("time [s]")
axs[0, 0].set_ylabel("Joint angle [rad]")
axs[0, 0].legend(["Policy", "Ideal", "FF 1^{st} order", "FF 2^{nd} order"])

# Plot the knee angle
axs[1, 0].plot(t, policy_knee, color='#0072BD')
axs[1, 0].plot(t, knee, color='#D95319', linestyle='--', linewidth=1.5)
axs[1, 0].plot(t, knee_hat, color='#77AC30', linewidth=1)
axs[1, 0].plot(t, knee_hat2, color='#7E2F8E', linestyle='-.', linewidth=0.75)
axs[1, 0].grid(True)
axs[1, 0].set_title("Knee Angle")
axs[1, 0].set_xlabel("time [s]")
axs[1, 0].set_ylabel("Joint angle [rad]")
axs[1, 0].legend(["Policy", "Ideal", "FF 1^{st} order", "FF 2^{nd} order"])

# Plot the hip velocity
axs[0, 1].plot(t, policy_dhip, color='#0072BD')
axs[0, 1].plot(t, dhip, color='#D95319', linestyle='--', linewidth=1.5)
axs[0, 1].plot(t, dhip_hat, color='#77AC30', linewidth=1)
axs[0, 1].plot(t, dhip_hat2, color='#7E2F8E', linestyle='-.', linewidth=0.75)
axs[0, 1].grid(True)
axs[0, 1].set_title("Hip Velocity")
axs[0, 1].set_xlabel("time [s]")
axs[0, 1].set_ylabel("Joint velocity [rad/s]")
axs[0, 1].legend(["Policy", "Ideal", "FF 1^{st} order", "FF 2^{nd} order"])

# Plot the knee velocity
axs[1, 1].plot(t, policy_dknee, color='#0072BD')
axs[1, 1].plot(t, dknee, color='#D95319', linestyle='--', linewidth=1.5)
axs[1, 1].plot(t, dknee_hat, color='#77AC30', linewidth=1)
axs[1, 1].plot(t, dknee_hat2, color='#7E2F8E', linestyle='-.', linewidth=0.75)
axs[1, 1].grid(True)
axs[1, 1].set_title("Knee Velocity")
axs[1, 1].set_xlabel("time [s]")
axs[1, 1].set_ylabel("Joint velocity [rad/s]")
axs[1, 1].legend(["Policy", "Ideal", "FF 1^{st} order", "FF 2^{nd} order"])

# Show the plot
plt.show()

# Create the figure for the estimation errors
fig2, axs2 = plt.subplots(2, 2, figsize=(10, 10))
fig2.suptitle("Estimation Errors")

# Plot the hip angle error
axs2[0, 0].plot(t, (hip_hat - hip) * 180/np.pi, color='#77AC30')
axs2[0, 0].plot(t, (hip_hat2 - hip) * 180/np.pi, color='#7E2F8E')
axs2[0, 0].grid(True)
axs2[0, 0].set_title("Hip Angle")
axs2[0, 0].set_xlabel("time [s]")
axs2[0, 0].set_ylabel("Joint angle error [deg]")
axs2[0, 0].legend(["FF 1^{st} order", "FF 2^{nd} order"])

# Plot the knee angle error
axs2[1, 0].plot(t, (knee_hat - knee) * 180/np.pi, color='#77AC30')
axs2[1, 0].plot(t, (knee_hat2 - knee) * 180/np.pi, color='#7E2F8E')
axs2[1, 0].grid(True)
axs2[1, 0].set_title("Knee Angle")
axs2[1, 0].set_xlabel("time [s]")
axs2[1, 0].set_ylabel("Joint angle error [deg]")
axs2[1, 0].legend(["FF 1^{st} order", "FF 2^{nd} order"])

# Plot the hip velocity error
axs2[0, 1].plot(t, (dhip_hat - dhip) * 180/np.pi, color='#77AC30')
axs2[0, 1].plot(t, (dhip_hat2 - dhip) * 180/np.pi, color='#7E2F8E')
axs2[0, 1].grid(True)
axs2[0, 1].set_title("Hip Velocity")
axs2[0, 1].set_xlabel("time [s]")
axs2[0, 1].set_ylabel("Joint velocity error [deg/s]")
axs2[0, 1].legend(["FF 1^{st} order", "FF 2^{nd} order"])

# Plot the knee velocity error
axs2[1, 1].plot(t, (dknee_hat - dknee) * 180/np.pi, color='#77AC30')
axs2[1, 1].plot(t, (dknee_hat2 - dknee) * 180/np.pi, color='#7E2F8E')
axs2[1, 1].grid(True)
axs2[1, 1].set_title("Knee Velocity")
axs2[1, 1].set_xlabel("time [s]")
axs2[1, 1].set_ylabel("Joint velocity error [deg/s]")
axs2[1, 1].legend(["FF 1^{st} order", "FF 2^{nd} order"])

# Show the plot
plt.show()

# # Create the figure and subplots
# fig, axs = plt.subplots(2, 2, figsize=(10, 10))
# fig.suptitle("Medie & Deviazioni Standard")
#
# # Plot the signal mean
# axs[0, 0].bar(['Hip Angle', 'Knee Angle', 'Hip Vel', 'Knee Vel'],
#               [torch.mean(hip)*(180/np.pi), torch.mean(policy_hip)*(180/np.pi), torch.mean(hip_hat)*(180/np.pi), torch.mean(hip_hat2)*(180/np.pi),
#                torch.mean(knee)*(180/np.pi), torch.mean(policy_knee)*(180/np.pi), torch.mean(knee_hat)*(180/np.pi), torch.mean(knee_hat2)*(180/np.pi),
#                torch.mean(dhip)*(180/np.pi), torch.mean(policy_dhip)*(180/np.pi), torch.mean(dhip_hat)*(180/np.pi), torch.mean(dhip_hat2)*(180/np.pi),
#                torch.mean(dknee)*(180/np.pi), torch.mean(policy_dknee)*(180/np.pi), torch.mean(dknee_hat)*(180/np.pi), torch.mean(dknee_hat2)*(180/np.pi)])
# axs[0, 0].set_title("Signal Mean")
# axs[0, 0].legend(['Ideale', 'Policy', 'FF 1^{st} ord','FF 2^{nd} ord'])
# axs[0, 0].grid(True)
#
# # Plot the signal standard deviation
# axs[0, 1].bar(['Hip Angle', 'Knee Angle', 'Hip Vel', 'Knee Vel'],
#               [torch.std(hip)*(180/np.pi), torch.std(policy_hip)*(180/np.pi), torch.std(hip_hat)*(180/np.pi), torch.std(hip_hat2)*(180/np.pi),
#                torch.std(knee)*(180/np.pi), torch.std(policy_knee)*(180/np.pi), torch.std(knee_hat)*(180/np.pi), torch.std(knee_hat2)*(180/np.pi),
#                torch.std(dhip)*(180/np.pi), torch.std(policy_dhip)*(180/np.pi), torch.std(dhip_hat)*(180/np.pi), torch.std(dhip_hat2)*(180/np.pi),
#                torch.std(dknee)*(180/np.pi), torch.std(policy_dknee)*(180/np.pi), torch.std(dknee_hat)*(180/np.pi), torch.std(dknee_hat2)*(180/np.pi)])
# axs[0, 1].set_title("Signal Standard Deviation")
# axs[0, 1].legend(['Ideale', 'Policy', 'FF 1^{st} ord','FF 2^{nd} ord'])
# axs[0, 1].grid(True)
#
# # Plot the error mean
# axs[1, 0].bar(['Hip Angle', 'Knee Angle', 'Hip Vel', 'Knee Vel'],
#               [0, torch.mean(policy_hip - hip)*(180/np.pi), torch.mean(hip_hat - hip)*(180/np.pi), torch.mean(hip_hat2 - hip)*(180/np.pi),
#                0, torch.mean(policy_knee - knee)*(180/np.pi), torch.mean(knee_hat - knee)*(180/np.pi), torch.mean(knee_hat2 - knee)*(180/np.pi),
#                0, torch.mean(policy_dhip - dhip)*(180/np.pi), torch.mean(dhip_hat - dhip)*(180/np.pi), torch.mean(dhip_hat2 - dhip)*(180/np.pi),
#                0, torch.mean(policy_dknee - dknee)*(180/np.pi), torch.mean(dknee_hat - dknee)*(180/np.pi), torch.mean(dknee_hat2 - dknee)*(180/np.pi)])
# axs[1, 0].set_title("Error Mean")
# axs[1, 0].legend(['Ideale', 'Policy', 'FF 1^{st} ord','FF 2^{nd} ord'])
# axs[1, 0].grid(True)
#
# # Plot the error standard deviation
# axs[1, 1].bar(['Hip Angle', 'Knee Angle', 'Hip Vel', 'Knee Vel'],
#               [0, torch.std(policy_hip - hip)*(180/np.pi), torch.std(hip_hat - hip)*(180/np.pi), torch.std(hip_hat2 - hip)*(180/np.pi),
#                0, torch.std(policy_knee - knee)*(180/np.pi), torch.std(knee_hat - knee)*(180/np.pi), torch.std(knee_hat2 - knee)*(180/np.pi),
#                0, torch.std(policy_dhip - dhip)*(180/np.pi), torch.std(dhip_hat - dhip)*(180/np.pi), torch.std(dhip_hat2 - dhip)*(180/np.pi),
#                0, torch.std(policy_dknee - dknee)*(180/np.pi), torch.std(dknee_hat - dknee)*(180/np.pi), torch.std(dknee_hat2 - dknee)*(180/np.pi)])
# axs[1, 1].set_title("Error Standard Deviation")
# axs[1, 1].legend(['Ideale', 'Policy', 'FF 1^{st} ord','FF 2^{nd} ord'])
# axs[1, 1].grid(True)
#
# # Show the plot
# plt.show()
