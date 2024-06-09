import numpy as np
import matplotlib.pyplot as plt

# Parametri di comando
np.random.seed(0)  # Imposta il seed per la riproducibilità
vxcomm = np.random.rand()                                 # [m/s] velocità longitudinale desiderata
yawcomm = (2 * np.random.rand() - 1) * np.pi              # [rad] yaw desiderato in [-pi; pi]
dt = 0.005                                                # [s] passo d'integrazione
yawratecomm = (2 * np.random.rand() - 1) * np.pi / 6 / dt  # [rad/s] yaw rate desiderato in [-pi/6; pi/6]

print("v_x desiderata:", vxcomm)
print("yaw desiderato:", yawcomm)
print("yaw_rate desiderato:", yawratecomm)

# Scale factor della reward
linearVelocityXYRewardScale = 1.0
varlinvel = 0.03
linearVelocityZRewardScale = -4.0
angularVelocityXYRewardScale = -0.05
angularVelocityZRewardScale = 0.5
varangvel = 2
orientationRewardScale = -0.
torqueRewardScale = -0.00002
jointAccRewardScale = -0.0005
baseHeightRewardScale = -500.0
varbase = -1 / baseHeightRewardScale
baseHeightRewardScale1 = 0.3
actionRateRewardScale = -0.01
fallenOverRewardScale = -1.0

# Parametri di sistema
vxlim = 1                   # [m/s] massima velocità longitudinale
yawmax = np.pi / 6          # [rad] massimo valore di yaw accettato
yawratelim = 1              # [rad/s] limite per lo yawrate
taumax = 5                  # [Nm] coppia massima generabile dai motori
initbaseheight = 0.35       # [m] altezza iniziale base
rollmax = np.pi / 6         # [rad] massimo rollio ammissibile
pitchmax = np.pi / 6        # [rad] massimo beccheggio ammissibile
kneeaccelmax = 2 * np.pi / 3 / dt  # [rad/s] massima velocità angolare relativa per knee tra 2 istanti successivi
hipaccelmax = np.pi / 3 / dt        # [rad/s] massima velocità angolare relativa per hip tra 2 istanti successivi
actionkneemax = 2 * np.pi / 3      # [rad] massima differenza tra action relative al knee tra 2 istanti successivi
actionhipmax = np.pi / 3           # [rad] massima differenza tra action relative al hip tra 2 istanti successivi
numrobots = 4096                    # [adim] numero robot simulati in contemporanea
desbaseheight = 0.25                # [m] altezza da terra desiderata per la base

# Variabili di stato / osservazioni
vx = np.linspace(-vxlim, vxlim)                            # [m/s] velocità attuale longitudinale robot
yawrate = np.linspace(-yawmax / dt, yawmax / dt, int(2 / dt))  # [rad/s] yaw rate attuale
tauknee = np.linspace(-taumax, taumax)                     # [Nm] coppia sul ginnocchio attuale
tauhip = np.linspace(-taumax, taumax)                      # [Nm] coppia sul hip attuale
vz = np.linspace(-initbaseheight / dt, initbaseheight / dt)  # [m/s] velocità di abbassamento della base
rollrate = np.linspace(-rollmax / dt, rollmax / dt)         # [rad/s] roll rate attuale
pitchrate = np.linspace(-pitchmax / dt, pitchmax / dt)      # [rad/s] pitch rate attuale
kneeaccel = np.linspace(-kneeaccelmax, kneeaccelmax)       # [rad/s] accel angolare knee vista come differenza di velocità angolari
hipaccel = np.linspace(-hipaccelmax, hipaccelmax)          # [rad/s] accel angolare hip vista come differenza di velocità angolari
deltactionknee = np.linspace(-actionkneemax, actionkneemax)  # [rad] velocità angolare knee vista come differenza di action
deltactionhip = np.linspace(-actionhipmax, actionhipmax)    # [rad] velocità angolare hip vista come differenza di action
fallen = np.linspace(0, numrobots)                         # [adim] numero di robot caduti
baseheight = np.linspace(-initbaseheight, initbaseheight)  # [m] altezza attuale della base

# Funzioni di reward (singolo robot)
rewlinvelxy = linearVelocityXYRewardScale * np.exp(-(vxcomm - vx) ** 2 / (2 * varlinvel))
rewangvelz = angularVelocityZRewardScale * np.exp(-(yawratecomm - yawrate) ** 2 / (2 * varangvel))
rewtorque = torqueRewardScale * (tauknee[:, None] ** 2 + tauhip ** 2)
rewlinvelz = linearVelocityZRewardScale * vz ** 2
rewbaseheight = baseHeightRewardScale * (baseheight - desbaseheight) ** 2

# Plot
plt.figure(figsize=(10, 6))

plt.subplot(2, 3, 1)
plt.plot(vx, rewlinvelxy, 'b')
plt.axvline(x=vxcomm, color='r', linestyle='--')
plt.grid(True)
plt.title('Reward velocità longitudinale')
plt.xlabel('vx [m/s]')
plt.ylabel('rew [adim]')
plt.ylim(0, 1.2 * linearVelocityXYRewardScale)

plt.subplot(2, 3, 2)
plt.plot(yawrate, rewangvelz, 'b')
plt.axvline(x=yawratecomm, color='r', linestyle='--')
plt.grid(True)
plt.title('Reward yaw rate')
plt.xlabel('yawrate [rad/s]')
plt.ylabel('rew [adim]')
plt.ylim(0, 1.1 * angularVelocityZRewardScale)

plt.subplot(2, 3, 3)
plt.plot(vz, rewlinvelz, 'b')
plt.axvline(x=0, color='r', linestyle='--')
plt.grid(True)
plt.title('Reward velocità verticale')
plt.xlabel('vz [m/s]')
plt.ylabel('rew [adim]')

plt.subplot(2, 3, 4)
plt.plot(fallen, fallenOverRewardScale * fallen)
plt.grid(True)
plt.title('Reward fallen over')
plt.xlabel('fallen robots [adim]')
plt.ylabel('rew [adim]')

plt.subplot(2, 3, 5)
plt.plot(baseheight, rewbaseheight, 'b')
plt.axvline(x=desbaseheight, color='r', linestyle='--')
plt.grid(True)
plt.title('Reward altezza base')
plt.xlabel('z [m]')
plt.ylabel('rew [adim]')

plt.tight_layout()
plt.show()
