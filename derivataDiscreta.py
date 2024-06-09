from scipy import signal
import numpy as np
import matplotlib.pyplot as plt


# Tempo di campionamento del segnale
Ts = 0.01
# Definisci la tua funzione di trasferimento continua
num = [1.0]  # Numeratore
den = [Ts/1000, 1.0]  # Denominatore
tc = signal.TransferFunction(num, den)

# Discretizza la funzione di trasferimento con il metodo di Tustin
dt = 0.01  # Intervallo di tempo di campionamento
num_d, den_d, dt = signal.cont2discrete((tc.num, tc.den), dt, method='bilinear')
print(f"Numeratore: ", num_d)
print(f"Denominatore: ", den_d)

#Verifica dei risultati
t_max = 10    # Tempo massimo di simulazione
time = np.linspace(0,t_max,1000)
y = np.zeros_like(time)
dy = np.zeros_like(time)
dy_diff = np.zeros_like(time)
dy_tust = np.zeros_like(time)
for i in range(1000):
    y[i] = 5*np.sin(Ts*i)
    dy[i] = 5*Ts*np.cos(Ts*i)
    if i != 0:
        dy_diff[i]=y[i]-y[i-1]
        dy_tust[i] = num_d[0,0]*(y[i]-y[i-1]) - den_d[1]*dy_tust[i-1]

plt.figure(1)
plt.plot(time, dy ,'r', time, dy_diff, 'b--', time, dy_tust, 'g--')
plt.grid()
plt.show()
plt.figure(2)
plt.plot(time, dy_diff, 'b--', time, dy_tust, 'g--')
plt.show()
