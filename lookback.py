import numpy as np
import matplotlib.pyplot as plt

# Constants
Mpc = 3.0856e22 # m
H_0 = 72e3 / Mpc # s^-1
t_H = 1 / H_0 

def Omega_k(O_m, O_l):
    return 1 - O_m - O_l

def E(z):
    return np.sqrt( (Omega_m * (1+z)**3) + (Omega_k(Omega_m, Omega_l) * (1+z)**2) + Omega_l)

def lookback_func(z):
    return t_H * (1 / ((1+z) * E(z)))
    
def lookback_integral(z):
    integral_sum = 0
    for i in range(1, len(z)):
        integral_sum += lookback_func(z[i]) * (z[i] - z[i-1])
    return integral_sum

def lookback_time(z):
    t_L = np.zeros_like(Z)
    for integral_value,i in enumerate(Z):
        t_L[integral_value] = lookback_integral(z)[integral_value]
    return np.array(t_L) / (3.1556926e7 * 1e9) # converting to Gyr

values = np.linspace(0,7,50) # upper limits of integral to give lookback time

b = [0.1, 1.0, 2.0, 4.0, 6.7]
total_arr = np.append(values,b)
Z = np.sort(total_arr)

x = np.linspace(0,Z,2000) # change accuracy via number of steps


Omega_m_values = [1, 2, 0.04, 0.27]
Omega_l_values = [0, 0, 0, 0.73]

universes = ('EDS', 'classical closed', 'baryon', 'currently popular')

all_times = []

plt.figure(figsize = (6,5))
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd']

for i in range(0,4):
    Omega_m = Omega_m_values[i]
    Omega_l = Omega_l_values[i]
    all_times.append(lookback_time(x))
    print(f'{universes[i]}: {lookback_time(x)}')
    print('\n')


    plt.plot(Z, all_times[i], color=colors[i], label=f'$\Omega_m={Omega_m}$, $\Omega_\Lambda={Omega_l}$')


plt.xlabel('z', fontsize = 14)
plt.ylabel('$t_L$', fontsize = 14)
plt.title('Lookback Times for Different Universes')
plt.legend()
# plt.savefig('Lookback Times Plot.png', dpi=300) 
plt.show()

Z_values = [0.1, 1.0, 2.0, 4.0, 6.7]

for l,k in enumerate(all_times):
    for z_val in Z_values:
        print(f'Time {l+1} for z = {z_val}')
        
        for i,j in enumerate(Z):
            if round(j,2) == z_val:
                print(f'{k[i]:.2f}')
        print('\n')
    print('-------------------------------')