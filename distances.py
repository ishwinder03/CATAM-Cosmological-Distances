import numpy as np
import matplotlib.pyplot as plt

# Constants
Mpc = 3.0856e22 # m
c = 2.998e8

H_0 = 72e3 / Mpc # s^-1
t_H = 1 / H_0 
D_H = c / H_0 # Hubble distance

def Omega_k(O_m, O_l):
    return 1 - O_m - O_l

def E(z):
    return np.sqrt( (Omega_m * (1+z)**3) + (Omega_k(Omega_m, Omega_l) * (1+z)**2) + Omega_l)

def distance_func(z):
    return D_H * (1 / E(z))

def midpoint_rule(func, a, b, n):
    h = (b - a) / n
    integral = 0
    for i in range(n):
        midpoint = a + (i + 0.5) * h
        integral += func(midpoint)
    integral *= h
    return integral
    
def distance_integral(z):  # this is D_C 
    return midpoint_rule(distance_func, 0, z, 1000)

def dimensionless_angular(z): 
    O_k =  Omega_k(Omega_m,Omega_l)
    if O_k > 0:
        da_over_dh = (1 / (np.sqrt(O_k) * (1+z))) * np.sinh(np.sqrt(O_k) * distance_integral(z) / D_H)
    
    elif O_k == 0:
        da_over_dh = distance_integral(z) / (D_H * (1+z))
    
    elif O_k < 0:
        da_over_dh = (1 / (np.sqrt(abs(O_k)) * (1+z))) * np.sin(np.sqrt(abs(O_k)) * distance_integral(z) / D_H)
    
    else:
        print('error')
    return da_over_dh

def dimensionless_luminosity(z):
    DL_over_DH = (1 + z)**2 * dimensionless_angular(z)
    return DL_over_DH

Z = np.linspace(0,7,10000) # upper limits of integral 'z', needs to be 0<z<7

Omega_m_values = [1, 0.04, 0.27]
Omega_l_values = [0, 0, 0.73]

ang_distances = []
lum_distances = []

fig, ax1 = plt.subplots(figsize=(6,5))
fig, ax2 = plt.subplots(figsize=(6,5))
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

for i in range(0,3):
    Omega_m = Omega_m_values[i]
    Omega_l = Omega_l_values[i]
    DA_over_DH = dimensionless_angular(Z)
    DL_over_DH = dimensionless_luminosity(Z)
    
    ang_distances.append(DA_over_DH)
    lum_distances.append(DL_over_DH)

    ax1.plot(Z, DA_over_DH,color = colors[i], linestyle = '-', label=f'$\Omega_m={Omega_m}$, $\Omega_\Lambda={Omega_l}$')
    ax2.plot(Z, DL_over_DH, color = colors[i], linestyle = '-', label=f'$\Omega_m={Omega_m}$, $\Omega_\Lambda={Omega_l}$')

ax1.set_xlabel('z', fontsize = 14)
ax2.set_xlabel('z', fontsize = 14)

ax1.set_ylabel('$D_A \, / \, D_H$', fontsize = 14)
ax2.set_ylabel('$D_L \, / \, D_H$', fontsize = 14)

ax1.set_title('Dimensionless Angular Diameter Distance for Different Universes')
ax2.set_title('Dimensionless Luminosity Distance for Different Universes')

ax1.legend()
ax2.legend()

# ax1.get_figure().savefig('Angular Distance Plot.png', dpi=300)
# ax2.get_figure().savefig('Luminosity Distance Plot.png', dpi=300)

plt.show()

ang_distances_1, ang_distances_2, ang_distances_3 = ang_distances
lum_distances_1, lum_distances_2, lum_distances_3 = lum_distances

all_distances = (ang_distances_1, ang_distances_2, ang_distances_3, lum_distances_1, lum_distances_2, lum_distances_3)

Z_values = [1, 1.25, 2.0, 4.0]

for l,k in enumerate(all_distances):
    for z_val in Z_values:
        print(f'Distance {l+1} for z = {z_val}')
        
        for i,j in enumerate(Z):
            if round(j,3) == z_val:
                print(f'{k[i]:.2f}')
        print('\n')
    print('-------------------------------')