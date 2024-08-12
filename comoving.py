import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import brentq

df = pd.read_csv("/Users/ishu/computing/CATAM/Cosmological Distances/quasar.csv")

df[['z', 'f/f0']] = df['z      f/f0'].str.split(n=1, expand=True) # .dat file had data all together so this separates into actual columns

df['z'] = pd.to_numeric(df['z'])  # convert from string values to float
df['f/f0'] = pd.to_numeric(df['f/f0'])

df.drop(columns=['z      f/f0'], inplace=True)

z_values = df['z'].values
f_over_f0 = df['f/f0'].values

# Constants
Mpc = 3.0856e22  # m
c = 2.998e8  # m/s

H_0 = 72e3 / Mpc  # s^-1
t_H = 1 / H_0 
D_H = c / H_0  # Hubble distance

def E(z):
    return np.sqrt(Omega_m * (1 + z)**3 + Omega_l)

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

def volume(z):
    return (4 * np.pi * distance_integral(z)**3) / 3

V = []
Z = z_values

Omega_m_values = [1, 0.27]
Omega_l_values = [0, 0.73]

for i in range(0, 2):
    Omega_m = Omega_m_values[i]
    Omega_l = Omega_l_values[i]
    V.append([volume(z) for z in Z])

V1,V2 = V
print(V1)
print(V2)

### FINDING Z_MAX
def calc_z_max(D_Lmax):
    def nested_func(z):
        return ((1 + z) * distance_integral(z)) - D_Lmax  #eqn 22 in project but set to 0
    
    z_max, _ = brentq(nested_func, 0, 10, full_output=True)  
    return z_max

Omega_m = 0.27
Omega_l = 0.73
z_maxes = np.zeros_like(z_values)

for i in range(0,114):
    D_L = (1 + z_values[i]) * distance_integral(z_values[i])
    D_Lmax = D_L * np.sqrt(f_over_f0[i])
    z_max = calc_z_max(D_Lmax)
    z_maxes[i] = z_max

V_max2 = volume(z_maxes)
a = V2 / V_max2
print(f'<V/V_max> = {np.mean(a)}')

small_z = []
small_f_f0 = []

for i in range(0,114):
    if z_values[i] < 0.3:
        small_z.append(z_values[i]) 
        small_f_f0.append(f_over_f0[i])

small_z = np.array(small_z)
small_f_f0 = np.array(small_f_f0)

small_z_maxes = np.zeros_like(small_z)

for i in range(0,len(small_z)):
    D_L = (1 + small_z[i]) * distance_integral(small_z[i])
    small_D_Lmax = D_L * np.sqrt(small_f_f0[i])
    small_z_max = calc_z_max(small_D_Lmax)
    small_z_maxes[i] = small_z_max


v2 = volume(small_z)
v_max2 = volume(small_z_maxes)

volume_ratio = v2 / v_max2
f_ratio = small_f_f0**(-3/2)

def theoretical_model(x,alpha, beta):   
    return (alpha * x) + beta

popt, pcov = curve_fit(theoretical_model, f_ratio, volume_ratio)
alpha_fit = popt[0]
beta_fit = popt[1]

plt.figure(figsize=(7,5))
plt.scatter(f_ratio, volume_ratio, alpha = 0.9, marker = 'x',s=70, label = '$V/V_{{max}} \, $ test data')
plt.plot(f_ratio, theoretical_model(f_ratio, alpha_fit, beta_fit), color = 'r', label = f'Linear fit, $V/V_{{max}}$ = ${alpha_fit:.4f} \, (f/f_0)^{{-3/2}} + {beta_fit:.4f}$')
plt.title('Euclidean limit')
plt.xlabel('$(f/f_0)^{{-\\frac{{3}}{{2}}}}$', fontsize = 12)
plt.ylabel('$\\frac{{V}}{{V_{{max}}}}$', fontsize = 14)
plt.legend(fontsize = 11)
# plt.savefig('Euclidean limit', dpi=300, bbox_inches='tight')
plt.show()

print(alpha_fit)
print(beta_fit)