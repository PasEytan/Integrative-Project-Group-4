#written by Nigel D'Mello and Yuhan Peter Bauer
import numpy as np
import matplotlib . pyplot as plt

authors = "Nigel D'Mello & Yuhan Peter Bauer"

# - - - - - - - - - - - - - - - -
# Physical constants (natural units)
# -------------------------------

Au_in_m = 149597870691
G = 4 * np.pi**2       # AU^3/(yr^2 * solar mass)

# -------------------------------
# Numerical parameters
# -------------------------------

Delta_t = 1.908e-6     # time step in years (1 min in years)
N_steps = 100000       # number of time steps

# -------------------------------
# Initial Condiditons
# -------------------------------
C_name = "Jupiter"                                #name of the central body
names = ["Io", "Europa", "Ganymede", "Callisto"]  #name of the orbiting object
e = [0.004, 0.009, 0.001, 0.009]                  #eccentricity
SMA = [2.82e-3, 4.49e-3, 7.16e-3, 1.26e-2]        #semi-major axis (AU)
M = 9.49e-4                                       # mass of central body (Solar Masses)

# DO NOT CHANGE ANYTHING BELOW
SMa = []
t = {}      #dictionary for the different time lists
x = {}      #dictionary for the x-pos lists
y = {}      #dictionary for the y-pos lists
z = {}      #dictionary for the z-pos lists
vx = {}     #dictionary for the x-velocity lists
vy = {}     #dictionary for the y-velocity lists
zy = {}     #dictionary for the z-velocity lists
x_list = {} #dictionary to store the updating x-positions of the orbits
y_list = {} #dictionary to store the updating y-positions of the orbits
z_list = {} #dictionary to store the updating z-positions of the orbits

for i in range(len((SMA))):
    SMa.append(SMA[i] * Au_in_m) #Converts semi-major axis from AU to meters
    
    def xpos(Xpos):             #Calculates the Rp (AU)
        Xpos = Xpos*(1-e[i])
        return Xpos
    
    def yvel(Xrad):             #Calculates the y-velocity at the Rp (AU)
        Yvel = (G*M*((2/(AU(Xrad)))-(1/(SMA[i]))))**0.5
        return Yvel
    
    def AU(meters):                 #Converts units from meters to AU
        AU = meters / Au_in_m
        return AU
    
    
    time, xaxis, yaxis, zaxis, xvelocity, yvelocity, zvelocity = f't{i}', f'x{i}', f'y{i}', f'z{i}', f'vx{i}', f'vy{i}', f'vz{i}'
    t[time], x[xaxis], y[yaxis], z[zaxis] = 0, xpos(SMa[i]), 0, 0
    vx[xvelocity], vy[yvelocity] = 0, yvel(int(x[f'x{i}']))
    x[xaxis] = AU(x[xaxis])
    xlistname, ylistname, xlistval, ylistval = f'x{i}_list', f'y{i}_list', [], []
    x_list[xlistname], y_list[ylistname] = xlistval, ylistval
    

for i in range(len(names)):                     #loops for each orbiting object
    for j in range(N_steps):                    #loops for the defined number of steps
        x_list[f'x{i}_list'].append(x[f'x{i}'])
        y_list[f'y{i}_list'].append(y[f'y{i}'])
        
        r = np.sqrt(x[f'x{i}']**2 + y[f'y{i}']**2)
        ax = -G*M*(x[f'x{i}'])/r**3
        ay = -G*M*(y[f'y{i}'])/r**3
        
        x[f'x{i}'] = x[f'x{i}'] + vx[f'vx{i}']*Delta_t + 0.5*ax*Delta_t**2
        y[f'y{i}'] = y[f'y{i}'] + vy[f'vy{i}']*Delta_t + 0.5*ay*Delta_t**2
        
        r_new = (x[f'x{i}']**2 + y[f'y{i}']**2)**0.5
        ax_new = -G*M*(x[f'x{i}'])/r_new**3
        ay_new = -G*M*(y[f'y{i}'])/r_new**3
        vx[f'vx{i}'] = vx[f'vx{i}'] + 0.5*(ax + ax_new)*Delta_t
        vy[f'vy{i}'] = vy[f'vy{i}'] + 0.5*(ay + ay_new)*Delta_t

for i in range(len(names)):
    plt.plot(x_list[f'x{i}_list'], y_list[f'y{i}_list'], label = f"{names[i]}'s Orbit")
plt.plot (0, 0, 'o', label = C_name)
plt.xlabel('x (AU)')
plt.ylabel('y (AU)')
plt.axis('equal')
string = f'Orbit of {names} around {C_name} - {authors}'
plt.title(string.replace('[','').replace(']','').replace("'",""))
plt.legend()
plt.show()


    