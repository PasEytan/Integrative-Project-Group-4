import datetime
import math
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import matplotlib.path as mpath
import solarsystem
import numpy as np

authors = ["Astra", "Ariel", "Eytan", "Kangming"]

# --- Constants ---
M = 1 #Solar Masses
G = 4 * np.pi**2
Au_in_m = 149597870691
N_steps = 100000       # number of time steps

# --- Functions ---
def xpos(Xpos):             #Calculates the Rp (AU)
    Xpos = Xpos*(1-ECCENTRICITIES[i])
    return Xpos

def yvel(Xrad):             #Calculates the y-velocity at the Rp (AU)
    Yvel = (G*M*((2/(AU(Xrad)))-(1/(SMA[i]))))**0.5
    return Yvel

def AU(meters):                 #Converts units from meters to AU
    AU = meters / Au_in_m
    return AU


# --- Configuration ---
TARGET_PLANETS = [
    'Mercury', 'Venus', 'Earth', 'Mars',
    'Jupiter', 'Saturn', 'Uranus', 'Neptune'
]

MASS_PLANETS = [
    3.3011e23, 4.86731e24, 5.97217e24, 6.4171e23,
    1.898125e27, 5.68317e26, 8.68099e25, 1.024092e26
]

ECCENTRICITIES = [
    0.205630, 0.006772, 0.0167086, 0.0934,
    0.0489, 0.0565, 0.04717, 0.008678
]

SMA = [
    0.387098, 0.723332, 1, 1.52368055,
    5.2038, 9.5826, 19.19126, 30.07
]

PLANET_COLORS = {
    'Mercury': 'silver',
    'Venus': 'orange',
    'Earth': 'deepskyblue',
    'Mars': 'orangered',
    'Jupiter': 'peru',
    'Saturn': 'khaki',
    'Uranus': 'paleturquoise',
    'Neptune': 'dodgerblue'
}




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


TIME_STEP = datetime.timedelta(days=10)
Delta_t = int(TIME_STEP.total_seconds()) / 3.514e7
# Radius in AU: 3.0 frames the inner planets (Mercury-Mars); set to ~32.0 to fit Neptune
VIEW_LIMIT = 22.0
LABEL_OFFSET = 0.08

# --- Setup Figure & Static Elements ---
plt.ion()
fig, ax = plt.subplots(figsize=(10, 10))
#Path = mpath.Path
ax.set_facecolor('#0b0e14')
ax.set_xlim((-VIEW_LIMIT, VIEW_LIMIT))
ax.set_ylim((-VIEW_LIMIT, VIEW_LIMIT))
ax.set_aspect('equal')
ax.grid(True, color='#21262d', linestyle='--', linewidth=0.5)
# Sun fixed at origin
ax.plot(0, 0, 'o', color='gold', markersize=14, label='Sun')

# Initialize artists once to avoid tearing down the canvas each frame
orbits = {}
planet_markers = {}
planet_labels = {}

for i in range(len(TARGET_PLANETS)):
    color = PLANET_COLORS.get(TARGET_PLANETS[i], 'white')
    SMa.append(SMA[i] * Au_in_m) #Converts semi-major axis from AU to meters



    time, xaxis, yaxis, zaxis, xvelocity, yvelocity, zvelocity = f't{i}', f'x{TARGET_PLANETS[i]}', f'y{TARGET_PLANETS[i]}', f'z{i}', f'vx{TARGET_PLANETS[i]}', f'vy{TARGET_PLANETS[i]}', f'vz{i}'
    t[time], x[xaxis], y[yaxis], z[zaxis] = 0, xpos(SMa[i]), 0, 0
    vx[xvelocity], vy[yvelocity] = 0, yvel(int(x[f'x{TARGET_PLANETS[i]}']))
    x[xaxis] = AU(x[xaxis])
    xlistname, ylistname, xlistval, ylistval = f'x{TARGET_PLANETS[i]}_list', f'y{TARGET_PLANETS[i]}_list', [], []
    x_list[xlistname], y_list[ylistname] = xlistval, ylistval

    orbit = pat.PathPatch((0, 0), label = f"{TARGET_PLANETS[i]}'s Orbit", pickradius=0)
    ax.add_patch(orbit)
    (marker,) = ax.plot([], [], 'o', color=color, markersize=8, label=TARGET_PLANETS[i])
    label = ax.text(0, 0, f" {TARGET_PLANETS[i]}", color='white', fontsize=9, alpha=0.85)

    orbits[TARGET_PLANETS[i]] = orbit
    planet_markers[TARGET_PLANETS[i]] = marker
    planet_labels[TARGET_PLANETS[i]] = label

title_artist = ax.set_title("", color='white', fontsize=12, pad=12)
ax.legend(loc='lower left', facecolor='#161b22', edgecolor='none', labelcolor='white', fontsize=8)



# --- Simulation Loop ---
sim_time = datetime.datetime.now(datetime.timezone.utc)

try:
    while True:
        h_rect = solarsystem.Heliocentric(
            year=sim_time.year,
            month=sim_time.month,
            day=sim_time.day,
            hour=sim_time.hour,
            minute=sim_time.minute,
            UT=0,
            dst=0,
            view='rectangular',
            precession=True
        )
        positions = h_rect.planets()

        for planet in TARGET_PLANETS:                    #loops for each orbiting object
            x_list[f'x{planet}_list'].append(x[f'x{planet}'])
            y_list[f'y{planet}_list'].append(y[f'y{planet}'])

            r = np.sqrt(x[f'x{planet}']**2 + y[f'y{planet}']**2)
            a_x = -G*M*(x[f'x{planet}'])/r**3
            ay = -G*M*(y[f'y{planet}'])/r**3

            x[f'x{planet}'] = x[f'x{planet}'] + vx[f'vx{planet}']*Delta_t + 0.5*a_x*Delta_t**2
            y[f'y{planet}'] = y[f'y{planet}'] + vy[f'vy{planet}']*Delta_t + 0.5*ay*Delta_t**2

            r_new = (x[f'x{planet}']**2 + y[f'y{planet}']**2)**0.5
            ax_new = -G*M*(x[f'x{planet}'])/r_new**3
            ay_new = -G*M*(y[f'y{planet}'])/r_new**3
            vx[f'vx{planet}'] = vx[f'vx{planet}'] + 0.5*(a_x + ax_new)*Delta_t
            vy[f'vy{planet}'] = vy[f'vy{planet}'] + 0.5*(ay + ay_new)*Delta_t

            orbits[planet].set_pickradius(r)
            planet_markers[planet].set_data([x[f'x{planet}']], [y[f'y{planet}']])
            planet_labels[planet].set_position((x[f'x{planet}'] + LABEL_OFFSET, y[f'y{planet}'] + LABEL_OFFSET))

except KeyboardInterrupt:
    print("\nSimulation stopped.")
    plt.close()