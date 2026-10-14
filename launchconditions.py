# attempting to create functions to obtain initial velocity vector and ideal launch window for a transfer orbit
import math

#values for some constants in SI units
M_sun = 2e30
M_earth = 5.97e24
G = 6.67e-11
mu_sun = M_sun * G
mu_earth = M_earth * G
radius_LEO = 6357 + 800 #800 being altitude of assumed LEO

def create_vector(point1, point2):
    vector = []
    for i in point1:
        vector.append(point2-point1)
    return vector 

def magnitude(vector):
    ans = 0
    for i in vector:
        ans += i^2
    ans = math.sqrt(ans)
    return ans

def calculate_delta_v_burn(sun_coords, planet_coords, planet_R, destination_R): #direction of delta v is pi - B degrees from radial heliocentric vector
    r_vector = create_vector(planet_coords, sun_coords)
    v_infinity = math.sqrt(mu_sun/planet_R) * (math.sqrt(2*destination_R/(planet_R+destination_R)) - 1) #magnitude of v infinity
    h_eccentricity = (1 + (radius_LEO * v_infinity^2)mu_earth)^-1
    beta = math.cos(1/h_eccentricity)