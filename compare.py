import csv
import datetime
import math

G = 6.67430e-11
M_SUN = 1.98847e30
MU_SUN = G * M_SUN
AU_TO_M = 1.495978707e11
DAY_TO_S = 86400.0
PROBE_MASS = 721.9  # kg

INPUT_CSV = "gravity_assists.csv"
OUTPUT_CSV = "compare.csv"


def sans_gam_velocity(r_i_m: float, r_f_m: float) -> float:
    """
    Strictly from sansGM.pdf:
    1/2 * v_i^2 - G*M/r_i = -G*M/r_f  (as v_f -> 0)
    v_i = sqrt(2 * G * M * (1/r_i - 1/r_f))
    """
    return math.sqrt(2.0 * MU_SUN * ((1.0 / r_i_m) - (1.0 / r_f_m)))


def calc_fuel_cost(v_i: float, mass: float) -> float:
    """
    Strictly from sansGM.pdf:
    Energy Cost = 1/2 * m * v_i^2
    """
    return 0.5 * mass * (v_i ** 2)


def simulate_sans_gam_time_days(init_x_au: float, init_y_au: float, final_x_au: float, final_y_au: float, v_i: float) -> tuple:
    """
    Simulates the travel time from (init_x_au, init_y_au) to the final distance
    r_f using ONLY the Sun's gravity and the kinematic step equations from sim.py,
    with no planets / gravity assists.
    """
    r_i_au = math.hypot(init_x_au, init_y_au)
    r_f_au = math.hypot(final_x_au, final_y_au)

    # Direction unit vector pointing radially outward so it climbs directly from r_i to r_f
    ux = init_x_au / r_i_au
    uy = init_y_au / r_i_au

    x_m = init_x_au * AU_TO_M
    y_m = init_y_au * AU_TO_M
    vx = v_i * ux
    vy = v_i * uy

    dt_s = DAY_TO_S  # 1-day step
    elapsed_seconds = 0.0

    # Run until probe reaches r_f 
    while True:
        r_m = math.hypot(x_m, y_m)
        r_au = r_m / AU_TO_M

        if r_au >= r_f_au:
            break

       
        r_cubed = r_m ** 3
        ax = -MU_SUN * (x_m / r_cubed)
        ay = -MU_SUN * (y_m / r_cubed)

        
        v_fx = vx + ax * dt_s
        v_fy = vy + ay * dt_s

        # Check if radial velocity reverses 
        if (v_fx * x_m + v_fy * y_m) <= 0:
            break

        dx = 0.5 * (vx + v_fx) * dt_s
        dy = 0.5 * (vy + v_fy) * dt_s

        x_m += dx
        y_m += dy
        vx, vy = v_fx, v_fy
        elapsed_seconds += dt_s

    elapsed_days = elapsed_seconds / DAY_TO_S
    final_speed = math.hypot(vx, vy)
    return elapsed_days, vx, vy, final_speed


with open(INPUT_CSV, mode='r', newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    rows = list(reader)


# Filter out any previously appended comparison rows
valid_rows = [
    r for r in rows
    if not str(r['trial_id']).startswith('SANS_GAM')
]

# Read coordinates from gravity assist 
ref_row = valid_rows[0]
init_x_au = float(ref_row['init_x_au'])
init_y_au = float(ref_row['init_y_au'])
final_x_au = float(ref_row['final_x_au'])
final_y_au = float(ref_row['final_y_au'])

r_i_au = math.hypot(init_x_au, init_y_au)
r_f_au = math.hypot(final_x_au, final_y_au)

r_i_m = r_i_au * AU_TO_M
r_f_m = r_f_au * AU_TO_M

# Velocity and Fuel 
v_i = sans_gam_velocity(r_i_m, r_f_m)
fuel_cost_j = calc_fuel_cost(v_i, PROBE_MASS)

# Flight time without any gravity assist maneuvers
flight_days, f_vx, f_vy, f_speed = simulate_sans_gam_time_days(
    init_x_au, init_y_au, final_x_au, final_y_au, v_i
)
flight_years = flight_days / 365.25

launch_dt = datetime.datetime.strptime(ref_row['launch_time_utc'], '%Y-%m-%d %H:%M:%S')
arrival_dt = launch_dt + datetime.timedelta(days=flight_days)

ux = init_x_au / r_i_au
uy = init_y_au / r_i_au

sans_gam_row = {
    'cumulative_assist_counter': len(valid_rows) + 1,
    'trial_id': 'SANS_GAM',
    'total_gravity_assists': 0,
    'assist_number_in_trial': 0,
    'assist_planet': 'None',
    'assist_time_utc': 'N/A',
    'closest_approach_au': 0.0,
    'assist_speed_change_m_s': 0.0,
    'launch_time_utc': ref_row['launch_time_utc'],
    'arrival_time_utc': arrival_dt.strftime('%Y-%m-%d %H:%M:%S'),
    'total_flight_time_days': round(flight_days, 2),
    'total_flight_time_years': round(flight_years, 3),
    'target_r2_au': round(r_f_au, 3),
    'v_inf_m_s': round(v_i, 2),
    'init_x_au': init_x_au,
    'init_y_au': init_y_au,
    'init_vx_m_s': round(v_i * ux, 2),
    'init_vy_m_s': round(v_i * uy, 2),
    'init_speed_m_s': round(v_i, 2),
    'final_vx_m_s': round(f_vx, 2),
    'final_vy_m_s': round(f_vy, 2),
    'final_speed_m_s': round(f_speed, 2),
    'final_heliocentric_dist_au': round(r_f_au, 3),
    'final_x_au': final_x_au,
    'final_y_au': final_y_au,
    'fuel_usage_gm': round(fuel_cost_j, 2)
}

valid_rows.append(sans_gam_row)
calc_fuel_cost
with open(OUTPUT_CSV, mode='w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(valid_rows)


print(f"Start Coords : ({init_x_au}, {init_y_au}) AU -> r_i = {r_i_au:.5f} AU")
print(f"Final Coords : ({final_x_au}, {final_y_au}) AU -> r_f = {r_f_au:.3f} AU")
print(f"Initial v_i  : {v_i:,.2f} m/s")
print(f"Energy Cost  : {fuel_cost_j:.3e} J ({fuel_cost_j:,.2f} J)")
print(f"Travel Time  : {flight_days:,.2f} days ({flight_years:.3f} years)")
print(f"\nSaved all rows + sans-GAM comparison row to '{OUTPUT_CSV}'.")

