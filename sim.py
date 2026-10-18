# Fuel Equation
# Camparison File 
# Python portion of paper
# * graphing file


import csv
import datetime
import math
import solarsystem

# Constants
G = 6.67430e-11                  
AU_TO_M = 1.495978707e11         
DAY_TO_S = 86400.0               
OUTER_BOUNDARY_AU = 121.0         

# Masses of Celestial Objects (kg)
MASSES = {
    'Sun':     1.98847e30,
    'Mercury': 3.3011e23,
    'Venus':   4.8675e24,
    'Earth':   5.9722e24,
    'Mars':    6.4171e23,
    'Jupiter': 1.8982e27,
    'Saturn':  5.6834e26,
    'Uranus':  8.6810e25,
    'Neptune': 1.02413e26
}

# 1/2 m v^2

MU_SUN = G * MASSES['Sun']

TARGET_PLANETS = [
    'Mercury', 'Venus', 'Earth', 'Mars',
    'Jupiter', 'Saturn', 'Uranus', 'Neptune'
]

# Gravity assist detection radius around each planet (in AU)
ASSIST_RADIUS_AU = {
    'Mercury': 0.015,
    'Venus':   0.040,
    'Earth':   0.040,
    'Mars':    0.030,
    'Jupiter': 0.650,
    'Saturn':  0.750,
    'Uranus':  0.800,
    'Neptune': 1.150
}

# Safety Net Radius
MIN_APPR_M = {
    'Mercury': 2.44e6 * 1.05, 
    'Venus':   6.05e6 * 1.05,
    'Earth':   6.37e6 * 1.12,
    'Mars':    3.39e6 * 1.12,
    'Jupiter': 7.15e7 * 1.12,
    'Saturn':  6.03e7 * 1.12,    
    'Uranus':  2.56e7 * 1.12,
    'Neptune': 2.48e7 * 1.12,
}

# Time step
DT_CLOSE = datetime.timedelta(minutes=1)        # Approaching planet Danger zone
DT_FLYBY = datetime.timedelta(hours=1)          # 0.04 days when near a planet
DT_INNER = datetime.timedelta(days=1)           # 1 day inside 3.5 AU
DT_OUTER_PLANETS = datetime.timedelta(days=5)   # 5 days between 3.5 AU and 32.5 AU
DT_DEEP_SPACE = datetime.timedelta(days=20)     # 20 days beyond 32.5 AU (past Neptune)

# Simulation Search Grid Parameters
START_TIME = datetime.datetime.now(datetime.timezone.utc).replace(
    hour=0, minute=0, second=0, microsecond=0
)
LAUNCH_WINDOW_DAYS = 1200         # Search launch dates over the next x years
LAUNCH_STEP_DAYS = 1             # Launch trial every 5 days
TARGET_R2_VALUES_AU = [121.0]  # Trial R2 values for V_inf

MAX_SIM_YEARS = 100               # Cutoff if probe never reaches 121 AU
OUTPUT_CSV = "gravity_assists.csv"


# 2. EPHEMERIS CACHE & ADAPTIVE DT SELECTION
planet_pos_cache = {}

def get_planet_positions_au(sim_time: datetime.datetime) -> dict:
    """
    Returns positions in AU for the Sun and all TARGET_PLANETS,
    caching timestamps so repeated trials reuse ephemeris calculations.
    """
    key = (sim_time.year, sim_time.month, sim_time.day, sim_time.hour)
    if key not in planet_pos_cache:
        h_rect = solarsystem.Heliocentric(
            year=sim_time.year,
            month=sim_time.month,
            day=sim_time.day,
            hour=sim_time.hour,
            minute=0,
            UT=0,
            dst=0,
            view='rectangular',
            precession=True
        )
        raw_pos = h_rect.planets()
        pos_2d = {'Sun': (0.0, 0.0)}
        for planet in TARGET_PLANETS:
            x, y, _ = raw_pos[planet]
            pos_2d[planet] = (x, y)
        planet_pos_cache[key] = pos_2d
    return planet_pos_cache[key]


def select_adaptive_dt(probe_x_au: float, probe_y_au: float, heliocentric_dist_au: float, positions_au: dict, elapsed_days: float) -> datetime.timedelta:
    """
    Switched between the different time steps available depending on how
    far away the probe is to the sun.
    """
    # 1. Beyond Neptune's orbit (> 32.5 AU), use large 20-day deep-space steps
    if heliocentric_dist_au > 32.5:
        return DT_DEEP_SPACE

    # 2. Check if probe is approaching any planet (within 3x assist radius)
    for planet in TARGET_PLANETS:
        if planet == 'Earth' and elapsed_days < 30.0:
            continue
        px_au, py_au = positions_au[planet]
        dist_au = math.hypot(probe_x_au - px_au, probe_y_au - py_au)
        if dist_au <= 3.0 * ASSIST_RADIUS_AU[planet]:
            return DT_FLYBY
        # Danger zone
        if dist_au * AU_TO_M <= 20.0 * MIN_APPR_M[planet]:
            return DT_CLOSE

    # 3. Inner solar system (< 3.5 AU) vs Outer planetary zone (3.5 - 32.5 AU)
    if heliocentric_dist_au < 3.5:
        return DT_INNER
    return DT_OUTER_PLANETS

def closest_approach_in_step(rel0: tuple, rel1: tuple) -> float:
    """
    rel0, rel1: probe position minus planet position (m) at the start and end
    of the step. Treats the relative motion as a straight line and returns the
    minimum distance reached on that segment.
    """
    dx = rel1[0] - rel0[0]
    dy = rel1[1] - rel0[1]
    seg_sq = dx * dx + dy * dy
    if seg_sq == 0.0:
        t = 0.0
    else:
        t = -(rel0[0] * dx + rel0[1] * dy) / seg_sq
        t = max(0.0, min(1.0, t))        # clamp to the segment
    return math.hypot(rel0[0] + t * dx, rel0[1] + t * dy)

# V inf as start velocity
def calculate_v_infinity(r1_m: float, r2_m: float) -> float:
    """
    Calculates hyperbolic excess speed v_inf from the image equation:
    v_inf = sqrt(mu_sun / R1) * (sqrt(2 * R2 / (R1 + R2)) - 1)
    """
    v_circular = math.sqrt(MU_SUN / r1_m)
    transfer_factor = math.sqrt((2.0 * r2_m) / (r1_m + r2_m)) - 1.0
    return v_circular * transfer_factor


def calculate_acceleration(probe_x_m: float, probe_y_m: float, positions_au: dict) -> tuple:
    """
    Calculates 2D gravitational acceleration (a_x, a_y) on the probe from the Sun
    and all planets using the summation equation in GM.pdf:
    a_x = -G * sum(M_i * x_i / (x_i^2 + y_i^2)^(3/2))
    a_y = -G * sum(M_i * y_i / (x_i^2 + y_i^2)^(3/2))
    """
    sum_x = 0.0
    sum_y = 0.0

    for body_name, mass in MASSES.items():
        body_x_m = positions_au[body_name][0] * AU_TO_M
        body_y_m = positions_au[body_name][1] * AU_TO_M

        # Displacement from celestial object i to space probe
        x_i = probe_x_m - body_x_m
        y_i = probe_y_m - body_y_m

        r_sq = x_i * x_i + y_i * y_i
        if r_sq < 1e6:
            continue

        r_cubed = r_sq ** 1.5
        sum_x += mass * (x_i / r_cubed)
        sum_y += mass * (y_i / r_cubed)

    return -G * sum_x, -G * sum_y


def kinematic_step(x_m: float, y_m: float, vx: float, vy: float, ax: float, ay: float, dt_s: float) -> tuple:
    """
    Updates velocity and position over time-step dt_s using GM.pdf kinematics:
        V_fx = V_ix + a_x * dt
        V_fy = V_iy + a_y * dt
        V_fx^2 = V_ix^2 + 2 * a_x * d_x  ==>  d_x = 0.5 * (V_ix + V_fx) * dt
    """
    v_fx = vx + ax * dt_s
    v_fy = vy + ay * dt_s

    dx = 0.5 * (vx + v_fx) * dt_s
    dy = 0.5 * (vy + v_fy) * dt_s

    return x_m + dx, y_m + dy, v_fx, v_fy


# Single Trial Sim
def run_trial(trial_id: int, launch_time: datetime.datetime, r2_target_au: float) -> list:
    """
    Simulates a single probe trajectory from launch_time until it reaches 121 AU
    or exceeds MAX_SIM_YEARS. Returns a list of detected gravity assist records.
    """
    init_positions = get_planet_positions_au(launch_time)
    earth_x_au, earth_y_au = init_positions['Earth']
    r1_au = math.hypot(earth_x_au, earth_y_au)
    r1_m = r1_au * AU_TO_M
    r2_m = r2_target_au * AU_TO_M

    # Unit vector perpendicular to Earth's position vector (prograde tangent)
    perp_ux = -earth_y_au / r1_au
    perp_uy = earth_x_au / r1_au

    # Initial speeds: Earth circular speed + V_inf from image
    v_earth = math.sqrt(MU_SUN / r1_m)
    v_inf = calculate_v_infinity(r1_m, r2_m)
    init_speed = v_earth + v_inf

    init_vx = init_speed * perp_ux
    init_vy = init_speed * perp_uy

    # Start just outside Earth's SOI (~0.0062 AU) along the perpendicular direction
    earth_soi_au = 0.0062
    init_x_au = earth_x_au + earth_soi_au * perp_ux
    init_y_au = earth_y_au + earth_soi_au * perp_uy

    probe_x_m = init_x_au * AU_TO_M
    probe_y_m = init_y_au * AU_TO_M
    probe_x_au, probe_y_au = init_x_au, init_y_au
    heliocentric_dist_au = math.hypot(probe_x_au, probe_y_au)
    vx, vy = init_vx, init_vy

    sim_time = launch_time
    max_duration_days = MAX_SIM_YEARS * 365.25
    elapsed_days = 0.0

    active_encounters = {}
    completed_assists = []

    while elapsed_days < max_duration_days:
        positions_au = get_planet_positions_au(sim_time)

        # Dynamically choose computationally reasonable dt
        dt_step = select_adaptive_dt(
            probe_x_au, probe_y_au, heliocentric_dist_au, positions_au, elapsed_days
        )
        dt_s = dt_step.total_seconds()

        prev_x_m, prev_y_m = probe_x_m, probe_y_m 

        # Compute acceleration & advance state
        ax, ay = calculate_acceleration(probe_x_m, probe_y_m, positions_au)
        probe_x_m, probe_y_m, vx, vy = kinematic_step(probe_x_m, probe_y_m, vx, vy, ax, ay, dt_s)

        sim_time += dt_step
        elapsed_days = (sim_time - launch_time).total_seconds() / DAY_TO_S

        probe_x_au = probe_x_m / AU_TO_M
        probe_y_au = probe_y_m / AU_TO_M
        heliocentric_dist_au = math.hypot(probe_x_au, probe_y_au)
        current_speed = math.hypot(vx, vy)

        # if probe gets too close to the planet
        if heliocentric_dist_au <= 32.5:
            positions_end_au = get_planet_positions_au(sim_time)
            for planet in TARGET_PLANETS:
                if planet == 'Earth' and elapsed_days < 30.0:
                    continue
                rel0 = (prev_x_m - positions_au[planet][0] * AU_TO_M,
                        prev_y_m - positions_au[planet][1] * AU_TO_M)
                rel1 = (probe_x_m - positions_end_au[planet][0] * AU_TO_M,
                        probe_y_m - positions_end_au[planet][1] * AU_TO_M)
                if closest_approach_in_step(rel0, rel1) < MIN_APPR_M[planet]:
                    return []      # crashed 

        # Check for Gravity Assist encounters (skipping initial Earth departure)
        if heliocentric_dist_au <= 32.5:
            for planet in TARGET_PLANETS:
                if planet == 'Earth' and elapsed_days < 30.0:
                    continue

                px_au, py_au = positions_au[planet]
                dist_to_planet_au = math.hypot(probe_x_au - px_au, probe_y_au - py_au)
                threshold_au = ASSIST_RADIUS_AU[planet]

                if dist_to_planet_au <= threshold_au:
                    if planet not in active_encounters:
                        active_encounters[planet] = {
                            'planet': planet,
                            'entry_time': sim_time,
                            'closest_time': sim_time,
                            'min_dist_au': dist_to_planet_au,
                            'entry_speed_m_s': current_speed
                        }
                    else:
                        if dist_to_planet_au < active_encounters[planet]['min_dist_au']:
                            active_encounters[planet]['min_dist_au'] = dist_to_planet_au
                            active_encounters[planet]['closest_time'] = sim_time
                else:
                    if planet in active_encounters:
                        enc = active_encounters.pop(planet)
                        enc['exit_time'] = sim_time
                        enc['exit_speed_m_s'] = current_speed
                        enc['delta_v_assist_m_s'] = current_speed - enc['entry_speed_m_s']
                        completed_assists.append(enc)

        # Check if probe has reached the outer solar system boundary (121 AU)
        if heliocentric_dist_au >= OUTER_BOUNDARY_AU:
            for planet, enc in list(active_encounters.items()):
                enc['exit_time'] = sim_time
                enc['exit_speed_m_s'] = current_speed
                enc['delta_v_assist_m_s'] = current_speed - enc['entry_speed_m_s']
                completed_assists.append(enc)
            active_encounters.clear()

            if not completed_assists:
                return []

            total_assists = len(completed_assists)
            rows = []
            for idx, assist in enumerate(completed_assists, start=1):
                rows.append({
                    'trial_id': trial_id,
                    'total_gravity_assists': total_assists,
                    'assist_number_in_trial': idx,
                    'assist_planet': assist['planet'],
                    'assist_time_utc': assist['closest_time'].strftime('%Y-%m-%d %H:%M:%S'),
                    'closest_approach_au': round(assist['min_dist_au'], 5),
                    'assist_speed_change_m_s': round(assist['delta_v_assist_m_s'], 2),
                    'launch_time_utc': launch_time.strftime('%Y-%m-%d %H:%M:%S'),
                    'arrival_time_utc': sim_time.strftime('%Y-%m-%d %H:%M:%S'),
                    'total_flight_time_days': round(elapsed_days, 2),
                    'total_flight_time_years': round(elapsed_days / 365.25, 3),
                    # Beginning parameters of the trial
                    'target_r2_au': r2_target_au,
                    'v_inf_m_s': round(v_inf, 2),
                    'init_x_au': round(init_x_au, 5),
                    'init_y_au': round(init_y_au, 5),
                    'init_vx_m_s': round(init_vx, 2),
                    'init_vy_m_s': round(init_vy, 2),
                    'init_speed_m_s': round(init_speed, 2),
                    # End parameters of the trial
                    'final_vx_m_s': round(vx, 2),
                    'final_vy_m_s': round(vy, 2),
                    'final_speed_m_s': round(current_speed, 2),
                    'final_heliocentric_dist_au': round(heliocentric_dist_au, 3),
                    # Final (x, y) coordinates where it finishes the trial
                    'final_x_au': round(probe_x_au, 5),
                    'final_y_au': round(probe_y_au, 5)
                })
            return rows

    return []


# 5. Running + csv file
def main():
    fieldnames = [
        'cumulative_assist_counter',
        'trial_id',
        'total_gravity_assists',
        'assist_number_in_trial',
        'assist_planet',
        'assist_time_utc',
        'closest_approach_au',
        'assist_speed_change_m_s',
        'launch_time_utc',
        'arrival_time_utc',
        'total_flight_time_days',
        'total_flight_time_years',
        'target_r2_au',
        'v_inf_m_s',
        'init_x_au',
        'init_y_au',
        'init_vx_m_s',
        'init_vy_m_s',
        'init_speed_m_s',
        'final_vx_m_s',
        'final_vy_m_s',
        'final_speed_m_s',
        'final_heliocentric_dist_au',
        'final_x_au',
        'final_y_au'
    ]

    all_assist_records = []
    trial_id = 1
    num_launch_steps = LAUNCH_WINDOW_DAYS // LAUNCH_STEP_DAYS

    print(f"Starting 2D Gravity Assist Simulation to {OUTER_BOUNDARY_AU} AU...")
    print(f"Running {num_launch_steps * len(TARGET_R2_VALUES_AU)} trials with adaptive time-stepping...")

    for day_offset in range(0, LAUNCH_WINDOW_DAYS, LAUNCH_STEP_DAYS):
        launch_time = START_TIME + datetime.timedelta(days=day_offset)
        for r2_au in TARGET_R2_VALUES_AU:
            trial_rows = run_trial(trial_id, launch_time, r2_au)
            if trial_rows:
                all_assist_records.extend(trial_rows)
                for row in trial_rows:
                    print(
                        f"[Trial {row['trial_id']:04d}] Assist {row['assist_number_in_trial']}/{row['total_gravity_assists']}: "
                        f"{row['assist_planet']:<7} on {row['assist_time_utc']} | "
                        f"dV: {row['assist_speed_change_m_s']:+7.1f} m/s | "
                        f"121 AU Arrival: {row['arrival_time_utc']} ({row['total_flight_time_years']} yrs)"
                    )
            trial_id += 1

    # Sort results by fastest flight time to 121 AU
    all_assist_records.sort(key=lambda r: (r['total_flight_time_days'], r['trial_id'], r['assist_number_in_trial']))

    # Populate the global cumulative counter after sorting
    for global_idx, row in enumerate(all_assist_records, start=1):
        row['cumulative_assist_counter'] = global_idx

    with open(OUTPUT_CSV, mode='w', newline='', encoding='utf-8') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_assist_records)

    print(f"\nSimulation complete. Total gravity assists found: {len(all_assist_records)}")
    print(f"Results saved to '{OUTPUT_CSV}'.")


if __name__ == '__main__':
    main()