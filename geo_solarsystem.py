import datetime
import math
from collections import defaultdict, deque
import matplotlib.pyplot as plt
import solarsystem

# --- Configuration ---
TARGET_PLANETS = [
    'Sun', 'Moon', 'Mercury', 'Venus',
    'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune'
]

PLANET_COLORS = {
    'Sun': 'gold',
    'Mercury': 'silver',
    'Venus': 'orange',
    'Moon': 'lightgray',
    'Earth': 'deepskyblue',
    'Mars': 'orangered',
    'Jupiter': 'peru',
    'Saturn': 'khaki',
    'Uranus': 'paleturquoise',
    'Neptune': 'dodgerblue'
}

MAX_TRAIL_POINTS = 500
TIME_STEP = datetime.timedelta(days=8)
MOON_VISUAL_SCALE = 50.0  # Scale Moon orbit (~0.00257 AU) to be visible on the plot


def get_moon_geocentric(sim_time: datetime.datetime) -> tuple[float, float, float]:
    """Approximate geocentric Moon coordinates in AU relative to ecliptic."""
    epoch = datetime.datetime(2000, 1, 1, 12, 0, tzinfo=datetime.timezone.utc)
    days = (sim_time - epoch).total_seconds() / 86400.0

    theta = math.radians((218.316 + 13.176396 * days) % 360.0)
    dist_au = 0.00257
    inc = math.radians(5.14)

    x = dist_au * math.cos(theta)
    y = dist_au * math.sin(theta) * math.cos(inc)
    z = dist_au * math.sin(theta) * math.sin(inc)
    return x, y, z


def to_cartesian(body: str, ra: float, dec: float, dist: float) -> tuple[float, float, float]:
    """Convert spherical equatorial coordinates (RA, Dec, Dist) to 3D Cartesian."""
    # Convert Earth radii to AU if Moon distance is returned in Earth radii
    if body == 'Moon' and dist > 1.0:
        dist *= 4.2635e-5

    ra_rad = math.radians(ra)
    dec_rad = math.radians(dec)

    x = dist * math.cos(dec_rad) * math.cos(ra_rad)
    y = dist * math.cos(dec_rad) * math.sin(ra_rad)
    z = dist * math.sin(dec_rad)
    return x, y, z


# --- Setup Figure & Static Elements ---
plt.ion()
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_facecolor('#0b0e14')
ax.set_xlim((-10, 10))
ax.set_ylim((-10, 10))
ax.set_aspect('equal')
ax.grid(True, color='#21262d', linestyle='--', linewidth=0.5)

# Earth fixed at center (Geocentric view)
ax.plot(0, 0, 'o', color=PLANET_COLORS['Earth'], markersize=10, label='Earth')

# Initialize artist references once to avoid redrawing canvas from scratch
trail_lines = {}
body_markers = {}
body_labels = {}
orbit_history = defaultdict(lambda: deque(maxlen=MAX_TRAIL_POINTS))

for body in TARGET_PLANETS:
    color = PLANET_COLORS.get(body, 'white')
    marker_size = 12 if body == 'Sun' else (4 if body == 'Moon' else 6)

    (trail_line,) = ax.plot([], [], '-', color=color, alpha=0.35, linewidth=1.0)
    (body_marker,) = ax.plot([], [], 'o', color=color, markersize=marker_size, label=body)
    label = ax.text(0, 0, f" {body}", color='white', fontsize=9, alpha=0.85)

    trail_lines[body] = trail_line
    body_markers[body] = body_marker
    body_labels[body] = label

title_artist = ax.set_title("", color='white', fontsize=12, pad=12)
ax.legend(loc='lower left', facecolor='#161b22', edgecolor='none', labelcolor='white', fontsize=8)

# --- Simulation Loop ---
sim_time = datetime.datetime.now(datetime.timezone.utc)

try:
    while True:
        geo_calc = solarsystem.Geocentric(
            year=sim_time.year,
            month=sim_time.month,
            day=sim_time.day,
            hour=sim_time.hour,
            minute=sim_time.minute,
            UT=0,
            dst=0,
            precession=True
        )
        raw_positions = geo_calc.position()

        # Compute Cartesian positions
        positions = {
            body: to_cartesian(body, *raw_positions[body])
            for body in TARGET_PLANETS if body != 'Moon'
        }
        positions['Moon'] = get_moon_geocentric(sim_time)

        # Print state to console
        print(f"\n--- Geocentric Positions at {sim_time.strftime('%Y-%m-%d %H:%M UTC')} (+{TIME_STEP.days}d step) ---")
        print(f"{'Planet':<10} {'X (AU)':>10} {'Y (AU)':>10} {'Z (AU)':>10} {'Dist (AU)':>10}")
        print("-" * 54)

        for planet in TARGET_PLANETS:
            x, y, z = positions[planet]
            dist = math.sqrt(x**2 + y**2 + z**2)
            print(f"{planet:<10} {x:10.3f} {y:10.3f} {z:10.3f} {dist:10.3f}")

        # Update plot artists in-place
        for body in TARGET_PLANETS:
            x, y, _ = positions[body]

            if body == 'Moon':
                x *= MOON_VISUAL_SCALE
                y *= MOON_VISUAL_SCALE

            history = orbit_history[body]
            history.append((x, y))

            # Update trail and positions
            tx, ty = zip(*history)
            trail_lines[body].set_data(tx, ty)
            body_markers[body].set_data([x], [y])
            body_labels[body].set_position((x + 0.08, y + 0.08))

        title_artist.set_text(
            f"Geocentric Solar System\nSimulated Time: {sim_time.strftime('%Y-%m-%d %H:%M UTC')}"
        )

        sim_time += TIME_STEP
        plt.pause(0.001)

except KeyboardInterrupt:
    print("\nSimulation stopped.")
    plt.close()