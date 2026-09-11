import datetime
import math
import matplotlib.pyplot as plt
import solarsystem

# --- Configuration ---
TARGET_PLANETS = [
    'Mercury', 'Venus', 'Earth', 'Mars',
    'Jupiter', 'Saturn', 'Uranus', 'Neptune'
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

TIME_STEP = datetime.timedelta(days=10)
VIEW_LIMIT = 22.0
LABEL_OFFSET = 0.08

# --- Setup Figure & Static Elements ---
plt.ion()
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_facecolor('#0b0e14')
ax.set_aspect('equal')
ax.grid(True, color='#21262d', linestyle='--', linewidth=0.5)

# Sun fixed at origin
sun_marker, = ax.plot(0, 0, 'o', color='gold', markersize=14, label='Sun')

# Initialize artists once
orbit_circles = {}
planet_markers = {}
planet_labels = {}

for planet in TARGET_PLANETS:
    color = PLANET_COLORS.get(planet, 'white')

    circle = plt.Circle((0, 0), radius=0, color=color, fill=False, linestyle=':', alpha=0.35, linewidth=1.2)
    ax.add_patch(circle)

    (marker,) = ax.plot([], [], 'o', color=color, markersize=8, label=planet)
    label = ax.text(0, 0, f" {planet}", color='white', fontsize=9, alpha=0.85)

    orbit_circles[planet] = circle
    planet_markers[planet] = marker
    planet_labels[planet] = label

title_artist = ax.set_title("", color='white', fontsize=12, pad=12)
ax.legend(loc='lower left', facecolor='#161b22', edgecolor='none', labelcolor='white', fontsize=8)

# --- Simulation Loop ---
sim_time = datetime.datetime.now(datetime.timezone.utc)
grow = True

try:
    while True:
        # Breathing zoom effect
        if grow:
            VIEW_LIMIT += 0.1
            if VIEW_LIMIT > 33:
                grow = False
        else:
            VIEW_LIMIT -= 0.1
            if VIEW_LIMIT < 0.5:
                grow = True

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

        # Center camera on Earth
        earth_x, earth_y, _ = positions['Earth']
        ax.set_xlim((earth_x - VIEW_LIMIT, earth_x + VIEW_LIMIT))
        ax.set_ylim((earth_y - VIEW_LIMIT, earth_y + VIEW_LIMIT))

        # Terminal readout
        print(f"\n--- Positions at {sim_time.strftime('%Y-%m-%d %H:%M UTC')} ---")
        for planet in TARGET_PLANETS:
            x, y, z = positions[planet]
            dist = math.hypot(x, y, z)
            print(f"{planet:<10} {x:10.3f} {y:10.3f} {z:10.3f} {dist:10.3f}")

        # Update persistent artists
        for planet in TARGET_PLANETS:
            x, y, _ = positions[planet]
            r = math.hypot(x, y)

            orbit_circles[planet].set_radius(r)
            planet_markers[planet].set_data([x], [y])
            planet_labels[planet].set_position((x + LABEL_OFFSET, y + LABEL_OFFSET))

        title_artist.set_text(
            f"Earth-Centered Camera (Heliocentric Coordinates)\nSimulated Time: {sim_time.strftime('%Y-%m-%d %H:%M UTC')}"
        )

        sim_time += TIME_STEP
        plt.pause(0.001)

except KeyboardInterrupt:
    print("\nSimulation stopped.")
    plt.close()