import datetime
import math
import matplotlib
matplotlib.use('TkAgg')  # Or 'Qt5Agg' / 'MacOSX' depending on your OS
import matplotlib.pyplot as plt
import solarsystem

TARGET_PLANETS = [
    'Mercury', 'Venus', 'Earth', 'Mars', 
    #'Jupiter', 'Saturn', 'Uranus', 'Neptune'
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

plt.ion()
fig, ax = plt.subplots(figsize=(10, 10))

# Start time: initialize once outside the loop
sim_time = datetime.datetime.now(datetime.timezone.utc)
time_step = datetime.timedelta(hours=12)

try:
    while True:
        # 1. Query heliocentric 3D rectangular coordinates for the simulated time
        H_rect = solarsystem.Heliocentric(
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
        positions = H_rect.planets()

        # 2. Print coordinates to terminal
        print(f"\n--- Solar System Positions at {sim_time.strftime('%Y-%m-%d %H:%M UTC')} (+12h step) ---")
        print(f"{'Planet':<10} {'X (AU)':>10} {'Y (AU)':>10} {'Z (AU)':>10} {'Dist (AU)':>10}")
        print("-" * 54)
        
        for planet in TARGET_PLANETS:
            x, y, z = positions[planet]
            distance = math.sqrt(x**2 + y**2 + z**2)
            print(f"{planet:<10} {x:10.3f} {y:10.3f} {z:10.3f} {distance:10.3f}")

        # 3. Update the orbital plot
        ax.cla()
        ax.set_facecolor('#0b0e14')
        ax.set_xlim((-3, 3))
        ax.set_ylim((-3, 3))
        ax.set_aspect('equal')
        ax.grid(True, color='#21262d', linestyle='--', linewidth=0.5)

        # Plot Sun
        ax.plot(0, 0, 'o', color='gold', markersize=14, label='Sun')

        # Plot orbits and planets
        for planet in TARGET_PLANETS:
            x, y, _ = positions[planet]
            r = math.sqrt(x**2 + y**2)
            color = PLANET_COLORS.get(planet, 'white')

            orbit = plt.Circle((0, 0), r, color=color, fill=False, linestyle=':', alpha=0.35, linewidth=1.2)
            ax.add_artist(orbit)

            ax.plot(x, y, 'o', color=color, markersize=8, label=planet)
            ax.text(x + 0.6, y + 0.6, planet, color='white', fontsize=9, alpha=0.85)

        ax.set_title(f"Heliocentric Solar System\nSimulated Time: {sim_time.strftime('%Y-%m-%d %H:%M UTC')}", 
                     color='white', fontsize=12, pad=12)
        ax.legend(loc='lower left', facecolor='#161b22', edgecolor='none', labelcolor='white', fontsize=8)

        # 4. Advance simulation by 12 hours
        sim_time += time_step

        # 5. Pause for 1 second in real time to handle the GUI redraw
        plt.pause(0.001)

except KeyboardInterrupt:
    print("\nSimulation stopped.")
    plt.close()