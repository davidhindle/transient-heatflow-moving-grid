import matplotlib.pyplot as plt
import numpy as np
import time
import os

# x y plot dimensions on page in cms-------------------------
xplot = 10  # cm
yplot = 16.6666
# -----------------------------------------------------------
xpoff = '2c'
ypoff1 = '2c'

# total time (yrs) and no. of saves, for plot annotation-------
time = 25000000
saves = 30
timeint = time / saves
print(timeint)
print(time, saves, timeint)
# --------------------------------------------------------------

# Setting up plot parameters
plt.figure(figsize=(xplot, yplot))
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

# x,y min and max plot values---------
ymin = -200000
ymax = 0
xmin = 0
xmax = 2000
# -------------------------------------
INT = 500

REGION = [xmin, xmax, ymin, ymax]  # region of plot

# total number of input files to plot -----------------
LIMIT = 30
# -----------------------------------------------------
pref1 = '00'
pref2 = '0'

# Loop for generating plots
for idn in range(1, LIMIT + 1):
    if idn < 10:
        name = f"{pref1}{idn}"
    elif idn == 10:
        name = f"{pref2}{idn}"
    elif idn < 100:
        name = f"{pref2}{idn}"
    elif idn == 100:
        name = f"{idn}"
    else:
        name = f"{idn}"

    suff = 'yr'  # suffix for plot numbering

    name1 = idn * timeint
    print(name, name1)

    while not os.path.isfile(f"sol{name}.dat"):
        print(f"user-file \"sol{name}.dat\" is not there")
        time.sleep(1)

    print(f"user-file \"sol{name}.dat\" is there")

    # Read data from sol{name}.dat and lay{name}.dat files
    sol_data = np.loadtxt(f"sol{name}.dat")
    lay_data = np.loadtxt(f"lay{name}.dat")

    # Reversing the x-coordinate (depth) for plotting
    sol_data[:, 0] = -sol_data[:, 0]
    lay_data[:, 0] = -lay_data[:, 0]

    # Extracting depths and temperatures from lay_data
    depths = lay_data[:, 0]
    temperatures = lay_data[:, 1]

    # Plotting sol_data (depth vs temperature)
    plt.clf()
    plt.subplot(111)
    plt.plot(sol_data[:, 1], sol_data[:, 0], label='sol')

    # Plotting layers as shaded regions
    for i in range(len(depths) - 1):
        plt.fill_betweenx([depths[i], depths[i + 1]], 1500, 2000, alpha=0.3, color=f"C{i}")

    plt.scatter(lay_data[:, 1], lay_data[:, 0], s=10, color='black', label='lay')

    plt.xlabel('temp, K')
    plt.ylabel('depth, metres')
    plt.title(name1)
    plt.legend()

    # Invert y-axis
    plt.gca().invert_yaxis()


    # Save the plot as a JPEG file
    plt.savefig(f"solu{name}.jpg", dpi=300, bbox_inches='tight')

plt.close()

