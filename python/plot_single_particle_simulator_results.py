'''
Copyright (c) 2009-2022 Kartik Kumar (me@kartikkumar.com)
Distributed under the MIT License.
See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
'''

# Set up modules and packages.
# Plotting
# import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
# from matplotlib import rcParams
# from matplotlib import cm
# from matplotlib.font_manager import FontProperties
# from mpl_toolkits.mplot3d import Axes3D
# from mpl_toolkits.mplot3d import axes3d

# I/O
import json
from jsmin import jsmin
from pprint import pprint

# Numerical
import math
import numpy as np
import pandas as pd

# System
import sys
import time

# Start timer.
start_time = time.time( )

print ("")
print ("------------------------------------------------------------------")
print ("                             dustsim                              ")
print ("      Copyright (c) 2009-2023, K. Kumar (me@kartikkumar.com)      ")
print ("------------------------------------------------------------------")
print ("")

print ("")
print ("******************************************************************")
print ("                          Input parameters                        ")
print ("******************************************************************")
print ("")

# Parse JSON configuration file.
# Raise exception if wrong number of inputs are provided to script.
if len(sys.argv) != 2:
    raise Exception("Only provide a JSON config file as input!")

json_input_file = open(sys.argv[1])
with open(sys.argv[1], 'r') as json_input_file:
  json_input_string = jsmin(json_input_file.read())
config  = json.loads(json_input_string)
pprint(config)

print ("")
print ("******************************************************************")
print ("                            Operations                            ")
print ("******************************************************************")
print ("")

print ("Input data files being read ...")

input_path_prefix = config["input_directory"] + "/"
output_path_prefix = config["output_directory"] + '/'

# Read and store data files.
state_history = pd.read_csv(input_path_prefix + config["state_history_file"])

metadata = pd.read_csv(input_path_prefix + config["metadata_file"], header=None)
metadata_table = []
metadata_table.append(["Gravitational parameter",metadata[1][0],'${0}$'.format(metadata[2][0])])
metadata_table.append(["J2 coefficient",metadata[1][1],'${0}$'.format(metadata[2][1])])
metadata_table.append(["Equatorial radius",metadata[1][2],'${0}$'.format(metadata[2][2])])
metadata_table.append(["Initial state (Kepler)",metadata[1][3],'${0}$'.format(metadata[2][3])])
metadata_table.append(["Start epoch",metadata[1][4],'${0}$'.format(metadata[2][4])])
metadata_table.append(["End epoch",metadata[1][5],'${0}$'.format(metadata[2][5])])
metadata_table.append(["Time step",metadata[1][6],'${0}$'.format(metadata[2][6])])

print ("Input data files successfully read!")

print ("figures being generated ...")

# Generate figure with 2D views.
fig = plt.figure()
ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4, frameon=False)

# Plot X-Y projection.
ax1.plot(state_history['x'],state_history['y'],color='k')
ax1.scatter(0.0,0.0,s=100,marker='o',color='b')
ax1.set_xlabel('x [km]')
ax1.set_ylabel('y [km]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax1.grid()

# Plot X-Z projection.
ax2.plot(state_history['x'],state_history['z'],color='k')
ax2.scatter(0.0,0.0,s=100,marker='o',color='b')
ax2.set_xlabel('x [km]')
ax2.set_ylabel('z [km]')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax2.grid()

# Plot Y-Z projection.
ax3.plot(state_history['y'],state_history['z'],color='k')
ax3.scatter(0.0,0.0,s=100,marker='o',color='b')
ax3.set_xlabel('y [km]')
ax3.set_ylabel('z [km]')
ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax3.grid()

# Plot metadata table.
# # @TODO: Fix font size.
# ax4.axis('off')
# the_table = ax4.table(cellText=metadata_table,colLabels=None,cellLoc='center',loc='center')
# table_props = the_table.properties()
# table_cells = table_props['child_artists']
# for cell in table_cells:
#     cell.set_height(0.15)
# cell_dict = the_table.get_celld()
# for row in range(0,len(metadata_table)):
#     cell_dict[(row,2)].set_width(0.1)

# Save figure.
plt.tight_layout()
plt.savefig(output_path_prefix + config["2D_figure"], dpi=config["figure_dpi"])

# Generate figure with time histories of Keplerian elements.
fig = plt.figure()
ax1 = fig.add_subplot(2, 3, 1)
ax2 = fig.add_subplot(2, 3, 2)
ax3 = fig.add_subplot(2, 3, 3)
ax4 = fig.add_subplot(2, 3, 4)
ax5 = fig.add_subplot(2, 3, 5)
ax6 = fig.add_subplot(2, 3, 6)

# Plot time-history of semi-major axis.
ax1.plot(state_history['t'],state_history['a'],color='k')
ax1.set_xlabel('t [s]')
ax1.set_ylabel('a [km]')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax1.grid()

# Plot time-history of eccentricity.
ax2.plot(state_history['t'],state_history['e'],color='k')
ax2.set_xlabel('t [s]')
ax2.set_ylabel('e [-]')
ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax2.grid()

# Plot time-history of inclination.
ax3.plot(state_history['t'],state_history['i'].apply(math.degrees),color='k')
ax3.set_xlabel('t [s]')
ax3.set_ylabel('i [deg]')
ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax3.grid()

# Plot time-history of argument of periapsis.
ax4.plot(state_history['t'],np.unwrap(np.degrees(state_history['aop'])),color='k')
ax4.set_xlabel('t [s]')
ax4.set_ylabel('$\omega$ [deg]')
ax4.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax4.grid()

# Plot time-history of right ascension of ascending node.
ax5.plot(state_history['t'],np.unwrap(np.degrees(state_history['raan'])),color='k')
ax5.set_xlabel('t [s]')
ax5.set_ylabel('$\Omega$ [deg]')
ax5.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax5.grid()

# Plot time-history of true anomaly.
# ax6.plot(state_history['t'],np.unwrap(np.degrees(state_history['ta'])),color='k')
ax6.plot(state_history['t'],np.degrees(state_history['ta']),color='k')
ax6.set_xlabel('t [s]')
ax6.set_ylabel(r'$\theta$ [deg]')
ax6.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax6.grid()

# Save figure.
plt.tight_layout()
plt.savefig(output_path_prefix + config["kepler_figure"], dpi=config["figure_dpi"])

# Generate 3D figure if requested.
if config["show_3D_figure"]:
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.set_xlabel('x [km]')
    ax.set_ylabel('y [km]')
    ax.set_zlabel('z [km]')
    ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

    # Plot sphere for the central body.
    radius_central_body = config[ "central_body_radius" ] # km
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = radius_central_body * np.outer(np.cos(u), np.sin(v))
    y = radius_central_body * np.outer(np.sin(u), np.sin(v))
    z = radius_central_body * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='b',edgecolors='b')

    # Plot dust particle trajectory.
    ax.plot3D(state_history['x'],state_history['y'],state_history['z'],'k')

    # Create cubic bounding box to simulate equal aspect ratio.
    X = state_history['x']
    Y = state_history['y']
    Z = state_history['z']
    max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0
    mean_x = X.mean()
    mean_y = Y.mean()
    mean_z = Z.mean()
    ax.set_xlim(mean_x - max_range, mean_x + max_range)
    ax.set_ylim(mean_y - max_range, mean_y + max_range)
    ax.set_zlim(mean_z - max_range, mean_z + max_range)

    plt.grid()
    plt.show()

# Generate figure with time history of orbital energy.
fig = plt.figure()
gravitational_parameter = float(metadata[1][0])
plt.plot(state_history['t'],-gravitational_parameter/(2*state_history['a']),color='k')
plt.tight_layout()
plt.savefig(output_path_prefix + "orbital-energy-time-history.png", dpi=config["figure_dpi"])

# Generate figure with time history of true anomaly.
fig = plt.figure()
plt.plot(state_history['t'],np.degrees(state_history['ta']),color='k')
plt.tight_layout()
plt.savefig(output_path_prefix + "true-anomaly-time-history.png", dpi=config["figure_dpi"])

# Generate figure with time history of x-position.
fig = plt.figure()
plt.plot(state_history['t'],state_history['x'],color='k')
plt.tight_layout()
plt.savefig(output_path_prefix + "x-position-time-history.png", dpi=config["figure_dpi"])

# Generate figure with time history of y-position.
fig = plt.figure()
plt.plot(state_history['t'],state_history['y'],color='k')
plt.tight_layout()
plt.savefig(output_path_prefix + "y-position-time-history.png", dpi=config["figure_dpi"])

# Generate figure with time history of z-position.
fig = plt.figure()
plt.plot(state_history['t'],state_history['z'],color='k')
plt.tight_layout()
plt.savefig(output_path_prefix + "z-position-time-history.png", dpi=config["figure_dpi"])

# Generate figure with time history of x-velocity.
fig = plt.figure()
plt.plot(state_history['t'],state_history['xdot'],color='k')
plt.tight_layout()
plt.savefig(output_path_prefix + "x-velocity-time-history.png", dpi=config["figure_dpi"])

# Generate figure with time history of y-velocity.
fig = plt.figure()
plt.plot(state_history['t'],state_history['ydot'],color='k')
plt.tight_layout()
plt.savefig(output_path_prefix + "y-velocity-time-history.png", dpi=config["figure_dpi"])

# Generate figure with time history of z-velocity.
fig = plt.figure()
plt.plot(state_history['t'],state_history['zdot'],color='k')
plt.tight_layout()
plt.savefig(output_path_prefix + "z-velocity-time-history.png", dpi=config["figure_dpi"])

print ("Figures generated successfully!")
print ("")

print ("")
print ("******************************************************************")
print ("                              Metrics                             ")
print ("******************************************************************")
print ("")

# Print percentage change in orbital energy
percentage_change_in_orbital_energy = max(abs(state_history['a']-state_history['a'][0]))/state_history['a'][0]*100.0
print ("Percentage change in orbital energy: {0} %".format(percentage_change_in_orbital_energy))

print ("")
print ("------------------------------------------------------------------")
print ("                         Exited successfully!                     ")
print ("------------------------------------------------------------------")
print ("")

# Stop timer
end_time = time.time( )

# Print elapsed time
print ("(Script time: " + str("{:,g}".format(end_time - start_time)) + "s)")
print ("")
