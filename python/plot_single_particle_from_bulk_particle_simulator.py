'''
Copyright (c) 2009-2018, K. Kumar (me@kartikkumar.com)
Distributed under the MIT License.
See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
'''

# Set up modules and packages.
# Plotting
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import cm
from matplotlib.font_manager import FontProperties
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d
import matplotlib.animation as animation

# I/O
import json
import jstyleson
from pprint import pprint
import sqlite3

# Numerical
import math
import numpy as np
import pandas as pd

# System
import sys
import time

print ("")
print ("------------------------------------------------------------------")
print ("                             dustsim                              ")
print ("      Copyright (c) 2009-2018, K. Kumar (me@kartikkumar.com)      ")
print ("------------------------------------------------------------------")
print ("")

# Start timer.
start_time = time.time( )

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
  json_input_string = json_input_file.read()
config = jstyleson.loads(json_input_string)
jstyleson.dumps(config)
pprint(config)

print ("")
print ("******************************************************************")
print ("                            Operations                            ")
print ("******************************************************************")
print ("")

print ("Fetching data from database ...")

# Connect to SQLite database.
try:
    database = sqlite3.connect(config['database'])

except sqlite3.Error as error:

    print ("Error %s:" % error.args[0])
    sys.exit(1)

metadata = pd.read_sql("SELECT * FROM " + config['metadata_table'], database)
initial_state = pd.read_sql("SELECT * FROM " + config['initial_states_table'] + " WHERE simulation_id = " + str(config['simulation_id']), database)
simulation_results = pd.read_sql("SELECT * FROM " + config['simulation_results_table'] + " WHERE simulation_id = " + str(config['simulation_id']), database)

print ("Data successfully fetched!")
print ("")

print ("Generating figures ...")

# Pre-compute useful variables.
output_path_prefix = config["output_directory"] + '/'
simulation_time_in_years = simulation_results['time'] / (365.25*24.0*3600.0)

# Generate semi-major axis change figure.
fig = plt.figure()
ax1 = fig.add_subplot(2, 3, 1)
ax2 = fig.add_subplot(2, 3, 2)
ax3 = fig.add_subplot(2, 3, 3)
ax4 = fig.add_subplot(2, 3, 4)
ax5 = fig.add_subplot(2, 3, 5)
ax6 = fig.add_subplot(2, 3, 6)

# Plot semi-major axis time history.
ax1.set_xlabel(r'$t$ [s]')
ax1.set_ylabel(r'$a$ [km]')
ax1.grid()
ax1.plot(simulation_time_in_years, simulation_results['semi_major_axis'],color='k')

# Plot eccentricity time history.
ax2.set_xlabel(r'$t$ [s]')
ax2.set_ylabel(r'$e$ [-]')
ax2.grid()
ax2.plot(simulation_time_in_years, simulation_results['eccentricity'],color='k')

# Plot inclination time history.
ax3.set_xlabel(r'$t$ [s]')
ax3.set_ylabel(r'$i$ [deg]')
ax3.grid()
ax3.plot(simulation_time_in_years, simulation_results['inclination'].apply(math.degrees),color='k')

# Plot argument of periapsis time history.
ax4.set_xlabel(r'$t$ [s]')
ax4.set_ylabel(r'$\omega$ [deg]')
ax4.grid()
ax4.plot(simulation_time_in_years, simulation_results['argument_of_periapsis'].apply(math.degrees),color='k')

# Plot longitude of ascending node time history.
ax5.set_xlabel(r'$t$ [s]')
ax5.set_ylabel(r'$\Omega$ [deg]')
ax5.grid()
ax5.plot(simulation_time_in_years, simulation_results['longitude_of_ascending_node'].apply(math.degrees),color='k')

# Plot true anomaly time history.
ax6.set_xlabel(r'$t$ [s]')
ax6.set_ylabel(r'$\theta$ [deg]')
ax6.grid()
ax6.plot(simulation_time_in_years, simulation_results['true_anomaly'].apply(math.degrees),color='k')

# Save figure.
plt.tight_layout()
plt.savefig(output_path_prefix + config["eccentricity_change_figure"], dpi=config["figure_dpi"])


# number_of_simulations = len(initial_states)
# times = simulation_results[simulation_results['simulation_id'] == 1 ]['time']
# zero_change = pd.DataFrame(np.zeros((number_of_simulations)))

# # Generate histograms for initial_states in Keplerian elements.
# fig = plt.figure()
# ax1 = fig.add_subplot(2, 3, 1)
# ax2 = fig.add_subplot(2, 3, 2)
# ax3 = fig.add_subplot(2, 3, 3)
# ax4 = fig.add_subplot(2, 3, 4)
# ax5 = fig.add_subplot(2, 3, 5)
# ax6 = fig.add_subplot(2, 3, 6)

# # Plot semi-major axis histogram.
# ax1.hist(initial_states['semi_major_axis']-config['semi_major_axis_reference'],color='k')
# ax1.set_xlabel('a [km]')
# ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
# ax1.grid()

# # Plot eccentricity histogram.
# ax2.hist(initial_states['eccentricity'],color='k')
# ax2.set_xlabel('e [-]')
# ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
# ax2.grid()

# # Plot inclination histogram.
# ax3.hist(initial_states['inclination'].apply(math.degrees),color='k')
# ax3.set_xlabel('i [deg]')
# ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
# ax3.grid()

# # Plot argument of periapsis histogram.
# ax4.hist(initial_states['argument_of_periapsis'].apply(math.degrees),color='k')
# ax4.set_xlabel(r'$\omega$ [deg]')
# ax4.grid()

# # Plot longitude of ascending node histogram.
# ax5.hist(initial_states['longitude_of_ascending_node'].apply(math.degrees),color='k')
# ax5.set_xlabel(r'$\Omega$ [deg]')
# ax5.grid()

# # Plot true anomaly histogram.
# ax6.hist(initial_states['true_anomaly'].apply(math.degrees),color='k')
# ax6.set_xlabel(r'$\theta$ [deg]')
# ax6.grid()

# # Save figure.
# plt.tight_layout()
# plt.savefig(output_path_prefix + config["initial_states_figure"], dpi=config["figure_dpi"])

# # Generate semi-major axis change figure.
# fig = plt.figure()
# plt.xlabel(r'$a_{0}$ [km]')
# plt.ylabel(r'$\Delta a$ [km]')
# plt.grid()

# plt.plot(initial_states['semi_major_axis'],zero_change,marker='.',color='k',linestyle='None')

# simulation_id_mask = simulation_results['simulation_id'] == 1
# simulation_times = simulation_results['time'][simulation_id_mask]
# semi_major_axis_change = simulation_results['semi_major_axis'][simulation_id_mask] - simulation_results['semi_major_axis'][simulation_id_mask][0]

# print(initial_states['semi_major_axis'][initial_states['simulation_id'] == 1][0])

# # plt.plot(initial_states['semi_major_axis'],zero_change,marker='.',color='k',linestyle='None')


# # for x in range(1,len(times)):
# #     # a_time = simulation_results[simulation_results['time'] == times[x-1]]['time']
#     # semi_major_axis = simulation_results[simulation_results['time'] == times[x-1]]['semi_major_axis']
# #     print(simulation_results['time'])
# #     # semi_major_axis_change = pd.DataFrame(semi_major_axis.values-initial_states['semi_major_axis'].values)
# #     # plt.plot(initial_states['semi_major_axis'],semi_major_axis_change,marker='.',color='k',linestyle='None')

# # Save figure.
# plt.tight_layout()
# plt.savefig(output_path_prefix + config["semi_major_axis_change_figure"], dpi=config["figure_dpi"])

# # # Generate eccentricity change figure.
# # fig = plt.figure()
# # plt.xlabel(r'$e_{0}$ [-]')
# # plt.ylabel(r'$\Delta e$ [-]')
# # plt.grid()

# # plt.plot(initial_states['eccentricity'],zero_change,marker='.',color='k',linestyle='None')

# # for x in range(1,len(times)):
# #     # a_time = simulation_results[simulation_results['time'] == times[x-1]]['time']
# #     eccentricity = simulation_results[simulation_results['time'] == times[x-1]]['eccentricity']
# #     eccentricity_change = pd.DataFrame(eccentricity.values-initial_states['eccentricity'].values)
# #     plt.plot(initial_states['eccentricity'],eccentricity_change,marker='.',color='k',linestyle='None')

# # # Save figure.
# # plt.tight_layout()
# # plt.savefig(output_path_prefix + config["eccentricity_change_figure"], dpi=config["figure_dpi"])

print ("Figures generated successfully!")
print ("")

# # print ("Generating animation ...")

# # # Generate animation of change in Keplerian elements.
# # fig = plt.figure()
# # # plt.tight_layout()
# # ax1 = fig.add_subplot(2, 3, 1)
# # ax2 = fig.add_subplot(2, 3, 2)
# # ax3 = fig.add_subplot(2, 3, 3)
# # ax4 = fig.add_subplot(2, 3, 4)
# # ax5 = fig.add_subplot(2, 3, 5)
# # ax6 = fig.add_subplot(2, 3, 6)

# # # Generate animation of semi-major axis change.
# # ax1.set_xlim(metadata['semi_major_axis_minimum'][0], metadata['semi_major_axis_maximum'][0])
# # ax1.set_ylim(-5.0,5.0)
# # ax1.set_xlabel(r'$a_{0}$ [km]')
# # ax1.set_ylabel(r'$\Delta a$ [km]')
# # line1, = ax1.plot([],[],marker='o',color='k',linestyle='None')

# # # ax2.set_xlim(0.0, 1.0e-2)
# # # ax2.set_ylim(-1, 1)
# # ax2.set_xlabel(r'$e_{0}$ [-]')
# # ax2.set_ylabel(r'$\Delta e$ [-]')
# # line2, = ax2.plot([],[],marker='o',color='k',linestyle='None')

# # # Set up animation functions.
# # def init():
# #     ax1.plot(initial_states['semi_major_axis'],zero_change,marker='o',color='k',linestyle='None')
# #     ax2.plot(initial_states['eccentricity'],zero_change,marker='o',color='k',linestyle='None')

# # def animate(i):
# #     a_time = simulation_results[simulation_results['time'] == times[i]]['time']

# #     semi_major_axis = simulation_results[simulation_results['time'] == times[i]]['semi_major_axis']
# #     semi_major_axis_change = pd.DataFrame(semi_major_axis.values-initial_states['semi_major_axis'].values)
# #     line1.set_data(initial_states['semi_major_axis'],semi_major_axis_change)

# #     eccentricity = simulation_results[simulation_results['time'] == times[i]]['eccentricity']
# #     eccentricity_change = pd.DataFrame(eccentricity.values-initial_states['eccentricity'].values)
# #     line2.set_data(initial_states['eccentricity'],eccentricity_change)

# #     return line1, line2

# # # Generate and save animation.
# # animation_data = animation.FuncAnimation(fig,animate,init_func=init,blit=False,frames=len(times))
# # animation_path = output_path_prefix + "test.mp4"
# # animation_data.save(animation_path,fps=50,bitrate=6000)

# # print ("Animation generated successfully!")
# # print ("")

# # if config['show_figures']:
# #     plt.show()

# Stop timer
end_time = time.time( )

print ("")
print ("------------------------------------------------------------------")
print ("                         Exited successfully!                     ")
print ("------------------------------------------------------------------")
print ("")

# Print (elapsed time)
print ("(Script time: " + str("{:,g}".format(end_time - start_time)) + "s)")
print ("")
