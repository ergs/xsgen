# This is a default RunControl file which should get you started. Edit as you
# see fit.

burn_times = [0, 3, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
group_structure = [0, 10, 100, 1000]

track_nucs = ["U235", "U238"]

# in g/cc
fuel_density = 19.1
clad_density = 6.56
cool_density = 1.0

# in cm
fuel_cell_radius = 0.7
void_cell_radius = 0.8
clad_cell_radius = 0.9
unit_cell_pitch = 1.5

burn_regions = 1

# in W/g
fuel_specific_power = 1 # maybe a gross underestimation

fuel_material = {'U238': 0.96, 'U235': 0.04}

# currently available formats: "brightlite"
formats = ["brightlite"]

# currently available solvers: "openmc+origen"
solver = "openmc+origen"

reactor = "lwr"
