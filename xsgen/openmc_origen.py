"""An engine that uses OpenMC and ORIGEN2.2 to do burnup-criticality computations.
"""

from __future__ import print_function
import os
import subprocess

import numpy as np

from xsgen.statepoint import StatePoint

from pyne import rxname
from pyne import nucname
from pyne import origen22
from pyne.material import Material
from pyne.xs import data_source
from pyne.xs.cache import XSCache

from xsgen.utils import indir, NotSpecified

# templates are from openmc/examples/lattice/simple

SETTINGS_TEMPLATE = """<?xml version="1.0"?>
<settings>
  <cross_sections>{openmc_cross_sections}</cross_sections>
  <eigenvalue>
    <batches>{k_cycles}</batches>
    <inactive>{k_cycles_skip}</inactive>
    <particles>{k_particles}</particles>
  </eigenvalue>
  <source>
    <space type="box">
      <parameters>0 0 0 {unit_cell_height} {unit_cell_height} {unit_cell_height}</parameters>
    </space>
  </source>
</settings>
"""

MATERIALS_TEMPLATE = """<?xml version="1.0"?>
<materials>
  <default_xs>71c</default_xs>
  <material id="1">
    <density value="{fuel_density}" units="g/cc" />
    {_fuel_nucs}
  </material>
  <material id="2">
    <density value="{cool_density}" units="g/cc" />
    {_cool_nucs}
    <sab name="HH2O" xs="71t" />
  </material>
  <material id="3">
    <density value="{clad_density}" units="g/cc" />
    {_clad_nucs}
  </material>
</materials>
"""

GEOMETRY_TEMPLATE = """<?xml version="1.0"?>
<geometry>
  <cell id="1" fill="5" surfaces="1 -2 3 -4" />
  <cell id="101" universe="1" material="1" surfaces="-5" />
  <cell id="102" universe="1" material="void" surfaces="5 -6" />
  <cell id="103" universe="1" material="3" surfaces="6 -7" />
  <cell id="104" universe="1" material="2" surfaces="7" />
  <cell id="201" universe="2" material="2" surfaces="8" />
  <cell id="302" universe="3" material="3" surfaces="9" />
  <lattice id="5">
    <type>rectangular</type>
    <dimension>{_latt_shape0} {_latt_shape1}</dimension>
    <lower_left>-{_latt_x_half_pitch} -{_latt_y_half_pitch}</lower_left>
    <width>{unit_cell_pitch} {unit_cell_pitch}</width>
    <universes>
      {lattice}
    </universes>
  </lattice>
  <surface id="1" type="x-plane" coeffs="-{_latt_x_half_pitch}" boundary="reflective" />
  <surface id="2" type="x-plane" coeffs="{_latt_x_half_pitch}" boundary="reflective" />
  <surface id="3" type="y-plane" coeffs="-{_latt_y_half_pitch}" boundary="reflective" />
  <surface id="4" type="y-plane" coeffs="{_latt_y_half_pitch}" boundary="reflective" />
  <surface id="5" type="z-cylinder" coeffs="0.0 0.0 {fuel_cell_radius}" />
  <surface id="6" type="z-cylinder" coeffs="0.0 0.0 {void_cell_radius}" />
  <surface id="7" type="z-cylinder" coeffs="0.0 0.0 {clad_cell_radius}" />
  <surface id="8" type="z-cylinder" coeffs="0.0 0.0 0.0" />
  <surface id="9" type="z-cylinder" coeffs="0.0 0.0 0.0" />
</geometry>
"""

TALLIES_TEMPLATE = """<?xml version="1.0"?>
<tallies>
  <tally id="1">
    <label>flux</label>
    <filter type="energy" bins="{_egrid}" />
    <filter type="material" bins="1" />
    <scores>flux</scores>
    <nuclides>total</nuclides>
  </tally>
  <tally id="2">
    <label>cinderflux</label>
    <filter type="energy" bins="{_cds_egrid}" />
    <filter type="material" bins="1" />
    <scores>flux</scores>
    <nuclides>total</nuclides>
  </tally>
  <tally id="3">
    <label>eafflux</label>
    <filter type="energy" bins="{_eafds_egrid}" />
    <filter type="material" bins="1" />
    <scores>flux</scores>
    <nuclides>total</nuclides>
  </tally>
  <tally id="4">
    <label>omcflux</label>
    <filter type="energy" bins="{_omcds_egrid}" />
    <filter type="material" bins="1" />
    <scores>flux</scores>
    <nuclides>total</nuclides>
  </tally>
  <tally id="5">
    <label>s_gh</label>
    <filter type="energy" bins="{_egrid}" />
    <filter type="energyout" bins="{_egrid}" />
    <filter type="material" bins="1" />
    <scores>scatter</scores>
    <nuclides>total</nuclides>
  </tally>
</tallies>
"""

PLOTS_TEMPLATE = """<?xml version="1.0"?>
<plots>
  <plot id="1" color="mat">
    <origin>0. 0. 0.</origin>
    <width>{_latt_x_pitch} {_latt_y_pitch}</width>
    <pixels>600 600</pixels>
  </plot>
</plots>
"""

class OpenMCOrigen(object):
    """An that combines OpenMC for k-code calculations and ORIGEN for
    transmutation.
    """

    reactions = {rxname.id(_) for _ in ('total', 'absorption', 'gamma', 'gamma_1',
                 'z_2n', 'z_2n_1', 'z_3n', 'proton', 'alpha', 'fission')}

    def __init__(self, rc):
        self.rc = rc
        self.statelibs = {}
        self.builddir = 'build-' + rc.reactor
        if not os.path.isdir(self.builddir):
            os.makedirs(self.builddir)
        self.cinderds = data_source.CinderDataSource()
        self.eafds = data_source.EAFDataSource()
        self.omcds = data_source.OpenMCDataSource(
                        cross_sections=rc.openmc_cross_sections,
                        #src_group_struct=self.eafds.src_group_struct)
                        src_group_struct=np.logspace(1, -9, 1001))
        data_sources = [self.omcds]
        if not rc.is_thermal:
            data_sources.append(self.eafds)
        data_sources += [self.cinderds, data_source.SimpleDataSource(),
                         data_source.NullDataSource()]
        for ds in data_sources[1:]:
            ds.load()
        self.xscache = XSCache(data_sources=data_sources)

    def pwd(self, state, directory):
        """Path to directory we will be running specific physics codes in.

        Parameters
        ----------
        state : namedtuple (State)
            The state we are running physics codes for.
        directory : string
            The name of the sub-directory we would like.

        Returns
        -------
        str
            The path to the desired directory.
        """
        return os.path.join(self.builddir, str(hash(state)), directory)

    def context(self, state):
        """Unite parameters in the run-control file and  the current state.

        Parameters
        ----------
        state : namedtuple (State)
            A State tuple that contains perturbation parameters specific to that
            state.

        Returns
        -------
        ctx : dict
            A dictionary with all relevant perturbation parameters.
        """
        rc = self.rc
        ctx = dict(rc._dict)
        ctx.update(zip(rc.perturbation_params, state))
        return ctx

    def generate_run(self, run):
        """Generate transmutation tables, neutron production/destruction rates, and
        burnup statistics for a sequence of states with the same initial
        conditions.

        Parameters
        ----------
        run : list of States
            A list of States that has the same initial conditions at increasing
            burnup times.

        Returns
        -------
        libs : list of dicts
            Libraries to write out - one for the full fuel and one for each tracked nuclide.
        """
        libs = {"fuel": {
            "TIME": [],
            "NEUT_PROD": [],
            "NEUT_DEST": [],
            "BUd":  [],
            "transmutation": []}
        }
        libs.update(dict(zip(self.rc.track_nucs,
                             [{"TIME": [],
                               "NEUT_PROD": [],
                               "BUd": [],
                               "transmutation": []} for _ in self.rc.track_nucs])))
        for i, state in enumerate(run):
            if i < len(run) - 1 :
                transmute_time = run[i+1].burn_times - state.burn_times
            else:
                transmute_time = 0
            results = self.generate(state, transmute_time)
            # do stuff with results here:
            self.rc.fuel_material = results["fuel"]["transmutation"]
            libs = self._update_libs_with_results(libs, results)
            # looks like {"fuel": {"TIME": 0, "NEUT_PROD": ...}, 92350000:
            # {"TIME": ...}}
        nucs = []
        for mat in mat_hist:
            nucs.extend(mat.comp.keys())
        nucs = set(nucs)
        nucs = [(nucname.name(nuc), [mat.comp[nuc] for mat in mat_hist]) for nuc in nucs]
        nucs = dict(nucs)

        return libs

    def _update_libs_with_results(libs, results):
        """Update a set of libraries with results from a single timestep in-place.

        Parameters
        ----------
        libs : dict
            A set of libraries with nuclides as keys.
        results: dict
            A set of libraries for a single timestep.

        Returns
        -------
        libs : dict
            The updated library.
        """
        for mat in libs.keys():
            for key in libs[mat].keys():
                libs[mat][key].append(results[mat][key])
        return libs

    def generate(self, state, transmute_time):
        """Runs physics codes on a specific state. First runs OpenMC for transport,
        then uses the results to run ORIGEN for a single timestep.

        Parameters
        ----------
        state : namedtuple (State)
            A namedtuple containing the state parameters.
        transmute_time : float
            The length of the time step we would like to run ORIGEN for.

        Returns
        -------
        results : dict
            Dict of physics code results. Keys are either nuclide ID's or "fuel"
            for the full fuel results.
        """
        results = {"fuel": {}}
        results.update(dict(zip(self.rc.track_nucs, [{} for _ in self.rc.track_nucs])))
        if state in self.statelibs:
            return self.statelibs[state]
        k, phi_g, xs = self.openmc(state)
        self.statelibs[state] = (k, phi_g, xs)
        for key in results.keys():
            results[key] = self.origen(state, xs, transmute_time, phi_g, mat)
        # store what has become of each tracked nuclide
        # run origen on 1 kg of each nuclide, write out
        return results

    def openmc(self, state):
        """Runs OpenMC for a given state.

        Parameters
        ----------
        state : namedtuple (State)
            A namedtuple containing the state parameters.

        Returns
        -------
        k : float
            Neutron multiplication factor.
        phi_g : list of floats
            Group flux.
        xstab : list of tuples
            A list of tuples of the format (nuc, rx, xs).
        """
        # make inputs
        pwd = self.pwd(state, "omc")
        if not os.path.isdir(pwd):
            os.makedirs(pwd)
        self._make_omc_input(state)
        statepoint = _find_statepoint(pwd)
        if statepoint is None:
            with indir(pwd):
                subprocess.check_call(['openmc'])
            statepoint = _find_statepoint(pwd)
        # parse & prepare results
        k, phi_g = self._parse_statepoint(statepoint)
        xstab = self._generate_xs(phi_g)
        return k, phi_g, xstab

    def _make_omc_input(self, state):
        """Make OpenMC input files for a given state.

        Parameters
        ----------
        state : namedtuple (State)
            A namedtuple containing the state parameters.

        Returns
        -------
        None
        """
        pwd = self.pwd(state, "omc")
        ctx = self.context(state)
        rc = self.rc
        # settings
        settings = SETTINGS_TEMPLATE.format(**ctx)
        with open(os.path.join(pwd, 'settings.xml'), 'w') as f:
            f.write(settings)
        # materials
        valid_nucs = self.nucs_in_cross_sections()
        # core_nucs = set(ctx['core_transmute'])
        ctx['_fuel_nucs'] = _mat_to_nucs(rc.fuel_material[valid_nucs])
        ctx['_clad_nucs'] = _mat_to_nucs(rc.clad_material[valid_nucs])
        ctx['_cool_nucs'] = _mat_to_nucs(rc.cool_material[valid_nucs])
        materials = MATERIALS_TEMPLATE.format(**ctx)
        with open(os.path.join(pwd, 'materials.xml'), 'w') as f:
            f.write(materials)
        # geometry
        ctx['lattice'] = ctx['lattice'].strip().replace('\n', '\n      ')
        ctx['_latt_shape0'] = ctx['lattice_shape'][0]
        ctx['_latt_shape1'] = ctx['lattice_shape'][1]
        ctx['_latt_x_pitch'] = ctx['unit_cell_pitch'] * ctx['lattice_shape'][0]
        ctx['_latt_y_pitch'] = ctx['unit_cell_pitch'] * ctx['lattice_shape'][1]
        ctx['_latt_x_half_pitch'] = ctx['_latt_x_pitch'] / 2.0
        ctx['_latt_y_half_pitch'] = ctx['_latt_y_pitch'] / 2.0
        geometry = GEOMETRY_TEMPLATE.format(**ctx)
        with open(os.path.join(pwd, 'geometry.xml'), 'w') as f:
            f.write(geometry)
        # tallies
        ctx['_egrid'] = " ".join(map(str, sorted(ctx['group_structure'])))
        # ctx['_cds_egrid'] = " ".join(map(str, sorted(self.cinderds.src_group_struct)))
        ctx['_cds_egrid'] = " 1 10 100 1000" # I don't have cinder
        ctx['_eafds_egrid'] = " ".join(map(str, sorted(self.eafds.src_group_struct)))
        ctx['_omcds_egrid'] = " ".join(map(str, sorted(self.omcds.src_group_struct)))
        # nucs = core_nucs & valid_nucs
        # ctx['_nucs'] = " ".join([nucname.serpent(nuc) for nuc in sorted(nucs)])
        tallies = TALLIES_TEMPLATE.format(**ctx)
        with open(os.path.join(pwd, 'tallies.xml'), 'w') as f:
            f.write(tallies)
        # plots
        plots = PLOTS_TEMPLATE.format(**ctx)
        with open(os.path.join(pwd, 'plots.xml'), 'w') as f:
            f.write(plots)

    def nucs_in_cross_sections(self):
        """Returns the set of nuclides present in the cross_sections.xml file.

        Returns
        -------
        nucs : set
            A set of all nuclides known to OpenMC, in ID form.
        """
        return {n.nucid for n in self.omcds.cross_sections.ace_tables \
                if n.nucid is not None}

    def _parse_statepoint(self, statepoint):
        """Parses a statepoint file and reads in the relevant fluxes, assigns them
        to the DataSources or the XSCache, and returns k and phi_g.

        Parameters
        ----------
        statepoint : xsgen.statepoint.StatePoint
            An OpenMC StatePoint.

        Returns
        -------
        k : float
            Neutron multiplication factor.
        phi_g : list of floats
            Group flux.
        """
        sp = StatePoint(statepoint)
        sp.read_results()
        # compute group fluxes for data sources
        for tally, ds in zip(sp.tallies[1:4], (self.cinderds, self.eafds, self.omcds)):
            ds.src_phi_g = tally.results[::-1, :, 0].flatten()
            ds.src_phi_g /= ds.src_phi_g.sum()
        # compute return values
        k, kerr = sp.k_combined
        t_flux = sp.tallies[0]
        phi_g = t_flux.results[::-1, :, 0].flatten()
        phi_g /= phi_g.sum()
        return k, phi_g

    def _generate_xs(self, phi_g):
        """Grab xs data from cache depending on group flux.

        Parameters
        ----------
        phi_g : list of floats
            Group flux.

        Returns
        -------
        data : list of tuples
            A list of tuples of the format (nuc, rx, xs).
        """
        rc = self.rc
        verbose = rc.verbose
        xscache = self.xscache
        xscache['E_g'] = rc.group_structure
        xscache['phi_g'] = phi_g
        G = len(phi_g)
        temp = rc.temperature
        rxs = self.reactions
        nucs = rc.track_nucs
        dt = np.dtype([('nuc', 'i4'), ('rx', np.uint32), ('xs', 'f8', G)])
        data = np.empty(len(nucs)*len(rxs), dtype=dt)
        i = 0
        for nuc in nucs:
            for rx in rxs:
                xs = xscache[nuc, rx, temp]
                if verbose:
                    print("OpenMC XS:", nucname.name(nuc), rxname.name(rx), xs)
                data[i] = nuc, rx, xs
                i += 1
        return data

    def origen(self, state, transmute_time, phi_g):
        """Run ORIGEN on a state.

        Parameters
        ----------
        state : namedtuple (State)
            A namedtuple containing the state parameters.
        transmute_time : float
            Length of transmutation timestep.
        phi_g : list of floats
            Group flux.
        mat : str or int
            The material to start transmuting. Either "fuel" or a nuclide in ID
            form.

        Returns
        -------
        results : dict
            Dictionary with neutron production and destruction rates, burnup, and
            transmutation results.
        """
        pwd = self.pwd(state, "origen" + str(mat))
        if mat == "fuel":
            mat = self.rc.fuel_material
        else:
            mat = Material({mat: 1}, 1000, attrs={"units": "g"});

        if not os.path.isdir(pwd):
            os.makedirs(pwd)
        self._make_origen_input(state, transmute_time, phi_g, mat)

        if self.rc.origen_call is NotSpecified:
            if self.rc.is_thermal:
                origen_call = "o2_therm_linux.exe"
            else:
                origen_call = "o2_fast_linux.exe"
        else:
            origen_call = self.rc.origen_call

        with indir(pwd):
            subprocess.check_call(origen_call)
            tape6 = origen22.parse_tape6("TAPE6.OUT")
            material = tape6.materials[-1]
            burnup = tape6.burnup_MWD
            neutron_prod = tape6.neutron_production_rate
            neutron_dest = tape6.neutron_destruction_rate

        results = {
            "TIME": transmute_time,
            "NEUT_PROD": neutron_prod,
            "NEUT_DEST": neutron_dest,
            "BUd": burnup,
            "transmutation": material
            }
        return results

    def _make_origen_input(self, state, transmute_time, phi_g, mat):
        """Make ORIGEN input files for a given state.

        Parameters
        ----------
        state : namedtuple (State)
            A namedtuple containing the state parameters.
        mat : pyne.material.Material
            The fuel material to transmute.

        Returns
        -------
        None
        """
        pwd = self.pwd(state, "origen")
        ctx = self.context(state)
        with indir(pwd):
            # may need to filter tape4 for Bad Nuclides
            origen22.write_tape4(mat)
            total_flux = phi_g.sum()
            origen22.write_tape5_irradiation("IRF", transmute_time, total_flux)
            tape9 = origen22.make_tape9(self.rc.track_nucs)
            origen22.write_tape9(tape9)


def _mat_to_nucs(mat):
    """Convert a ``pyne.material.Material`` into OpenMC ``materials.xml`` format.

    Parameters
    ----------
    mat : ``pyne.material.Material``
        Material to convert.

    Returns
    -------
    nucs : string
        OpenMC-friendly XML tag with material composition.
    """
    nucs = []
    template = '<nuclide name="{nuc}" wo="{mass}" />'
    for nuc, mass in mat.comp.items():
        nucs.append(template.format(nuc=nucname.serpent(nuc), mass=mass*100))
    nucs = "\n    ".join(nucs)
    return nucs


def _find_statepoint(pwd):
    """Find a statepoint in a directory. Returns None if none found.

    Parameters
    ----------
    pwd : str
        Directory to search.

    Returns
    -------
    path : str or None
        The path of the statepoint directory, or None.
    """
    for f in os.listdir(pwd):
        if f.startswith('statepoint'):
            return os.path.join(pwd, f)
    return None
