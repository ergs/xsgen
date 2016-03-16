from __future__ import print_function
import os
import json
import subprocess
from pprint import pformat
from multiprocessing import Pool

import numpy as np

from openmc import statepoint

from pyne import rxname
from pyne import nucname
from pyne import origen22
from pyne.material import Material
from pyne.xs import data_source
from pyne.xs.cache import XSCache
from pyne.bins import stair_step

from matplotlib import pyplot as plt

from xsgen.utils import indir, NotSpecified
from xsgen.tape9 import brightlitetape9

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
  <energy_grid>{energy_grid}</energy_grid>
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
    <sab name="{sab}" xs="{sab_xs}" />
  </material>
  <material id="3">
    <density value="{clad_density}" units="g/cc" />
    {_clad_nucs}
  </material>
</materials>
"""  # TODO we need to make sab a thing in the rc...

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
    <label>eafflux</label>
    <filter type="energy" bins="{_eafds_egrid}" />
    <filter type="material" bins="1" />
    <scores>flux</scores>
    <nuclides>total</nuclides>
  </tally>
  <tally id="3">
    <label>omcflux</label>
    <filter type="energy" bins="{_omcds_egrid}" />
    <filter type="material" bins="1" />
    <scores>flux</scores>
    <nuclides>total</nuclides>
  </tally>
  <tally id="4">
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
        self.eafds = data_source.EAFDataSource()
        self.omcds = data_source.OpenMCDataSource(
                        cross_sections=rc.openmc_cross_sections,
                        src_group_struct=rc.openmc_group_struct)
        data_sources = [self.omcds]
        if not rc.is_thermal:
            data_sources.append(self.eafds)
        data_sources += [data_source.SimpleDataSource(),
                         data_source.NullDataSource()]
        for ds in data_sources[1:]:
            ds.load()
        self.xscache = XSCache(data_sources=data_sources)
        self.tape9 = None

        if self.rc.origen_call is NotSpecified:
            if self.rc.is_thermal:
                self.origen_call = "o2_therm_linux.exe"
            else:
                self.origen_call = "o2_fast_linux.exe"
        else:
            self.origen_call = self.rc.origen_call

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
        self.libs = {'xs': [], 'phi_g': {
            'E_g': {'EAF': self.eafds.src_group_struct.tolist(),
                    'OpenMC': self.omcds.src_group_struct.tolist()},
            'phi_g': []}, 
          "fuel": {
            "TIME": [0],
            "NEUT_PROD": [0],
            "NEUT_DEST": [0],
            "BUd":  [0],
            "material": [self.rc.fuel_material],
            "tracked_nucs": {nucname.name(n): [self.rc.fuel_material.comp.get(n, 0) * 1000]
                             for n in self.rc.track_nucs},
            "phi_tot": [0]
            }}

        for nuc in self.rc.track_nucs:
            self.libs[nuc] = {
                "TIME": [0],
                "NEUT_PROD": [0],
                "NEUT_DEST": [0],
                "BUd": [0],
                "material": [Material({nuc: 1}, 1000)],
                "tracked_nucs": {nucname.name(n): [0]
                                 for n in self.rc.track_nucs},
                "phi_tot": [0]
                }
            self.libs[nuc]["tracked_nucs"][nucname.name(nuc)] = [1000]

        print([state.burn_times for state in run])
        for i, state in enumerate(run):
            if i > 0:
                transmute_time = state.burn_times - run[i-1].burn_times
                results = self.generate(state, transmute_time)
                self.libs = self._update_libs_with_results(self.libs, results)
        return self.libs

    def _update_libs_with_results(self, matlibs, newlibs):
        """Update a set of libraries with results from a single timestep in-place.

        Parameters
        ----------
        matlibs : dict
            A set of libraries with nuclides as keys.
        newlibs: dict
            A set of libraries for a single timestep, with nuclides as keys.

        Returns
        -------
        libs : dict
            The updated library.
        """
        for mat, newlib in newlibs.items():
            if mat == 'xs':
                matlibs[mat].append(newlib)
                continue
            elif mat == 'phi_g':
                matlibs[mat][mat].append(newlib)
                continue
            oldlib = matlibs[mat]
            for nuc in self.rc.track_nucs:
                name = nucname.name(nuc)
                nuc_frac = newlib["material"].comp.get(nuc, 0)
                mass = newlib["material"].mass
                oldlib["tracked_nucs"][name].append(nuc_frac * mass)
            for key, value in newlib.items():
                if isinstance(key, int):
                    continue
                else:
                    oldlib[key].append(value)
        return matlibs

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
        print("generating for a state with transmute_time {}".format(transmute_time))
        if state in self.statelibs:
            return self.statelibs[state]
        rc = self.rc
        k, phi_g, xstab = self.openmc(state)
        results = {"fuel": {}}
        results.update(dict(zip(rc.track_nucs, [{} for _ in rc.track_nucs])))
        if 'flux' in rc:
            phi_tot = state.flux
        elif 'fuel_specific_power' in rc:
            raise RuntimeError('needs refactor for state')
            G = len(phi_g)
            fission_id = rxname.id("fission")
            if G == 1:
                fission_xs = {xs[0]: xs[2] * 1e-24 for xs in xstab  # xs is in barns not cm2
                              if xs[1] == fission_id}
            else:
                fission_xs = {xs[0]: np.sum(xs[2]) * 1e-24 / len(xs[2]) 
                              for xs in xstab  # xs is in barns not cm2
                              if xs[1] == fission_id}
            fuel_material = self.libs["fuel"]["material"][-1]
            fuel_material.atoms_per_molecule = sum([self.rc.fuel_chemical_form[m]
                                                    for m in self.rc.fuel_chemical_form])
            fuel_atom_frac = fuel_material.to_atom_frac()
            fuel_number_density = 6.022e23 * self.rc.fuel_density / \
                (fuel_material.molecular_mass() * fuel_material.atoms_per_molecule)
            number_densities = {nuc: fuel_atom_frac[nuc] * fuel_number_density
                                for nuc in fuel_material.comp}
            sum_N_i_sig_fi = sum([number_densities[nuc] * fission_xs.get(nuc, 0)
                                  for nuc in fuel_material.comp])
            sum_N_i_sig_fi = sum_N_i_sig_fi[sum_N_i_sig_fi != 0]
            fuel_specific_power_mwcc = self.rc.fuel_density * 1e-6 * state.fuel_specific_power
            # see http://iriaxp.iri.tudelft.nl/~leege/SCALE44/origens.PDF for formula
            # (search for "the specific power due to fission", on p. 22 of the PDF)
            phi_tot = sum(3.125e16*fuel_specific_power_mwcc/sum_N_i_sig_fi)
        results = self.run_all_the_origens(state, transmute_time, phi_tot, results)
        results['xs'] = xstab
        results['phi_g'] = {'EAF': self.eafds.src_phi_g.tolist(),
                            'OpenMC': self.omcds.src_phi_g.tolist()}
        self.statelibs[state] = results
        return results

    def run_all_the_origens(self, state, transmute_time, phi_tot, results):
        """Call ORIGEN as much as necessary and unite the results.

        Parameters
        ----------
        state : namedtuple (State)
            A namedtuple containing the state parameters.
        transmute_time : float
            The length of the transmutation timestep. Has units of [days].
        phi_tot : float
            The total neutron flux.
        results : dict
            A dict with material identifiers as keys, and dictionaries as
            values. The basic data structure to fill.

        Returns
        -------
        dict
           A dict of all the ORIGEN results.
        """
        if self.rc.verbose:
            print("making tape9 for {0} with phi={1}".format(state, phi_tot))
        self.tape9 = origen22.make_tape9(self.rc.track_nucs, self.xscache, nlb=(219, 220, 221))
        self.tape9 = origen22.merge_tape9((self.tape9,
                                          origen22.loads_tape9(brightlitetape9)))
        origen22.write_tape9(self.tape9)
        for mat_id in results.keys():
            pwd = self.pwd(state, "origen{}".format(mat_id))
            mat = self.libs[mat_id]["material"][-1]
            if not os.path.isdir(pwd):
                os.makedirs(pwd)
            with indir(pwd):
                if not os.path.isfile("TAPE6.OUT"):
                    self._make_origen_input(transmute_time, phi_tot, mat)
        origen_results = []
        if self.rc.threads == 1:
            for mat_id in results.keys():
                pwd = self.pwd(state, "origen{}".format(mat_id))
                origen_params = (state.burn_times,
                                 transmute_time,
                                 phi_tot,
                                 mat_id,
                                 self.libs[mat_id]["material"][-1],
                                 pwd,
                                 self.origen_call)
                origen_results.append(_origen(origen_params))
        else:
            origen_params_ls = [(state.burn_times,
                                 transmute_time,
                                 phi_tot,
                                 mat_id,
                                 self.libs[mat_id]["material"][-1],
                                 self.pwd(state, "origen{}".format(mat_id)),
                                 self.origen_call)
                                for mat_id in results]
            pool = Pool(self.rc.threads)
            origen_results = pool.map(_origen, origen_params_ls)
            pool.close()
            pool.join()
        for result in origen_results:
            result[1]["material"] = Material(dict(result[1]["material"]),
                                             1000,
                                             attrs={"units": "g"})
        return dict(origen_results)

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
                subprocess.check_call(['openmc', '-s', '{}'.format(self.rc.threads)])
            statepoint = _find_statepoint(pwd)
        # parse & prepare results
        k, phi_g, e_g = self._parse_statepoint(statepoint)
        if self.rc.plot_group_flux:
            with indir(pwd):
                self._plot_group_flux(e_g, phi_g)
        xstab = self._generate_xs(e_g, phi_g)
        return k, phi_g, xstab

    def _plot_group_flux(self, e_g, phi_g):
        """Plot the group flux output by OpenMC and save plot to file.

        Parameters
        ----------
        e_g : array
            Energy bins
        phi_g : array
            Flux values
        """
        fig = plt.figure(figsize=(12, 8))
        plt.loglog(*stair_step(e_g, phi_g), figure=fig)
        plt.title("Flux vs Energy in")
        plt.xlabel('E [MeV]')
        plt.ylabel('Flux [N/cm$^2\cdot$s]')
        plt.savefig("flux")
        plt.close()

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
        # discard Cd-119m1 as a valid nuc
        valid_nucs.discard(481190001)
        # core_nucs = set(ctx['core_transmute'])
        #ctx['_fuel_nucs'] = _mat_to_nucs(rc.fuel_material[valid_nucs])
        curr_fuel = self.libs['fuel']['material'][-1][valid_nucs]
        ctx['_fuel_nucs'] = _mat_to_nucs(curr_fuel[valid_nucs])        
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
        return {n.nucid for n in self.omcds.cross_sections.ace_tables
                if n.nucid is not None}

    def _parse_statepoint(self, statepoint_path, tally_id=1):
        """Parses a statepoint file and reads in the relevant fluxes, assigns them
        to the DataSources or the XSCache, and returns k, phi_g, and E_g.

        Parameters
        ----------
        statepoint : xsgen.statepoint.StatePoint
            An OpenMC StatePoint.
        tally_id : int
            The tally id we wish to read group flux from.

        Returns
        -------
        k : float
            Neutron multiplication factor.
        phi_g : list of floats
            Group flux.
        e_g : list of floats
            Group structure.
        """
        sp = statepoint.StatePoint(statepoint_path)
        print(sp.tallies[0])
        # compute group fluxes for data sources
        for tally, ds in zip(sp.tallies[1:3], (self.eafds, self.omcds)):
            ds.src_phi_g = tally.results[::-1, :, 0].flatten()
            ds.src_phi_g /= ds.src_phi_g.sum()
        # compute return values
        k, kerr = sp.k_combined
        tally = sp.tallies[tally_id]
        phi_g = tally.results[::-1, :, 0].flatten()
        phi_g /= phi_g.sum()
        e_g = tally.filters["energyin"].bins
        print("test " + k.to_string())
        return k, phi_g, e_g

    def _generate_xs(self, e_g, phi_g):
        """Grab xs data from cache depending on group flux.

        Parameters
        ----------
        e_g : list of floats
            Group structure.
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
        xscache.clear()
        xscache['E_g'] = e_g
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
                    print("OpenMC XS:", nucname.name(nuc), rxname.name(rx), xs, temp)
                data[i] = nuc, rx, xs
                i += 1
        return data

    def _make_origen_input(self, transmute_time, phi_tot, mat):
        """Make ORIGEN input files for a given state.

        Parameters
        ----------
        transmute_time : float
            The time, in days, to run the transmutation for.
        phi_tot: float
            Total neutron flux.
        mat : pyne.material.Material
            The fuel material to transmute.
        Returns
        -------
        None
        """
        # may need to filter tape4 for Bad Nuclides
        # if sum(mat.comp.values()) > 1:
        origen22.write_tape4(mat)
        origen22.write_tape5_irradiation("IRF",
                                         transmute_time,
                                         phi_tot,
                                         xsfpy_nlb=(219, 220, 221),
                                         cut_off=1e-300)
        origen22.write_tape9(self.tape9)


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


def _origen(origen_params):
    """Run ORIGEN on a state.

    Parameters
    ----------
    origen_params : dict
        A dictionary containing the following parameters:

        abs_time : float
            The absolute time (relative to 0 days) that the transmutation occurs at.
        transmute_time : float
            Length of transmutation timestep - time relative to absolute time of
            last timestep.
        phi_tot : float
            Total neutron flux.
        mat_id : str or int
            The identifier of a material to start transmuting. Either "fuel" or a
            nuclide in ID form.
        mat : pyne.material.Material
            The material to be transmuted.
        pwd : str
            The directory to run ORIGEN in.
        origen_call : str
            The shell command for running ORIGEN.

    Returns
    -------
    results : dict
        Dictionary with neutron production and destruction rates, burnup, and
        transmutation results.

    """
    abs_time, transmute_time, phi_tot, mat_id, mat, pwd, origen_call = origen_params
    mat.mass = 1000
    mat.attrs = {"units": "g"}

    with indir(pwd):
        if not os.path.isfile("TAPE6.OUT"):
            times_called = 0
            while times_called < 3:
                times_called += 1
                try:
                    subprocess.check_call(origen_call)
                    break
                except subprocess.CalledProcessError:
                    print("Warning: ORIGEN2.2 in " + pwd + "failed. Retrying.")

        print("Parsing " + pwd + "/TAPE6.OUT...")
        tape6 = origen22.parse_tape6("TAPE6.OUT")

    out_mat = tape6["materials"][-1]
    out_mat.mass = 1000
    out_mat.comp = {n: frac for n, frac in out_mat.comp.items() if frac != 0}
    if mat_id == "fuel":
        out_mat.atoms_per_molecule = 3
    burnup = tape6["burnup_MWD"][-1]
    neutron_prod = tape6["neutron_production_rate"][-1]
    neutron_dest = tape6["neutron_destruction_rate"][-1]

    results = (mat_id, {
        "TIME": abs_time,
        "NEUT_PROD": neutron_prod,
        "NEUT_DEST": neutron_dest,
        "BUd": burnup,
        "material": list(out_mat.comp.items()),
        "phi_tot": phi_tot
        })
    #if burnup < 0.0:
    #    msg = 'Negative burnup found for {0}:\n{1}'
    #    msg = msg.format(mat_id, pformat(results[1]))
    #    raise ValueError(msg)
    return results
