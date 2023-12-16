"""
========================================================================================================================
SimulationConfig
========================================================================================================================
"""

from qiskit_metal.analyses.quantization import EPRanalysis
from utils import *
from sweeper_helperfunctions import *
from qiskit_metal.analyses.quantization import LOManalysis

class SimulationConfig:
    def __init__(self, design_name="CavitySweep", renderer_type="hfss", sim_type="eigenmode",
                 setup_name="Setup", max_passes=50, max_delta_f=0.05, min_converged_passes=2, Lj=0, Cj=0):
        self.design_name = design_name
        self.renderer_type = renderer_type
        self.sim_type = sim_type
        self.setup_name = setup_name
        self.max_passes = max_passes
        self.max_delta_f = max_delta_f
        self.min_converged_passes = min_converged_passes
        self.Lj = Lj
        self.Cj = Cj

def get_sim_results(emode_df, lom_df):
    data_emode = emode_df["sim_results"]
    data_lom = lom_df["data"]

    data = {}

    cross2cpw = abs(lom_df["data"]["cross_to_claw"]) * 1e-15
    cross2ground = abs(lom_df["data"]["cross_to_ground"]) * 1e-15
    f_r = emode_df["sim_results"]["cavity_frequency"]
    Lj = lom_df["design_options"]["aedt_q3d_inductance"] * (1 if lom_df["design_options"]["aedt_q3d_inductance"] > 1e-9 else 1e-9)
    print(Lj)
    gg, aa, ff_q = find_g_a(cross2cpw, cross2ground, f_r, Lj, N=4)
    data = dict(
        f_cavity = f_r,
        Q = emode_df["sim_results"]["Q"],
        kappa = emode_df["sim_results"]["kappa"],
        g = gg,
        a = aa,
        f_qubit = ff_q
    )

    return data


def run_eigenmode(filename):
    design = metal.designs.design_planar.DesignPlanar()
    gui = metal.MetalGUI(design)
    design.overwrite_enabled = True

    sim_json = open(filename)
    sim_data = json.load(sim_json)
    geometry_dict = sim_data["design_options"]["geometry_dict"]

    cpw_length = int("".join(filter(str.isdigit, geometry_dict["cpw_opts"]["total_length"])))
    claw = create_claw(geometry_dict["claw_opts"], cpw_length, design)
    coupler = create_coupler(geometry_dict["cplr_opts"], design)
    cpw = create_cpw(geometry_dict["cpw_opts"], coupler, design)
    config = SimulationConfig(min_converged_passes=1)

    epra, hfss = start_simulation(design, config)
    setup = set_simulation_hyperparameters(epra, config)
    epra.sim.renderer.options.max_mesh_length_port = '7um'

    render_simulation_with_ports(epra, config.design_name, setup.vars, coupler)
    modeler = hfss.pinfo.design.modeler

    mesh_lengths = {'mesh1': {"objects": [f"prime_cpw_{coupler.name}", f"second_cpw_{coupler.name}", f"trace_{cpw.name}", f"readout_connector_arm_{claw.name}"], "MaxLength": '7um'}}
    #add_ground_strip_and_mesh(modeler, coupler, mesh_lengths=mesh_lengths)
    print(mesh_lengths)
    mesh_objects(modeler, mesh_lengths)
    f_rough, Q, kappa = get_freq_Q_kappa(epra, hfss)

    data = epra.get_data()

    data_df = {
        "design_options": {
            "coupling_type": "CLT",
            "geometry_dict": geometry_dict
        },
        "sim_options": {
            "sim_type": "epr",
            "setup": setup,
        },
        "sim_results": {
            "cavity_frequency": f_rough,
            "Q": Q,
            "kappa": kappa
        },
        "misc": data
    }

    return data_df

def run_LOM(filename):
    design = metal.designs.design_planar.DesignPlanar()
    gui = metal.MetalGUI(design)
    design.overwrite_enabled = True
    
    c1 = LOManalysis(design, "q3d")

    c1.sim.setup.reuse_selected_design = False
    c1.sim.setup.reuse_setup = False

    # example: update single setting
    c1.sim.setup.max_passes = 30
    c1.sim.setup.min_converged_passes = 1
    c1.sim.setup.percent_error = 0.1
    c1.sim.setup.name = 'sweep_setup'

    sim_json = open(filename)
    sim_data = json.load(sim_json)
    cross_dict = sim_data["design_options"]["geometry_dict"]["claw_opts"]

    # return cross_dict["connection_pads"].keys()

    qname = 'xmon'
    cnames = cross_dict["connection_pads"].keys()
    cname = list(cnames)[0]

    temp_arr = np.repeat(qname, len(cnames))
    ports_zip = zip(temp_arr, cnames)

    # xmon1_options["cross_length"] = f"{cross_length}um"
    # xmon1_options["connection_pads"]["readout"]["claw_length"] = f"{claw_length}um"
    # xmon1_options["connection_pads"]["readout"]["ground_spacing"] = f"{ground_spacing}um"
    # xmon1_options["aedt_hfss_inductance"]=Lj,
    # xmon1_options["hfss_inductance"]=Lj, 
    # xmon1_options["q3d_inductance"]=Lj,
    # xmon1_options["aedt_q3d_inductance"]=Lj
    q = TransmonCross(design, qname, options=cross_dict)
    design.rebuild()
    selection = [qname]
    open_pins = ports_zip #[(qname, cname)]
    print(q.options)
    c1.sim.renderer.clean_active_design()
    c1.sim.run(name = 'LOMv2.0', components=selection,
               open_terminations=open_pins)
    cap_df = c1.sim.capacitance_matrix

    # print(f'For Qubit_{i}, the cap matrix is...')
    print(cap_df)

    data = {
        "design_options": design.components[qname].options,
        "data": {
            "cross_to_ground": 0 if 'ground_main_plane' not in cap_df.loc[f'cross_{qname}'] else cap_df.loc[f'cross_{qname}']['ground_main_plane'],
            "claw_to_ground": 0 if 'ground_main_plane' not in cap_df.loc[f'{cname}_connector_arm_{qname}'] else cap_df.loc[f'{cname}_connector_arm_{qname}']['ground_main_plane'],
            "cross_to_claw": cap_df.loc[f'cross_{qname}'][f'{cname}_connector_arm_{qname}'],
            "cross_to_cross": cap_df.loc[f'cross_{qname}'][f'cross_{qname}'],
            "claw_to_claw": cap_df.loc[f'{cname}_connector_arm_{qname}'][f'{cname}_connector_arm_{qname}'],
            "ground_to_ground": 0 if 'ground_main_plane' not in cap_df.loc[f'cross_{qname}'] else cap_df.loc['ground_main_plane']['ground_main_plane']
        },
        "sim_info": {
            "setup": c1.sim.setup,
            "renderer_options": c1.sim.renderer.options
        }
    }
    # save_simulation_data_to_json(data, filename = f"qubitonly_num{i}_{comp_id}_v{version}")
    return data

def CLT_epr_sweep(design, sweep_opts, filename):    
    for param in extract_QSweep_parameters(sweep_opts):
        # if int("".join(filter(str.isdigit, param["cpw_opts"]["total_length"]))) < 2000:
        # param["claw_opts"].update({"pos_x": ("-1000um" if int("".join(filter(str.isdigit, param["cpw_opts"]["total_length"]))) < 2000 else "-1500um") })
        cpw_length = int("".join(filter(str.isdigit, param["cpw_opts"]["total_length"])))
        claw = create_claw(param["claw_opts"], cpw_length, design)
        coupler = create_coupler(param["cplr_opts"], design)
        cpw = create_cpw(param["cpw_opts"], coupler, design)
        # gui.rebuild()
        # gui.autoscale()

        config = SimulationConfig(min_converged_passes=3)

        epra, hfss = start_simulation(design, config)
        setup = set_simulation_hyperparameters(epra, config)
        epra.sim.renderer.options.max_mesh_length_port = '7um'

        render_simulation_with_ports(epra, config.design_name, setup.vars, coupler)
        modeler = hfss.pinfo.design.modeler

        mesh_lengths = {'mesh1': {"objects": [f"prime_cpw_{coupler.name}", f"second_cpw_{coupler.name}", f"trace_{cpw.name}", f"readout_connector_arm_{claw.name}"], "MaxLength": '7um'}}
        #add_ground_strip_and_mesh(modeler, coupler, mesh_lengths=mesh_lengths)
        print(mesh_lengths)
        mesh_objects(modeler, mesh_lengths)
        f_rough, Q, kappa = get_freq_Q_kappa(epra, hfss)

        data = epra.get_data()

        data_df = {
            "design_options": {
                "coupling_type": "CLT",
                "geometry_dict": param
            },
            "sim_options": {
                "sim_type": "epr",
                "setup": setup,
            },
            "sim_results": {
                "cavity_frequency": f_rough,
                "Q": Q,
                "kappa": kappa
            },
            "misc": data
        }
        
        # filename = f"CLT_cpw{cpw.options.total_length}_claw{claw.options.connection_pads.readout.claw_width}_clength{coupler.options.coupling_length}"
        save_simulation_data_to_json(data_df, filename)

def NCap_epr_sweep(design, sweep_opts):    
    for param in extract_QSweep_parameters(sweep_opts):
        claw = create_claw(param["claw_opts"], design)
        coupler = create_coupler(param["cplr_opts"], design)
        cpw = create_cpw(param["cpw_opts"], design)
        # gui.rebuild()
        # gui.autoscale()
        
        config = SimulationConfig()

        epra, hfss = start_simulation(design, config)
        setup = set_simulation_hyperparameters(epra, config)
        
        render_simulation_no_ports(epra, [cpw,claw], [(cpw.name, "start")], config.design_name, setup.vars)
        modeler = hfss.pinfo.design.modeler

        mesh_lengths = {'mesh1': {"objects": [f"trace_{cpw.name}", f"readout_connector_arm_{claw.name}"], "MaxLength": '4um'}}
        mesh_objects(modeler,  mesh_lengths)
        f_rough = get_freq(epra, hfss)

        data = epra.get_data()

        data_df = {
            "design_options": {
                "coupling_type": "NCap",
                "geometry_dict": param
            },
            "sim_options": {
                "sim_type": "epr",
                "setup": setup,
            },
            "sim_results": {
                "cavity_frequency": f_rough
            },
            "misc": data
        }
        
        filename = f"CLT_cpw{cpw.options.total_length}_claw{claw.options.connection_pads.readout.claw_width}_clength{coupler.options.coupling_length}"
        save_simulation_data_to_json(data_df, filename)

def NCap_LOM_sweep(design, sweep_opts):
    for param in extract_QSweep_parameters(sweep_opts):
        # claw = create_claw(param["claw_opts"], design)
        coupler = create_coupler(param, design)
        # coupler.options[""]
        # cpw = create_cpw(param["cpw_opts"], design)
        # gui.rebuild()
        # gui.autoscale()

        loma = LOManalysis(design, "q3d")
        loma.sim.setup.reuse_selected_design = False
        loma.sim.setup.reuse_setup = False

        # example: update single setting
        loma.sim.setup.max_passes = 30
        loma.sim.setup.min_converged_passes = 5
        loma.sim.setup.percent_error = 0.1
        loma.sim.setup.auto_increase_solution_order = 'False'
        loma.sim.setup.solution_order = 'Medium'

        loma.sim.setup.name = 'lom_setup'

        loma.sim.run(name = 'LOMv2.01', components=[coupler.name],
        open_terminations=[(coupler.name, pin_name) for pin_name in coupler.pin_names])
        cap_df = loma.sim.capacitance_matrix
        data = loma.get_data()
        setup = loma.sim.setup

        data_df = {
            "design_options": {
                "coupling_type": "NCap",
                "geometry_dict": param
            },
            "sim_options": {
                "sim_type": "lom",
                "setup": setup,
            },
            "sim_results": {
                "C_top2top" : abs(cap_df[f"cap_body_0_{coupler.name}"].values[0]),
                "C_top2bottom" : abs(cap_df[f"cap_body_0_{coupler.name}"].values[1]),
                "C_top2ground" : abs(cap_df[f"cap_body_0_{coupler.name}"].values[2]),
                "C_bottom2bottom" : abs(cap_df[f"cap_body_1_{coupler.name}"].values[1]),
                "C_bottom2ground" : abs(cap_df[f"cap_body_1_{coupler.name}"].values[2]),
                "C_ground2ground" : abs(cap_df[f"ground_main_plane"].values[2]),
            },
            "misc": data
        }

        filename = f"NCap_LOM_fingerwidth{coupler.options.cap_width}_fingercount{coupler.options.finger_count}_fingerlength{coupler.options.finger_length}_fingergap{coupler.options.cap_gap}"
        save_simulation_data_to_json(data_df, filename)

def start_simulation(design, config):
    """
    Starts the simulation with the specified design and configuration.

    :param design: The design to be simulated.
    :param config: The configuration settings for the simulation.
    :return: A tuple containing the EPR analysis object and the HFSS object.
    """
    epra = EPRanalysis(design, config.renderer_type)
    hfss = epra.sim.renderer
    print("Starting the Simulation")
    hfss.start()
    hfss.new_ansys_design(config.design_name, config.sim_type)
    return epra, hfss


def set_simulation_hyperparameters(epra, config):
    """
    Sets the simulation hyperparameters based on the provided configuration.

    :param epra: The EPR analysis object.
    :param config: The configuration settings for the simulation.
    :return: The setup object with the updated settings.
    """
    setup = epra.sim.setup
    setup.name = config.setup_name
    setup.max_passes = config.max_passes
    setup.max_delta_f = config.max_delta_f
    setup.min_converged = config.min_converged_passes
    setup.n_modes = 1
    setup.vars = {'Lj': f'{config.Lj}nH', 'Cj': f'{config.Cj}fF'}
    return setup


def render_simulation_with_ports(epra, ansys_design_name, setup_vars, coupler):
    """
    Renders the simulation into HFSS.

    :param epra: The EPR analysis object.
    :param ansys_design_name: The name of the Ansys design.
    :param setup_vars: The setup variables for the rendering.
    :param coupler: The coupler object.
    """
    print(epra.sim)
    epra.sim.renderer.clean_active_design()
    epra.sim._render(name=ansys_design_name,
                     solution_type='eigenmode',
                     vars_to_initialize=setup_vars,
                     open_pins=[(coupler.name, "prime_start"), (coupler.name, "prime_end")],
                     port_list=[(coupler.name, 'prime_start', 50), (coupler.name, "prime_end", 50)],
                     box_plus_buffer=True)
    print("Sim rendered into HFSS!")

def render_simulation_no_ports(epra, components, open_pins, ansys_design_name, setup_vars):
    """
    Renders the simulation into HFSS.

    :param epra: The EPR analysis object.
    :param ansys_design_name: The name of the Ansys design.
    :param setup_vars: The setup variables for the rendering.
    :param components: List of QComponent object.
    """
    epra.sim._render(name=ansys_design_name,
                     selection=[qcomp.name for qcomp in components],
                     open_pins=open_pins,
                     solution_type='eigenmode',
                     vars_to_initialize=setup_vars,
                     box_plus_buffer=True)
    print("Sim rendered into HFSS!")


if __name__ == "__main__":
    # Usage
    config = SimulationConfig()
    bbox = generate_bbox(coupler)
    epra, hfss = start_simulation(design, config)
    setup = set_simulation_hyperparameters(epra, config)
    render_simulation(epra, config.design_name, setup.vars, coupler)
    modeler = hfss.pinfo.design.modeler
