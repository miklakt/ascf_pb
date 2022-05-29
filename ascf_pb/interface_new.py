#%%
import ascf_pb.particle_geometry
import ascf_pb.free_energy
import ascf_pb.solver
import ascf_pb.topology
import importlib
import inspect
import numpy as np
from functools import lru_cache, partial, wraps
from collections import OrderedDict

#global variables
__topology_module = None
__particle_geometry_module = None
__callbacks = {}


__all__ = [
    "D",
    "phi",
    "Pi",
    #"critical_pore_radius",
    "osmotic_free_energy",
    "surface_free_energy",
    "total_free_energy",
    "set_config"
]


__arg_description = OrderedDict(
    N='polymer chain length',
    sigma='polymer brush grafting density',
    chi='Flory-Huggins parameter polymer-solvent',
    pore_Radius='grafted pore radius',
    R='distance from to grafting surface to the restriction',
    z='distance from grafting surface',
    ph="particle height",
    pw="particle width",
    pc="particle center",
    chi_PC='Flory-Huggins parameter polymer-particle',
    expansion_coefs="polynomial coefficients of gamma - phi dependency",
    eta="topological parameter",
)

__return_description = dict(
    D="polymer brush's thickness",
    osmotic_free_energy="osmotic part of the insertion free energy change",
    surface_free_energy="surface interaction contribution of the insertion free energy change",
    total_free_energy="total free energy of an inserted particle",
    phi='local polymer volume fraction',
    phi_D='local polymer volume fraction at the brush\'s end',
    Pi='local osmotic pressure',
    critical_pore_radius='minimal radius of an open pore'
)

__config = dict(
    vectorize=False,
    two_brackets=True,
    topology="plain",
    particle_geometry="cylinder"
)


def __load_modules():
    # loads specific modules for a given topology and particle geometry
    global __topology_module
    global __particle_geometry_module
    __topology_module = importlib.import_module(
        'ascf_pb.topology.' + __config["topology"])
    __particle_geometry_module = importlib.import_module(
        'ascf_pb.particle_geometry.' + __config["particle_geometry"])


def __ignore_extra_kwargs(func):
    """Decorator that allows to call functions with redundant kwargs, 
    so that the function only uses required ones.
    """
    parameters = inspect.signature(func).parameters

    def wrapped(**kwargs):
        new_kwargs = {k: kwargs[k] for k in kwargs if k in parameters}
        return func(**new_kwargs)
    wrapped.__annotations__ = func.__annotations__
    wrapped.__doc__ = func.__doc__
    wrapped.__name__ = func.__name__
    return wrapped


def __load_callbacks():
    # load specific for topology and particle geometry modules
    global __callbacks
    cb = dict(
        D=__topology_module.D_universal,
        phi_D=__topology_module.phi_D_universal,
        kappa=ascf_pb.topology.kappa.kappa,
        phi_solver=ascf_pb.solver.Phi,
        phi_to_Pi=ascf_pb.solver.Pi,
        integration_interval=__particle_geometry_module.get_integration_interval,
        volume=__particle_geometry_module.volume,
        surface=__particle_geometry_module.surface,
        volume_integrand=__particle_geometry_module.volume_integrand,
        surface_integrand=__particle_geometry_module.surface_integrand,
        #gamma_phi = ascf_pb.free_energy.gamma_phi()
    )
    if __config["topology"] == "pore":
        cb.update(dict(
            critical_pore_radius=ascf_pb.topology.pore.opening_pore_Radius
        ))
    __callbacks = dict((k, __ignore_extra_kwargs(v)) for k, v in cb.items())


def __generate_annotation(used_callbacks, exclude_keys=[]):
    # generate annotation based on used callbacks
    cb = [__callbacks[k] for k in used_callbacks]
    annotation = {param: T for func in cb
                  for param, T in func.__annotations__.items() if param in __arg_description}

    [annotation.pop(key, None) for key in exclude_keys]

    annotation_ordered = OrderedDict(
        ((k, annotation[k]) for k in __arg_description.keys() if k in annotation.keys()))

    return annotation_ordered


def __generate_docstring(annotation, return_name, return_type):
    # generate docstring for a given annotation
    docstring = f'\nCalculates {__return_description[return_name]} for given arguments' +\
        '\nArgs:\n' +\
        '\n'.join([f'{k} ({v.__name__}) : {__arg_description[k]}' for k, v in annotation.items()]) +\
        '\n\nReturns:' +\
        f'\n({return_type.__name__}): '+f'{__return_description[return_name]}'
    return docstring

# sort arguments in function signature consistent for all functions
# def __sort_args(func):
#    annotation = func.__annotations__
#
#    def wrapped(*args, **kwargs):
#        kwargs.update({k: arg for k, arg in zip(annotation.keys(), args)})
#        return func(**kwargs)
#    wrapped.__annotations__ = func.__annotations__
#    wrapped.__doc__ = func.__doc__
#    wrapped.__name__ = func.__name__
#    return wrapped


def __annotate(used_callbacks, return_name, return_type, add_annotation_from=None, exclude_keys=[]):
    # annotate function based on used callbacks and other used functions
    annotation = __generate_annotation(used_callbacks, exclude_keys)
    if add_annotation_from is not None:
        if isinstance(add_annotation_from, list):
            for ann in add_annotation_from:
                annotation.update(ann.__annotations__.items())
        else:
            annotation.update(add_annotation_from.__annotations__.items())

    [annotation.pop(key, None) for key in exclude_keys]

    annotation_ordered = OrderedDict(
        ((k, annotation[k]) for k in __arg_description.keys() if k in annotation.keys()))

    def decorator(func):
        func.__annotations__ = annotation_ordered
        func.__doc__ = __generate_docstring(
            annotation_ordered, return_name, return_type)
        return func
    return decorator


def _two_brackets(func):
    # modify function to use two bracket pairs f(a,b,c) -> f(a)(b,c) or f(a,b,c) -> f(b,c)(a)
    required_keys = func.__annotations__.keys()

    def first_brackets(**kwargs1):
        unused_keys = [k for k in required_keys if k not in kwargs1]

        def second_brackets(*args2, **kwargs2):
            # kwargs from args
            all_args = {k: v for k, v in zip(unused_keys, args2)}
            all_args.update(kwargs1)
            all_args.update(kwargs2)
            return func(**all_args)
        second_brackets.__doc__ = func.__doc__ + f"\n fixed:{kwargs1}"
        second_brackets.__annotations__ = {
            k: v for k, v in func.__annotations__.items() if k not in kwargs1}
        second_brackets.__name__ = func.__name__ + "_partial"
        return second_brackets
    first_brackets.__annotations__ = func.__annotations__
    first_brackets.__doc__ = func.__doc__
    first_brackets.__name__ = func.__name__
    if __config["vectorize"]:
            first_brackets =  np.vectorize(first_brackets)
    return first_brackets


def set_config(**kwargs):
    # set new config and reload all modules
    __config.update(kwargs)
    __load_modules()
    __load_callbacks()


set_config(
    topology="plain",
    particle_geometry="cylinder",
    two_brackets=True,
    vectorize=False
)


@__annotate(
    used_callbacks=["D", "kappa"],
    return_name="D",
    return_type=float
)
def D(**kwargs):
    kappa = __callbacks["kappa"](**kwargs)
    D = __callbacks["D"](kappa=kappa, **kwargs)
    return D


@__annotate(
    used_callbacks=["phi_D", "kappa"],
    return_name="phi_D",
    return_type=float,
    add_annotation_from=D
)
def phi_D(**kwargs):
    kappa = __callbacks["kappa"](**kwargs)
    phi_D = __callbacks["phi_D"](kappa=kappa, **kwargs)
    return phi_D


@__annotate(
    used_callbacks=["D", "phi_D", "kappa", "phi_solver"],
    return_name="phi",
    return_type=float
)
def phi(**kwargs):
    kappa = __callbacks["kappa"](**kwargs)
    D = __callbacks["D"](kappa=kappa, **kwargs)
    phi_D = __callbacks["phi_D"](kappa=kappa, **kwargs)
    phi_value = __callbacks["phi_solver"](
        kappa=kappa, d=D, phi_D=phi_D, **kwargs)

    return phi_value


@__annotate(
    used_callbacks=["phi_to_Pi"],
    return_name="Pi",
    return_type=float,
    add_annotation_from=phi
)
def Pi(**kwargs):
    kappa = __callbacks["kappa"](**kwargs)
    D = __callbacks["D"](kappa=kappa, **kwargs)
    phi_D = __callbacks["phi_D"](kappa=kappa, **kwargs)
    phi_value = __callbacks["phi_solver"](
        kappa=kappa, d=D, phi_D=phi_D, **kwargs)
    Pi_value = ascf_pb.solver.Pi(phi = phi_value, chi=kwargs['chi'])
    return Pi_value


if __config["topology"] == "pore":
    @__annotate(
        used_callbacks=["kappa", "critical_pore_radius"],
        return_name="critical_pore_radius",
        return_type=float,
        # add_annotation_from=phi
    )
    def critical_pore_radius(**kwargs):
        kappa = __callbacks["kappa"](**kwargs)
        pore_R = __callbacks["critical_pore_radius"](kappa=kappa, **kwargs)
        return pore_R


@__annotate(
    used_callbacks=["volume_integrand", "integration_interval"],
    return_name="osmotic_free_energy",
    return_type=float,
    add_annotation_from=[Pi, ascf_pb.free_energy.osmotic_free_energy],
    exclude_keys=["z"]
)
def osmotic_free_energy(**kwargs):
    Pi_integrand = partial(Pi, **kwargs)
    volume_integrand = __callbacks["volume_integrand"](**kwargs)
    z0, z1 = __callbacks["integration_interval"](**kwargs)
    fe = ascf_pb.free_energy.osmotic_free_energy(
        z0, z1, volume_integrand, Pi_integrand)
    return fe


@__annotate(
    used_callbacks=["surface_integrand", "integration_interval"],
    return_name="surface_free_energy",
    return_type=float,
    add_annotation_from=[
        phi,
        ascf_pb.free_energy.surface_free_energy,
    ],
    exclude_keys=["z"]
)
def surface_free_energy(**kwargs):
    phi_integrand = partial(phi, **kwargs)
    surface_integrand = __callbacks["surface_integrand"](**kwargs)
    z0, z1 = __callbacks["integration_interval"](**kwargs)
    fe = ascf_pb.free_energy.surface_free_energy(
        z0, z1,
        surface_integrand,
        phi_integrand,
    )
    return fe


@__annotate(
    used_callbacks=[],
    return_name="total_free_energy",
    return_type=float,
    add_annotation_from=[
        osmotic_free_energy,
        surface_free_energy,
    ]
)
def total_free_energy(**kwargs):
    sur = surface_free_energy(**kwargs)
    osm = surface_free_energy(**kwargs)
    return sur+osm

if __config["two_brackets"]:
    for k in __all__:
        globals()[k] = _two_brackets(globals()[k])
# %%
