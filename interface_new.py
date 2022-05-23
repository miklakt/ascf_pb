# %%
import ascf_pb.particle_geometry
import ascf_pb.free_energy
import ascf_pb.solver
import ascf_pb.topology
import importlib
import inspect
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
    "pore_radius",
    "osmotic_free_energy",
    "surface_free_energy",
    "set_config"
]

#docstring and annotation generation
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
    eta="topological parameter"
)

__return_description = dict(
    D="polymer brush's thickness",
    osmotic_free_energy="osmotic part of the insertion free energy change",
    phi='local polymer volume fraction',
    phi_D = 'local polymer volume fraction at the brush\'s end',
    Pi='local osmotic pressure',
    critical_pore_radius='minimal radius of an open pore'
)

__config = dict(
    vectorize=False,
    cache=True,
    topology="plain",
    particle_geometry="cylinder"
)

#loads specific for topology and particle geometry modules
def __load_modules():
    global __topology_module
    global __particle_geometry_module
    __topology_module = importlib.import_module(
        'ascf_pb.topology.' + __config["topology"])
    __particle_geometry_module = importlib.import_module(
        'ascf_pb.particle_geometry.' + __config["particle_geometry"])

#make function ignore redundant kwargs
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
    return wrapped

#load specific for topology and particle geometry modules
def __load_callbacks():
    global __callbacks
    cb = dict(
        D=__topology_module.D_universal,
        phi_D=__topology_module.phi_D_universal,
        kappa=ascf_pb.topology.kappa.kappa,
        phi_solver=ascf_pb.solver.Phi,
        phi_to_Pi=ascf_pb.solver.Pi,
    )
    if __config["topology"] == "pore":
        cb.update(dict(
            critical_pore_radius
        ))
    __callbacks = dict((k, __ignore_extra_kwargs(v)) for k, v in cb.items())


def __generate_annotation(used_callbacks, exclude_keys=[]):
    cb = [__callbacks[k] for k in used_callbacks]
    annotation = {param: T for func in cb
                  for param, T in func.__annotations__.items() if param in __arg_description}

    [annotation.pop(key, None) for key in exclude_keys]

    annotation_ordered = OrderedDict(
        ((k, annotation[k]) for k in __arg_description.keys() if k in annotation.keys()))

    return annotation_ordered


def __generate_docstring(annotation, return_name, return_type):
    docstring = f'\nCalculates {__return_description[return_name]} for given arguments' +\
        '\nArgs:\n' +\
        '\n'.join([f'{k} ({v.__name__}) : {__arg_description[k]}' for k, v in annotation.items()]) +\
        '\n\nReturns:' +\
        f'\n({return_type.__name__}): '+f'{__return_description[return_name]}'
    return docstring


def __annotate(used_callbacks, return_name, return_type, add_annotation_from=None, exclude_keys=[]):
    annotation = __generate_annotation(used_callbacks, exclude_keys)
    if add_annotation_from is not None:
        annotation.update(add_annotation_from.__annotations__.items())
    
    annotation_ordered = OrderedDict(
        ((k, annotation[k]) for k in __arg_description.keys() if k in annotation.keys()))
    
    def decorator(func):
        func.__annotations__ = annotation_ordered
        func.__doc__ = __generate_docstring(
            annotation_ordered, return_name, return_type)
        return func
    return decorator


def __sorted_args(func):
    annotation = func.__annotations__
    def wrapped(*args, **kwargs):
        kwargs.update({k: arg for k, arg in zip(annotation.keys(), args)})
        return func(**kwargs)
    wrapped.__annotations__ = func.__annotations__
    wrapped.__doc__ = func.__doc__
    return wrapped
    

def set_config(**kwargs):
    __config.update(kwargs)
    __load_modules()
    __load_callbacks()


set_config(topology="pore",
    particle_geometry="cylinder")

# %%
@__sorted_args
@__annotate(
    used_callbacks=["D", "kappa"],
    return_name="D",
    return_type=float
)
def D(**kwargs):
    kappa = __callbacks["kappa"](**kwargs)
    D = __callbacks["D"](kappa=kappa, **kwargs)
    return D


@__sorted_args
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


@__sorted_args
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


@__sorted_args
@__annotate(
    used_callbacks=["phi_to_Pi"],
    return_name="Pi",
    return_type=float,
    add_annotation_from=phi
)
def Pi(**kwargs):
    return ascf_pb.solver.Pi(**kwargs), kwargs['chi']


@__sorted_args
@__annotate(
    used_callbacks=["phi_to_Pi"],
    return_name="critical_pore_radius",
    return_type=float,
    add_annotation_from=phi
)


def critical_pore_radius(**kwargs):
    return ascf_pb.solver.Pi(**kwargs), kwargs['chi']

# %%
