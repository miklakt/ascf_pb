import importlib
import pathlib
import functools
import inspect
from typing import Callable, List
from collections import OrderedDict
import ascf_pb.topology
from ascf_pb import solver

def combine_signatures(
        callables: List[Callable],
        exclude=[],
        return_annotation=inspect.Signature.empty,
        positional_or_keyword_order : List[str] = None
):

    signatures = [inspect.signature(c) for c in callables]
    params = {key: param for sig in signatures for key,
              param in sig.parameters.items() if key not in exclude}


    if positional_or_keyword_order:
        #if any(params[key].default != inspect.Parameter.empty for key in positional_or_keyword_order)
        ordered_params = OrderedDict(params[key] for key in positional_or_keyword_order)
        
    
    signature = inspect.Signature(
        params.values(), return_annotation=return_annotation)

    return signature


def generate_docstring(signature: inspect.Signature, return_name):
    header = f'Calculates {__return_description[return_name]} for given arguments'

    args_desc = []
    for k, v in signature.parameters.items():
        if v.default is inspect.Parameter.empty:
            args_desc.append(
                f'{k} ({v.annotation.__name__}) : {__arg_description[k]}')
        else:
            args_desc.append(
                f'{k} ({v.annotation.__name__}, optional) : {__arg_description[k]}. Defaults to {v.default}')
    args_section = '\nArgs:\n' +\
        '\n'.join(args_desc)

    return_section = '\nReturns:' + \
        f'\n({signature.return_annotation.__name__}): ' + \
        f'{__return_description[return_name]}'

    return header+args_section+return_section


def generate_annotation(signature: inspect.Signature):
    annotations = {k: v.annotation for k, v in signature.parameters.items()}
    annotations["return"] = signature.return_annotation
    return annotations


def ignore_extra_kwargs(func : Callable) -> Callable:
    """Decorator that allows to call functions with redundant kwargs, 
    so that the function only uses required ones.
    """
    parameters = inspect.signature(func).parameters

    #if args_order is None:
    #    args_order = parameters.keys()
    #else:
    #    args_order = {k : parameters[k] for k in args_order if k in parameters}
    
    def wrapped(**kwargs):
        #new_args = {k : arg for k, arg in zip(args_order, args)}
        new_kwargs = {k: kwargs[k] for k in kwargs if k in parameters}
        print(func.__name__, new_kwargs)
        return func(**new_kwargs)
    #wrapped.__annotations__ = func.__annotations__
    wrapped.__doc__ = func.__doc__
    wrapped.__name__ = func.__name__
    wrapped.__module__ = func.__module__
    return wrapped


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
    g="dendron generations",
    q="dendron functionality"
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

# %%
class BrushPhiSolver:
    def __init__(self, geometry: str = 'plain') -> None:
        self.geometry = geometry
        self.__geometry_module = importlib.import_module(
            f"brush_geometry.{geometry}")
        self.phi_D_callback = self.__geometry_module.phi_D_universal
        self.D_callback = self.__geometry_module.D_universal
        self.kappa_callback = ascf_pb.topology.kappa
        
        self.phi = self.__build_phi_function()

    def __build_phi_function(self, positional_args_order = None) -> Callable:
        signature = combine_signatures(
            [self.phi_D_callback, self.D_callback, self.kappa_callback],
            exclude=["kappa"], 
            return_annotation=float
        )
        if positional_args_order is None:
            positional_args = signature.parameters
        else:
            positional_args = {k: signature.parameters[k] for k in positional_args_order}
        def phi(*args, **kwargs):
            print(kwargs)
            kappa = ignore_extra_kwargs(self.kappa_callback)(**kwargs)
            phi_D = ignore_extra_kwargs(self.phi_D_callback)(**kwargs, kappa=kappa)
            D = ignore_extra_kwargs(self.D_callback)(**kwargs, kappa=kappa)
            phi = ignore_extra_kwargs(solver.Phi)(**kwargs, d=D, phi_D = phi_D, kappa=kappa)
            return phi
        phi.__signature__ = signature
        phi.__annotations__ = generate_annotation(signature)
        phi.__doc__ = generate_docstring(signature, return_name="phi")
        return phi

b = BrushPhiSolver()
#%%
b.phi( z=10, chi = 0.5, N=1000, sigma = 0.02)
# %%
sig = combine_signatures([b.phi_D_callback, b.D_callback, ascf_pb.topology.kappa], exclude=[
                         "kappa"], return_annotation=float)
print(generate_docstring(sig, "phi"))
ann =generate_annotation(sig)
# %%
b.phi_D_callback.__annotations__ = ann
b.phi_D_callback.__signature__ = sig
# %%
inspect.signature(b.phi_D_callback)
# %%
order_args_ignore_extra_kwargs(b.D_callback)(0, 3.14/2000, 1000, 0.02, sigma =100)
# %%
