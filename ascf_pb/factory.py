"""
Provides interface to calculate brush thickness, 
polymer density and osmotic pressure profiles.

The script allows users to create closures(callbacks) for such calculations.
One has to provide system geometry (default - 'plain') 
and then can create a closure object by calling coresponding factory-function 
to study dependencies of brush thickness, 
polymer density and osmotic pressure profiles from varying arguments.

Recommended steps is to create closure, 
then call it with the rest of the arguments.
Here are possible ways to get a result:
    -----------------------------------------------------
    closure = factory(**kwargs_to_bind)
    result = closure(**unused_kwargs)
    -----------------------------------------------------
    result = factory(**kwargs_to_bind)(**unused_kwargs)
    -----------------------------------------------------
    result = factory(**all_required_kwargs)()
    -----------------------------------------------------

e.g to calculate local polymer density (phi) the next keys has are required: 
['N', 'sigma', 'chi', 'z'].
Thus, calling factory.phi gives us a closure object that can be used to study 
polymer density for a given distance (z) at varying grafting densities (sigma).
    -----------------------------------------------------
    #generate closure
    phi_from_sigma = factory.phi(N=1000, chi=0.3, z=10)
    #grafting density we want to make calculations with
    sigma = [0.01, 0.02, 0.03, 0.04]
    #list of resulting local polymer densities
    phi_varying_sigma = [phi_from_sigma(s) for s in sigma]
    -----------------------------------------------------

@author: Mikhail Laktionov
miklakt@gmail.com
"""
import importlib
import inspect
from functools import lru_cache

import ascf_pb.solver
from ascf_pb.topology import kappa

#used to generate descriptions and docstrings
__keys_description = dict(
        N = 'polymer chain length',
        sigma = 'polymer brush grafting density',
        chi = 'Flory-Huggins parameter polymer-solvent',
        z = 'distance from grafting surface',
        pore_Radius = 'grafted pore radius',
        R = 'distance from to grafting surface to the restriction',
        D = "polymer brush's thickness",
        phi = 'local polymer volume fraction',
        Pi = 'local osmotic pressure'
    )

def __get_required_keys(topology : str) -> list:
    """Returns keys required to calculate a local polymer density
    or local osmotic pressure, by inspecting module. An importing module 
    contains required_keys list with the keys.

    Args:
        topology (str): topology module name, geometry of the system (plain, pore)

    Returns:
        list: required keys
    """    
    topology_module = importlib.import_module('ascf_pb.topology.' + topology)
    return topology_module.required_keys[:]

def __generate_docstring(keys, return_):
    docstring = f'\nCalculates {__keys_description[return_]} for given args'+\
        '\nArgs:\n'+\
        '\n'.join([f'{k} : {__keys_description[k]}' for k in keys])+\
        '\n\nReturns:'+\
        '\n(float): '+f'{__keys_description[return_]}'
    return docstring

def __ignore_extra_kwargs(func):
    """Decorator that allows to call functions with redundant kwargs, 
    so that the function only uses required ones.
    """    
    parameters = inspect.signature(func).parameters
    def wrapped(**kwargs):
        new_kwargs = {k:kwargs[k] for k in kwargs if k in parameters}
        return func(**new_kwargs)
    return wrapped


################################################################################
def _D(kappa_cb, topology : str, **kwargs) -> float:
    """Calculates polymer brush thickness, imports necessary functions from
    coresponding topology module.

    Args:
        kappa_cb ([type]): kappa parameter callback
        topology (str): topology module name, geometry of the system (plain, pore)

    Returns:
        float : polymer brush's thickness
    """
    #import module with routines    
    topology_module = importlib.import_module('ascf_pb.topology.' + topology)
    #import callbacks to calculate brush thickness, 
    #ensure it ignores redundant kwargs
    D_cb = __ignore_extra_kwargs(topology_module.D_universal)
    #call kappa parameter callback immediately to get kappa parameter value
    kappa = __ignore_extra_kwargs(kappa_cb)(**kwargs)
    #calculate polymer brush's thickness
    D = D_cb(kappa = kappa,**kwargs)
    return D

def _phi(kappa_cb, topology : str, **kwargs) -> float:
    """Calculates local polymer density, imports necessary functions from
    coresponding topology module.

    Args:
        kappa_cb (Callable): kappa parameter callback
        topology (str): topology module name, geometry of the system (plain, pore)

    Returns:
        float: local polymer density
    """
    #import module with routines  
    topology_module = importlib.import_module('ascf_pb.topology.' + topology)
    #import callbacks to calculate brush thickness, local polymer density 
    #ensure they ignore redundant kwargs
    D_cb = __ignore_extra_kwargs(topology_module.D_universal)
    phi_D_cb = __ignore_extra_kwargs(topology_module.phi_D_universal)
    #call kappa parameter callback immediately to get kappa parameter value
    kappa = __ignore_extra_kwargs(kappa_cb)(**kwargs)
    #calculate polymer brush's thickness
    D = D_cb(kappa = kappa,**kwargs)
    #calculate polymer density at the grafting surface
    phi_D = phi_D_cb(kappa = kappa, **kwargs)

    chi = kwargs['chi']
    z = kwargs['z']
    #calculate phi with obtained D and phi_D
    return ascf_pb.solver.Phi(z, chi, kappa, D, phi_D)

def _Pi(kappa_cb, topology : str, **kwargs) -> float:
    """Calculates local osmotic pressure, imports necessary functions from
    coresponding topology module.

    Args:
        kappa_cb (Callable): kappa parameter callback
        topology (str): topology module name, geometry of the system (plain, pore)

    Returns:
        float: local polymer density
    """ 
    return ascf_pb.solver.Pi(_phi(kappa_cb, topology, **kwargs), kwargs['chi'])


def D(kappa_cb = kappa.kappa, topology : str = 'plain', **kwargs):
    """Factory function produces closure with bound parameters provided in
    keyword arguments. Unused keys has to be provided in sequential calls.
    
    Examples:
        ---------------------------------------------------
        closure = factory(**kwargs_to_bind)
        result = closure(**unused_kwargs)
        ---------------------------------------------------
        result = factory(**kwargs_to_bind)(**unused_kwargs)
        ---------------------------------------------------
        result = factory(**all_required_kwargs)()

    Args:
        kappa_cb (Callable, optional): kappa parameter callback. Defaults to kappa.kappa.
        topology (str, optional): topology module name, geometry of the system (plain, pore). Defaults to 'plain'.

    Returns:
        Callable: closure to calculate brush thickness
    """
    #check the comments in function phi, to know what's going on
    required_keys = __get_required_keys(topology)
    required_keys.remove('z') 
    unused_keys = [k for k in required_keys if k not in kwargs]
    print (f'Closure created, bound keys: {list(kwargs.keys())}')
    print (f'Keys unused: {unused_keys}')
    @lru_cache()
    def wrapped(*new_args, **new_kwargs):
        unused_args = {k:v for k,v in zip(unused_keys, new_args)}
        unused_args.update(new_kwargs)
        unused_args.update(kwargs)
        return _D(kappa_cb = kappa_cb, topology = topology, **unused_args)
    wrapped.__doc__=__generate_docstring(unused_keys, 'D')
    return wrapped

def phi(kappa_cb = kappa.kappa, topology : str = 'plain', **kwargs):
    """Factory function produces closure with bound parameters provided in
    keyword arguments. Unused keys has to be provided in sequential calls.
    
    Examples:
        ---------------------------------------------------
        closure = factory(**kwargs_to_bind)
        result = closure(**unused_kwargs)
        ---------------------------------------------------
        result = factory(**kwargs_to_bind)(**unused_kwargs)
        ---------------------------------------------------
        result = factory(**all_required_kwargs)()

    Args:
        kappa_cb (Callable, optional): kappa parameter callback. Defaults to kappa.kappa.
        topology (str, optional): topology module name, geometry of the system (plain, pore). Defaults to 'plain'.

    Returns:
        Callable: closure to calculate local polymer density
    """
    #get all keys required to finnish calculation
    required_keys = __get_required_keys(topology)
    #compare it with provided in the function call
    unused_keys = [k for k in required_keys if k not in kwargs]
    print (f'Closure created, bound keys: {list(kwargs.keys())}')
    print (f'Keys unused: {unused_keys}')
    #create cached inner function to use as closure
    @lru_cache()
    def wrapped(*new_args, **new_kwargs):
        #turn args to kwargs in an order of required_keys
        #and collected all args/kwargs 
        #from all the calls (**kwargs)(*new_args, **new_kwargs)
        unused_args = {k:v for k,v in zip(unused_keys, new_args)}
        unused_args.update(new_kwargs)
        unused_args.update(kwargs)
        #get final result
        return _phi(kappa_cb = kappa_cb, topology = topology, **unused_args)
    #automated docstring for closure, to make use of built-in help()
    wrapped.__doc__=__generate_docstring(unused_keys, 'phi')
    #note the absence of (), we return a closure, not the final result
    return wrapped

def Pi(kappa_cb = kappa.kappa, topology : str = 'plain', **kwargs): 
    """Factory function produces closure with bound parameters provided in
    keyword arguments. Unused keys has to be provided in sequential calls.
    
    Examples:
        ---------------------------------------------------
        closure = factory(**kwargs_to_bind)
        result = closure(**unused_kwargs)
        ---------------------------------------------------
        result = factory(**kwargs_to_bind)(**unused_kwargs)
        ---------------------------------------------------
        result = factory(**all_required_kwargs)()

    Args:
        kappa_cb (Callable, optional): kappa parameter callback. Defaults to kappa.kappa.
        topology (str, optional): topology module name, geometry of the system (plain, pore). Defaults to 'plain'.

    Returns:
        Callable: closure to calculate local osmotic pressure
    """
    #check the comments in function phi, to know what's going on
    required_keys = __get_required_keys(topology)  
    unused_keys = [k for k in required_keys if k not in kwargs]
    print (f'Closure created, bound keys: {list(kwargs.keys())}')
    print (f'Keys unused: {unused_keys}')
    @lru_cache()
    def wrapped(*new_args, **new_kwargs):
        unused_args = {k:v for k,v in zip(unused_keys, new_args)}
        unused_args.update(new_kwargs)
        unused_args.update(kwargs)
        return _Pi(kappa_cb = kappa_cb, topology = topology, **unused_args)
    wrapped.__doc__=__generate_docstring(unused_keys, 'Pi')
    return wrapped


################################################################################
def pore_radius(kappa_cb = kappa.kappa, **kwargs):
    """Factory function produces closure with bound parameters provided in
    keyword arguments. Unused keys has to be provided in sequential calls.
    
    Examples:
        ---------------------------------------------------
        closure = factory(**kwargs_to_bind)
        result = closure(**unused_kwargs)
        ---------------------------------------------------
        result = factory(**kwargs_to_bind)(**unused_kwargs)
        ---------------------------------------------------
        result = factory(**all_required_kwargs)()

    Args:
        kappa_cb (Callable, optional): kappa parameter callback. Defaults to kappa.kappa.
        topology (str, optional): topology module name, geometry of the system (plain, pore). Defaults to 'plain'.

    Returns:
        Callable: closure to calculate critical radius of a pore, when it starts to open.
    """
    #check the comments in function phi, to know what's going on
    from ascf_pb.topology.pore import opening_pore_Radius
    required_keys = inspect.signature(opening_pore_Radius).parameters
    unused_keys = [k for k in required_keys if k not in kwargs]
    unused_keys.remove('kappa')
    print (f'Closure created, bound keys: {kwargs.keys()}')
    print (f'Keys unused: {unused_keys}')
    def wrapped(*new_args, **new_kwargs):
        unused_args = {k:v for k,v in zip(unused_keys, new_args)}
        unused_args.update(new_kwargs)
        unused_args.update(kwargs)
        unused_args['kappa'] = __ignore_extra_kwargs(kappa_cb)(**unused_args)
        return opening_pore_Radius(**unused_args)
    wrapped.__doc__=__generate_docstring(unused_keys, 'Pi')
    return wrapped
