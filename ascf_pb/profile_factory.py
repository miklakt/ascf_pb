import importlib
import ascf_pb.solver
from ascf_pb.topology.kappa import kappa_plain

__keys_description = dict(
        N = 'chain length',
        sigma = 'grafting density',
        chi = 'Flory-Huggins parameter polymer-solvent',
        z = 'distance from grafting surface',
        pore_Radius = 'pore radius',
        R = 'distance from grafting surface to restriction',
        D = "polymer brush's thickness",
        phi = 'polymer concentration',
        Pi = 'osmotic pressure'
    )

def __get_required_keys(topology : str):
    topology_module = importlib.import_module('ascf_pb.topology.' + topology)
    return topology_module.required_keys[:]

def __generate_docstring(keys, return_):
    docstring = '\nCalculates a polymer brush density for given args'+\
        '\nArgs:\n'+\
        '\n'.join([f'{k} : {__keys_description[k]}' for k in keys])+\
        '\n\nReturns:'+\
        '\n(float): '.join(f'{__keys_description[return_]}')
    return docstring

def _D(kappa_cb, topology : str, **kwargs):
    topology_module = importlib.import_module('ascf_pb.topology.' + topology)
    D_cb = topology_module.D_universal
    
    kappa = kappa_cb(**kwargs)
    D = D_cb(kappa = kappa,**kwargs)
    return D

def _phi(kappa_cb, topology : str, **kwargs):
    topology_module = importlib.import_module('ascf_pb.topology.' + topology)
    D_cb = topology_module.D_universal
    phi_D_cb = topology_module.phi_D_universal
    
    kappa = kappa_cb(**kwargs)
    D = D_cb(kappa = kappa,**kwargs)
    phi_D = phi_D_cb(kappa = kappa, **kwargs)
    chi = kwargs['chi']
    z = kwargs['z']
    return ascf_pb.solver.Phi(z, chi, kappa, D, phi_D)

def _Pi(kappa_cb, topology : str, **kwargs):
    return ascf_pb.solver.Pi(_phi(kappa_cb, topology, **kwargs), kwargs['chi'])


def D(kappa_cb = kappa_plain, topology : str = 'plain', **kwargs): 
    required_keys = __get_required_keys(topology)
    required_keys.remove('z') 
    unused_keys = [k for k in required_keys if k not in kwargs]
    print ('Keys unused:', unused_keys)
    def wrapped(*new_args, **new_kwargs):
        unused_args = {k:v for k,v in zip(unused_keys, new_args)}
        unused_args.update(new_kwargs)
        unused_args.update(kwargs)
        return _D(kappa_cb = kappa_cb, topology = topology, **unused_args)
    wrapped.__doc__=__generate_docstring(unused_keys, 'D')
    return wrapped

def phi(kappa_cb = kappa_plain, topology : str = 'plain', **kwargs): 
    required_keys = __get_required_keys(topology)  
    unused_keys = [k for k in required_keys if k not in kwargs]
    print ('Keys unused:', unused_keys)
    def wrapped(*new_args, **new_kwargs):
        unused_args = {k:v for k,v in zip(unused_keys, new_args)}
        unused_args.update(new_kwargs)
        unused_args.update(kwargs)
        return _phi(kappa_cb = kappa_cb, topology = topology, **unused_args)
    wrapped.__doc__=__generate_docstring(unused_keys, 'phi')
    return wrapped

def Pi(kappa_cb = kappa_plain, topology : str = 'plain', **kwargs): 
    required_keys = __get_required_keys(topology)  
    unused_keys = [k for k in required_keys if k not in kwargs]
    print ('Keys unused:', unused_keys)
    def wrapped(*new_args, **new_kwargs):
        unused_args = {k:v for k,v in zip(unused_keys, new_args)}
        unused_args.update(new_kwargs)
        unused_args.update(kwargs)
        return _Pi(kappa_cb = kappa_cb, topology = topology, **unused_args)
    wrapped.__doc__=__generate_docstring(unused_keys, 'Pi')
    return wrapped