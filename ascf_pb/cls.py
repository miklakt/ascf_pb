#%%
import importlib
from typing import Callable
import ascf_pb.topology
from ascf_pb import solver
import inspect
import numpy as np
import itertools
import functools

from ascf_pb.decorators import doc_decorator, ignore_extra_kwargs
import ascf_pb.decorators

brush_solver_default_key_order = ["chi", "N", "sigma"]



class BrushSolver:
    def __init__(self, geometry: str = 'plain', vectorize = True, v_product = True) -> None:
        self.geometry = geometry
        self.__geometry_module = importlib.import_module(
            f"ascf_pb.brush_geometry.{geometry}")
        self.phi_D_callback = functools.lru_cache()(self.__geometry_module.phi_D_universal)
        self.D_callback = functools.lru_cache()(self.__geometry_module.D_universal)
        self.kappa_callback = functools.lru_cache()(ascf_pb.topology.kappa)
        #self.__vectorize = vectorize
        #self.__product = v_product
        if vectorize:
            self.__deco = ascf_pb.decorators.vectorize(v_product)
        else:
            self.__deco = ascf_pb.decorators.dummy
        self.build_functions()
        
    
    def build_functions(self):
        self.D = self.__deco(self._D())
        self.phi_D = self.__deco(self._phi_D())
        self.phi = self.__deco(self._phi())
        self.Pi = self.__deco(self._Pi())

    
    def _D(self) -> Callable:
        # generate signature combining all signatures of all functions 
        # called in the body of phi function
        @doc_decorator(
            #return_name="D",
            callables=[self.D_callback, self.kappa_callback],
            exclude=["kappa"],
            return_annotation=float,
            positional_or_keyword_order=brush_solver_default_key_order
        )
        def D(**kwargs):
            kappa = ignore_extra_kwargs(self.kappa_callback)(**kwargs)
            D_ = ignore_extra_kwargs(self.D_callback)(**kwargs, kappa=kappa)
            return D_
        return D


    def _phi_D(self) -> Callable:
        # generate signature combining all signatures of all functions 
        # called in the body of phi function
        @doc_decorator(
            callables=[self.D_callback, self.kappa_callback],
            exclude=["kappa"],
            return_annotation=float,
            positional_or_keyword_order=brush_solver_default_key_order
        )
        def phi_D(**kwargs):
            kappa = ignore_extra_kwargs(self.kappa_callback)(**kwargs)
            _phi_D = ignore_extra_kwargs(self.phi_D_callback)(**kwargs, kappa=kappa)
            return _phi_D
        return phi_D


    def _phi(self) -> Callable:
        @doc_decorator(
            callables=[self.D, self.kappa_callback, self.phi_D, solver.Phi],
            exclude=["kappa", "d", "phi_D"],
            return_annotation=float,
            positional_or_keyword_order=brush_solver_default_key_order
        )
        def phi(**kwargs):
            # geometry and topology dependent calls
            kappa = ignore_extra_kwargs(self.kappa_callback)(**kwargs)
            phi_D = ignore_extra_kwargs(self.phi_D_callback)(**kwargs, kappa=kappa)
            D = ignore_extra_kwargs(self.D_callback)(**kwargs, kappa=kappa)

            phi_ = ignore_extra_kwargs(solver.Phi)(**kwargs, d=D, phi_D = phi_D, kappa=kappa)
            return phi_

        return phi


    def _Pi(self) -> Callable:
        @doc_decorator(
            callables=[self.D_callback, self.kappa_callback, self.phi_D_callback, solver.Phi, solver.Pi],
            exclude=["kappa", "d", "phi_D", "phi"],
            return_annotation=float,
            positional_or_keyword_order=brush_solver_default_key_order
        )
        def Pi(**kwargs):
            # geometry and topology dependent calls
            kappa = ignore_extra_kwargs(self.kappa_callback)(**kwargs)
            phi_D = ignore_extra_kwargs(self.phi_D_callback)(**kwargs, kappa=kappa)
            D = ignore_extra_kwargs(self.D_callback)(**kwargs, kappa=kappa)
            phi = ignore_extra_kwargs(solver.Phi)(**kwargs, d=D, phi_D = phi_D, kappa=kappa)
            Pi_ = solver.Pi(phi, kwargs["chi"])
            return Pi_

        return Pi
    

class Particle:
    def __init__(self, particle = "cylinder") -> None:
        self.__particle_module = importlib.import_module(
            f"ascf_pb.particle_geometry.{particle}")
        self.surface_integrand = functools.lru_cache()(self.__particle_module.surface_integrand)
        self.volume_integrand = functools.lru_cache()(self.__particle_module.volume_integrand)
        self.surface = functools.lru_cache()(self.__particle_module.surface)
        self.volume = functools.lru_cache()(self.__particle_module.volume)
        self.integration_interval = functools.lru_cache()(self.__particle_module.integration_interval)


class BrushInsertionFreeEnergy:
    def __init__(self, geometry: str = 'plain', vectorize=True, v_product=True, particle = "cylinder") -> None:
        self.Particle = Particle(particle)
        self.Brush = BrushSolver(geometry, vectorize=False)
        
        if vectorize:
            self.__deco = ascf_pb.decorators.vectorize(v_product)
        else:
            self.__deco = ascf_pb.decorators.dummy

        self.__free_energy_module = importlib.import_module("ascf_pb.free_energy")

        self.osmotic_free_energy_callback = functools.lru_cache()(self.__free_energy_module.osmotic_free_energy)
        self.surface_free_energy_callback = functools.lru_cache()(self.__free_energy_module.surface_free_energy)

        self.build_functions()
        
    
    def build_functions(self):
        self.phi_D = self.__deco(self.Brush._phi_D())
        self._D = self.__deco(self.Brush._D())
        self.phi = self.__deco(self.Brush._phi())
        self.Pi = self.__deco(self.Brush._Pi())

        self.osmotic_free_energy = self.__deco(self._osmotic_free_energy())
        #self.surface_free_energy = self.__deco(self._surface_free_energy())
    
    
    def _osmotic_free_energy(self):
        @doc_decorator(
            callables=[
                self.Particle.integration_interval,
                self.Particle.volume_integrand,
                self.Brush.Pi,
                self.osmotic_free_energy_callback
                ],
            exclude=["kappa", "d", "phi_D", "phi", "z0", "z1", "z"],
            return_annotation=float,
            positional_or_keyword_order=brush_solver_default_key_order+["ph", "pw", "pc"]
        )
        def osmotic_free_energy(**kwargs):
            z0, z1 = ignore_extra_kwargs(self.Particle.integration_interval)(**kwargs)
            volume_integrand = ignore_extra_kwargs(self.Particle.volume_integrand)(**kwargs)
            Pi_z = ascf_pb.decorators.make_partial(**kwargs, ignore_keyerror =True)(self.Brush.Pi)
            fe = self.osmotic_free_energy_callback(
                z0=z0, 
                z1=z1, 
                volume_integrand = volume_integrand, 
                Pi_z = Pi_z
                )
            return fe

        return osmotic_free_energy

    def _surface_free_energy(self):
        @doc_decorator(
            callables=[
                self.Particle.integration_interval,
                self.Particle.surface_integrand,
                self.Brush.phi,
                self.surface_free_energy_callback
                ],
            exclude=["kappa", "d", "phi_D", "phi", "z0", "z1", "z"],
            return_annotation=float,
            positional_or_keyword_order=brush_solver_default_key_order+["ph", "pw", "pc"]
        )
        def surface_free_energy(**kwargs):
            kwargs["z0"], kwargs["z1"] = ignore_extra_kwargs(self.Particle.integration_interval)(**kwargs)
            kwargs["surface_integrand"] = ignore_extra_kwargs(self.Particle.surface_integrand)(**kwargs)
            kwargs["phi_z"] = ascf_pb.decorators.make_partial(**kwargs, ignore_keyerror =True)(self.Brush.phi)
            
            fe = ignore_extra_kwargs(self.surface_free_energy_callback)(**kwargs)
            return fe

        return surface_free_energy

#%%
if __name__ == "__main__":
    b = BrushInsertionFreeEnergy(vectorize=False)
# %%
    b.osmotic_free_energy(chi = 0, N=1000, sigma = 0.02, ph =4, pw =4, pc = 12)
# %%
    b.surface_free_energy(chi = 0, N=1000, sigma = 0.02, ph =4, pw =4, pc = 12)
# %%
    b.surface_free_energy
# %%
