#%%
from ascf_pb.topology import pore
from ascf_pb.topology import kappa
chi = 0.0
sigma = 0.02
N=1000
kappa_cb = kappa.kappa_plain
pore.opening_pore_Radius(chi,kappa_cb(N),N,sigma)