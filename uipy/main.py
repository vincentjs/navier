from ctypes import *
from decimal import *

cfdlib = cdll.LoadLibrary("./cfd.so")

nx = 41
ny = 41
nt = 1000
nit = 100
c = 1

dx = c_double(2.0 / (nx-1))
dy = c_double(2.0 / (ny-1))

rho = 1
nu = c_double(0.1)
dt = c_double(0.001)

sz = nx*ny
py_u = [0] * sz
py_v = [0] * sz
py_p = [0] * sz

u = (c_int * len(py_u))(*py_u)
v = (c_int * len(py_v))(*py_v)
p = (c_int * len(py_p))(*py_p)

cavityflow = cfdlib.CavityFlow_init()

cfdlib.CavityFlow_cavityFlow(cavityflow, nt, nit, u, sz, v, sz, dt, dx, dy, p, sz, rho, nu, nx, ny)
