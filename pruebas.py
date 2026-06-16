from esquemas_clasicos import Simplificado, Shamir, ShamirRampa, McElieceSarwate
from esquemas_cuanticos import CGL, Ogawa, ZhangMatsumoto

from qiskit.quantum_info import Statevector
from numpy.linalg import norm
from galois import GF

secreto = Statevector([1,0,0,0,0,0,0,1])
secreto = secreto/norm(secreto)

cgl = CGL(GF(2**3), 4, ['a','b','c','d','e'])
part = cgl.codificacion(secreto)
sv = cgl.decodificacion(part)
print(sv.draw('latex_source'))