from esquemas_clasicos import Simplificado, Shamir, ShamirRampa, McElieceSarwate
from esquemas_cuanticos import CGL, Ogawa, ZhangMatsumoto

from qiskit.quantum_info import Statevector
from galois import GF, Poly, lagrange_poly

from utils import bytes_a_int, int_a_bytes

import numpy as np
import random
import itertools


nombres_esq = {
    'Simplificado': 'Esquema Simplificado',
    'Shamir': 'Esquema de Shamir',
    'ShamirRampa': 'Esquema de Shamir en rampa',
    'McElieceSarwate': 'Esquema de McEliece-Sarwate',
    'CGL': 'Esquema de Cleve-Gottesman-Lo',
    'Ogawa': 'Esquema de Ogawa et al.',
    'ZhangMatsumoto': 'Esquema de Zhang-Matsumoto',
}
clasico_umbral = [Shamir]
clasico_rampa = [ShamirRampa, McElieceSarwate]
cuantico_umbral = [CGL]
cuantico_rampa = [Ogawa, ZhangMatsumoto]


def construir_polinomios(gf, s, r):  # Genera el conjunto D_r(s)
    s = gf(s)
    for extra in itertools.product(range(gf.order), repeat=r - len(s)):
        coef = np.concatenate([s, extra])
        yield Poly(coef, field=gf, order='asc')

def construir_polinomios_E(gf, A, s, r):  # Genera el conjunto E_r(A,s)
    A = gf(A)
    s = gf(s)
    lagrange = lagrange_poly(A, s)
    zeros = Poly.Roots(A, field=gf)
    for extra in itertools.product(range(gf.order), repeat=r - len(s)):
        yield lagrange + zeros * Poly(extra, field=gf)

def funcionamiento_clasico(Esquema, cuerpo, r=None, l=None, participantes=None):
    print('=' * 80)
    print('[Test]', nombres_esq[Esquema.__name__])
    print('=' * 80)
    if Esquema is Simplificado:
        esquema = Esquema(cuerpo, participantes)
        secreto = random.randint(0, cuerpo.order - 1)
        r = len(participantes)
    elif Esquema in clasico_umbral:
        esquema = Esquema(cuerpo, r, participantes)
        secreto = random.randint(0, cuerpo.order-1)
    else:
        esquema = Esquema(cuerpo, r, l, participantes)
        secreto = [random.randint(0, cuerpo.order-1) for _ in range(l)]
    print('Secreto =', secreto)
    secreto = int_a_bytes(secreto)
    participaciones = esquema.codificacion(secreto)
    part = random.sample(participaciones, r)
    sec = esquema.decodificacion(part)
    print('Secreto decodificado =', bytes_a_int(sec))
    print()

def comparticion_anticipada_clasico(Esquema, cuerpo, r=None, l=None, participantes=None):
    print('=' * 80)
    print(f'[Test] {nombres_esq[Esquema.__name__]}: compartición anticipada')
    print('=' * 80)
    if Esquema is Simplificado:
        esquema = Esquema(cuerpo, participantes)
        secreto = random.randint(0, cuerpo.order - 1)
        r = len(participantes)
        anticipados = random.sample(participantes, r - 1)
    elif Esquema in clasico_umbral:
        esquema = Esquema(cuerpo, r, participantes)
        secreto = random.randint(0, cuerpo.order-1)
        anticipados = random.sample(participantes, r - 1)
    else:
        esquema = Esquema(cuerpo, r, l, participantes)
        secreto = [random.randint(0, cuerpo.order-1) for _ in range(l)]
        anticipados = random.sample(participantes, r - l)
    participaciones_anticipadas = esquema.comparticion_anticipada(anticipados)
    print('Participaciones anticipadas:', participaciones_anticipadas)
    print('Secreto =', secreto)
    secreto = int_a_bytes(secreto)
    participaciones = esquema.codificacion(secreto)
    print('Resto de participaciones:', participaciones)
    part = participaciones + participaciones_anticipadas
    random.shuffle(part)
    sec = esquema.decodificacion(part)
    print('Secreto decodificado =', bytes_a_int(sec))

def funcionamiento_cuantico(Esquema, cuerpo, r=None, l=None, participantes=None):
    print('=' * 80)
    print('[Test]', nombres_esq[Esquema.__name__])
    print('=' * 80)
    if Esquema in cuantico_umbral:
        dims = (cuerpo.order,)
        esquema = Esquema(cuerpo, r, participantes)
        n1, n2 = random.sample(range(cuerpo.order), 2)
        secreto = (1+1j)*Statevector.from_int(n1,dims) + 5*Statevector.from_int(n2,dims) # (1+1j)|n1> + 5|n2>
        secreto = secreto/np.linalg.norm(secreto) # normalizar el secreto
        m = np.arange(2*r-1,0,-1)
        funcion = construir_polinomios
        secreto_vector = [[n1], [n2]]
    else:
        dims = tuple(cuerpo.order for _ in range(l))
        esquema = Esquema(cuerpo, r, l, participantes)
        n1, n2 = np.random.randint(cuerpo.order, size=(2,l))
        base = cuerpo.order**np.arange(l)
        secreto = (1+3j)*Statevector.from_int((n1*base).sum(),dims) + 2j*Statevector.from_int((n2*base).sum(),dims)  # (1+3j)|n1l,..., n11> + 2j|n2l,..., n21>
        secreto = secreto / np.linalg.norm(secreto)  # normalizar el secreto
        secreto_vector = np.array([n1,n2])
        if Esquema is Ogawa:
            m = np.arange(2 * r - l, 0, -1)
            funcion = construir_polinomios
        else:
            m = np.arange(2 * r - 1, l-1, -1)
            funcion = lambda gf, s, r: construir_polinomios_E(gf, np.arange(l),s, r)
    print('Secreto =', secreto.to_dict())
    participaciones = esquema.codificacion(secreto)
    estado_codificado = Statevector(getattr(esquema,f'_{Esquema.__name__}__circuito'))
    estados_codificado = Statevector(estado_codificado.data, dims=tuple(cuerpo.order for _ in range(len(m)))).to_dict().keys()
    print('Estados base codificado:', estados_codificado)
    estados_matematicos = [''.join(f(m).view(np.ndarray).astype(str)) for s in secreto_vector for f in funcion(cuerpo, s, r)]
    print('Estados base matematica;', estados_matematicos)
    print('Estados matemáticos == estados codificados:', set(estados_matematicos) == estados_codificado)
    print('Participaciones:', participaciones)
    part = random.sample(participaciones, r)
    sec = esquema.decodificacion(part)
    sec = Statevector(sec.data, dims)
    print('Secreto decodificado =', sec.to_dict())
    print('secreto == sec:', secreto.equiv(sec))
    print()

def comparticion_anticipada_cuantico(Esquema, cuerpo, r=None, l=None, participantes=None):
    print('=' * 80)
    print('[Test]', nombres_esq[Esquema.__name__])
    print('=' * 80)
    if Esquema in cuantico_umbral:
        dims = (cuerpo.order,)
        esquema = Esquema(cuerpo, r, participantes)
        n1, n2 = random.sample(range(cuerpo.order), 2)
        secreto = (1 + 1j) * Statevector.from_int(n1, dims) + 5 * Statevector.from_int(n2, dims)  # (1+1j)|n1> + 5|n2>
        secreto = secreto / np.linalg.norm(secreto)  # normalizar el secreto
        m = np.arange(2 * r - 1, 0, -1)
        funcion = construir_polinomios
        secreto_vector = [[n1], [n2]]
        anticipados = random.sample(participantes, r - 1)
    else:
        dims = tuple(cuerpo.order for _ in range(l))
        esquema = Esquema(cuerpo, r, l, participantes)
        n1, n2 = np.random.randint(cuerpo.order, size=(2, l))
        base = cuerpo.order ** np.arange(l)
        secreto = (1 + 3j) * Statevector.from_int((n1 * base).sum(), dims) + 2j * Statevector.from_int((n2 * base).sum(), dims)  # (1+3j)|n1l,..., n11> + 2j|n2l,..., n21>
        secreto = secreto / np.linalg.norm(secreto)  # normalizar el secreto
        secreto_vector = np.array([n1, n2])
        anticipados = random.sample(participantes, r - l)
        if Esquema is Ogawa:
            m = np.arange(2 * r - l, 0, -1)
            funcion = construir_polinomios
        else:
            m = np.arange(2 * r - 1, l - 1, -1)
            funcion = lambda gf, s, r: construir_polinomios_E(gf, np.arange(l), s, r)
    participaciones_anticipadas = esquema.comparticion_anticipada(anticipados)
    print('Participaciones anticipadas:', participaciones_anticipadas)
    print('Secreto =', secreto.to_dict())
    participaciones = esquema.codificacion(secreto)
    estado_codificado = Statevector(getattr(esquema,f'_{Esquema.__name__}__circuito'))
    estados_codificado = Statevector(estado_codificado.data, dims=tuple(cuerpo.order for _ in range(len(m)))).to_dict().keys()
    print('Estados base codificado:', estados_codificado)
    estados_matematicos = [''.join(f(m).view(np.ndarray).astype(str)) for s in secreto_vector for f in funcion(cuerpo, s, r)]
    print('Estados base matematica;', estados_matematicos)
    print('Estados matemáticos == estados codificados:', set(estados_matematicos) == estados_codificado)
    print('Participaciones:', participaciones)
    part = random.sample(participaciones + participaciones_anticipadas, r)
    sec = esquema.decodificacion(part)
    sec = Statevector(sec.data, dims)
    print('Secreto decodificado =', sec.to_dict())
    print('secreto == sec:', secreto.equiv(sec))
    print()


print('#' * 80)
print(' ' * 20 + 'FUNCIONAMIENTO DE LOS ESQUEMAS CLÁSICOS')
print('#' * 80)
for Esquema in [Simplificado] + clasico_umbral + clasico_rampa:
    funcionamiento_clasico(Esquema, GF(2, 8*10), r=5, l=3, participantes=['a', 'b', 'c', 'd', 'e', 'f'])

print('#' * 80)
print(' ' * 15 + 'COMPARTICIÓN ANTICIPADA DE LOS ESQUEMAS CLÁSICOS')
print('#' * 80)
for Esquema in [Simplificado] + clasico_umbral + clasico_rampa:
    comparticion_anticipada_clasico(Esquema, GF(2, 8 * 6), r=5, l=2, participantes=['a', 'b', 'c', 'd', 'e', 'f'])

print('#' * 80)
print(' ' * 18 + 'FUNCIONAMIENTO DE LOS ESQUEMAS CUÁNTICOS')
print('#' * 80)
for Esquema in cuantico_umbral + cuantico_rampa:
    funcionamiento_cuantico(Esquema, GF(2,3), r=3, l=2, participantes=['a', 'b', 'c', 'd'])

print('#' * 80)
print(' ' * 15 + 'COMPARTICIÓN ANTICIPADA DE LOS ESQUEMAS CUÁNTICOS')
print('#' * 80)
for Esquema in cuantico_umbral + cuantico_rampa:
    comparticion_anticipada_cuantico(Esquema, GF(2, 3), r=3, l=2, participantes=['a', 'b', 'c', 'd'])