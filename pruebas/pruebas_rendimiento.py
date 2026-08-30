from classic_schemes import Simplificado, Shamir, ShamirRampa, McElieceSarwate
from quantum_schemes import CGL, Ogawa, ZhangMatsumoto

from qiskit.quantum_info import Statevector
from qiskit_aer import StatevectorSimulator
from galois import GF

import numpy as np
import time
import secrets
import random

N_BYTES = 51
N_QUBITS_1 = 6
N_QUBITS_2 = 4

num_cores = 9
sim = StatevectorSimulator(max_parallel_threads=num_cores, max_parallel_experiments=num_cores)

def tiempos_clasico(esquema, r, l=1, anticipada=None):
    # Compartición anticipada
    if anticipada:
        inicio_anticipada = time.time()
        part_anticipadas = esquema.advance_sharing(anticipada)
        fin_anticipada = time.time()
        print('Duración de la compartición anticipada:', fin_anticipada - inicio_anticipada)
    if l==1: # Elección aleatoria del secreto
        secreto = secrets.token_bytes(N_BYTES)
    else:
        secreto = [secrets.token_bytes(N_BYTES) for _ in range(l)]
    # Codificación
    inicio_codificacion = time.time()
    part_resto = esquema.codificacion(secreto)
    fin_codificacion = time.time()
    print('Duración de la codificación:', fin_codificacion - inicio_codificacion)
    if anticipada: # Eleccion aleatoria de los participantes que reconstruirán el secreto
        part_reconstruccion = random.sample(part_anticipadas + part_resto, r)
    else:
        part_reconstruccion = random.sample(part_resto, r)
    # Decodificación
    inicio_decodificacion = time.time()
    secreto_dec = esquema.decodificacion(part_reconstruccion)
    fin_decodificacion = time.time()
    print('Duración de la decodificación:', fin_decodificacion - inicio_decodificacion)

def tiempos_cuantico(esquema, r, l=1, anticipada=None, tipo=1):
    # Compartición anticipada
    tiempo_anticipada = 0
    if anticipada:
        inicio_anticipada = time.time()
        part_anticipadas = esquema.advance_sharing(anticipada)
        sim.run(getattr(esquema, f'_{type(esquema).__name__}__circuito'))
        fin_anticipada = time.time()
        tiempo_anticipada = fin_anticipada - inicio_anticipada
        print('Duración de la compartición anticipada:', tiempo_anticipada)
    if tipo == 1: # Elección aleatoria del secreto
        secreto = Statevector(np.random.normal(size=2**(N_QUBITS_1*l)) + 1j*np.random.normal(size=2**(N_QUBITS_1*l)))
    else:
        secreto = Statevector(np.random.normal(size=2**(N_QUBITS_2*l)) + 1j*np.random.normal(size=2**(N_QUBITS_2*l)))
    secreto = secreto/np.linalg.norm(secreto)
    # Codificación
    inicio_codificacion = time.time()
    part_resto = esquema.codificacion(secreto)
    sim.run(getattr(esquema, f'_{type(esquema).__name__}__circuito'))
    fin_codificacion = time.time()
    tiempo_codificacion = fin_codificacion - inicio_codificacion - tiempo_anticipada
    print('Duración de la codificación:', tiempo_codificacion)
    if anticipada: # Eleccion aleatoria de los participantes que reconstruirán el secreto
        part_reconstruccion = random.sample(part_anticipadas + part_resto, r)
    else:
        part_reconstruccion = random.sample(part_resto, r)
    # Decodificación
    inicio_decodificacion = time.time()
    secreto_dec = esquema.decodificacion(part_reconstruccion)
    fin_decodificacion = time.time()
    print('Duración de la decodificación:', fin_decodificacion - inicio_decodificacion - tiempo_codificacion - tiempo_anticipada)


def rendimiento_clasico():
    cuerpo = GF(2, N_BYTES * 8 + 1)
    participantes = list('abcdefghijklmnopqrstuvwxyz1234567890')
    print('=' * 80)
    print('[Test] Esquema Simplificado')
    print('=' * 80)
    simplificado = Simplificado(cuerpo, participantes)
    tiempos_clasico(simplificado, len(participantes))
    print('=' * 80)
    print('[Test] Esquema Simplificado: compartición anticipada')
    print('=' * 80)
    simplificado = Simplificado(cuerpo, participantes)
    anticipadas = random.sample(participantes, len(participantes)-1)
    tiempos_clasico(simplificado, len(participantes), anticipada=anticipadas)
    print('=' * 80)
    print('[Test] Esquema de Shamir')
    print('=' * 80)
    shamir = Shamir(cuerpo, 20, participantes)
    tiempos_clasico(shamir, 20)
    print('=' * 80)
    print('[Test] Esquema de Shamir: compartición anticipada')
    print('=' * 80)
    shamir = Shamir(cuerpo, 20, participantes)
    anticipadas = random.sample(participantes, 19)
    tiempos_clasico(shamir, 20, anticipada=anticipadas)
    print('=' * 80)
    print('[Test] Esquema de Shamir en rampa')
    print('=' * 80)
    shamir_rampa = ShamirRampa(cuerpo, 20, 8, participantes)
    tiempos_clasico(shamir_rampa, 20, l=8)
    print('=' * 80)
    print('[Test] Esquema de Shamir en rampa: compartición anticipada')
    print('=' * 80)
    shamir_rampa = ShamirRampa(cuerpo, 20, 8, participantes)
    anticipadas = random.sample(participantes, 12)
    tiempos_clasico(shamir_rampa, 20, l=8, anticipada=anticipadas)
    print('=' * 80)
    print('[Test] Esquema de McEliece-Sarwate')
    print('=' * 80)
    mceliece_sarwate = McElieceSarwate(cuerpo, 20, 8, participantes)
    tiempos_clasico(mceliece_sarwate, 20, l=8)
    print('=' * 80)
    print('[Test] Esquema de McEliece-Sarwate: compartición anticipada')
    print('=' * 80)
    mceliece_sarwate = McElieceSarwate(cuerpo, 20, 8, participantes)
    anticipadas = random.sample(participantes, 12)
    tiempos_clasico(mceliece_sarwate, 20, l=8, anticipada=anticipadas)

def rendimiento_cuantico_1():
    cuerpo = GF(2, N_QUBITS_1)
    participantes = list('1234')
    print('=' * 80)
    print('[Test] Esquema de Cleve-Gottesman-Lo')
    print('=' * 80)
    cgl = CGL(cuerpo, 3, participantes)
    tiempos_cuantico(cgl, 3)
    print('=' * 80)
    print('[Test] Esquema de Cleve-Gottesman-Lo: compartición anticipada')
    print('=' * 80)
    cgl = CGL(cuerpo, 3, participantes)
    anticipadas = random.sample(participantes, 2)
    tiempos_cuantico(cgl, 3, anticipada=anticipadas)
    print('=' * 80)
    print('[Test] Esquema de Ogawa et al.')
    print('=' * 80)
    ogawa = Ogawa(cuerpo, 3, 2, participantes)
    tiempos_cuantico(ogawa, 3, l=2)
    print('=' * 80)
    print('[Test] Esquema de Ogawa et al: compartición anticipada')
    print('=' * 80)
    ogawa = Ogawa(cuerpo, 3, 2, participantes)
    anticipadas = random.sample(participantes, 1)
    tiempos_cuantico(ogawa, 3, l=2, anticipada=anticipadas)
    print('=' * 80)
    print('[Test] Esquema de Zhang-Matsumoto')
    print('=' * 80)
    zhang_matsumoto = ZhangMatsumoto(cuerpo, 3, 2, participantes)
    tiempos_cuantico(zhang_matsumoto, 3, l=2)
    print('=' * 80)
    print('[Test] Esquema de Zhang-Matsumoto: compartición anticipada')
    print('=' * 80)
    zhang_matsumoto = ZhangMatsumoto(cuerpo, 3, 2, participantes)
    anticipadas = random.sample(participantes, 1)
    tiempos_cuantico(zhang_matsumoto, 3, l=2, anticipada=anticipadas)

def rendimiento_cuantico_2():
    cuerpo = GF(2, N_QUBITS_2)
    participantes = list('1234567')
    print('=' * 80)
    print('[Test] Esquema de Cleve-Gottesman-Lo')
    print('=' * 80)
    cgl = CGL(cuerpo, 4, participantes)
    tiempos_cuantico(cgl, 4, tipo=2)
    print('=' * 80)
    print('[Test] Esquema de Cleve-GottesmaLo: compartición anticipada')
    print('=' * 80)
    cgl = CGL(cuerpo, 4, participantes)
    anticipadas = random.sample(participantes, 3)
    tiempos_cuantico(cgl, 4, anticipada=anticipadas, tipo=2)
    print('=' * 80)
    print('[Test] Esquema de Ogawa et al.')
    print('=' * 80)
    ogawa = Ogawa(cuerpo, 5, 2, participantes)
    tiempos_cuantico(ogawa, 5, l=2, tipo=2)
    print('=' * 80)
    print('[Test] Esquema de Ogawa et al: compartición anticipada')
    print('=' * 80)
    ogawa = Ogawa(cuerpo, 5, 2, participantes)
    anticipadas = random.sample(participantes, 3)
    tiempos_cuantico(ogawa, 5, l=2, anticipada=anticipadas, tipo=3)
    print('=' * 80)
    print('[Test] Esquema de Zhang-Matsumoto')
    print('=' * 80)
    zhang_matsumoto = ZhangMatsumoto(cuerpo, 5, 2, participantes)
    tiempos_cuantico(zhang_matsumoto, 5, l=2, tipo=3)
    print('=' * 80)
    print('[Test] Esquema de Zhang-Matsumoto: compartición anticipada')
    print('=' * 80)
    zhang_matsumoto = ZhangMatsumoto(cuerpo, 5, 2, participantes)
    anticipadas = random.sample(participantes, 3)
    tiempos_cuantico(zhang_matsumoto, 5, l=2, anticipada=anticipadas, tipo=3)

if __name__ == '__main__':
    rendimiento_clasico()