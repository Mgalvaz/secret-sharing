import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, transpile, ClassicalRegister
from qiskit.quantum_info import Statevector, partial_trace
import matplotlib.pyplot as plt
from qiskit.circuit.library import LinearFunction
from qiskit.visualization import plot_bloch_multivector
from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram
import matplotlib.pyplot as plt

sim = AerSimulator(method='statevector')
d = 4
l = 1
n = 3
r = 2
t = r-l
dim_participantes = (d-1).bit_length()

# Funciones utilidad
def run_statevector(qc):
    circ = qc.copy()
    circ.save_statevector()
    job = sim.run(circ)
    return job.result().get_statevector(circ)

def simulacion(qc):
    obj = sim.run(qc, shots=1024)
    return obj, qc

def suma(qc, origen, destino):
    """Suma bit a bit (XOR) en F4: destino = destino + origen"""
    qc.cx(origen[0], destino[0])
    qc.cx(origen[1], destino[1])

def mult_omega_y_suma(qc, origen, destino):
    """Multiplica por omega y suma: destino = destino + (omega * origen)"""
    # Regla: omega * (a1*omega + a0) = (a1+a0)*omega + a1
    # Por tanto: dest1 = dest1 + orig1 + orig0, dest0 = dest0 + orig1
    qc.cx(origen[1], destino[1])
    qc.cx(origen[0], destino[1])
    qc.cx(origen[1], destino[0])

def mult_omega2_y_suma(qc, origen, destino):
    """Multiplica por omega^2 y suma: destino = destino + (omega^2 * origen)"""
    # Regla: omega^2 * (a1*omega + a0) = a0*omega + (a1+a0)
    # Por tanto: dest1 = dest1 + orig0, dest0 = dest0 + orig1 + orig0
    qc.cx(origen[0], destino[1])
    qc.cx(origen[1], destino[0])
    qc.cx(origen[0], destino[0])


if __name__ == "__main__":
    participantes = [QuantumRegister(dim_participantes, f'p{i}') for i in range(n)]
    resultado = [ClassicalRegister(dim_participantes, f's{i}') for i in range(l)]
    qc = QuantumCircuit(*participantes, *resultado)
    ## PASO 1: inicializar los coeficientes de todos los posibles polinomios
    # Inicializar el secreto en los primeros l registros
    svect = [0, 1, 0, 1]
    print('secreto: ',f'{'00 + ' * svect[0]}{'01 + ' * svect[1]}{'10 + ' * svect[2]}{'11' * svect[3]}')
    # svect = np.random.random(4)
    svect = svect / np.linalg.norm(svect)
    svsec = Statevector(svect)
    qc.initialize(svsec, participantes[0])
    # Superponer todos los posibles valores de los coeficientes del polinomio
    for participante in participantes[l:r]:
        qc.h(participante)

    # PASO 2: Aplicar la isometría
    # --- Participante 3: x = omega^2 ---
    # y3 = c1*(1) + c2*(omega^2)
    suma(qc, participantes[0], participantes[2])
    mult_omega2_y_suma(qc, participantes[1], participantes[2])
    # Como los participantes 1 y 2 tienen en sus qudits información, no se puede hacer lo mismo que en el caso anterior.
    # --- Participante 1: x = 1 y Participante 2: x = omega ---
    # Hacemos el siguiente mapeo lineal: y1 = c1*(1) + c2*(1), y2 = c1*(1) + c2*(omega)
    matriz = [
        [True, False, True, False],  # y1_1 = c1_1 + c2_1
        [False, True, False, True],  # y1_0 = c1_0 + c2_0
        [True, False, True, True],  # y2_1 = c1_1 + c2_0 + c2_1 (Recordar regla de c2*omega)
        [False, True, True, False]  # y2_0 = c1_0 + c2_1 (Recordar regla de c2*omega)
    ]
    puerta_v1_v2 = LinearFunction(matriz)
    qc.append(puerta_v1_v2, [part[i] for part in participantes[:r] for i in range(r - 1, -1, -1)])
    qc = qc.decompose('Linear_function')
    sv = run_statevector(qc)
    print('Estado tras la codificación:')
    sv.draw('qsphere')
    plt.show(block=False)

    print('traza parcial del primer (último) participante del estado codificado: ', partial_trace(sv, range(2, 6)))

    ### DECODIFICACION DEL SECRETO A PARTIR DE LOS PARTICIPANTES P2 Y P3
    # PASO 1: Obtener c2 en el registro p2
    # Como c2 = y2 + y3, simplemente sumamos p3 en p2
    suma(qc, participantes[2], participantes[1])
    # Ahora p2 almacena c2

    # PASO 2: Limpiar p3 para que solo quede c1
    # De la ec: y3 = c1 + omega^2 * c2  =>  c1 = y3 + omega^2 * c2
    # Como ahora tenemos c2 en p2, multiplicamos p2 por omega^2 y lo sumamos a p3
    mult_omega2_y_suma(qc, participantes[1], participantes[2])

    # PASO 3: separar el secreto del resto
    suma(qc, participantes[2], participantes[1])
    sv = run_statevector(qc)
    print('Estado tras paso 2:')
    print('traza parcial del tercer (primero) participante del estado codificado: ', partial_trace(sv, range(4)))
    sv.draw('bloch')
    plt.show()

