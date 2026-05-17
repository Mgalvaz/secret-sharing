import numpy as np
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.quantum_info import partial_trace
from qiskit.circuit.library import LinearFunction
from qiskit_aer import AerSimulator

from utils import extender_matriz

sim = AerSimulator(method='statevector')
# Funciones utilidad
def run_statevector(quantum_circuit):
    circ_copy = quantum_circuit.copy()
    circ_copy.save_statevector()
    job = sim.run(circ_copy)
    return job.result().get_statevector(circ_copy)


class CGL:
    r"""
    Esquema cuántico de compartición de secretos de Cleve-Gottesman-Lo sobre el espacio de Hilbert complejo $\mathcal{H}_{2^m}$.

    Ejemplo:
        Crea un esquema de Cleve-Gottesman-Lo con parámetro r = 4 sobre $\mathcal{H}_{2^5}$ para los participantes ['a', 'b', 'c', 'd', 'e', 'f'] con .

        .. ipython:: python

            cuerpo = galois.GF(2**5)
            qss = CGL(cuerpo, 4, ['a', 'b', 'c', 'd', 'e', 'f'])
    """
    def __init__(self, cuerpo, r, participantes):
        r"""
        Crea un esquema de compartición de secretos de Cleve-Gottesman-Lo sobre el el espacio de Hilbert complejo $\mathcal{H}_{2^m}$.
        :param cuerpo: El cuerpo finito que actúa como base del espacio de Hilbert sobre el que el esquema esta construido.
        :param r: Parámetro de reconstrucción del esquema (número mínimo de participantes necesarios para reconstruir el secreto).
        :param participantes: Lista de los identificadores únicos de cada participante del esquema.
        """
        # Verificación de condiciones
        if cuerpo.characteristic != 2:
            raise ValueError(f'El cuerpo introducido debe tener como elemento base el 2.')
        if 2*r-1 < len(participantes):
            raise ValueError(f'El numero de participantes ({len(participantes)}) debe ser menor o igual a 2*r-1 ({2*r-1}).')
        if cuerpo.order <= 2*r-1:
            raise ValueError(f'El numero de participantes totales ({2*r-1}) debe ser menor que el orden del cuerpo de trabajo ({cuerpo.order}).')
        if len(participantes) != len(set(participantes)):
            raise ValueError(f'Se han encontrado participantes duplicados.')
        if len(participantes) < r:
            raise ValueError(f'El parámetro de reconstrucción ({r}) debe ser menor o igual que el número de participantes ({len(participantes)}).')
        if r < 2:
            raise ValueError(f'El parámetro de reconstrucción ({r}) debe ser mayor que 1.')

        self.cuerpo = cuerpo
        self.reconstruccion = r
        self.__participaciones_anticipadas = []
        self.elementos_no_anticipados = None
        self.__resto_anticipadas = []
        participaciones_reales = [QuantumRegister(cuerpo.degree, participante) for participante in participantes]
        farticipaciones_extra = [QuantumRegister(cuerpo.degree, f'p{i}') for i in range(len(participantes) + 1, 2 * r)]
        self.__participaciones = participaciones_reales + farticipaciones_extra
        self.__circuito = QuantumCircuit(*self.__participaciones)
        self.participantes_numero = {nombre: i for i, nombre in enumerate(participantes, 1)}

    def comparticion_anticipada(self, participantes_anticipados):
        """
        Crea participaciones anticipadas para cada participante especificado.
        El formato de las participaciones es: (nombre, participación), donde participación es un qudit del circuito cuántico.
        :param participantes_anticipados: Listado de los participantes a entregar participaciones anticipadas.
        :return: Una lista que contienene las participaciones anticipadas asignadas a cada participante especificado.
        """
        # Verificación de condiciones
        qc = self.__circuito
        r = self.reconstruccion
        if self.__participaciones_anticipadas is None:
            raise AttributeError(f'Ya se han repartido todas las participaciones.')
        if self.elementos_no_anticipados is not None:
            raise AttributeError(f'Ya se han repartido todas las participaciones anticipadas.')
        if len(participantes_anticipados) != len(set(participantes_anticipados)):
            raise ValueError(f'Se han encontrado participantes duplicados.')
        for nombre in participantes_anticipados:
            if nombre not in self.participantes_numero:
                raise ValueError(f'El participante {nombre} no está registrado.')
        if self.reconstruccion <= len(participantes_anticipados):
            raise ValueError(f'El numero de participaciones anticipadas ({len(participantes_anticipados)}) debe ser menor o igual que el parámetro de privacidad ({self.reconstruccion - 1})')

        elem_extra = [self.participantes_numero[nombre] for nombre in participantes_anticipados]
        extras = [self.__participaciones[idx-1] for idx in elem_extra]
        elem_resto = self.cuerpo(np.setdiff1d(range(1, 2 * r), elem_extra))
        self.elementos_no_anticipados = elem_resto[:len(self.participantes_numero)-len(participantes_anticipados)].tolist()
        resto = [self.__participaciones[idx - 1] for idx in elem_resto.tolist()]
        self.__participaciones_anticipadas = extras + resto[:r-len(extras)-1] # Se amplian las participaciones anticipadas hasta ser r-1
        self.__resto_anticipadas = resto[r-len(extras)-1:]
        # Crear el estado completamente mezcado (sum |y>|y>)
        for qudit_e, qudit_r in zip(self.__participaciones_anticipadas, self.__resto_anticipadas):
            qc.h(qudit_e)
            qc.cx(qudit_e, qudit_r)
        return extras

    def codificacion(self, secreto):
        """
        Crea las participaciones de todos los participantes de acuerdo al secreto recibido.
        El formato de las participaciones es: (nombre, participación), donde participación es un qudit del circuito cuántico.
        Si se han distribuido participaciones anticipadas, las participaciones serán coherentes con las mismas.
        :param secreto: Secreto que se quiere codificar entre todos los participantes.
        :return: Una lista que contienene las participaciones de cada participante que no ha participado en la distribución avanzada.
        """
        # Verificación de condiciones
        if not secreto.is_valid():
            raise ValueError(f'No se ha introducido un estado cuántico válido.')
        if len(secreto) != self.cuerpo.order:
            raise ValueError(f'Se esperaba un vector de estado de dimensión {self.cuerpo.order}, pero se ha recibido uno de dimensión {len(secreto)}.')
        if self.__participaciones_anticipadas is None:
            raise AttributeError(f'Ya se han repartido todas las participaciones.')

        qc = self.__circuito
        r = self.reconstruccion
        m = 2*r-1 # Numero de participantes totales
        # Procedimiento estandar
        if len(self.__participaciones_anticipadas) == 0:
            x = list(range(1, m+1))
            qc.initialize(secreto, [self.__participaciones[0]])  # Inicializar el secreto
            for participacion in self.__participaciones[1:r]:  # Superponer todos los posibles valores de los coeficientes del polinomio
                qc.h(participacion)
            # Evaluar en cada registro los polinomios en los elementos de los participantes
            vandermonde = self.cuerpo(x)[:, None] ** np.arange(r)
            matriz = self.cuerpo(np.column_stack([vandermonde, np.zeros((m, r - 1))])) # Matriz de evaluación
            matriz = extender_matriz(matriz) # Extender la matriz de numeros de F_q a vectores de F_2
            orden_participantes = [qubit for participacion in self.__participaciones for qubit in reversed(participacion)]  # Como qiskit es Little Endian, pero la matriz extendida está en Big Endian, hay que invertir el orden de los qubits de los participantes
            x = x[:len(self.participantes_numero)]

        # Compartición anticipada
        else:
            qc.initialize(secreto, self.__resto_anticipadas[-1])
            x = self.elementos_no_anticipados
            elem_extras = self.cuerpo([self.participantes_numero[participacion.name] for participacion in self.__participaciones_anticipadas] + [0])
            elem_resto = self.cuerpo(np.setdiff1d(range(1, 2 * r), elem_extras))
            matriz_p1 = np.linalg.inv(elem_extras[:, None] ** np.arange(r))  # Matriz del primer paso del procediemiento de decodificacion
            matriz_p2 = elem_resto[:, None] ** np.arange(r)  # Matriz del segundo paso del procediemiento de decodificacion
            matriz = extender_matriz(matriz_p2 @ matriz_p1)
            orden_participantes = [qubit for participacion in self.__resto_anticipadas for qubit in reversed(participacion)]

        qc.append(LinearFunction(matriz), orden_participantes)  # Aplicar matriz de evaluación
        self.__circuito = qc.decompose('Linear_function')
        self.__participaciones_anticipadas =  None
        self.__resto_anticipadas = None
        # Generar el resto de las participaciones
        return list(self.__participaciones[i-1] for i in x)

    def decodificacion(self, participaciones):
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (nombre, participación), donde participación es un qudit del circuito cuántico.
        :param participaciones: Secuencia con las participaciones de los participantes que desean obtener el secreto.
        :return: El secreto.
        """
        # Verificación de condiciones
        r = self.reconstruccion
        qc = self.__circuito
        if len(participaciones) < r:
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto')

        # Obetener los elementos asociados a cada participante
        elementos = self.cuerpo([self.participantes_numero[participacion.name] for participacion in participaciones[:r]])
        elementos_resto = self.cuerpo(np.setdiff1d(np.arange(2*r), elementos))
        orden_participantes = [qubit for participacion in participaciones[:r] for qubit in reversed(participacion)]
        matriz_p1 = np.linalg.inv(elementos[:,None]**np.arange(r)) # Matriz del primer paso del procediemiento de decodificacion
        matriz_p2 = elementos_resto[:,None]**np.arange(r) # Matriz del segundo paso del procediemiento de decodificacion
        matriz = extender_matriz(matriz_p2@matriz_p1)
        qc.append(LinearFunction(matriz), orden_participantes) # Realizar los dos pasos en uno
        sv = run_statevector(qc.decompose('Linear_function')) # Silmualar el circuito
        elementos_traza = list(range((int(elementos[0])-1)*self.cuerpo.degree)) + list(range(int(elementos[0])*self.cuerpo.degree, (2*self.reconstruccion-1)*self.cuerpo.degree))  # Posicion de los qubits a trazar
        return partial_trace(sv, elementos_traza).to_statevector()