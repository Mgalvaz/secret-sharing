import numpy as np
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.quantum_info import partial_trace
from qiskit.circuit.library import LinearFunction
from qiskit_aer import StatevectorSimulator

from utils import extender_matriz


sim = StatevectorSimulator()
class CGL:
    r"""
    Esquema cuántico de compartición de secretos de Cleve-Gottesman-Lo sobre el espacio de Hilbert complejo $\mathcal{H}_{2^m}$.

    Ejemplo:
        Crea un esquema de Cleve-Gottesman-Lo con parámetro r = 4 sobre $\mathcal{H}_{2^5}$ para los participantes ['a', 'b', 'c', 'd', 'e', 'f'].

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
        if r < 2:
            raise ValueError(f'El parámetro de reconstrucción ({r}) debe ser mayor que 1.')
        if len(participantes) < r:
            raise ValueError(f'El parámetro de reconstrucción ({r}) debe ser menor o igual que el número de participantes ({len(participantes)}).')

        self.cuerpo = cuerpo
        self.reconstruccion = r
        self.__participaciones_anticipadas = []
        self.elem_anticipadas_no_repartidas = []
        self.__participacion_secreto_anticipado = None
        participaciones_reales = [QuantumRegister(cuerpo.degree, participante) for participante in participantes]
        participaciones_ficticias = [QuantumRegister(cuerpo.degree, f'p{i}') for i in range(len(participantes) + 1, 2 * r)]
        self.__participaciones = participaciones_reales + participaciones_ficticias
        self.__circuito = QuantumCircuit(*self.__participaciones)
        self.participantes_numero = {nombre: i for i, nombre in enumerate(participantes, 1)}

    def comparticion_anticipada(self, participantes_anticipados):
        """
        Crea participaciones anticipadas para cada participante especificado.
        Cada participación es un registro cuántico (clase QuantumRegister) con el identificador único del participante al que le corresponde.
        :param participantes_anticipados: Listado de los participantes a entregar participaciones anticipadas.
        :return: Una lista que contienene las participaciones anticipadas asignadas a cada participante especificado.
        """
        # Verificación de condiciones
        if self.__participaciones_anticipadas is None:
            raise AttributeError(f'Ya se han repartido todas las participaciones.')
        if self.__participacion_secreto_anticipado is not None:
            raise AttributeError(f'Ya se han repartido las participaciones anticipadas.')
        if len(participantes_anticipados) != len(set(participantes_anticipados)):
            raise ValueError(f'Se han encontrado participantes duplicados.')
        for nombre in participantes_anticipados:
            if nombre not in self.participantes_numero:
                raise ValueError(f"El participante '{nombre}' no está registrado.")
        if self.reconstruccion <= len(participantes_anticipados):
            raise ValueError(f'El numero de participaciones anticipadas ({len(participantes_anticipados)}) debe ser menor o igual que el parámetro de privacidad ({self.reconstruccion - 1})')

        qc = self.__circuito
        r = self.reconstruccion
        elem_anticipadas = [self.participantes_numero[nombre] for nombre in participantes_anticipados]
        part_anticipadas = [self.__participaciones[idx - 1] for idx in elem_anticipadas]
        elem_resto = np.setdiff1d(range(1, 2 * r), elem_anticipadas) # Elementos de todos los participantes no anticipados
        resto = [self.__participaciones[idx - 1] for idx in elem_resto] # Participaciones de todos los participantes no anticipados
        no_entregadas = resto[:r - len(part_anticipadas) - 1]
        psi1 = part_anticipadas + no_entregadas # Se amplian las participaciones anticipadas hasta ser r-1
        psi2 = resto[r - len(part_anticipadas) - 1:-1] # Participaciones que no son anticipadas
        self.__participacion_secreto_anticipado = resto[-1] # Participación en la que se inicializará el secreto
        # Crear el estado máximamente entrelazado (sum |y>|y>)
        for qudit_1, qudit_2 in zip(psi1, psi2):
            qc.h(qudit_1)
            qc.cx(qudit_1, qudit_2)
        self.__participaciones_anticipadas = psi1
        self.elem_anticipadas_no_repartidas = elem_resto[:r - len(part_anticipadas) - 1].astype(int)
        return part_anticipadas

    def codificacion(self, secreto):
        """
        Crea las participaciones de todos los participantes de acuerdo al secreto recibido.
        Cada participación es un registro cuántico (clase QuantumRegister) con el identificador único del participante al que le corresponde.
        Si se han distribuido participaciones anticipadas, las participaciones serán coherentes con las mismas.
        :param secreto: Secreto que se quiere codificar entre todos los participantes.
        :return: Una lista que contienene las participaciones de cada participante que no ha participado en la distribución anticipada.
        """
        # Verificación de condiciones
        if self.__participaciones_anticipadas is None:
            raise AttributeError(f'Ya se han repartido todas las participaciones.')
        if not secreto.is_valid():
            raise ValueError(f'No se ha introducido un estado cuántico válido.')
        if secreto.dim != self.cuerpo.order:
            raise ValueError(f'Se esperaba un vector de estado de dimensión {self.cuerpo.order}, pero se ha recibido uno de dimensión {secreto.dim}.')

        qc = self.__circuito
        r = self.reconstruccion
        # Procedimiento estandar
        if len(self.__participaciones_anticipadas) == 0:
            x = np.arange(1, 2*r)
            qc.initialize(secreto, [self.__participaciones[0]])  # Inicializar el secreto
            for participacion in self.__participaciones[1:r]:  # Superponer todos los posibles valores de los coeficientes del polinomio
                qc.h(participacion)
            # Evaluar en cada registro los polinomios en los elementos de los participantes
            vandermonde = self.cuerpo(x)[:, None] ** np.arange(r)
            matriz = self.cuerpo(np.column_stack([vandermonde, np.vstack([np.eye(r-1), np.zeros((r, r-1))])])) # Matriz de evaluación
            matriz = extender_matriz(matriz) # Extender la matriz de numeros de F_q a vectores de F_2
            orden_participantes = [qubit for participacion in self.__participaciones for qubit in reversed(participacion)]  # Como qiskit es Little Endian, pero la matriz extendida está en Big Endian, hay que invertir el orden de los qubits de los participantes
        # Compartición anticipada
        else:
            qc.initialize(secreto, self.__participacion_secreto_anticipado)
            elem_anticipados = [self.participantes_numero[participacion.name] for participacion in self.__participaciones_anticipadas]
            elem_resto = np.setdiff1d(range(1, 2*r), elem_anticipados) # Elementos de todos los participantes no anticipados
            x = np.concat([elem_resto, self.elem_anticipadas_no_repartidas])
            elem_anticipados = self.cuerpo(elem_anticipados + [0])
            part_resto = [self.__participaciones[idx - 1] for idx in elem_resto] # Participaciones de todos los participantes no anticipados
            matriz_p1 = np.linalg.inv(elem_anticipados[:, None] ** np.arange(r))
            matriz_p2 = self.cuerpo(elem_resto)[:, None] ** np.arange(r)
            matriz = extender_matriz(matriz_p2 @ matriz_p1)
            orden_participantes = [qubit for participacion in part_resto for qubit in reversed(participacion)]
        qc.append(LinearFunction(matriz), orden_participantes)  # Aplicar matriz de evaluación
        self.__circuito = qc.decompose('Linear_function')
        self.__participaciones_anticipadas = None  # Eliminación de las participaciones anticipadas para mayor seguridad
        # Generar el resto de las participaciones reales
        return list(self.__participaciones[i-1] for i in x[x <= len(self.participantes_numero)])

    def decodificacion(self, participaciones):
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        Cada participación es un registro cuántico (clase QuantumRegister) con el identificador único del participante al que le corresponde.
        :param participaciones: Secuencia con las participaciones de los participantes que desean obtener el secreto.
        :return: El secreto hasta una fase global.
        """
        # Verificación de condiciones
        if self.__circuito is None:
            raise AttributeError(f'Ya se ha realizado el procedimiento de decodificación.')
        if self.__participaciones_anticipadas is not None:
            raise AttributeError(f'Todavía no se ha realizado el procedimiento de codificación')
        if len(participaciones) < self.reconstruccion:
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto')
        if len(participaciones) != len(set(participaciones)):
            raise ValueError(f'Se han encontrado participantes duplicados.')
        conjunto_participaciones = set(self.__participaciones)
        for participacion in participaciones:
            if participacion not in conjunto_participaciones:
                raise ValueError(f"El participante '{participacion.name}' no está registrado o ha entregado un qudit incorrecto.")

        r = self.reconstruccion
        qc = self.__circuito
        # Obetener los elementos asociados a cada participante
        elementos = self.cuerpo([self.participantes_numero[participacion.name] for participacion in participaciones[:r]])
        elementos_resto = self.cuerpo(np.setdiff1d(np.arange(2*r), elementos))
        orden_participantes = [qubit for participacion in participaciones[:r] for qubit in reversed(participacion)]
        matriz_p1 = np.linalg.inv(elementos[:,None]**np.arange(r)) # Matriz del primer paso del procediemiento de decodificacion
        matriz_p2 = elementos_resto[:,None]**np.arange(r) # Matriz del segundo paso del procediemiento de decodificacion
        matriz = extender_matriz(matriz_p2@matriz_p1)
        qc.append(LinearFunction(matriz), orden_participantes) # Realizar los dos pasos en uno
        sv = sim.run(qc.decompose('Linear_function')).result().get_statevector()
        elementos_traza = list(range((int(elementos[0])-1)*self.cuerpo.degree)) + list(range(int(elementos[0])*self.cuerpo.degree, (2*r-1)*self.cuerpo.degree))  # Posicion de los qubits a trazar
        self.__circuito = None  # Indicar que ya se ha realizado el procedimiento de decodificación
        return partial_trace(sv, elementos_traza).to_statevector()