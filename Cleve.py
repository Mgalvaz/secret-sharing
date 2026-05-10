import numpy as np
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.quantum_info import Statevector, partial_trace
from qiskit.circuit.library import LinearFunction
from qiskit_aer import AerSimulator
import re

from galois import GF

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
            raise ValueError(f'El cuerpo introducido debe tener como elemento base  el 2.')
        if cuerpo.order <= len(participantes):
            raise ValueError(f'El numero de participantes ({len(participantes)}) debe ser menor que el orden del cuerpo de trabajo ({cuerpo.order}).')
        if 2*r-l <= len(participantes):
            raise ValueError(f'El numero de participantes ({len(participantes)}) debe ser menor que 2*r-1 ({2*r-1}).')
        if len(participantes) < 2:
            raise ValueError(f'El numero de participantes ({len(participantes)}) debe ser mayor que 1.')
        if len(participantes) != len(set(participantes)):
            raise ValueError(f'Se han encontrado participantes duplicados.')
        if len(participantes) < r:
            raise ValueError(f'El parámetro de reconstrucción ({r}) debe ser menor o igual que el número de participantes ({len(participantes)}).')
        if r < 2:
            raise ValueError(f'El parámetro de reconstrucción ({r}) debe ser mayor que 1.')


        elementos_participantes = gf.Range(1, m + 1)
        qc = QuantumCircuit(*participantes)

        self._cuerpo = cuerpo
        self._reconstruccion = r
        self.__participaciones_anticipadas = None
        participantes_q = [QuantumRegister(dim_participantes, participante) for participante in participantes]
        QuantumRegister(2,'p').
        resto_q = [QuantumRegister(dim_participantes, f'p{i}') for i in range(len(participantes)+1, 2*r)]
        self.__participantes = participantes_q + resto_q
        self.circuito = QuantumCircuit(*self.__participantes)

        self._participantes_nombre = [None] # Array para pasar de numero -> nombre
        self._participantes_numero = {} # Diccionario para pasar nombre -> numero
        for i, nombre in enumerate(participantes, 1):
            self._participantes_nombre.append(nombre)
            self._participantes_numero[nombre] = i

    def crear_anticipadas(self, participantes_anticipados):
        """
        Crea participaciones participantes_anticipados para cada participante especificado.
        El formato de las participaciones es: (Identificador, Participación).
        :param participantes_anticipados: Listado de los participantes a entregar participaciones participantes_anticipados.
        :return: Una lista que contienene las participaciones participantes_anticipados asignadas a cada participante especificado.
        """

        # Verificación de condiciones
        if len(participantes_anticipados) < 1:
            raise ValueError(f'El numero de participaciones participantes_anticipados ({len(participantes_anticipados)}) debe ser al menos 1.')
        if self._reconstruccion - 1 < len(participantes_anticipados):
            raise ValueError(f'El numero de participaciones participantes_anticipados ({len(participantes_anticipados)}) debe ser menor o igual que el parámetro de privacidad ({self._reconstruccion - 1})')
        self._verificar_nombres(participantes_anticipados)

        # Generar las participaciones participantes_anticipados, que son elementos aleatorios del cuerpo
        aleatoriedad = array_aleatorio(self._cuerpo.order, len(participantes_anticipados))
        aleatoriedad_b64 = int_a_b64str(aleatoriedad, self._longitud_bytes)
        self.__participaciones_anticipadas = list(zip(participantes_anticipados, aleatoriedad_b64))
        return self.__participaciones_anticipadas

    def crear_participaciones(self, secreto):
        """
        Crea las participaciones de todos los participantes de acuerdo al secreto recibido.
        El formato de las participaciones es: (Identificador, Participación).
        Si se han distribuido participaciones participantes_anticipados, las participaciones serán coherentes con las mismas.
        :param secreto: Secreto que se quiere codificar entre todos los participantes.
        :return: Una lista que contienene las participaciones de cada participante que no ha participado en la distribución avanzada.
        """

        # Si se han repartido participaciones participantes_anticipados, se realiza compartición avanzada
        if self.__participaciones_anticipadas is not None:
            # Obtener el elemento asociado a cada participante y decodificar su participación
            nombres, valores_b64 = zip(*self.__participaciones_anticipadas)
            puntos_anticipados = self._cuerpo(list(self._participantes_numero[nombre] for nombre in nombres) + [0])
            valores_anticipados = self._cuerpo(b64str_a_int(valores_b64) + [secreto_i])
            x = np.setdiff1d(list(self._participantes_numero.values()), puntos_anticipados)

            # Se determina un polinomio de grado r-1 compatible con las participaciones participantes_anticipados
            lagrange = lagrange_poly(puntos_anticipados, valores_anticipados)
            if len(puntos_anticipados) < self._reconstruccion - 1: # Si el número de participaciones participantes_anticipados es menor que r-1, hay que completar el polinomio con aleatoriedad
                polinomio = lagrange + Poly.Roots(puntos_anticipados, field=self._cuerpo) * polinomio_aleatorio(self._cuerpo, self._reconstruccion - len(puntos_anticipados) - 1)
            else: # Si no, el único polinomio disponible es el de Lagrange
                polinomio = lagrange
            del self.__participaciones_anticipadas # Eliminación de las participaciones participantes_anticipados para mayor seguridad

        # Si no, se sigue el proceso estándar
        else:
            polinomio = Poly(array_aleatorio(self._cuerpo.order, self._reconstruccion - 1) + [secreto_i], field=self._cuerpo)
            x = list(self._participantes_numero.values())

        # Generar el resto de las participaciones
        participaciones_b64 = int_a_b64str(polinomio(x), self._longitud_bytes)
        return list(zip((self._participantes_nombre[p] for p in x), participaciones_b64))

    def recuperar_secreto_v1(self, participaciones):
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (Identificador, Participación).
        Esta versión reconstruye primero el polinomio generador y a partir de él, devuelve el secreto.
        :param participaciones: Secuencia con las participaciones de los participantes que desean obtener el secreto.
        :return: El secreto.
        """

        # Verificación de condiciones
        if len(participaciones) < self._reconstruccion:
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto')
        nombres, valores_b64 = zip(*participaciones[:self._reconstruccion])
        self._verificar_nombres(nombres)

        # Obtener el elemento asociado a cada participante y decodificar su participación
        puntos = self._cuerpo(list(self._participantes_numero[nombre] for nombre in nombres))
        valores = self._cuerpo(b64str_a_int(valores_b64))

        # Reconstruir el polinomio generador y el secreto como su coeficiente independiente
        polinomio = lagrange_poly(puntos, valores)
        return int_a_bytes(polinomio.coefficients(order="asc")[0])

    def recuperar_secreto_v2(self, participaciones):
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (Identificador, Participación).
        Esta versión reconstruye el secreto a partir de la fórmula del polinomio interpolador de Lagrange evaluado en 0.
        :param participaciones: Secuencia con las participaciones de los participantes que desean obtener el secreto.
        :return: El secreto.
        """

        # Verificación de condiciones
        r = self._reconstruccion
        if len(participaciones) < r:
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto')
        nombres, valores_b64 = zip(*participaciones[:r])
        self._verificar_nombres(nombres)

        # Obtener el elemento asociado a cada participante y decodificar su participación
        puntos = self._cuerpo(list(self._participantes_numero[nombre] for nombre in nombres))
        valores = self._cuerpo(b64str_a_int(valores_b64))

        # Calcular el valor del polinomio generador en 0 sin reconstruirlo
        mascara = ~np.eye(r, dtype=bool) # Máscara de los elementos x_h de la fórmula
        puntos_matriz = np.broadcast_to(puntos, (r, r)) # Se crea una matriz que cada fila es el array puntos
        puntos_matriz = puntos_matriz[mascara].reshape(r - 1, r) # Al usar la mascara, la matriz se aplana por lo que hay que usar reshape (trabajaremos por columnas)
        denominador = np.prod(puntos_matriz - puntos, axis=0) # Productorio del denominador
        numerador = np.prod(puntos_matriz, axis=0) # Productorio del numerador
        coef = numerador / denominador # Cálculo de l_j
        return int_a_bytes(np.sum(valores * coef)) # Se devuelve la suma y_j * l_j

    recuperar_secreto = recuperar_secreto_v2 # Alias para recuperar secreto version 2

    def _verificar_nombres(self, nombres):
        """
        Verifica que los participantes sean válidos, es decir, que no haya nombres duplicados y todos los nombres estén registrados como participantes.
        :param nombres: La secuencia de nombres que se quiere comprobar
        """
        # Comprobar no elementos duplicados
        if len(nombres) != len(set(nombres)):
            raise ValueError(f'Se han encontrado participantes duplicados.')
        # Comprobar que los participantes existen
        conjunto_nombres = self._participantes_numero
        for nombre in nombres:
            if nombre not in conjunto_nombres:
                raise ValueError(f'El participante {nombre} no está registrado.')



if __name__ == "__main__":
    n = int(input('numero participantes: '))
    r = int(input('parámetro reconstruccion: '))
    m = 2*r-1
    dim_participantes = (m - 1).bit_length()
    todos_participantes = np.arange(dim_participantes*m)
    gf = GF(2**dim_participantes)
    l = 1
    participantes = [QuantumRegister(dim_participantes, f'p{i}') for i in range(m)]
    elementos_participantes = gf.Range(1, m+1)
    qc = QuantumCircuit(*participantes)


    ## PASO 1: inicializar los coeficientes de todos los posibles polinomios
    # Inicializar el secreto en los primeros l registros

    texto = input('escriba el secreto: ')
    numeros = re.findall(r'\d+', texto)
    svect = [int(n) for n in numeros]
    print('secreto: ', ' + '.join([format(i, f'0{dim_participantes}b') * svect[i] for i in range(2**dim_participantes)]))
    svect = svect / np.linalg.norm(svect)
    svect = Statevector(svect)
    qc.initialize(svect, [participantes[0]]) # Inicializar el secreto
    for participante in participantes[l:r]: # Superponer todos los posibles valores de los coeficientes del polinomio
        qc.h(participante)
    # Evaluar en cada registro los polinomios en los elementos de los participantes
    vandermonde = elementos_participantes[:, None] ** np.arange(r)
    matriz = extender_matriz(gf(np.column_stack([vandermonde, np.zeros((m, r-1))])))
    orden_participantes = [part[i] for part in participantes for i in range(dim_participantes-1, -1, -1)] # Como qiskit es Little Endian, pero la matriz extendida está en Big Endian, hay que invertir el orden de los qubits de los participantes
    qc.append(LinearFunction(matriz), orden_participantes)
    qc = qc.decompose('Linear_function')

    #sv = run_statevector(qc)
    #print('Estado tras la codificación:')
    #print('traza parcial del primer (último) participante del estado codificado: ', partial_trace(sv, np.setdiff1d(todos_participantes, np.arange(dim_participantes)).tolist()))
    #print('traza parcial del primer y segundo (último) participante del estado codificado: ', partial_trace(sv, np.setdiff1d(todos_participantes, np.arange(dim_participantes*2)).tolist()))

    ### DECODIFICACION DEL SECRETO A PARTIR DE LOS PARTICIPANTES P2 Y P3
    numeros = [2,3]
    resto = np.setdiff1d(elementos_participantes, numeros)
    participantes_reconstruccion = [participantes[i-1] for i in numeros]
    vandermonde = gf(numeros)[:, None] ** np.arange(r)
    matriz_r1 = np.linalg.inv(vandermonde)
    matriz_r2 = np.vstack([[1] + [0]*(r-1) , resto[:, None] ** np.arange(r)])
    matriz = extender_matriz(matriz_r2@matriz_r1)
    orden_participantes_reconstruccion = [part[i] for part in participantes_reconstruccion for i in range(dim_participantes - 1, -1, -1)]  # Como qiskit es Little Endian, pero la matriz extendida está en Big Endian, hay que invertir el orden de los qubits de los participantes
    qc.append(LinearFunction(matriz), orden_participantes_reconstruccion)
    qc = qc.decompose('Linear_function')

    sv = run_statevector(qc)
    print('Estado tras paso 2:')
    print('traza parcial del segundo (primero) participante del estado descodificado: ', partial_trace(sv, np.setdiff1d(todos_participantes, np.arange(dim_participantes, dim_participantes*2)).tolist()))
    print('vector de estado del segundo (primero) participante del estado descodificado: ', partial_trace(sv, np.setdiff1d(todos_participantes, np.arange(dim_participantes, dim_participantes*2)).tolist()).to_statevector())

