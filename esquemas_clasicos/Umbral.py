import numpy as np
from galois import Poly, lagrange_poly

from utils import array_aleatorio, polinomio_aleatorio, bytes_a_int, int_a_bytes, int_a_b64str, b64str_a_int

class Shamir:
    r"""
    Esquema de compartición de secretos de Shamir sobre el cuerpo $\mathbb{F}_{p^m}$.

    Ejemplo:
        Crea un esquema de Shamir de umbral-(4,6) sobre el cuerpo $\mathbb{F}_{3^5}$ para los participantes ['a', 'b', 'c', 'd', 'e', 'f'].

        .. ipython:: python

            cuerpo = galois.GF(3, 5)
            sh = Shamir(cuerpo, 4, ['a', 'b', 'c', 'd', 'e', 'f'])
    """
    def __init__(self, cuerpo, r, participantes):
        r"""
        Crea un esquema de compartición de secretos de Shamir sobre el cuerpo $\mathbb{F}_{p^m}$.
        :param cuerpo: El cuerpo finito sobre sobre el que el esquema está construido.
        :param r: Parámetro de reconstrucción del esquema (número mínimo de participantes necesarios para reconstruir el secreto).
        :param participantes: Lista de los nombres únicos de cada participante del esquema.
        """
        # Verificación de condiciones
        if cuerpo.order <= len(participantes):
            raise ValueError(f'El numero de participantes ({len(participantes)}) debe ser menor que el orden del cuerpo de trabajo ({cuerpo.order}).')
        if len(participantes) != len(set(participantes)):
            raise ValueError(f'Se han encontrado participantes duplicados.')
        if len(participantes) < r:
            raise ValueError(f'El parámetro de reconstrucción ({r}) debe ser menor o igual que el número de participantes ({len(participantes)}).')
        if r < 2:
            raise ValueError(f'El parámetro de reconstrucción ({r}) debe ser mayor que 1.')

        self.cuerpo = cuerpo
        self.reconstruccion = r
        self.__participaciones_anticipadas = []
        self.longitud_bytes = ((cuerpo.order - 1).bit_length() + 7) // 8
        self.participantes_nombre = np.array([None] + participantes)
        self.participantes_numero = {nombre: i for i, nombre in enumerate(participantes, 1)}

    def comparticion_anticipada(self, participantes_anticipados):
        """
        Crea participaciones anticipadas para cada participante especificado.
        El formato de las participaciones es: (nombre, participación).
        :param participantes_anticipados: Listado de los participantes a entregar participaciones anticipadas.
        :return: Una lista que contienene las participaciones anticipadas asignadas a cada participante especificado.
        """
        # Verificación de condiciones
        if self.__participaciones_anticipadas is None:
            raise AttributeError(f'Ya se han repartido todas las participaciones.')
        if len(self.__participaciones_anticipadas) > 0:
            conjunto_nombres_anticipados = set(list(zip(*self.__participaciones_anticipadas))[0])
            for nombre in participantes_anticipados:
                if nombre in conjunto_nombres_anticipados:
                    raise ValueError(f"El participante '{nombre}' ya ha recibido una participación anticipada.")
        self._verificar_nombres(participantes_anticipados)
        if self.reconstruccion - len(self.__participaciones_anticipadas) <= len(participantes_anticipados):
            raise ValueError(f'El numero de participaciones anticipadas ({len(participantes_anticipados) + len(self.__participaciones_anticipadas)}) debe ser menor o igual que el parámetro de privacidad ({self.reconstruccion - 1}).')

        # Generar las participaciones anticipadas, que son elementos aleatorios del cuerpo
        aleatoriedad = array_aleatorio(self.cuerpo.order, len(participantes_anticipados))
        aleatoriedad_b64 = int_a_b64str(aleatoriedad, self.longitud_bytes)
        extras = list(zip(participantes_anticipados, aleatoriedad_b64))
        self.__participaciones_anticipadas.extend(extras)
        return extras

    def codificacion(self, secreto):
        """
        Crea las participaciones de todos los participantes de acuerdo al secreto recibido.
        El formato de las participaciones es: (nombre, participación).
        Si se han distribuido participaciones anticipadas, las participaciones serán coherentes con las mismas.
        :param secreto: Secreto que se quiere codificar entre todos los participantes.
        :return: Una lista que contienene las participaciones de cada participante que no ha participado en la distribución anticipada.
        """
        if self.__participaciones_anticipadas is None:
            raise AttributeError(f'Ya se han repartido todas las participaciones.')
        secreto_i = bytes_a_int(secreto)
        if secreto_i >= self.cuerpo.order:
            raise ValueError(f'El secreto proporcionado debe ser menor que el orden del cuerpo de trabajo ({self.cuerpo.order}).')

        # Procedimiento estandar
        if len(self.__participaciones_anticipadas) == 0:
            polinomio = Poly(array_aleatorio(self.cuerpo.order, self.reconstruccion - 1) + [secreto_i], field=self.cuerpo)
            x = np.arange(1, len(self.participantes_nombre))
        # Compartición anticipada
        else:
            # Obtener el elemento asociado a cada participante y decodificar su participación
            nombres, valores_b64 = zip(*self.__participaciones_anticipadas)
            puntos_anticipados = self.cuerpo(list(self.participantes_numero[nombre] for nombre in nombres) + [0])
            valores_anticipados = self.cuerpo(b64str_a_int(valores_b64) + [secreto_i])
            x = np.setdiff1d(np.arange(1, len(self.participantes_nombre)), puntos_anticipados)

            # Se determina un polinomio de grado r-1 compatible con las participaciones anticipadas
            lagrange = lagrange_poly(puntos_anticipados, valores_anticipados)
            if len(puntos_anticipados) < self.reconstruccion - 1: # Si el número de participaciones anticipadas es menor que r-1, hay que completar el polinomio con aleatoriedad
                polinomio = lagrange + Poly.Roots(puntos_anticipados, field=self.cuerpo) * polinomio_aleatorio(self.cuerpo, self.reconstruccion - len(puntos_anticipados) - 2)
            else: # Si no, el único polinomio disponible es el de Lagrange
                polinomio = lagrange

        # Generar el resto de las participaciones
        participaciones_b64 = int_a_b64str(polinomio(x), self.longitud_bytes)
        self.__participaciones_anticipadas = None  # Se eliminan las participaciones anticipadas almacenadas para mayor seguridad
        return list(zip(self.participantes_nombre[x], participaciones_b64))

    def _decodificacion_alternativa(self, participaciones):
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (nombre, participación).
        Esta versión reconstruye primero el polinomio generador y a partir de él, devuelve el secreto.
        :param participaciones: Secuencia con las participaciones de los participantes que desean obtener el secreto.
        :return: El secreto.
        """
        # Verificación de condiciones
        if len(participaciones) < self.reconstruccion:
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto')
        nombres, valores_b64 = zip(*participaciones[:self.reconstruccion])
        self._verificar_nombres(nombres)

        # Obtener el elemento asociado a cada participante y decodificar su participación
        puntos = self.cuerpo(list(self.participantes_numero[nombre] for nombre in nombres))
        valores = self.cuerpo(b64str_a_int(valores_b64))
        # Reconstruir el polinomio generador y el secreto como su coeficiente independiente
        polinomio = lagrange_poly(puntos, valores)
        return int_a_bytes(polinomio.coefficients(order="asc")[0])

    def decodificacion(self, participaciones):
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (nombre, participación).
        Esta versión reconstruye el secreto a partir de la fórmula del polinomio interpolador de Lagrange evaluado en 0.
        :param participaciones: Secuencia con las participaciones de los participantes que desean obtener el secreto.
        :return: El secreto.
        """
        # Verificación de condiciones
        r = self.reconstruccion
        if len(participaciones) < r:
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto.')
        nombres, valores_b64 = zip(*participaciones[:r])
        self._verificar_nombres(nombres)

        # Obtener el elemento asociado a cada participante y decodificar su participación
        puntos = self.cuerpo(list(self.participantes_numero[nombre] for nombre in nombres))
        valores = self.cuerpo(b64str_a_int(valores_b64))
        # Calcular el valor del polinomio generador en 0 sin reconstruirlo
        mascara = ~np.eye(r, dtype=bool) # Máscara de los elementos x_h de la fórmula
        puntos_matriz = np.broadcast_to(puntos, (r, r)) # Matriz en la que cada fila es el array puntos
        puntos_matriz = puntos_matriz[mascara].reshape((r, r-1)).T # Al usar la mascara, la matriz se aplana, por lo que se usa reshape
        numerador = np.prod(puntos_matriz, axis=0)  # Productorio del numerador
        denominador = np.prod(puntos_matriz - puntos, axis=0) # Productorio del denominador
        coef = numerador / denominador # Cálculo de l_i
        return int_a_bytes(np.sum(valores * coef)) # Se devuelve la suma y_i * l_i

    def _verificar_nombres(self, nombres):
        """
        Verifica que los participantes sean válidos, es decir, que no haya nombres duplicados y todos los nombres estén registrados como participantes.
        :param nombres: La secuencia de nombres que se quiere comprobar
        """
        # Comprobar no elementos duplicados
        if len(nombres) != len(set(nombres)):
            raise ValueError(f'Se han encontrado participantes duplicados.')
        # Comprobar que los participantes existen
        conjunto_nombres = self.participantes_numero
        for nombre in nombres:
            if nombre not in conjunto_nombres:
                raise ValueError(f"El participante '{nombre}' no está registrado.")


class Simplificado:
    r"""
    Esquema de compartición de secretos simplificado sobre cuerpo $\mathbb{F}_{p^m}$.

    Ejemplo:
        Crea un esquema de umbral-(5,5) sobre el cuerpo $\mathbb{F}_{5^4}$ para los participantes ['a', 'b', 'c', 'd', 'e'].

        .. ipython:: python

            cuerpo = galois.GF(5, 4)
            sh = Simplificado(cuerpo, ['a', 'b', 'c', 'd', 'e'])
    """
    def __init__(self, cuerpo, participantes):
        r"""
        Crea un esquema de compartición de secretos simlificado sobre el cuerpo $\mathbb{F}_{p^m}$.
        :param cuerpo: El cuerpo finito sobre sobre el que el esquema está construido.
        :param participantes: Lista de los nombres únicos de cada participante del esquema.
        """
        # Verificación de condiciones
        if len(participantes) < 2:
            raise ValueError(f'El numero de participantes ({len(participantes)}) debe ser mayor que 1.')
        if len(participantes) != len(set(participantes)):
            raise ValueError(f'Se han encontrado participantes duplicados.')

        self.participantes = participantes
        self.cuerpo = cuerpo
        self.longitud_bytes = ((cuerpo.order - 1).bit_length() + 7) // 8
        self.__participaciones_anticipadas = []

    def comparticion_anticipada(self, participantes_anticipados):
        """
        Crea participaciones anticipadas para cada participante especificado.
        El formato de las participaciones es: (nombre, participación).
        :param participantes_anticipados: Listado de los participantes a entregar participaciones anticipadas.
        :return: Una lista que contienene las participaciones anticipadas asignadas a cada participante especificado.
        """
        # Verificación de condiciones
        if self.__participaciones_anticipadas is None:
            raise AttributeError(f'Ya se han repartido todas las participaciones.')
        if len(self.__participaciones_anticipadas) > 0:
            conjunto_nombres_anticipados = set(list(zip(*self.__participaciones_anticipadas))[0])
            for nombre in participantes_anticipados:
                if nombre in conjunto_nombres_anticipados:
                    raise ValueError(f"El participante '{nombre}' ya ha recibido una participación anticipada.")
        self._verificar_nombres(participantes_anticipados)
        if len(self.participantes) - len(self.__participaciones_anticipadas) <= len(participantes_anticipados):
            raise ValueError(f'El numero de participaciones anticipadas ({len(participantes_anticipados) + len(self.__participaciones_anticipadas)}) debe ser menor o igual que el parámetro de privacidad ({len(self.participantes) - 1}).')

        # Generar las participaciones anticipadas, que son elementos aleatorios del cuerpo
        aleatoriedad = array_aleatorio(self.cuerpo.order, len(participantes_anticipados))
        aleatoriedad_b64 = int_a_b64str(aleatoriedad, self.longitud_bytes)
        extras = list(zip(participantes_anticipados, aleatoriedad_b64))
        self.__participaciones_anticipadas.extend(extras)
        return extras

    def codificacion(self, secreto):
        """
        Crea las participaciones de todos los participantes de acuerdo al secreto recibido.
        El formato de las participaciones es: (nombre, participación).
        Si se han distribuido participaciones anticipadas, las participaciones serán coherentes con las mismas.
        :param secreto: Secreto que se quiere codificar entre todos los participantes.
        :return: Una lista que contienene las participaciones de cada participante que no ha participado en la distribución anticipada.
        """
        if self.__participaciones_anticipadas is None:
            raise AttributeError(f'Ya se han repartido todas las participaciones.')
        secreto_i = bytes_a_int(secreto)
        if secreto_i >= self.cuerpo.order:
            raise ValueError(f'El secreto proporcionado debe ser menor que el orden del cuerpo de trabajo ({self.cuerpo.order}).')

        # Proceso estándar
        if len(self.__participaciones_anticipadas) == 0:
            x = self.participantes
            suma = self.cuerpo(0)

        # Compartición anticipada
        else:
            # Obtener el elemento asociado a cada participante y decodificar su participación
            nombres, valores_b64 = zip(*self.__participaciones_anticipadas)
            valores = self.cuerpo(b64str_a_int(valores_b64))
            x = np.setdiff1d(self.participantes, nombres).tolist()
            suma = valores.sum()

        # Se generan el resto de participaciones
        participaciones = self.cuerpo(array_aleatorio(self.cuerpo.order, len(x) - 1))
        participaciones = np.append(participaciones, self.cuerpo(secreto_i) - participaciones.sum() - suma) # La última debe ser igual al secreto menos la suma de todas las anteriores
        participaciones_b64 = int_a_b64str(participaciones, self.longitud_bytes)
        self.__participaciones_anticipadas = None  # Se eliminan las participaciones anticipadas almacenadas para mayor seguridad
        return list(zip(x, participaciones_b64))

    def decodificacion(self, participaciones):
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (nombre, participación).
        :param participaciones: Secuencia con las participaciones de los participantes que desean obtener el secreto.
        :return: El secreto.
        """
        # Verificación de condiciones
        if len(participaciones) < len(self.participantes):
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto.')
        nombres, valores_b64 = zip(*participaciones[:len(self.participantes)])
        self._verificar_nombres(nombres)

        # Obtener las participaciones
        valores = self.cuerpo(b64str_a_int(valores_b64))
        # El secreto es la suma de todas las participaciones
        return int_a_bytes(valores.sum())

    def _verificar_nombres(self, nombres):
        """
        Verifica que los participantes sean válidos, es decir, que no haya nombres duplicados y todos los nombres estén registrados como participantes.
        :param nombres: La secuencia de nombres que se quiere comprobar
        """
        # Comprobar no elementos duplicados
        if len(nombres) != len(set(nombres)):
            raise ValueError(f'Se han encontrado participantes duplicados.')
        # Comprobar que los participantes existen
        conjunto_nombres = set(self.participantes)
        for nombre in nombres:
            if nombre not in conjunto_nombres:
                raise ValueError(f"El participante '{nombre}' no está registrado.")