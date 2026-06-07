import numpy as np
from galois import Poly, lagrange_poly

from utils import *


class ShamirRampa:
    r"""
        Esquema de compartición de secretos de Shamir en rampa sobre el cuerpo $\mathbb{F}_{p^m}$.

        Ejemplo:
            Crea un esquema de Shamir en rampa con parámetros r = 4 y l = 3 sobre el cuerpo $\mathbb{F}_{3^5}$ para los participantes ['a', 'b', 'c', 'd', 'e', 'f'].

            .. ipython:: python

                cuerpo = galois.GF(3**5)
                rsh = ShamirRampa(cuerpo, 4, 3, ['a', 'b', 'c', 'd', 'e', 'f'])
        """

    def __init__(self, cuerpo, r, l,  participantes):
        r"""
        Crea un esquema de compartición de secretos de Shamir en rampa sobre el cuerpo $\mathbb{F}_{p^m}$.
        :param cuerpo: El cuerpo finito sobre sobre el que el esquema está construido.
        :param r: Parámetro de reconstrucción del esquema (número mínimo de participantes necesarios para reconstruir el secreto).
        :param l: Longitud del secreto que se quiere repartir.
        :param participantes: Lista de los identificadores únicos de cada participante del esquema.
        """
        # Verificación de condiciones
        if cuerpo.order <= len(participantes):
            raise ValueError(f'El numero de participantes ({len(participantes)}) debe ser menor que el orden del cuerpo de trabajo ({cuerpo.order}).')
        if len(participantes) != len(set(participantes)):
            raise ValueError(f'Se han encontrado participantes duplicados.')
        if len(participantes) < r:
            raise ValueError(f'El parámetro de reconstrucción ({r}) debe ser menor o igual que el número de participantes ({len(participantes)}).')
        if r <= l:
            raise ValueError(f'La longitud del secreto ({l}) debe ser menor que el parámetro de reconstrucción ({r}).')
        if  l < 2:
            raise ValueError(f'La longitud del secreto ({l}) debe ser mayor que 1.')

        self.cuerpo = cuerpo
        self.reconstruccion = r
        self.longitud_secreto = l
        self.__participaciones_anticipadas = []
        self.longitud_bytes = ((cuerpo.order - 1).bit_length() + 7) // 8
        self.participantes_nombre = np.array([None] + participantes)
        self.participantes_numero = {nombre: i for i, nombre in enumerate(participantes, 1)}

    def comparticion_anticipada(self, participantes_anticipados):
        """
        Crea participaciones anticipadas para cada participante especificado.
        El formato de las participaciones es: (Identificador, Participación).
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
                    raise ValueError(f'El participante {nombre} ya ha recibido una participación anticipada.')
        self._verificar_nombres(participantes_anticipados)
        if self.reconstruccion - self.longitud_secreto < len(participantes_anticipados) + len(self.__participaciones_anticipadas):
            raise ValueError(f'El numero de participaciones anticipadas ({len(participantes_anticipados) + len(self.__participaciones_anticipadas)}) debe ser menor o igual que el parámetro de privacidad ({self.reconstruccion - self.longitud_secreto})')

        # Generar las participaciones anticipadas, que son elementos aleatorios del cuerpo
        aleatoriedad = array_aleatorio(self.cuerpo.order, len(participantes_anticipados))
        aleatoriedad_b64 = int_a_b64str(aleatoriedad, self.longitud_bytes)
        extras = list(zip(participantes_anticipados, aleatoriedad_b64))
        self.__participaciones_anticipadas.extend(extras)
        return extras

    def codificacion(self, secreto):
        """
        Crea las participaciones de todos los participantes de acuerdo al secreto recibido.
        El formato de las participaciones es: (Identificador, Participación).
        Si se han distribuido participaciones anticipadas, las participaciones serán coherentes con las mismas.
        :param secreto: Secreto que se quiere codificar entre todos los participantes.
        :return: Una lista que contienene las participaciones de cada participante que no ha participado en la distribución anticipada.
        """
        # Verificación de condiciones
        if self.__participaciones_anticipadas is None:
            raise AttributeError(f'Ya se han repartido todas las participaciones.')
        if len(secreto) != self.longitud_secreto:
            raise ValueError(f'Se esperaba un secreto de longitud {self.longitud_secreto}, pero se ha recibido uno de longitud {len(secreto)}.')

        secreto_i = bytes_a_int(secreto)
        # Procedimiento estandar
        if len(self.__participaciones_anticipadas) == 0:
            polinomio = Poly(secreto_i + array_aleatorio(self.cuerpo.order, self.reconstruccion - self.longitud_secreto), field=self.cuerpo, order='asc')
            x = np.arange(1, len(self.participantes_nombre))
            # Generar el resto de las participaciones
            participaciones_b64 = int_a_b64str(polinomio(x), self.longitud_bytes)
        # Compartición anticipada
        else:
            # Obtener el elemento asociado a cada participante y decodificar su participación
            nombres, valores_b64 = zip(*self.__participaciones_anticipadas)
            puntos_anticipados = self.cuerpo(list(self.participantes_numero[nombre] for nombre in nombres))
            valores_anticipados = self.cuerpo(b64str_a_int(valores_b64))
            x = np.setdiff1d(np.arange(1, len(self.participantes_nombre)), puntos_anticipados)
            # Se determina un polinomio de grado r-1 compatible con las participaciones anticipadas
            polinomio_secreto = Poly(secreto_i, field=self.cuerpo, order='asc') # f_s
            lagrange = lagrange_poly(puntos_anticipados, (valores_anticipados - polinomio_secreto(puntos_anticipados)) / puntos_anticipados ** self.longitud_secreto)
            if len(puntos_anticipados) < self.reconstruccion - self.longitud_secreto:  # Si el número de participaciones anticipadas es menor que r-l, hay que completar el polinomio con aleatoriedad
                polinomio = lagrange + Poly.Roots(puntos_anticipados, field=self.cuerpo) * polinomio_aleatorio(self.cuerpo, self.reconstruccion - self.longitud_secreto - len(puntos_anticipados) - 1)
            else:  # Si no, el único polinomio disponible es el de Lagrange
                polinomio = lagrange
            # Generar el resto de las participaciones
            participaciones_b64 = int_a_b64str(polinomio_secreto(x) + self.cuerpo(x) ** self.longitud_secreto * polinomio(x), self.longitud_bytes)

        self.__participaciones_anticipadas = None  # Se eliminan las participaciones anticipadas almacenadas para mayor seguridad
        return list(zip(self.participantes_nombre[x], participaciones_b64))

    def decodificacion(self, participaciones):
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (Identificador, Participación).
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
        # Reconstruir el polinomio generador y el secreto como sus l primeros coeficientes
        polinomio = lagrange_poly(puntos, valores)
        return int_a_bytes(polinomio.coefficients(order="asc")[:self.longitud_secreto])

    def _verificar_nombres(self, nombres):
        """
        Verifica que los participantes sean válidos, es decir, que no haya nombres duplicados y todos los nombres estén registrados como participantes.
        :param nombres: La secuencia de nombres que se quiere comprobar
        """
        # Comprobar no elementos duplicados
        if len(nombres) != len(set(nombres)):
            raise ValueError(f"Se han encontrado participantes duplicados.")
        # Comprobar que los participantes existen
        conjunto_nombres = self.participantes_numero
        for nombre in nombres:
            if nombre not in conjunto_nombres:
                raise ValueError(f"El participante {nombre} no está registrado.")

class McElieceSarwate:
    r"""
        Esquema de compartición de secretos de McEliece-Sarwate sobre el cuerpo $\mathbb{F}_{p^m}$.

        Ejemplo:
            Crea un esquema de McEliece-Sarwate con parámetros r = 4 y l = 3 sobre el cuerpo $\mathbb{F}_{3^5}$ para los participantes ['a', 'b', 'c', 'd', 'e', 'f'].

            .. ipython:: python

                cuerpo = galois.GF(3**5)
                rsh = McElieceSarwate(cuerpo, 4, 3, ['a', 'b', 'c', 'd', 'e', 'f'])
        """

    def __init__(self, cuerpo, r, l,  participantes):
        r"""
        Crea un esquema de compartición de secretos de McEliece-Sarwate sobre el cuerpo $\mathbb{F}_{p^m}$.
        :param cuerpo: El cuerpo finito sobre sobre el que el esquema está construido.
        :param r: Parámetro de reconstrucción del esquema (número mínimo de participantes necesarios para reconstruir el secreto).
        :param l: Longitud del secreto que se quiere repartir.
        :param participantes: Lista de los identificadores únicos de cada participante del esquema.
        """
        # Verificación de condiciones
        if cuerpo.order - l < len(participantes):
            raise ValueError(f'El numero de participantes ({len(participantes)}) debe ser menor o igual que el orden del cuerpo de trabajo menos la longitud del secreto ({cuerpo.order - l}).')
        if len(participantes) != len(set(participantes)):
            raise ValueError(f'Se han encontrado participantes duplicados.')
        if len(participantes) < r:
            raise ValueError(f'El parámetro de reconstrucción ({r}) debe ser menor o igual que el número de participantes ({len(participantes)}).')
        if r <= l:
            raise ValueError(f'La longitud del secreto ({l}) debe ser menor que el parámetro de reconstrucción ({r}).')
        if  l < 2:
            raise ValueError(f'La longitud del secreto ({l}) debe ser mayor que 1.')

        self.cuerpo = cuerpo
        self.reconstruccion = r
        self.longitud_secreto = l
        self.__participaciones_anticipadas = []
        self.longitud_bytes = ((cuerpo.order - 1).bit_length() + 7) // 8
        self.participantes_nombre = np.array([None] * l + participantes)
        self.participantes_numero = {nombre: i for i, nombre in enumerate(participantes, l)}

    def comparticion_anticipada(self, participantes_anticipados):
        """
        Crea participaciones anticipadas para cada participante especificado.
        El formato de las participaciones es: (Identificador, Participación).
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
                    raise ValueError(f'El participante {nombre} ya ha recibido una participación anticipada.')
        self._verificar_nombres(participantes_anticipados)
        if self.reconstruccion - self.longitud_secreto < len(participantes_anticipados) + len(self.__participaciones_anticipadas):
            raise ValueError(f'El numero de participaciones anticipadas ({len(participantes_anticipados) + len(self.__participaciones_anticipadas)}) debe ser menor o igual que el parámetro de privacidad ({self.reconstruccion - self.longitud_secreto})')

        # Generar las participaciones anticipadas, que son elementos aleatorios del cuerpo
        aleatoriedad = array_aleatorio(self.cuerpo.order, len(participantes_anticipados))
        aleatoriedad_b64 = int_a_b64str(aleatoriedad, self.longitud_bytes)
        extras = list(zip(participantes_anticipados, aleatoriedad_b64))
        self.__participaciones_anticipadas.extend(extras)
        return extras

    def codificacion(self, secreto):
        """
        Crea las participaciones de todos los participantes de acuerdo al secreto recibido.
        El formato de las participaciones es: (Identificador, Participación).
        Si se han distribuido participaciones anticipadas, las participaciones serán coherentes con las mismas.
        :param secreto: Secreto que se quiere codificar entre todos los participantes.
        :return: Una lista que contienene las participaciones de cada participante que no ha participado en la distribución anticipada.
        """
        # Verificación de condiciones
        if self.__participaciones_anticipadas is None:
            raise AttributeError(f'Ya se han repartido todas las participaciones.')
        if len(secreto) != self.longitud_secreto:
            raise ValueError(f'Se esperaba un secreto de longitud {self.longitud_secreto}, pero se ha recibido uno de longitud {len(secreto)}.')

        secreto_i = self.cuerpo(bytes_a_int(secreto))
        alpha = self.cuerpo.Range(0, self.longitud_secreto)
        # Procedimiento estandar
        if len(self.__participaciones_anticipadas) == 0:
            x = np.arange(self.longitud_secreto, len(self.participantes_nombre))
            # Hay que construir un polinomio que interpole al secreto en sus respectivos puntos y que sea de grado r - 1
            lagrange = lagrange_poly(alpha, secreto_i)
            polinomio = lagrange + Poly.Roots(alpha, field=self.cuerpo) * polinomio_aleatorio(self.cuerpo,self.reconstruccion - self.longitud_secreto - 1)
        # Compartición anticipada
        else:
            # Obtener el elemento asociado a cada participante y decodificar su participación
            nombres_anticipados, valores_anticipados_b64 = zip(*self.__participaciones_anticipadas)
            puntos_anticipados = self.cuerpo(list(self.participantes_numero[nombre] for nombre in nombres_anticipados))
            valores_anticipados = self.cuerpo(b64str_a_int(valores_anticipados_b64))
            x = np.setdiff1d(list(self.participantes_numero.values()), puntos_anticipados)
            # Se determina un polinomio de grado r-1 compatible con las participaciones anticipadas
            puntos_lagrange = self.cuerpo(np.concatenate([alpha, puntos_anticipados]))
            valores_lagrange = self.cuerpo(np.concatenate([secreto_i, valores_anticipados]))
            lagrange = lagrange_poly(puntos_lagrange, valores_lagrange)
            if len(puntos_anticipados) < self.reconstruccion - self.longitud_secreto:  # Si el número de participaciones anticipadas es menor que r-l, hay que completar el polinomio con aleatoriedad
                polinomio = lagrange + Poly.Roots(np.concatenate([alpha, puntos_anticipados]), field=self.cuerpo) * polinomio_aleatorio(self.cuerpo, self.reconstruccion - self.longitud_secreto - len(puntos_anticipados) - 1)
            else:  # Si no, el único polinomio disponible es el de Lagrange
                polinomio = lagrange

        # Generar el resto de las participaciones
        participaciones_b64 = int_a_b64str(polinomio(x), self.longitud_bytes)
        self.__participaciones_anticipadas = None  # Eliminación de las participaciones anticipadas para mayor seguridad
        return list(zip(self.participantes_nombre[x], participaciones_b64))

    def _decodificacion_alternativa(self, participaciones):
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (Identificador, Participación).
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
        # Reconstruir el polinomio generador y el secreto como su evaluacion en los elementos alpha_j
        polinomio = lagrange_poly(puntos, valores)
        return int_a_bytes(polinomio(np.arange(self.longitud_secreto)))

    def decodificacion(self, participaciones):
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (Identificador, Participación).
        Esta versión reconstruye el secreto a partir de la fórmula del polinomio interpolador de Lagrange evaluado en 0, ..., l-1.
        :param participaciones: Secuencia con las participaciones de los participantes que desean obtener el secreto.
        :return: El secreto.
        """
        # Verificación de condiciones
        r = self.reconstruccion
        if len(participaciones) < r:
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto')
        nombres, valores_b64 = zip(*participaciones[:r])
        self._verificar_nombres(nombres)

        # Obtener el elemento asociado a cada participante y decodificar su participación
        puntos = self.cuerpo(list(self.participantes_numero[nombre] for nombre in nombres))
        valores = self.cuerpo(b64str_a_int(valores_b64))
        # Calcular el valor del polinomio generador en 0, ..., l-1 sin reconstruirlo
        mascara = ~np.eye(self.reconstruccion, dtype=bool)  # Máscara de los elementos x_h de la fórmula
        coef = self.cuerpo.Zeros((self.longitud_secreto, self.reconstruccion))
        alphas = self.cuerpo.Range(0, self.longitud_secreto)[:, None]
        for i in range(self.reconstruccion):
            numerador = np.prod(alphas - puntos[mascara[i]], axis=1)  # Productorio del numerador aj - xh
            denominador = np.prod(puntos[i] - puntos[mascara[i]])  # Productorio del denominador xi - xh
            coef[:, i] = numerador / denominador  # Cálculo de l_i
        return int_a_bytes(np.sum(valores * coef, axis=1))  # Se devuelve la suma y_i * l_i(a_j)

    def _verificar_nombres(self, nombres):
        """
        Verifica que los participantes sean válidos, es decir, que no haya nombres duplicados y todos los nombres estén registrados como participantes.
        :param nombres: La secuencia de nombres que se quiere comprobar
        """
        # Comprobar no elementos duplicados
        if len(nombres) != len(set(nombres)):
            raise ValueError(f"Se han encontrado participantes duplicados.")
        # Comprobar que los participantes existen
        conjunto_nombres = self.participantes_numero
        for nombre in nombres:
            if nombre not in conjunto_nombres:
                raise ValueError(f"El participante {nombre} no está registrado.")