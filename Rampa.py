import secrets
from galois import GF, Poly, lagrange_poly

from Shamir import Shamir
from utils import *


class ShamirRampa:
    r"""
        Esquema de compartición de secretos de Shamir en rampa sobre el cuerpo $\mathbb{F}_{p^m}$.

        Ejemplo:
            Crea un esquema de Shamir en rampa con parámetros r = 4 y l = 3 sobre el cuerpo $\mathbb{F}_{3^5}$ para los participantes ['a', 'b', 'c', 'd', 'e', 'f'].

            .. ipython:: python

                gf = galois.GF(3**5)
                rsh = ShamirRampa(gf, 4, 3, ['a', 'b', 'c', 'd', 'e', 'f'])
        """

    def __init__(self, cuerpo: GF, r: int, l: int,  participantes: Sequence[str]):
        r"""
        Crea un esquema de compartición de secretos de Shamir en rampa sobre el cuerpo $\mathbb{F}_{p^m}$.
        :param cuerpo: El cuerpo finito sobre sobre el que el esquema está construido.
        :param r: Parámetro de reconstrucción del esquema (número mínimo de participantes necesarios para reconstruir el secreto).
        :param l: Longitud del secreto que se quiere repartir.
        :param participantes: Lista de los identificadores únicos de cada participante del esquema.
        """

        # Verificación de condiciones
        if cuerpo.order <= len(participantes) < 2:
            raise ValueError(f'El numero de participantes ({len(participantes)}) debe ser mayor que 1 y menor que el orden del cuerpo de trabajo ({cuerpo.order}).')
        if len(participantes) != len(set(participantes)):
            raise ValueError(f'Se han encontrado participantes duplicados.')
        if len(participantes) < r < 2:
            raise ValueError(f'El parámetro de reconstrucción ({r}) debe ser mayor que 1 y menor o igual que el número de participantes ({len(participantes)}).')
        if r <= l < 2:
            raise ValueError(f'La longitud del secreto ({l}) debe ser mayor que 1 y menor que el parámetro de reconstrucción ({r}).')

        self._cuerpo = cuerpo
        self._reconstruccion = r
        self._longitud_secreto = l
        self.__participaciones_anticipadas = None
        self._longitud_bytes = ((cuerpo.order - 1).bit_length() + 7) // 8
        self._participantes_nombre = [None]  # Array para pasar de numero -> nombre
        self._participantes_numero = {}  # Diccionario para pasar nombre -> numero
        for i, nombre in enumerate(participantes, 1):
            self._participantes_nombre.append(nombre)
            self._participantes_numero[nombre] = i

    def crear_anticipadas(self, participantes: Sequence[str]) -> list[tuple[str, str]]:
        """
        Crea participaciones anticipadas para cada participante especificado.
        El formato de las participaciones es: (Identificador, Participación).
        :param participantes: Listado de los participantes a entregar participaciones anticipadas.
        :return: Una lista que contienene las participaciones anticipadas asignadas a cada participante especificado.
        """

        # Verificación de condiciones
        if not 0 <= len(participantes) < self._reconstruccion - self._longitud_secreto:
            raise ValueError(f'El numero de participaciones anticipadas ({len(participantes)}) debe estar entre 0 y el parámetro de privacidad ({self._reconstruccion - self._longitud_secreto}).')
        self._verificar_nombres(participantes)

        # Generar las participaciones anticipadas, que son elementos aleatorios del cuerpo
        aleatoriedad = [secrets.randbelow(self._cuerpo.order) for _ in range(len(participantes))]
        aleatoriedad_b64 = int_a_b64str(aleatoriedad, self._longitud_bytes)
        self.__participaciones_anticipadas = list(zip(participantes, aleatoriedad_b64))
        return self.__participaciones_anticipadas

    def crear_participaciones(self, secreto: Sequence[Buffer]) -> list[tuple[str, str]]:
        """
        Crea las participaciones de todos los participantes de acuerdo al secreto recibido.
        El formato de las participaciones es: (Identificador, Participación).
        Si se han distribuido participaciones anticipadas, las participaciones serán coherentes con las mismas.
        :param secreto: Secreto que se quiere codificar entre todos los participantes.
        :return: Una lista que contienene las participaciones de cada participante que no ha participado en la distribución avanzada.
        """

        # Verificación de condiciones
        if len(secreto) != self._longitud_secreto:
            raise ValueError(f'Se esperaba un secreto de longitud {self._longitud_secreto}, pero se ha recibido {len(secreto)}.')

        secreto_i = bytes_a_int(secreto)

        # Si se han repartido participaciones anticipadas, se realiza compartición avanzada
        if self.__participaciones_anticipadas is not None:
            # Obtener el elemento asociado a cada participante y decodificar su participación
            nombres, valores_b64 = zip(*self.__participaciones_anticipadas)
            puntos = self._cuerpo(list(self._participantes_numero[nombre] for nombre in nombres) + [0])
            valores = self._cuerpo(b64str_a_int(valores_b64) + [secreto_i])
            x = np.setdiff1d(list(self._participantes_numero.values()), puntos)

            # Se determina un polinomio de grado r-1 compatible con las participaciones anticipadas
            lagrange = lagrange_poly(puntos, valores)
            if len(puntos) < self._reconstruccion - 1:  # Si el número de participaciones anticipadas es menor que r-1, hay que completar el polinomio con aleatoriedad
                polinomio = lagrange + Poly.Roots(puntos, field=self._cuerpo) * Poly.Random(
                    self._reconstruccion - len(puntos) - 1, field=self._cuerpo)
            else:  # Si no, el único polinomio disponible es el de Lagrange
                polinomio = lagrange
            del self.__participaciones_anticipadas  # Eliminación de las participaciones anticipadas para mayor seguridad

        # Si no, se sigue el proceso estándar
        else:
            polinomio = Poly([*self._cuerpo.Random(self._reconstruccion - 1), secreto_i], field=self._cuerpo)
            x = list(self._participantes_numero.values())

        # Generar el resto de las participaciones
        participaciones_b64 = int_a_b64str(polinomio(x), self._longitud_bytes)
        return list(zip((self._participantes_nombre[p] for p in x), participaciones_b64))

    def recuperar_secreto_v1(self, participaciones: Sequence[tuple[str, str]]) -> bytes:
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

    def recuperar_secreto_v2(self, participaciones: Sequence[tuple[str, str]]) -> bytes:
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (Identificador, Participación).
        Esta versión reconstruye el secreto a partir de la fórmula del polinomio interpolador de Lagrange evaluado en 0.
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

        # Calcular el valor del polinomio generador el el 0 sin reconstruirlo
        mascara = ~np.eye(self._reconstruccion, dtype=bool)  # Máscara de los elementos x_h de la fórmula
        coef = self._cuerpo.Zeros(self._reconstruccion)
        for j in range(self._reconstruccion):
            denominador = np.prod(puntos[mascara[j]] - puntos[j])  # Productorio del denominador
            numerador = np.prod(puntos[mascara[j]])  # Productorio del numinador
            coef[j] = numerador / denominador  # Cálculo de l_j
        return int_a_bytes(np.sum(valores * coef))  # Se devuelve la suma y_j * l_j

    recuperar_secreto = recuperar_secreto_v2  # Alias para recuperar secreto version 2

    def _verificar_nombres(self, nombres: Sequence[str]) -> None:
        """
        Verifica que los participantes sean válidos, es decir, que no haya nombres duplicados y todos los nombres estén registrados como participantes.
        :param nombres: La secuencia de nombres que se quiere comprobar
        """
        # Comprobar no elementos duplicados
        if len(nombres) != len(set(nombres)):
            raise ValueError(f"Se han encontrado participantes duplicados.")
        # Comprobar que los participantes existen
        conjunto_nombres = self._participantes_numero
        for nombre in nombres:
            if nombre not in conjunto_nombres:
                raise ValueError(f"El participante {nombre} no está registrado.")

