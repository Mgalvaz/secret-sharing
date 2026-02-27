import secrets
import galois
import numpy as np
from galois import GF, Poly
from collections.abc import Buffer

from utils import *


class Shamir:
    r"""
    Esquema de compartición de secretos de Shamir sobre el cuerpo $\mathbb{F}_{p^m}$.

    Ejemplo:
        Crea un esquema de Shamir de umbral-(4,6) sobre el cuerpo $\mathbb{F}_{3^5}$ para los participantes ['a', 'b', 'c', 'd', 'e', 'f'].

        .. ipython:: python

            gf = galois.GF(3**5)
            sh = Shamir(fg, 4, ['a', 'b', 'c', 'd', 'e', 'f'])
    """
    def __init__(self, cuerpo: GF, r: int, participantes: Sequence[str]):
        r"""
        Crea un esquema de compartición de secretos de Shamir sobre el cuerpo $\mathbb{F}_{p^m}$.
        :param cuerpo: El cuerpo finito sobre sobre el que el esquema está construido.
        :param r: Parámetro de reconstrucción del esquema (número mínimo de participantes necesarios para reconstruir el secreto).
        :param participantes: Lista de los identificadores únicos de cada participante del esquema.
        """
        if len(participantes) >= cuerpo.order:
            raise ValueError(f'El numero de participantes ({len(participantes)}) debe ser menor que el orden del cuerpo de trabajo ({cuerpo.order}).')
        if r > len(participantes):
            raise ValueError(f'El parámetro de reconstrucción ({r}) debe ser menor o igual que el número de participantes ({len(participantes)}).')
        self._participantes = participantes
        self._cuerpo = cuerpo
        self._reconstruccion = r
        self.__participaciones_anticipadas = None
        self._longitud_bytes = ((cuerpo.order - 1).bit_length() + 7) // 8

    def crear_anticipadas(self, participantes: Sequence[str]) -> list[tuple[str, str]]:
        """
        Crea participaciones anticipadas para cada participante especificado.
        El formato de las participaciones es: (Identificador, Participación).
        :param participantes: Listado de los participantes a entregar participaciones anticipadas.
        :return: Una lista que contienene las participaciones anticipadas asignadas a cada participante especificado.
        """
        if  not  0 <= len(participantes) < self._reconstruccion:
            raise ValueError(f'El numero de participaciones anticipadas ({len(participantes)}) debe estar entre 0 y el parámetro de privacidad ({self._reconstruccion - 1}).')
        aleatoriedad = []
        for _ in range(len(participantes)):
            aleatoriedad.append(secrets.randbelow(self._cuerpo.order))
        aleatoriedad_b64 = int_a_b64str(aleatoriedad, self._longitud_bytes)
        self.__participaciones_anticipadas = list(zip(participantes, aleatoriedad_b64))
        return self.__participaciones_anticipadas

    def crear_participaciones(self, secreto: Buffer) -> list[tuple[str, str]]:
        """
        Crea las participaciones de todos los participantes de acuerdo al secreto recibido.
        El formato de las participaciones es: (Identificador, Participación).
        Si se han distribuido participaciones anticipadas, las participaciones serán coherentes con las mismas.
        :param secreto: Secreto que se quiere codificar entre todos los participantes.
        :return: Una lista que contienene las participaciones de cada participante que no ha participado en la distribución avanzada.
        """
        secreto_i = bytes_a_int(secreto)
        if self.__participaciones_anticipadas is not None:
            nombres, valores_b64 = zip(*self.__participaciones_anticipadas)
            puntos = self._cuerpo(list(self._participantes.index(nombre) + 1 for nombre in nombres) + [0])
            valores = self._cuerpo(b64str_a_int(valores_b64) + [secreto_i])
            lagrange = galois.lagrange_poly(puntos, valores)
            if lagrange.degree < self._reconstruccion - 1:
                polinomio = lagrange + Poly.Roots(puntos, field=self._cuerpo)*Poly.Random(self._reconstruccion - len(puntos) - 1, field=self._cuerpo)
            else:
                polinomio = lagrange
            x = np.setdiff1d(np.arange(1, len(self._participantes) + 1), puntos)
            del self.__participaciones_anticipadas
        else:
            polinomio = Poly([*self._cuerpo.Random(self._reconstruccion-1), secreto_i], field=self._cuerpo)
            x = np.arange(1, len(self._participantes) + 1)
        participaciones_b64 = int_a_b64str(polinomio(x), self._longitud_bytes)
        return list(zip((self._participantes[p-1] for p in x), participaciones_b64))

    def recuperar_secreto_v1(self, participaciones: Sequence[tuple[str, str]]) -> bytes:
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (Identificador, Participación).
        Esta versión reconstruye primero el polinomio generador y a partir de él, devuelve el secreto.
        :param participaciones: Secuencia con las participaciones de los participantes que desean obtener el secreto.
        :return: El secreto.
        """
        if len(participaciones) < self._reconstruccion:
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto')
        nombres, valores_b64 = zip(*participaciones[:self._reconstruccion])
        puntos = self._cuerpo(list(self._participantes.index(nombre) + 1 for nombre in nombres))
        valores = self._cuerpo(b64str_a_int(valores_b64))
        polinomio = galois.lagrange_poly(puntos, valores)
        return int_a_bytes(polinomio.coefficients(order="asc")[0])

    def recuperar_secreto_v2(self, participaciones: Sequence[tuple[str, str]]) -> bytes:
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (Identificador, Participación).
        Esta versión reconstruye el secreto a partir de la fórmula del polinomio interpolador de Lagrange evaluado en 0.
        :param participaciones: Secuencia con las participaciones de los participantes que desean obtener el secreto.
        :return: El secreto.
        """
        if len(participaciones) < self._reconstruccion:
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto')
        nombres, valores_b64 = zip(*participaciones[:self._reconstruccion])
        puntos = self._cuerpo(list(self._participantes.index(nombre) + 1 for nombre in nombres))
        valores = self._cuerpo(b64str_a_int(valores_b64))
        mascara = ~np.eye(self._reconstruccion, dtype=bool)
        coef = self._cuerpo.Zeros(self._reconstruccion)
        for j in range(self._reconstruccion):
            denominador = np.prod(puntos[mascara[j]] - puntos[j])
            numerador = np.prod(puntos[mascara[j]])
            coef[j] = numerador / denominador
        return int_a_bytes(np.sum(valores * coef))

    recuperar_secreto = recuperar_secreto_v2 # Alias para recuperar secreto version 2

class ShamirSimplificado:
    r"""
    Esquema de compartición de secretos de Shamir simplificado sobre cuerpo $\mathbb{F}_{p^m}$.

    Ejemplo:
        Crea un esquema de umbral-(5,5) sobre el cuerpo $\mathbb{F}_{5^4}$ para los participantes ['a', 'b', 'c', 'd', 'e'].

        .. ipython:: python

            gf = galois.GF(5**4)
            sh = ShamirSimplificado(fg, ['a', 'b', 'c', 'd', 'e'])
    """
    def __init__(self, cuerpo: GF, participantes: Sequence[str]):
        r"""
        Crea un esquema de compartición de secretos de Shamir simlificado sobre el cuerpo $\mathbb{F}_{p^m}$.
        :param cuerpo: El cuerpo finito sobre sobre el que el esquema está construido.
        :param participantes: Lista de los identificadores únicos de cada participante del esquema.
        """
        if len(participantes) >= cuerpo.order:
            raise ValueError(f'El numero de participantes ({len(participantes)}) debe ser menor que el orden del cuerpo de trabajo ({cuerpo.order}).')
        self._participantes = participantes
        self._cuerpo = cuerpo
        self.__participaciones_anticipadas = None
        self._longitud_bytes = ((cuerpo.order - 1).bit_length() + 7) // 8

    def crear_anticipadas(self, participantes: Sequence[str]) -> list[tuple[str, str]]:
        """
        Crea participaciones anticipadas para cada participante especificado.
        El formato de las participaciones es: (Identificador, Participación).
        :param participantes: Listado de los participantes a entregar participaciones anticipadas.
        :return: Una lista que contienene las participaciones anticipadas asignadas a cada participante especificado.
        """
        if  not 0 <= len(participantes) < len(self._participantes):
            raise ValueError(f'El numero de participaciones anticipadas ({len(participantes)}) debe estar entre 0 y el parámetro de privacidad ({len(self._participantes) - 1}).')
        aleatoriedad = []
        for _ in range(len(participantes)):
            aleatoriedad.append(secrets.randbelow(self._cuerpo.order))
        aleatoriedad_b64 = int_a_b64str(aleatoriedad, self._longitud_bytes)
        self.__participaciones_anticipadas = list(zip(participantes, aleatoriedad_b64))
        return self.__participaciones_anticipadas

    def crear_participaciones(self, secreto: Buffer) -> list[tuple[str, str]]:
        """
        Crea las participaciones de todos los participantes de acuerdo al secreto recibido.
        El formato de las participaciones es: (Identificador, Participación).
        Si se han distribuido participaciones anticipadas, las participaciones serán coherentes con las mismas.
        :param secreto: Secreto que se quiere codificar entre todos los participantes.
        :return: Una lista que contienene las participaciones de cada participante que no ha participado en la distribución avanzada.
        """
        secreto_i = bytes_a_int(secreto)
        if self.__participaciones_anticipadas is not None:
            nombres, valores_b64 = zip(*self.__participaciones_anticipadas)
            valores = self._cuerpo(b64str_a_int(valores_b64))
            x = np.setdiff1d(self._participantes, nombres)
            participaciones = self._cuerpo.Zeros(len(x))
            for i in range(len(x)):
                participaciones[i] = secrets.randbelow(self._cuerpo.order)
            participaciones[-1] = self._cuerpo(secreto_i) - participaciones[:-1].sum() - valores.sum()
            del self.__participaciones_anticipadas
        else:
            x = np.arange(len(self._participantes))
            participaciones = self._cuerpo.Zeros(len(x))
            for i, punto in enumerate(x[:-1]):
                participaciones[i] = secrets.randbelow(self._cuerpo.order)
            participaciones[-1] = self._cuerpo(secreto_i) - participaciones[:-1].sum()
        participaciones_b64 = int_a_b64str(participaciones, self._longitud_bytes)
        return list(zip((self._participantes[p] for p in x), participaciones_b64))

    def recuperar_secreto(self, participaciones: Sequence[tuple[str, str]]) -> bytes:
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (Identificador, Participación).
        :param participaciones: Secuencia con las participaciones de los participantes que desean obtener el secreto.
        :return: El secreto.
        """
        if len(participaciones) < len(self._participantes):
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto')
        _, valores_b64 = zip(*participaciones[:len(self._participantes)])
        valores = self._cuerpo(b64str_a_int(valores_b64))
        return int_a_str(valores.sum(), self._longitud_bytes)

if __name__ == '__main__':
    gf = GF(2 ** 64)


    a = str_a_int([input(i) for i in range(5)])
    print(int_a_str(np.array(a)))


    exit()


    participantes = []
    r = int(input('Escriba el parámetro de reconstrucción (nº minimo de participantes para recuperar el secreto): '))
    n = int(input('Escriba el número de particiantes: '))
    for i in range(n):
        part = input(f'Escriba el nombre del participante nº{i+1}: ')
        participantes.append(part)
    if r == n:
        sh = ShamirSimplificado(gf, participantes)
    else:
        sh = Shamir(gf, r, participantes)

    yn = input('Desea repartir participación anticipada? (y/n): ')
    if yn.lower() in ('si', 's', 'y', 'yes'):
        n_anticipados = int(input('Introduzca el número de participaciones anticipadas: '))
        anticipadas = []
        for i in range(n_anticipados):
            part_ant = input(f'Escriba el nombre del participante anticipado nº{i+1}: ')
            anticipadas.append(part_ant)
        participaciones_anticipadas = sh.crear_anticipadas(anticipadas)
        for nombre, ant_b64 in participaciones_anticipadas:
            print(f'{nombre}: {ant_b64}')


    secreto = input('Escriba el secreto: ')
    participantes = sh.crear_participaciones(secreto)
    for nombre, part_b64 in participantes:
        print(f'{nombre}: {part_b64}')

    print('Escriba el nombre de los participantes que busca reconstruir el secreto y su participacion.')
    conjunto = []
    for _ in range(r):
        nombre = input('Nombre: ')
        participacion = input('Participacion: ')
        conjunto.append((nombre, participacion))

    secreto = sh.recuperar_secreto(conjunto)
    print('secreto:', secreto)


