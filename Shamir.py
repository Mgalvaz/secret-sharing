from typing import Sequence
import time

import galois
import numpy as np
from galois import GF, Poly
from galois.typing import ArrayLike


class Shamir:
    def __init__(self, cuerpo: GF, r: int) -> None:
        self.__participantes = None
        self.__participaciones_anticipadas = None
        self._cuerpo = cuerpo
        self._reconstruccion = r


    def crear_participantes(self, n: int, *, anticipadas: int = 0):
        if self.__participantes is not None:
            raise AttributeError('Ya se han creado los participantes.')
        if n >= self._cuerpo.order:
            raise ValueError(f'El numero de participantes ({n}) debe ser menor que el orden del cuerpo de trabajo ({self._cuerpo.order}).')
        if  not  0 <= anticipadas < self._reconstruccion:
            raise ValueError(f'El numero de participaciones anticipadas ({anticipadas}) debe estar entre 0 y el parámetro de privacidad ({self._reconstruccion - 1}).')
        idx = np.arange(1, n+1) # idx = np.random.default_rng().choice(range(1, self._cuerpo.order), size=n, replace=False)
        self.__participantes = self._cuerpo(idx)
        if anticipadas:
            self.__participaciones_anticipadas = zip(self.__participantes[0:anticipadas], self._cuerpo.Random(anticipadas))
            self.__participantes = self.__participantes[anticipadas:]

    """def anadir_participantes(self, n: int):
        if self.__participantes is None:
            self.crear_participantes(n)
        else:
            if len(self.__participantes) + n >= self._cuerpo.order:
                raise ValueError(f'El numero de participantes ({len(self.__participantes) + n}) debe ser menor que el orden del cuerpo de trabajo ({self._cuerpo.order})')
            self.__participantes = np.append(self.__participantes, self._cuerpo.Random(n))"""


    def crear_participaciones(self, secreto: ArrayLike):
        if self.__participantes is None:
            raise AttributeError('No se han definido ningun participante')
        if self.__participaciones_anticipadas:
            puntos, valores = zip(*self.__participaciones_anticipadas)
            puntos = np.append(puntos, self._cuerpo(0))
            valores = np.append(valores, secreto)
            polinomio = galois.lagrange_poly(puntos, valores)
            del self.__participaciones_anticipadas

        else:
            polinomio = Poly([*self._cuerpo.Random(self._reconstruccion-1), secreto], field=self._cuerpo)
        return zip(self.__participantes, polinomio(self.__participantes))

    def recuperar_secreto_v1(self, participaciones: Sequence[tuple[ArrayLike, ArrayLike]]) -> ArrayLike:
        if len(participaciones) < self._reconstruccion:
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto')
        puntos, valores = zip(*participaciones[:self._reconstruccion])
        puntos = self._cuerpo(puntos)
        valores = self._cuerpo(valores)
        polinomio = galois.lagrange_poly(puntos, valores)
        return polinomio(0)

    def recuperar_secreto_v2(self, participaciones: Sequence[tuple[ArrayLike, ArrayLike]]) -> ArrayLike:
        if len(participaciones) < self._reconstruccion:
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto')
        puntos, valores = zip(*participaciones[:self._reconstruccion])
        puntos = self._cuerpo(puntos)
        valores = self._cuerpo(valores)
        mascara = ~np.eye(self._reconstruccion, dtype=bool)
        coef = np.empty(self._reconstruccion, dtype=int)
        for j in range(self._reconstruccion):
            denominador = np.prod(puntos[mascara[j]] - puntos[j])
            numerador = np.prod(puntos[mascara[j]])
            coef[j] = numerador / denominador
        return np.sum(valores * self._cuerpo(coef))

if __name__ == '__main__':
    secreto = bytes(input('Escriba el secreto: '), 'utf-8')
    if len(secreto) > 8:
        print(f'Warning: se han mandado más de 8 bytes, se procede a truncar')
        secreto = secreto[:8]

    gf = GF(2**64)
    sh = Shamir(gf ,3)
    sh.crear_participantes(15)
    participantes = list(sh.crear_participaciones(gf(4782391)))
    sec2 = sh.recuperar_secreto_v2(participantes)
    print('secreto', sec2)


