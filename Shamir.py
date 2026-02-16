from dataclasses import field
from typing import Sequence
import time
import secrets

import galois
import numpy as np
from galois import GF, Poly
from galois.typing import ArrayLike


class Shamir:
    def __init__(self, cuerpo: GF, r: int, n: int) -> None:
        if n >= cuerpo.order:
            raise ValueError(f'El numero de participantes ({n}) debe ser menor que el orden del cuerpo de trabajo ({cuerpo.order}).')
        if r > n:
            raise ValueError(f'El parámetro de reconstrucción ({r}) debe ser menor o igual que el número de participantes ({n}).')
        self._num_participantes = n
        self._cuerpo = cuerpo
        self._reconstruccion = r
        self.__participaciones_anticipadas = None

        # Elegir elementos aleatorios para cada participante
        self.__participantes = []
        for _ in range(n):
            participante = cuerpo(secrets.randbelow(cuerpo.order))
            while participante in self.__participantes or participante == 0:
                participante = cuerpo(secrets.randbelow(cuerpo.order))
            self.__participantes.append(participante)



    def crear_anticipadas(self, anticipadas: int):
        if  not  0 <= anticipadas < self._reconstruccion:
            raise ValueError(f'El numero de participaciones anticipadas ({anticipadas}) debe estar entre 0 y el parámetro de privacidad ({self._reconstruccion - 1}).')
        #self.__participaciones_anticipadas = list(zip(gf(np.arange(1, anticipadas + 1)), self._cuerpo.Random(anticipadas)))
        self.__participaciones_anticipadas = list(zip(self.__participantes[:anticipadas], self._cuerpo.Random(anticipadas)))
        return self.__participaciones_anticipadas

    def crear_participaciones(self, secreto: ArrayLike):
        if self._num_participantes is None:
            raise AttributeError('No se han definido ningun participante')
        if self.__participaciones_anticipadas is not None:
            puntos, valores = zip(*self.__participaciones_anticipadas)
            puntos = np.append(puntos, self._cuerpo(0))
            valores = np.append(valores, secreto)
            lagrange = galois.lagrange_poly(puntos, valores)
            if lagrange.degree < self._reconstruccion - 1:
                polinomio = lagrange + Poly.Roots(puntos, field=self._cuerpo)*Poly.Random(self._reconstruccion - len(puntos) - 1, field=self._cuerpo)
            else:
                polinomio = lagrange
            #x = np.arange(len(puntos), self._num_participantes + 1)
            x = self.__participantes[len(puntos) - 1:]
            del self.__participaciones_anticipadas
        else:
            polinomio = Poly([*self._cuerpo.Random(self._reconstruccion-1), secreto], field=self._cuerpo)
            #x = np.arange(1, self._num_participantes + 1)
            x = self.__participantes
        return list(zip(self._cuerpo(x), polinomio(x)))

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
        coef = self._cuerpo.Zeros(self._reconstruccion)
        for j in range(self._reconstruccion):
            denominador = np.prod(puntos[mascara[j]] - puntos[j])
            numerador = np.prod(puntos[mascara[j]])
            coef[j] = numerador / denominador
        return np.sum(valores * self._cuerpo(coef))

if __name__ == '__main__':
    secreto_b = input('Escriba el secreto: ').encode()
    if len(secreto_b) > 8:
        print(f'Warning: se han mandado más de 8 bytes, se procede a truncar los primeros 8 bytes')
        secreto_b = secreto_b[:8]
    secreto = int.from_bytes(secreto_b)

    gf = GF(2**64)
    sh = Shamir(gf ,10, 30)
    anticipadas = sh.crear_anticipadas(7)
    participantes = sh.crear_participaciones(gf(secreto))

    conjunto = anticipadas[2:6] + participantes[15:23]
    sec2 = sh.recuperar_secreto_v2(conjunto)

    secreto_b = int(sec2).to_bytes(8).lstrip(b'\x00')
    print('secreto:', secreto_b.decode())


