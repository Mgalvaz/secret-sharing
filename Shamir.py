from typing import Sequence
from warnings import deprecated

import galois
import numpy as np
from galois import GF, Poly
from galois.typing import ArrayLike


class Shamir:
    def __init__(self, cuerpo: GF, r: int) -> None:
        self.__participantes = None
        self._cuerpo = cuerpo
        self._reconstruccion = r
        self.__generador = np.random.default_rng()


    def crear_participantes(self, n: int):
        if self.__participantes is not None:
            raise AttributeError('Ya se han creado los participantes')
        if n >= self._cuerpo.order:
            raise ValueError(f'El numero de participantes ({n}) debe ser menor que el orden del cuerpo de trabajo ({self._cuerpo.order})')
        idx = np.random.default_rng().choice(range(1, self._cuerpo.order), size=n, replace=False)
        self.__participantes = self._cuerpo(idx)

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
        polinomio = Poly([*self._cuerpo.Random(self._reconstruccion-1), secreto], field=self._cuerpo)
        return zip(self.__participantes, polinomio(self.__participantes))

    @deprecated('Usar recuperar_cecreto_v2 en su lugar')
    def recuperar_secreto_v1(self, participaciones: Sequence[tuple[ArrayLike, ArrayLike]]) -> ArrayLike:
        if len(participaciones) < self._reconstruccion:
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto')
        puntos, valores = zip(*participaciones)
        puntos = self._cuerpo(puntos)
        valores = self._cuerpo(valores)
        polinomio = galois.lagrange_poly(puntos, valores)
        return polinomio(0)








if __name__ == '__main__':
    gf = GF(2**10)
    sh = Shamir(gf ,7)
    sh.crear_participantes(15)
    participantes = list(sh.crear_participaciones(gf(950)))
    print(sh.recuperar_secreto_v1(participantes))


