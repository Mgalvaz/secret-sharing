from typing import Sequence
import secrets
from base64 import b64encode, b64decode

import galois
import numpy as np
from galois import GF, Poly
from galois.typing import ArrayLike

from utils import int_a_b64str, b64str_a_int, gf_a_b64str,int_a_str


class Shamir:
    def __init__(self, cuerpo: GF, r: int, participantes: ArrayLike):
        if n >= cuerpo.order:
            raise ValueError(f'El numero de participantes ({len(participantes)}) debe ser menor que el orden del cuerpo de trabajo ({cuerpo.order}).')
        if r > n:
            raise ValueError(f'El parámetro de reconstrucción ({r}) debe ser menor o igual que el número de participantes ({len(participantes)}).')
        self._participantes = participantes
        self._cuerpo = cuerpo
        self._reconstruccion = r
        self.__participaciones_anticipadas = None

        """# Elegir elementos aleatorios para cada participante
        self.__participantes = []
        for _ in range(n):
            participante = cuerpo(secrets.randbelow(cuerpo.order))
            while participante in self.__participantes or participante == 0:
                participante = cuerpo(secrets.randbelow(cuerpo.order))
            self.__participantes.append(participante)"""



    def crear_anticipadas(self, anticipadas: Sequence[str]) -> list[tuple[str, str]]:
        if  not  0 <= len(anticipadas) < self._reconstruccion:
            raise ValueError(f'El numero de participaciones anticipadas ({len(anticipadas)}) debe estar entre 0 y el parámetro de privacidad ({self._reconstruccion - 1}).')
        aleatoriedad = []
        for _ in range(len(anticipadas)):
            aleatoriedad.append(secrets.randbelow(self._cuerpo.order))
        aleatoriedad_b64 = int_a_b64str(aleatoriedad, 8)
        self.__participaciones_anticipadas = list(zip(anticipadas, aleatoriedad_b64))
        #self.__participaciones_anticipadas = list(zip(self.__participantes[:anticipadas], self._cuerpo.Random(anticipadas)))
        return self.__participaciones_anticipadas

    def crear_participaciones(self, secreto: str) -> list[tuple[Sequence[str], str]]:
        secreto_i = int.from_bytes(secreto.encode())
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
            #x = self.__participantes[len(puntos) - 1:]
            del self.__participaciones_anticipadas
        else:
            polinomio = Poly([*self._cuerpo.Random(self._reconstruccion-1), secreto_i], field=self._cuerpo)
            x = np.arange(1, len(self._participantes) + 1)
            #x = self.__participantes
        participaciones_b64 = gf_a_b64str(polinomio(x), 8)
        return list(zip(self._cuerpo(x), participaciones_b64))

    def recuperar_secreto_v1(self, participaciones: Sequence[tuple[str, str]]) -> str:
        if len(participaciones) < self._reconstruccion:
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto')
        nombres, valores_b64 = zip(*participaciones[:self._reconstruccion])
        puntos = self._cuerpo(list(self._participantes.index(nombre) + 1 for nombre in nombres))
        valores = self._cuerpo(b64str_a_int(valores_b64))
        polinomio = galois.lagrange_poly(puntos, valores)
        return int_a_str(polinomio(0), 8)

    def recuperar_secreto_v2(self, participaciones: Sequence[tuple[str, str]]) -> str:
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
    gf = GF(2 ** 64)

    print(type(gf([4238194])))
    print(int_a_str(gf([461278346283174, 1]), 8))

    participantes = []
    r = int(input('Escriba el parámetro de reconstrucción (nº minimo de participantes para recuperar el secreto): '))
    n = int(input('Escriba el número de particiantes: '))
    for i in range(n):
        part = input(f'Escriba el nombre del participante nº{i+1}: ')
        participantes.append(part)

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

    secreto = sh.recuperar_secreto_v1(conjunto)
    print('secreto:', secreto)


