from base64 import b64decode, b64encode
from typing import Sequence

import numpy as np
from galois.typing import ArrayLike

def str_a_int(cadena: str | Sequence[str]) -> list[int] | int:
    """
    Obtiene la representacion en entero de un string.
    Si se proporciona una secuencia de strings, devuelve una lista de las representaciones de cada cadena.
    Primero codifica la cadena en 'utf-8' y posteriormente obtiene la representacion en entero de los bytes.
    :param cadena: La cadena de caracteres cuya representacion en enetros se quiere obtener.
    :return: El numero entero que representa a la cadena proporcionada.
    """
    if isinstance(cadena, str):
        return int.from_bytes(cadena.encode())
    return list(int.from_bytes(n.encode()) for n in cadena)

def int_a_str(numero: int | ArrayLike) -> list[str] | str:
    """
    Obtiene la representacion en string de un entero.
    Si se proporciona una secuancia de enteros, devuelve una lista de las representaciones de cada número.
    Primero obtiene la representacon en bytes del número y después decodifica los bytes en 'utf-8'.

    No siempre existe una representacion en string de un numero debido a la codificación en 'utf-8'.

    :param numero: El número entero cuya representación en string se quiere obtener.
    :return: La representacion en string del número entero.
    """
    if isinstance(numero, int):
        longitud = (numero.bit_length() + 7) // 8 if numero > 0 else 1
        return numero.to_bytes(longitud).decode()
    if isinstance(numero, np.ndarray):
        if numero.ndim == 0:
            n_int = int(numero)
            longitud = (n_int.bit_length() + 7) // 8 if n_int > 0 else 1
            return n_int.to_bytes(longitud).decode()
        representacion_str = []
        for n_int in numero.tolist():
            longitud = (n_int.bit_length() + 7) // 8 if n_int > 0 else 1
            representacion_str.append(n_int.to_bytes(longitud).decode())
        return representacion_str
    representacion_str = []
    for n in numero:
        n_int = int(n)
        longitud = (n_int.bit_length() + 7) // 8 if n_int > 0 else 1
        representacion_str.append(n_int.to_bytes(longitud).decode())
    return representacion_str


def gf_a_b64str(gf_array: ArrayLike, longitud: int) -> list[str]:
    return list(b64encode(int(numero).to_bytes(longitud)).decode() for numero in gf_array)

def int_a_b64str(lista_int: Sequence[int], longitud: int) -> list[str]:
    return list(b64encode(numero.to_bytes(longitud)).decode() for numero in lista_int)

def b64str_a_int(lista_str_b64: list[str]) -> list[int]:
    return list(int.from_bytes(b64decode(str_b64)) for str_b64 in lista_str_b64)

# def varificar_instancias(funcion):
#     nombres_argumentos = funcion.__code__.co_varnames
#
#     @wraps(funcion)
#     def verificado(*args, **kwargs):
#         for nombre, valor in zip(nombres_argumentos, args):
#             if not isinstance(argumento, tipos):
#                 raise TypeError(f'El argumento {nombre!r} debe ser una instancia de {tipos}, no {type(argumento)}.")