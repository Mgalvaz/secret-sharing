from base64 import b64decode, b64encode
from collections.abc import Buffer
from typing import Sequence

import numpy as np
from galois.typing import ArrayLike

def bytes_a_int(cadena: Buffer | list[Buffer]) -> int | list[int]:
    """
    Obtiene la representacion en entero de un string.
    Si se proporciona una secuencia de strings, devuelve una lista de las representaciones de cada cadena.
    :param cadena: La cadena de caracteres cuya representacion en enetros se quiere obtener.
    :return: El numero entero que representa a la cadena proporcionada.
    """
    if isinstance(cadena, Buffer):
        return int.from_bytes(cadena)
    return list(int.from_bytes(n) for n in cadena)

def int_a_bytes(numero: int | ArrayLike) -> bytes | list[bytes]:
    """
    Obtiene la representacion en bytes de un entero.
    Si se proporciona una secuancia de enteros, devuelve una lista de las representaciones de cada número.
    :param numero: El número entero cuya representación en string se quiere obtener.
    :return: La representacion en string del número entero.
    """
    if isinstance(numero, int):
        longitud = (numero.bit_length() + 7) // 8 if numero > 0 else 1
        return numero.to_bytes(longitud)
    if isinstance(numero, np.ndarray):
        if numero.ndim == 0:
            n_int = int(numero)
            longitud = (n_int.bit_length() + 7) // 8 if n_int > 0 else 1
            return n_int.to_bytes(longitud)
        representacion_str = []
        for n_int in numero.tolist():
            longitud = (n_int.bit_length() + 7) // 8 if n_int > 0 else 1
            representacion_str.append(n_int.to_bytes(longitud))
        return representacion_str
    representacion_str = []
    for n in numero:
        n_int = int(n)
        longitud = (n_int.bit_length() + 7) // 8 if n_int > 0 else 1
        representacion_str.append(n_int.to_bytes(longitud))
    return representacion_str

def int_a_b64str(lista_int: ArrayLike, longitud: int) -> list[str]:
    """
    Obtiene la codificación en base64 de cada número de un array de enteros.
    :param lista_int: El array de enteros.
    :param longitud: El número de bytes que se desean usar para representar cada número.
        Se crea un OverflowError si el número no se puede representar con el la longitud dada.
    :return: La lista con la codificación de cada entero.
    """
    if isinstance(lista_int, np.ndarray):
        return list(b64encode(numero.to_bytes(longitud)).decode() for numero in lista_int.tolist())
    return list(b64encode(numero.to_bytes(longitud)).decode() for numero in lista_int)

def b64str_a_int(lista_b64: list[str]) -> list[int]:
    """
    Obtiene la decodificación en base64 de cada string de una lista.
    :param lista_b64: La lista de strings codificados.
    :return: La lista con la decodificación de cada string.
    """
    return list(int.from_bytes(b64decode(str_b64)) for str_b64 in lista_b64)