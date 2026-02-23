from base64 import b64decode, b64encode
from typing import Sequence
from galois.typing import ArrayLike

def str_a_int(numero_str: str | list[str]) -> list[int] | int:
    if isinstance(numero_str, list):
        return list(int.from_bytes(n.encode()) for n in numero_str)
    return int.from_bytes(numero_str.encode())

def int_a_str(numero: int | ArrayLike | list[int] , longitud: int) -> list[str] | str:
    if isinstance(numero, list):
        return list(n.to_bytes(longitud).lstrip(b'\x00').decode() for n in numero)
    elif isinstance(numero, int):
        return numero.to_bytes(longitud).lstrip(b'\x00').decode()
    else:
        if numero.ndim == 0:
            return int(numero).to_bytes(longitud).lstrip(b'\x00').decode()
    return list(n.to_bytes(longitud).lstrip(b'\x00').decode() for n in numero.tolist())

def gf_a_b64str(gf_array: ArrayLike, longitud: int) -> list[str]:
    return list(b64encode(int(numero).to_bytes(longitud)).decode() for numero in gf_array)

def int_a_b64str(lista_int: Sequence[int], longitud: int) -> list[str]:
    return list(b64encode(numero.to_bytes(longitud)).decode() for numero in lista_int)

def b64str_a_int(lista_str_b64: list[str]) -> list[int]:
    return list(int.from_bytes(b64decode(str_b64.encode())) for str_b64 in lista_str_b64)

# def varificar_instancias(funcion):
#     nombres_argumentos = funcion.__code__.co_varnames
#
#     @wraps(funcion)
#     def verificado(*args, **kwargs):
#         for nombre, valor in zip(nombres_argumentos, args):
#             if not isinstance(argumento, tipos):
#                 raise TypeError(f'El argumento {nombre!r} debe ser una instancia de {tipos}, no {type(argumento)}.")