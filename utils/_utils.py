from functools import wraps

# def varificar_instancias(funcion):
#     nombres_argumentos = funcion.__code__.co_varnames
#
#     @wraps(funcion)
#     def verificado(*args, **kwargs):
#         for nombre, valor in zip(nombres_argumentos, args):
#             if not isinstance(argumento, tipos):
#                 raise TypeError(f'El argumento {nombre!r} debe ser una instancia de {tipos}, no {type(argumento)}.")