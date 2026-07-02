from esquemas_cuanticos import CGL, Ogawa, ZhangMatsumoto
from galois import GF
from qiskit.quantum_info import Statevector
from numpy.linalg import norm
import ast

from utils import pedir_entero

def pedir_secreto(dimension, n=-1):
    while True:
        if n == -1:
            secreto_str = input(f'Escriba el secreto en forma de lista: ')
        else:
            secreto_str = input(f'Escriba el secreto nº{n} en forma de lista: ')
        try:
            secreto = Statevector(ast.literal_eval(secreto_str))
        except (ValueError, SyntaxError):
            print('No se ha introducido una lista válida de números, introduzca otra.')
            continue
        if len(secreto) != dimension:
            print(
                f'Se esperaba que el secreto tuviera dimensión {dimension}, pero se ha recibido {len(secreto)}, introduzca otro.')
        elif not secreto.is_valid():
            print('El secreto obtenido no es válido, ¿quiere introducir otro o normalizarlo?\n1.- Introducir otro.\n2.- Normalizar.')
            op = pedir_entero(f'Respuesta: ', f'No se ha introducido un numero válido.', lambda x: 1 <= x <= 2)
            if op == 2:
                norma = norm(secreto)
                if norma == 0:
                    print('Se ha introducido un vector de ceros, introduzca otra lista.')
                else:
                    secreto = secreto / norma
                    return secreto
        else:
            return secreto

def programa_cuantico():
    # Pedir datos del esquema
    cuerpo =  GF(2, 3)

    l = pedir_entero('Escriba el número de secretos que se desea compartir: ',
                     f'El número de secretos debe ser al menos 1.', lambda x: 1 <= x)

    r = pedir_entero('Escriba número de participantes necesarios para recuperar el secreto: ',
                     f'El número de participantes necesarios para recuperar el secreto debe ser mayor que la longitud del secreto y menor que el orden del cuerpo de trabajo menos la longitud del secreto entre 2 ({(cuerpo.order + l + 1) // 2})', lambda x: l < x < (cuerpo.order + l + 1) // 2)

    # Pedir participantes
    n = pedir_entero('Escriba el número de particiantes: ',
                     f'El número de participantes debe ser como mínimo el parámetro de recostrucción ({r}) y no puede superar ({2 * r - l}).', lambda x: r <= x <= 2 * r - l)
    participantes = []
    for i in range(1, n+1):
        participante = input(f'Escriba el nombre del participante nº{i}: ')
        while participante in participantes:
            participante = input(f'El participante introducido ya existe, por favor introduzca otro: ')
        participantes.append(participante)

    # Preguntar por el esquema que se desea usar
    if l == 1:
        ss = CGL(cuerpo, r, participantes)
    else:
        print('¿Cual de los dos siguientes esquemas desea realizar?\n1.- Esquema de Ogawa et al.\n2.- Esquema de Zhang-Matsumoto.')
        esq = pedir_entero('Respuesta: ',f'No se ha introducido un numero válido.', lambda x: 1 <= x <= 2)
        if esq == 1:
            ss = Ogawa(cuerpo, r, l, participantes)
        else:
            if 2*r > cuerpo.order:
                print('Debido a que el número de participantes totales es mayor que el orden del cuerpo menos la longitud del secreto, no se puede realizar el esquema de Zhang-Matsumoto, se procede con el esquema de Ogawa et al.')
                print()
                ss = Ogawa(cuerpo, r, l, participantes)
            else:
                ss = ZhangMatsumoto(cuerpo, r, l, participantes)

    # Preguntar por participaciones anticipadas
    yn = input('¿Desea repartir participaciones anticipadas? (y/n): ')
    diccionario_participaciones = {}
    if yn.lower() in ('si', 's', 'y', 'yes'):
        n_anticipados = pedir_entero(f'Introduzca el número de participaciones anticipadas (1 - {r-l}): ',
                     f'El número de participaciones anticipadas debe ser al menos 1 y menor que el parámetro de privacidad {r-l}', lambda x: 1 <= x <= r-l)
        participantes_anticipados = []
        # Añadir participantes anticipados
        for i in range(1, n_anticipados+1):
            part_anticipado = input(f'Escriba el nombre del participante anticipado nº{i}: ')
            # Comporbar validez de los participantes anticipados
            while True:
                if part_anticipado in participantes_anticipados:
                    print('Se ha introducido un participante duplicado.')
                elif part_anticipado not in participantes:
                    print('Se ha introducido un participante no registrado.')
                else:
                    break
                part_anticipado = input(f'Escriba el nombre del participante anticipado nº{i}: ')
            participantes_anticipados.append(part_anticipado)
        # Crear particpaciones anticipadas
        participaciones_anticipadas = ss.comparticion_anticipada(participantes_anticipados)
        for registro in participaciones_anticipadas:
            print(f'{registro.name}: {registro}')
            diccionario_participaciones[registro.name] = registro

    # Pedir el secreto y codificarlo
    if l == 1:
        secreto = pedir_secreto(cuerpo.order)
    else:
        print('¿Desea introducir el vector de estado compuesto o cada uno de los subsistemas?\n1.- Vector de estado global.\n2.- Separar por subsistemas.')
        elec = pedir_entero(f'Respuesta: ', f'No se ha introducido un numero válido.', lambda x: 1 <= x <= 2)
        if elec == 1:
            secreto = pedir_secreto(cuerpo.order**l)
        else:
            secreto = Statevector([1]) # Inicializamos el secreto
            for i in range(1,l+1):
                sec_i = pedir_secreto(cuerpo.order, i)
                secreto = secreto.tensor(sec_i)

    # Crear participaciones
    participaciones = ss.codificacion(secreto)
    for registro in participaciones:
        print(f'{registro.name}: {registro}')
        diccionario_participaciones[registro.name] = registro

    # Reconstrucción del secreto
    print('Escriba el nombre de los participantes que buscan reconstruir el secreto.')
    registros = []
    nombres = []
    # Pedir las participaciones
    for _ in range(r):
        nombre = input('Nombre: ')
        # Comprobar validez de participantes
        while True:
            if nombre in nombres:
                print('Se ha introducido un participante duplicado.')
            elif nombre not in participantes:
                print('Se ha introducido un participante no registrado.')
            else:
                break
            nombre = input('Nombre: ')
        registros.append(diccionario_participaciones[nombre])
        nombres.append(nombre)

    # Reconstruir el secreto
    secreto = ss.decodificacion(registros)
    secreto = Statevector(secreto, dims=tuple(cuerpo.order for _ in range(l)))
    print('Secreto:', secreto.to_dict())

if __name__ == '__main__':
    programa_cuantico()