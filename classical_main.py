from classic_schemes import Shamir, Simplificado, ShamirRampa, McElieceSarwate
from galois import GF

from utils import ask_int

def programa_clasico():
    cuerpo = GF(2, 64)

    # Pedir los participantes y los datos del esquema
    participantes = []
    n = ask_int('Escriba el número de particiantes: ',
                f'El número de participantes debe ser al menos 2 y menor que el orden del cuerpo de trabajo ({cuerpo.characteristic}^{cuerpo.degree})',
                lambda x: 2 <= x < cuerpo.order)
    for i in range(1, n+1):
        participante = input(f'Escriba el nombre del participante nº{i}: ')
        while participante in participantes:
            participante = input(f'El participante introducido ya existe, por favor introduzca otro: ')
        participantes.append(participante)
    print() # Imprimir espacio en blanco para mejor claridad visual

    r = ask_int('Escriba número de participantes necesarios para recuperar el secreto: ',
                f'El número de participantes necesarios para recuperar el secreto debe ser al menos 2 y menor o igual que el número de participantes ({n})',
                lambda x: 2 <= x <= n)

    l = ask_int('Escriba el número de secretos que se desea compartir: ',
                f'El número de secretos debe ser al menos 1 y menor que el parámetro de reconstrucción ({r}).',
                lambda x: 1 <= x < r)

    # Preguntar por el esquema que se desea usar
    if l == 1:
        if r == n:
            yn = input('Se ha detectado que el número de participantes necesarios para reconstruir el esquema coincide con el número de participantes, ¿desea realizar el esquema simplificado? (y/n): ')
            if yn.lower() in ('si', 's', 'y', 'yes'):
                ss = Simplificado(cuerpo, participantes)
            else:
                ss = Shamir(cuerpo, r, participantes)
        else:
            ss = Shamir(cuerpo, r, participantes)
    else:
        print('¿Cual de los dos siguientes esquemas desea realizar?\n1.- Esquema de Shamir en rampa.\n2.- Esquema de McEliece-Sarwate.')
        esq = ask_int('Respuesta: ', f'No se ha introducido un numero válido.', lambda x: 1 <= x <= 2)
        if esq == 1:
            ss = ShamirRampa(cuerpo, r, l, participantes)
        else:
            if n+l > cuerpo.order:
                print('Debido a que el número de participantes es mayor que el orden del cuerpo menos la longitud del secreto, no se puede realizar el esquema de McEliece-Sarwate, se procede con el esquema de Shamir en rampa.')
                print()
                ss = ShamirRampa(cuerpo, r, l, participantes)
            else:
                ss = McElieceSarwate(cuerpo, r, l, participantes)

    # Preguntar por participaciones anticipadas
    yn = input('¿Desea repartir participaciones anticipadas? (y/n): ')
    if yn.lower() in ('si', 's', 'y', 'yes'):
        n_anticipados = ask_int(f'Introduzca el número de participaciones anticipadas (1 - {r - l}): ',
                                f'El número de participaciones anticipadas debe ser al menos 1 y menor que el parámetro de privacidad {r - l}',
                                lambda x: 1 <= x <= r - l)
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
        participaciones_anticipadas = ss.advance_sharing(participantes_anticipados)
        for nombre, ant_b64 in participaciones_anticipadas:
            print(f'{nombre}: {ant_b64}')

    # Pedir el secreto y codificarlo
    if l == 1:
        secreto = input(f'Escriba el secreto: ').encode()
        while len(secreto) > 8:
            secreto = input('Se ha introducido un secreto con longitud de bytes mayor que 8, introduzca otro: ').encode()
    else:
        secreto = []
        for i in range(1, l+1):
            sec = input(f'Escriba el secreto nº{i}: ').encode()
            while len(sec) > 8:
                sec = input(
                    'Se ha introducido un secreto con longitud de bytes mayor que 8, introduzca otro: ').encode()
            secreto.append(sec)

    # Crear participaciones
    participaciones = ss.codificacion(secreto)
    for nombre, part_b64 in participaciones:
        print(f'{nombre}: {part_b64}')

    # Reconstrucción del secreto
    print('Escriba el nombre de los participantes que buscan reconstruir el secreto y su participacion.')
    conjunto = []
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
        participacion = input('Participacion: ')
        conjunto.append((nombre, participacion))
        nombres.append(nombre)

    # Reconstruir el secreto
    secreto = ss.decodificacion(conjunto)
    try:
        if l == 1:
            print('Secreto:', secreto.decode())
        else:
            print('Secreto:', [s.decode() for s in secreto])
    except UnicodeDecodeError:
        print('El secreto obtenido no es válido')

if __name__ == '__main__':
    programa_clasico()