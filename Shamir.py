import secrets
from galois import GF, Poly, lagrange_poly

from utils import *


class Shamir:
    r"""
    Esquema de compartición de secretos de Shamir sobre el cuerpo $\mathbb{F}_{p^m}$.

    Ejemplo:
        Crea un esquema de Shamir de umbral-(4,6) sobre el cuerpo $\mathbb{F}_{3^5}$ para los participantes ['a', 'b', 'c', 'd', 'e', 'f'].

        .. ipython:: python

            gf = galois.GF(3**5)
            sh = Shamir(gf, 4, ['a', 'b', 'c', 'd', 'e', 'f'])
    """
    def __init__(self, cuerpo: GF, r: int, participantes: Sequence[str]):
        r"""
        Crea un esquema de compartición de secretos de Shamir sobre el cuerpo $\mathbb{F}_{p^m}$.
        :param cuerpo: El cuerpo finito sobre sobre el que el esquema está construido.
        :param r: Parámetro de reconstrucción del esquema (número mínimo de participantes necesarios para reconstruir el secreto).
        :param participantes: Lista de los identificadores únicos de cada participante del esquema.
        """

        # Verificación de condiciones
        if cuerpo.order <= len(participantes):
            raise ValueError(f'El numero de participantes ({len(participantes)}) debe ser menor que el orden del cuerpo de trabajo ({cuerpo.order}).')
        if len(participantes) < 2:
            raise ValueError(f'El numero de participantes ({len(participantes)}) debe ser mayor que 1.')
        if len(participantes) != len(set(participantes)):
            raise ValueError(f'Se han encontrado participantes duplicados.')
        if len(participantes) < r:
            raise ValueError(f'El parámetro de reconstrucción ({r}) debe ser menor o igual que el número de participantes ({len(participantes)}).')
        if r < 2:
            raise ValueError(f'El parámetro de reconstrucción ({r}) debe ser mayor que 1.')

        self._cuerpo = cuerpo
        self._reconstruccion = r
        self.__participaciones_anticipadas = None
        self._longitud_bytes = ((cuerpo.order - 1).bit_length() + 7) // 8
        self._participantes_nombre = [None] # Array para pasar de numero -> nombre
        self._participantes_numero = {} # Diccionario para pasar nombre -> numero
        for i, nombre in enumerate(participantes, 1):
            self._participantes_nombre.append(nombre)
            self._participantes_numero[nombre] = i

    def crear_anticipadas(self, participantes: Sequence[str]) -> list[tuple[str, str]]:
        """
        Crea participaciones anticipadas para cada participante especificado.
        El formato de las participaciones es: (Identificador, Participación).
        :param participantes: Listado de los participantes a entregar participaciones anticipadas.
        :return: Una lista que contienene las participaciones anticipadas asignadas a cada participante especificado.
        """

        # Verificación de condiciones
        if len(participantes) < 1:
            raise ValueError(f'El numero de participaciones anticipadas ({len(participantes)}) debe ser al menos 1.')
        if self._reconstruccion - 1 < len(participantes):
            raise ValueError(f'El numero de participaciones anticipadas ({len(participantes)}) debe ser menor o igual que el parámetro de privacidad ({self._reconstruccion - 1})')
        self._verificar_nombres(participantes)

        # Generar las participaciones anticipadas, que son elementos aleatorios del cuerpo
        aleatoriedad = [secrets.randbelow(self._cuerpo.order) for _ in range(len(participantes))]
        aleatoriedad_b64 = int_a_b64str(aleatoriedad, self._longitud_bytes)
        self.__participaciones_anticipadas = list(zip(participantes, aleatoriedad_b64))
        return self.__participaciones_anticipadas

    def crear_participaciones(self, secreto: Buffer) -> list[tuple[str, str]]:
        """
        Crea las participaciones de todos los participantes de acuerdo al secreto recibido.
        El formato de las participaciones es: (Identificador, Participación).
        Si se han distribuido participaciones anticipadas, las participaciones serán coherentes con las mismas.
        :param secreto: Secreto que se quiere codificar entre todos los participantes.
        :return: Una lista que contienene las participaciones de cada participante que no ha participado en la distribución avanzada.
        """
        secreto_i = bytes_a_int(secreto)

        # Si se han repartido participaciones anticipadas, se realiza compartición avanzada
        if self.__participaciones_anticipadas is not None:
            # Obtener el elemento asociado a cada participante y decodificar su participación
            nombres, valores_b64 = zip(*self.__participaciones_anticipadas)
            puntos = self._cuerpo(list(self._participantes_numero[nombre] for nombre in nombres) + [0])
            valores = self._cuerpo(b64str_a_int(valores_b64) + [secreto_i])
            x = np.setdiff1d(list(self._participantes_numero.values()), puntos)

            # Se determina un polinomio de grado r-1 compatible con las participaciones anticipadas
            lagrange = lagrange_poly(puntos, valores)
            if len(puntos) < self._reconstruccion - 1: # Si el número de participaciones anticipadas es menor que r-1, hay que completar el polinomio con aleatoriedad
                polinomio = lagrange + Poly.Roots(puntos, field=self._cuerpo)*Poly([secrets.randbelow(self._cuerpo.order) for _ in range(self._reconstruccion - len(puntos) - 1)], field=self._cuerpo)
            else: # Si no, el único polinomio disponible es el de Lagrange
                polinomio = lagrange
            del self.__participaciones_anticipadas # Eliminación de las participaciones anticipadas para mayor seguridad

        # Si no, se sigue el proceso estándar
        else:
            polinomio = Poly([secreto_i, *(secrets.randbelow(self._cuerpo.order) for _ in range(self._reconstruccion - 1))], field=self._cuerpo, order='asc')
            x = list(self._participantes_numero.values())

        # Generar el resto de las participaciones
        participaciones_b64 = int_a_b64str(polinomio(x), self._longitud_bytes)
        return list(zip((self._participantes_nombre[p] for p in x), participaciones_b64))

    def recuperar_secreto_v1(self, participaciones: Sequence[tuple[str, str]]) -> bytes:
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (Identificador, Participación).
        Esta versión reconstruye primero el polinomio generador y a partir de él, devuelve el secreto.
        :param participaciones: Secuencia con las participaciones de los participantes que desean obtener el secreto.
        :return: El secreto.
        """

        # Verificación de condiciones
        if len(participaciones) < self._reconstruccion:
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto')
        nombres, valores_b64 = zip(*participaciones[:self._reconstruccion])
        self._verificar_nombres(nombres)

        # Obtener el elemento asociado a cada participante y decodificar su participación
        puntos = self._cuerpo(list(self._participantes_numero[nombre] for nombre in nombres))
        valores = self._cuerpo(b64str_a_int(valores_b64))

        # Reconstruir el polinomio generador y el secreto como su coeficiente independiente
        polinomio = lagrange_poly(puntos, valores)
        return int_a_bytes(polinomio.coefficients(order="asc")[0])

    def recuperar_secreto_v2(self, participaciones: Sequence[tuple[str, str]]) -> bytes:
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (Identificador, Participación).
        Esta versión reconstruye el secreto a partir de la fórmula del polinomio interpolador de Lagrange evaluado en 0.
        :param participaciones: Secuencia con las participaciones de los participantes que desean obtener el secreto.
        :return: El secreto.
        """

        # Verificación de condiciones
        r = self._reconstruccion
        if len(participaciones) < r:
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto')
        nombres, valores_b64 = zip(*participaciones[:r])
        self._verificar_nombres(nombres)

        # Obtener el elemento asociado a cada participante y decodificar su participación
        puntos = self._cuerpo(list(self._participantes_numero[nombre] for nombre in nombres))
        valores = self._cuerpo(b64str_a_int(valores_b64))

        # Calcular el valor del polinomio generador en 0 sin reconstruirlo
        mascara = ~np.eye(r, dtype=bool) # Máscara de los elementos x_h de la fórmula
        puntos_matriz = np.broadcast_to(puntos, (r, r)) # Se crea una matriz que cada fila es el array puntos
        puntos_matriz = puntos_matriz[mascara].reshape(r - 1, r) # Al usar la mascara, la matriz se aplana por lo que hay que usar reshape (trabajaremos por columnas)
        denominador = np.prod(puntos_matriz - puntos, axis=0) # Productorio del denominador
        numerador = np.prod(puntos_matriz, axis=0) # Productorio del numerador
        coef = numerador / denominador # Cálculo de l_j
        return int_a_bytes(np.sum(valores * coef)) # Se devuelve la suma y_j * l_j

    recuperar_secreto = recuperar_secreto_v2 # Alias para recuperar secreto version 2

    def _verificar_nombres(self, nombres: Sequence[str]) -> None:
        """
        Verifica que los participantes sean válidos, es decir, que no haya nombres duplicados y todos los nombres estén registrados como participantes.
        :param nombres: La secuencia de nombres que se quiere comprobar
        """
        # Comprobar no elementos duplicados
        if len(nombres) != len(set(nombres)):
            raise ValueError(f'Se han encontrado participantes duplicados.')
        # Comprobar que los participantes existen
        conjunto_nombres = self._participantes_numero
        for nombre in nombres:
            if nombre not in conjunto_nombres:
                raise ValueError(f'El participante {nombre} no está registrado.')


class ShamirSimplificado:
    r"""
    Esquema de compartición de secretos de Shamir simplificado sobre cuerpo $\mathbb{F}_{p^m}$.

    Ejemplo:
        Crea un esquema de umbral-(5,5) sobre el cuerpo $\mathbb{F}_{5^4}$ para los participantes ['a', 'b', 'c', 'd', 'e'].

        .. ipython:: python

            gf = galois.GF(5**4)
            sh = ShamirSimplificado(gf, ['a', 'b', 'c', 'd', 'e'])
    """
    def __init__(self, cuerpo: GF, participantes: Sequence[str]):
        r"""
        Crea un esquema de compartición de secretos de Shamir simlificado sobre el cuerpo $\mathbb{F}_{p^m}$.
        :param cuerpo: El cuerpo finito sobre sobre el que el esquema está construido.
        :param participantes: Lista de los identificadores únicos de cada participante del esquema.
        """

        # Verificación de condiciones
        if cuerpo.order <= len(participantes):
            raise ValueError(f'El numero de participantes ({len(participantes)}) debe ser menor que el orden del cuerpo de trabajo ({cuerpo.order}).')
        if len(participantes) < 2:
            raise ValueError(f'El numero de participantes ({len(participantes)}) debe ser mayor que 1.')
        if len(participantes) != len(set(participantes)):
            raise ValueError(f'Se han encontrado participantes duplicados.')

        self._participantes = participantes
        self._cuerpo = cuerpo
        self.__participaciones_anticipadas = None
        self._longitud_bytes = ((cuerpo.order - 1).bit_length() + 7) // 8

    def crear_anticipadas(self, participantes: Sequence[str]) -> list[tuple[str, str]]:
        """
        Crea participaciones anticipadas para cada participante especificado.
        El formato de las participaciones es: (Identificador, Participación).
        :param participantes: Listado de los participantes a entregar participaciones anticipadas.
        :return: Una lista que contienene las participaciones anticipadas asignadas a cada participante especificado.
        """
        # Verificación de condiciones
        # Verificación de condiciones
        if len(participantes) < 1:
            raise ValueError(f'El numero de participaciones anticipadas ({len(participantes)}) debe ser al menos 1.')
        if len(self._participantes) < len(participantes):
            raise ValueError(f'El numero de participaciones anticipadas ({len(participantes)}) debe ser menor o igual que el parámetro de privacidad ({len(self._participantes) - 1})')
        self._verificar_nombres(participantes)

        # Generar las participaciones anticipadas, que son elementos aleatorios del cuerpo
        aleatoriedad = [secrets.randbelow(self._cuerpo.order) for _ in range(len(participantes))]
        aleatoriedad_b64 = int_a_b64str(aleatoriedad, self._longitud_bytes)
        self.__participaciones_anticipadas = list(zip(participantes, aleatoriedad_b64))
        return self.__participaciones_anticipadas

    def crear_participaciones(self, secreto: Buffer) -> list[tuple[str, str]]:
        """
        Crea las participaciones de todos los participantes de acuerdo al secreto recibido.
        El formato de las participaciones es: (Identificador, Participación).
        Si se han distribuido participaciones anticipadas, las participaciones serán coherentes con las mismas.
        :param secreto: Secreto que se quiere codificar entre todos los participantes.
        :return: Una lista que contienene las participaciones de cada participante que no ha participado en la distribución avanzada.
        """
        secreto_i = bytes_a_int(secreto)

        # Si se han repartido participaciones anticipadas, se realiza compartición avanzada
        if self.__participaciones_anticipadas is not None:
            # Obtener el elemento asociado a cada participante y decodificar su participación
            nombres, valores_b64 = zip(*self.__participaciones_anticipadas)
            valores = self._cuerpo(b64str_a_int(valores_b64))
            x = np.setdiff1d(self._participantes, nombres)
            suma = valores.sum()
            del self.__participaciones_anticipadas # Se eliminan las participaciones anticipadas almacenadas para mayor seguridad

        # Si no, se sigue el proceso estándar
        else:
            x = self._participantes
            suma = self._cuerpo(0)

        # Se generan el resto de participaciones
        participaciones = self._cuerpo.Zeros(len(x))
        for i in range(len(x) - 1):
            participaciones[i] = secrets.randbelow(self._cuerpo.order)
        participaciones[-1] = self._cuerpo(secreto_i) - participaciones.sum() - suma # La última debe ser igual al secreo menos la suma de todas las anteriores

        # Generar el resto de las participaciones
        participaciones_b64 = int_a_b64str(participaciones, self._longitud_bytes)
        return list(zip(x, participaciones_b64))

    def recuperar_secreto(self, participaciones: Sequence[tuple[str, str]]) -> bytes:
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (Identificador, Participación).
        :param participaciones: Secuencia con las participaciones de los participantes que desean obtener el secreto.
        :return: El secreto.
        """

        # Verificación de condiciones
        if len(participaciones) < len(self._participantes):
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto')
        nombres, valores_b64 = zip(*participaciones[:len(self._participantes)])
        self._verificar_nombres(nombres)

        # Obtener las participaciones
        valores = self._cuerpo(b64str_a_int(valores_b64))

        # El secreto es la suma de todas las participaciones
        return int_a_bytes(valores.sum())

    def _verificar_nombres(self, nombres: Sequence[str]) -> None:
        """
        Verifica que los participantes sean válidos, es decir, que no haya nombres duplicados y todos los nombres estén registrados como participantes.
        :param nombres: La secuencia de nombres que se quiere comprobar
        """
        # Comprobar no elementos duplicados
        if len(nombres) != len(set(nombres)):
            raise ValueError(f'Se han encontrado participantes duplicados.')
        # Comprobar que los participantes existen
        conjunto_nombres = set(self._participantes)
        for nombre in nombres:
            if nombre not in conjunto_nombres:
                raise ValueError(f'El participante {nombre} no está registrado.')

if __name__ == '__main__':
    gf = GF(2 ** 64)

    participantes = []
    r = int(input('Escriba el parámetro de reconstrucción (nº minimo de participantes para recuperar el secreto): '))
    n = int(input('Escriba el número de particiantes: '))
    for i in range(n):
        part = input(f'Escriba el nombre del participante nº{i+1}: ')
        participantes.append(part)
    if r == n:
        sh = ShamirSimplificado(gf, participantes)
    else:
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


    secreto = input('Escriba el secreto: ').encode()
    participantes = sh.crear_participaciones(secreto)
    for nombre, part_b64 in participantes:
        print(f'{nombre}: {part_b64}')

    print('Escriba el nombre de los participantes que busca reconstruir el secreto y su participacion.')
    conjunto = []
    for _ in range(r):
        nombre = input('Nombre: ')
        participacion = input('Participacion: ')
        conjunto.append((nombre, participacion))

    secreto = sh.recuperar_secreto_v2(conjunto)
    print('secreto:', secreto.decode())


