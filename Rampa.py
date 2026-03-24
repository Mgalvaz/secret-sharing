import secrets

from galois import GF, Poly, lagrange_poly

from utils import *


class ShamirRampa:
    r"""
        Esquema de compartición de secretos de Shamir en rampa sobre el cuerpo $\mathbb{F}_{p^m}$.

        Ejemplo:
            Crea un esquema de Shamir en rampa con parámetros r = 4 y l = 3 sobre el cuerpo $\mathbb{F}_{3^5}$ para los participantes ['a', 'b', 'c', 'd', 'e', 'f'].

            .. ipython:: python

                gf = galois.GF(3**5)
                rsh = ShamirRampa(gf, 4, 3, ['a', 'b', 'c', 'd', 'e', 'f'])
        """

    def __init__(self, cuerpo: GF, r: int, l: int,  participantes: Sequence[str]):
        r"""
        Crea un esquema de compartición de secretos de Shamir en rampa sobre el cuerpo $\mathbb{F}_{p^m}$.
        :param cuerpo: El cuerpo finito sobre sobre el que el esquema está construido.
        :param r: Parámetro de reconstrucción del esquema (número mínimo de participantes necesarios para reconstruir el secreto).
        :param l: Longitud del secreto que se quiere repartir.
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
        if r <= l:
            raise ValueError(f'La longitud del secreto ({l}) debe ser menor que el parámetro de reconstrucción ({r}).')
        if  l < 2:
            raise ValueError(f'La longitud del secreto ({l}) debe ser mayor que 1.')

        self._cuerpo = cuerpo
        self._reconstruccion = r
        self._longitud_secreto = l
        self.__participaciones_anticipadas = None
        self._longitud_bytes = ((cuerpo.order - 1).bit_length() + 7) // 8
        self._participantes_nombre = [None]  # Array para pasar de numero -> nombre
        self._participantes_numero = {}  # Diccionario para pasar nombre -> numero
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
        if self._reconstruccion - self._longitud_secreto < len(participantes):
            raise ValueError(f'El numero de participaciones anticipadas ({len(participantes)}) debe ser menor o igual que el parámetro de privacidad ({self._reconstruccion - self._longitud_secreto})')
        self._verificar_nombres(participantes)

        # Generar las participaciones anticipadas, que son elementos aleatorios del cuerpo
        aleatoriedad = [secrets.randbelow(self._cuerpo.order) for _ in range(len(participantes))]
        aleatoriedad_b64 = int_a_b64str(aleatoriedad, self._longitud_bytes)
        self.__participaciones_anticipadas = list(zip(participantes, aleatoriedad_b64))
        return self.__participaciones_anticipadas

    def crear_participaciones(self, secreto: Sequence[Buffer]) -> list[tuple[str, str]]:
        """
        Crea las participaciones de todos los participantes de acuerdo al secreto recibido.
        El formato de las participaciones es: (Identificador, Participación).
        Si se han distribuido participaciones anticipadas, las participaciones serán coherentes con las mismas.
        :param secreto: Secreto que se quiere codificar entre todos los participantes.
        :return: Una lista que contienene las participaciones de cada participante que no ha participado en la distribución avanzada.
        """

        # Verificación de condiciones
        if len(secreto) != self._longitud_secreto:
            raise ValueError(f'Se esperaba un secreto de longitud {self._longitud_secreto}, pero se ha recibido uno de longitud {len(secreto)}.')

        secreto_i = bytes_a_int(secreto)

        # Si se han repartido participaciones anticipadas, se realiza compartición avanzada
        if self.__participaciones_anticipadas is not None:
            # Obtener el elemento asociado a cada participante y decodificar su participación
            nombres, valores_b64 = zip(*self.__participaciones_anticipadas)
            puntos = self._cuerpo(list(self._participantes_numero[nombre] for nombre in nombres))
            valores = self._cuerpo(b64str_a_int(valores_b64))
            x = self._cuerpo(np.setdiff1d(np.arange(1, len(self._participantes_nombre)), puntos))

            # Se determina un polinomio de grado r-1 compatible con las participaciones anticipadas
            polinomio_secreto = Poly(secreto_i, field=self._cuerpo, order='asc')
            lagrange = lagrange_poly(puntos, (valores - polinomio_secreto(puntos))/puntos**self._longitud_secreto)
            if len(puntos) < self._reconstruccion - self._longitud_secreto - 1:  # Si el número de participaciones anticipadas es menor que r-l-1, hay que completar el polinomio con aleatoriedad
                polinomio = lagrange + Poly.Roots(puntos, field=self._cuerpo) * Poly([secrets.randbelow(self._cuerpo.order) for _ in range(self._reconstruccion - len(puntos) - 1)], field=self._cuerpo)
            else:  # Si no, el único polinomio disponible es el de Lagrange
                polinomio = lagrange
            del self.__participaciones_anticipadas  # Eliminación de las participaciones anticipadas para mayor seguridad

            # Generar el resto de las participaciones
            participaciones_b64 = int_a_b64str(polinomio_secreto(x) + x**self._longitud_secreto * polinomio(x), self._longitud_bytes)

        # Si no, se sigue el proceso estándar
        else:
            polinomio = Poly([*secreto_i, *(secrets.randbelow(self._cuerpo.order) for _ in range(self._reconstruccion - self._longitud_secreto))], field=self._cuerpo, order='asc')
            x = np.arange(1, len(self._participantes_nombre))

            # Generar el resto de las participaciones
            participaciones_b64 = int_a_b64str(polinomio(x), self._longitud_bytes)

        return list(zip((self._participantes_nombre[p] for p in x.tolist()), participaciones_b64))

    def recuperar_secreto(self, participaciones: Sequence[tuple[str, str]]) -> list[bytes]:
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
        return int_a_bytes(polinomio.coefficients(order="asc")[:self._longitud_secreto])

    def _verificar_nombres(self, nombres: Sequence[str]) -> None:
        """
        Verifica que los participantes sean válidos, es decir, que no haya nombres duplicados y todos los nombres estén registrados como participantes.
        :param nombres: La secuencia de nombres que se quiere comprobar
        """
        # Comprobar no elementos duplicados
        if len(nombres) != len(set(nombres)):
            raise ValueError(f"Se han encontrado participantes duplicados.")
        # Comprobar que los participantes existen
        conjunto_nombres = self._participantes_numero
        for nombre in nombres:
            if nombre not in conjunto_nombres:
                raise ValueError(f"El participante {nombre} no está registrado.")

class McElieceSarwate:
    r"""
        Esquema de compartición de secretos de McEliece-Sarwate sobre el cuerpo $\mathbb{F}_{p^m}$.

        Ejemplo:
            Crea un esquema de McEliece-Sarwate con parámetros r = 4 y l = 3 sobre el cuerpo $\mathbb{F}_{3^5}$ para los participantes ['a', 'b', 'c', 'd', 'e', 'f'].

            .. ipython:: python

                gf = galois.GF(3**5)
                rsh = McElieceSarwate(gf, 4, 3, ['a', 'b', 'c', 'd', 'e', 'f'])
        """

    def __init__(self, cuerpo: GF, r: int, l: int,  participantes: Sequence[str]):
        r"""
        Crea un esquema de compartición de secretos de McEliece-Sarwate sobre el cuerpo $\mathbb{F}_{p^m}$.
        :param cuerpo: El cuerpo finito sobre sobre el que el esquema está construido.
        :param r: Parámetro de reconstrucción del esquema (número mínimo de participantes necesarios para reconstruir el secreto).
        :param l: Longitud del secreto que se quiere repartir.
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
        if r <= l:
            raise ValueError(f'La longitud del secreto ({l}) debe ser menor que el parámetro de reconstrucción ({r}).')
        if  l < 2:
            raise ValueError(f'La longitud del secreto ({l}) debe ser mayor que 1.')

        self._cuerpo = cuerpo
        self._reconstruccion = r
        self._longitud_secreto = l
        self.__participaciones_anticipadas = None
        self._longitud_bytes = ((cuerpo.order - 1).bit_length() + 7) // 8
        self._participantes_nombre: list[str | None] = [None] * l  # Array para pasar de numero -> nombre
        self._participantes_numero = {}  # Diccionario para pasar nombre -> numero
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
        if self._reconstruccion - self._longitud_secreto < len(participantes):
            raise ValueError(f'El numero de participaciones anticipadas ({len(participantes)}) debe ser menor o igual que el parámetro de privacidad ({self._reconstruccion - self._longitud_secreto})')
        self._verificar_nombres(participantes)

        # Generar las participaciones anticipadas, que son elementos aleatorios del cuerpo
        aleatoriedad = [secrets.randbelow(self._cuerpo.order) for _ in range(len(participantes))]
        aleatoriedad_b64 = int_a_b64str(aleatoriedad, self._longitud_bytes)
        self.__participaciones_anticipadas = list(zip(participantes, aleatoriedad_b64))
        return self.__participaciones_anticipadas

    def crear_participaciones(self, secreto: Sequence[Buffer]) -> list[tuple[str, str]]:
        """
        Crea las participaciones de todos los participantes de acuerdo al secreto recibido.
        El formato de las participaciones es: (Identificador, Participación).
        Si se han distribuido participaciones anticipadas, las participaciones serán coherentes con las mismas.
        :param secreto: Secreto que se quiere codificar entre todos los participantes.
        :return: Una lista que contienene las participaciones de cada participante que no ha participado en la distribución avanzada.
        """

        # Verificación de condiciones
        if len(secreto) != self._longitud_secreto:
            raise ValueError(f'Se esperaba un secreto de longitud {self._longitud_secreto}, pero se ha recibido uno de longitud {len(secreto)}.')
        secreto_i = self._cuerpo(bytes_a_int(secreto))

        # Si se han repartido participaciones anticipadas, se realiza compartición avanzada
        if self.__participaciones_anticipadas is not None:
            # Obtener el elemento asociado a cada participante y decodificar su participación
            nombres, valores_b64 = zip(*self.__participaciones_anticipadas)
            puntos = self._cuerpo(list(self._participantes_numero[nombre] for nombre in nombres))
            valores = self._cuerpo(b64str_a_int(valores_b64))
            x = self._cuerpo(np.setdiff1d(list(self._participantes_numero.values()), puntos))

            # Se determina un polinomio de grado r-1 compatible con las participaciones anticipadas
            polinomio_secreto = Poly(secreto_i, field=self._cuerpo, order='asc')
            lagrange = lagrange_poly(puntos, (valores - polinomio_secreto(puntos))/puntos**self._longitud_secreto)
            if len(puntos) < self._reconstruccion - self._longitud_secreto - 1:  # Si el número de participaciones anticipadas es menor que r-l-1, hay que completar el polinomio con aleatoriedad
                polinomio = lagrange + Poly.Roots(puntos, field=self._cuerpo) * Poly([secrets.randbelow(self._cuerpo.order) for _ in range(self._reconstruccion - len(puntos) - 1)], field=self._cuerpo)
            else:  # Si no, el único polinomio disponible es el de Lagrange
                polinomio = lagrange
            del self.__participaciones_anticipadas  # Eliminación de las participaciones anticipadas para mayor seguridad

            # Generar el resto de las participaciones
            participaciones_b64 = int_a_b64str(polinomio_secreto(x) + x**self._longitud_secreto * polinomio(x), self._longitud_bytes)

        # Si no, se sigue el proceso estándar
        else:
            x = np.arange(self._longitud_secreto, len(self._participantes_nombre))
            # Hay que construir un polinomio que interpole al secreto en sus respectivos puntos y que sea de grado r - 1
            puntos_secreto = self._cuerpo(np.arange(self._longitud_secreto))
            lagrange = lagrange_poly(puntos_secreto, secreto_i)
            polinomio = lagrange + Poly.Roots(puntos_secreto, field=self._cuerpo) * Poly([secrets.randbelow(self._cuerpo.order) for _ in range(self._reconstruccion - len(puntos_secreto) - 1)], field=self._cuerpo)

            # Generar el resto de las participaciones
            participaciones_b64 = int_a_b64str(polinomio(x), self._longitud_bytes)

        return list(zip((self._participantes_nombre[p] for p in x.tolist()), participaciones_b64))

    def recuperar_secreto(self, participaciones: Sequence[tuple[str, str]]) -> list[bytes]:
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
        return int_a_bytes(polinomio.coefficients(order="asc")[:self._longitud_secreto])

    def _verificar_nombres(self, nombres: Sequence[str]) -> None:
        """
        Verifica que los participantes sean válidos, es decir, que no haya nombres duplicados y todos los nombres estén registrados como participantes.
        :param nombres: La secuencia de nombres que se quiere comprobar
        """
        # Comprobar no elementos duplicados
        if len(nombres) != len(set(nombres)):
            raise ValueError(f"Se han encontrado participantes duplicados.")
        # Comprobar que los participantes existen
        conjunto_nombres = self._participantes_numero
        for nombre in nombres:
            if nombre not in conjunto_nombres:
                raise ValueError(f"El participante {nombre} no está registrado.")

if __name__ == '__main__':
    gf = GF(2 ** 64)

    participantes = []
    r = int(input('Escriba el parámetro de reconstrucción (nº minimo de participantes para recuperar el secreto): '))
    l = int(input('Escriba la longitud del secreto: '))
    n = int(input('Escriba el número de particiantes: '))
    for i in range(n):
        part = input(f'Escriba el nombre del participante nº{i+1}: ')
        participantes.append(part)

    sh = ShamirRampa(gf, r, l, participantes)

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


    secreto = []
    for i in range(l):
        part = input(f'Escriba el secreto nº{i+1}: ').encode()
        secreto.append(part)
    participantes = sh.crear_participaciones(secreto)
    for nombre, part_b64 in participantes:
        print(f'{nombre}: {part_b64}')

    print('Escriba el nombre de los participantes que busca reconstruir el secreto y su participacion.')
    conjunto = []
    for _ in range(r):
        nombre = input('Nombre: ')
        participacion = input('Participacion: ')
        conjunto.append((nombre, participacion))

    secreto = sh.recuperar_secreto(conjunto)
    print('secreto:', [s.decode() for s in secreto])