import numpy as np
from galois import Poly, lagrange_poly, GF

from utils import random_array, random_polynomial, bytes_to_int, int_to_bytes, int_to_b64str, b64str_to_int


class RampShamir:
    r"""
        Shamir's ramp secret sharing scheme over the finite field $\mathbb{F}_{p^m}$.

        Example:
            Creates a Shamir ramp scheme with parameters r = 4 and l = 3 over the field $\mathbb{F}_{3^5}$ for participants ['a', 'b', 'c', 'd', 'e', 'f'].

            .. ipython:: python

                rsh = ShamirRampa(3**5, 4, 3, ['a', 'b', 'c', 'd', 'e', 'f'])
        """

    def __init__(self, order, r, l,  participants):
        r"""
        Creates a Shamir ramp secret sharing scheme over the finite field $\mathbb{F}_{p^m}$.
        :param order: The order of the finite field over which the scheme is constructed.
        :param r: The reconstruction threshold of the scheme, i.e., the minimum number of participants required to reconstruct the secret.
        :param l: Length of the secret to be shared.
        :param participants: A list containing the unique names of all participants in the scheme.
        """
        # Condition checks
        if order <= len(participants):
            raise ValueError(f'The number of participants ({len(participants)}) must be less than the order of the underlying field ({order}).')
        if len(participants) != len(set(participants)):
            raise ValueError('Duplicate participants were found.')
        if  l < 2:
            raise ValueError(f'The length of the secret ({l}) must be greater than 1.')
        if r <= l:
            raise ValueError(f'The length of the secret ({l}) must be less than the reconstruction threshold ({r}).')
        if len(participants) < r:
            raise ValueError(f'The reconstruction threshold ({r}) must be less than or equal to the number of participants ({len(participants)}).')

        self.field = GF(order)
        self.reconstruction = r
        self.secret_length = l
        self.__advance_shares = []
        self.byte_length = ((order - 1).bit_length() + 7) // 8
        self.participants_name = np.array([None] + participants)
        self.participants_number = {name: i for i, name in enumerate(participants, 1)}

    def advance_sharing(self, advance_participants):
        """
        Creates advance shares for each specified participant.
        The shares are represented as tuples of the form (name, share).
        :param advance_participants: List of participants to receive advance shares.
        :return: A list containing the advance shares assigned to each specified participant.
        """
        # Condition checks
        if self.__advance_shares is None:
            raise AttributeError('All shares have already been distributed.')
        if len(self.__advance_shares) > 0:
            advance_names_set = set(list(zip(*self.__advance_shares))[0])
            for name in advance_participants:
                if name in advance_names_set:
                    raise ValueError(f"Participant '{name}' has already received an advance share.")
        self._validate_names(advance_participants)
        if self.reconstruction - self.secret_length < len(advance_participants) + len(self.__advance_shares):
            raise ValueError(f'The number of advance shares ({len(advance_participants) + len(self.__advance_shares)}) must be less than or equal to the privacy threshold ({self.reconstruction - self.secret_length}).')

        # Generar las participaciones anticipadas, que son elementos aleatorios del cuerpo
        randomness = random_array(self.field.order, len(advance_participants))
        shares_b64 = int_to_b64str(randomness, self.byte_length)
        extras = list(zip(advance_participants, shares_b64))
        self.__advance_shares.extend(extras)
        return extras

    def distribute(self, secret):
        """
        Creates the shares for all participants according to the given secret.
        The shares are represented as tuples of the form (name, share).
        If advance shares have been assigned, the generated shares will be consistent with them.
        :param secret: The secret to be shared among the participants.
        :return: A list containing the shares of all participants who did not receive an advance share.
        """
        # Condition checks
        if self.__advance_shares is None:
            raise AttributeError('All shares have already been distributed')
        if len(secret) != self.secret_length:
            raise ValueError(f'A secret with length {self.secret_length} was expected, but one with length {len(secret)} was provided.')
        secret_int = bytes_to_int(secret)
        for i, sec in enumerate(secret_int, 1):
            if sec >= self.field.order:
                raise ValueError(f'Secret nº{i} must be smaller than the order of the underlying field ({self.field.order}).')

        # Standard procedure
        if len(self.__advance_shares) == 0:
            polynomial = Poly(secret_int + random_array(self.field.order, self.reconstruction - self.secret_length), field=self.field, order='asc')
            x = np.arange(1, len(self.participants_name))
            # Generate all shares
            shares_b64 = int_to_b64str(polynomial(x), self.byte_length)
        # Advance sharing
        else:
            # Obtain the element associated with each participant and decode their advance share
            names, values_b64 = zip(*self.__advance_shares)
            x_advance = self.field(list(self.participants_number[nombre] for nombre in names))
            y_advance = self.field(b64str_to_int(values_b64))
            x = np.setdiff1d(np.arange(1, len(self.participants_name)), x_advance)

            # Determine a polynomial of degree r-1 consistent with the advance shares
            polynomial_s = Poly(secret_int, field=self.field, order='asc') # f_s
            lagrange = lagrange_poly(x_advance, (y_advance - polynomial_s(x_advance)) / x_advance ** self.secret_length)
            if len(x_advance) < self.reconstruction - self.secret_length:  # If the number of advance shares is less than r-l, the plynomial must be inserted with randomness
                polynomial = lagrange + Poly.Roots(x_advance, field=self.field) * random_polynomial(self.field, self.reconstruction - self.secret_length - len(x_advance) - 1)
            else:  # Else, the only possible polynomial is Lagrange's
                polynomial = lagrange
            # Generate remaining shares
            shares_b64 = int_to_b64str(polynomial_s(x) + self.field(x) ** self.secret_length * polynomial(x), self.byte_length)

        self.__advance_shares = None  # Delete the stored advance shares for further security
        return list(zip(self.participants_name[x], shares_b64))

    def reconstruct(self, participaciones):
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (Identificador, Participación).
        :param participaciones: Secuencia con las participaciones de los participantes que desean obtener el secreto.
        :return: El secreto.
        """
        # Condition checks
        if len(participaciones) < self.reconstruccion:
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto.')
        nombres, valores_b64 = zip(*participaciones[:self.reconstruccion])
        self._verificar_nombres(nombres)

        # Obtener el elemento asociado a cada participante y decodificar su participación
        puntos = self.field(list(self.participants_number[nombre] for nombre in nombres))
        valores = self.field(b64str_to_int(valores_b64))
        # Reconstruir el polinomio generador y el secreto como sus l primeros coeficientes
        polinomio = lagrange_poly(puntos, valores)
        return int_to_bytes(polinomio.coefficients(order="asc")[:self.secret_length])

    def _verificar_nombres(self, nombres):
        """
        Verifica que los participantes sean válidos, es decir, que no haya nombres duplicados y todos los nombres estén registrados como participantes.
        :param nombres: La secuencia de nombres que se quiere comprobar
        """
        # Comprobar no elementos duplicados
        if len(nombres) != len(set(nombres)):
            raise ValueError(f"Se han encontrado participantes duplicados.")
        # Comprobar que los participantes existen
        conjunto_nombres = self.participants_number
        for nombre in nombres:
            if nombre not in conjunto_nombres:
                raise ValueError(f"El participante '{nombre}' no está registrado.")

class McElieceSarwate:
    r"""
        Esquema de compartición de secretos de McEliece-Sarwate sobre el cuerpo $\mathbb{F}_{p^m}$.

        Ejemplo:
            Crea un esquema de McEliece-Sarwate con parámetros r = 4 y l = 3 sobre el cuerpo $\mathbb{F}_{3^5}$ para los participantes ['a', 'b', 'c', 'd', 'e', 'f'].

            .. ipython:: python

                cuerpo = galois.GF(3, 5)
                rsh = McElieceSarwate(cuerpo, 4, 3, ['a', 'b', 'c', 'd', 'e', 'f'])
        """

    def __init__(self, cuerpo, r, l,  participantes):
        r"""
        Crea un esquema de compartición de secretos de McEliece-Sarwate sobre el cuerpo $\mathbb{F}_{p^m}$.
        :param cuerpo: El cuerpo finito sobre sobre el que el esquema está construido.
        :param r: Parámetro de reconstrucción del esquema (número mínimo de participantes necesarios para reconstruir el secreto).
        :param l: Longitud del secreto que se quiere repartir.
        :param participantes: Lista de los identificadores únicos de cada participante del esquema.
        """
        # Condition checks
        if cuerpo.order - l < len(participantes):
            raise ValueError(f'El numero de participantes ({len(participantes)}) debe ser menor o igual que el orden del cuerpo de trabajo menos la longitud del secreto ({cuerpo.order - l}).')
        if len(participantes) != len(set(participantes)):
            raise ValueError(f'Se han encontrado participantes duplicados.')
        if  l < 2:
            raise ValueError(f'La longitud del secreto ({l}) debe ser mayor que 1.')
        if r <= l:
            raise ValueError(f'La longitud del secreto ({l}) debe ser menor que el parámetro de reconstrucción ({r}).')
        if len(participantes) < r:
            raise ValueError(f'El parámetro de reconstrucción ({r}) debe ser menor o igual que el número de participantes ({len(participantes)}).')

        self.cuerpo = cuerpo
        self.reconstruccion = r
        self.longitud_secreto = l
        self.__participaciones_anticipadas = []
        self.longitud_bytes = ((cuerpo.order - 1).bit_length() + 7) // 8
        self.participantes_nombre = np.array([None] * l + participantes)
        self.participantes_numero = {nombre: i for i, nombre in enumerate(participantes, l)}

    def comparticion_anticipada(self, participantes_anticipados):
        """
        Crea participaciones anticipadas para cada participante especificado.
        El formato de las participaciones es: (Identificador, Participación).
        :param participantes_anticipados: Listado de los participantes a entregar participaciones anticipadas.
        :return: Una lista que contienene las participaciones anticipadas asignadas a cada participante especificado.
        """
        # Condition checks
        if self.__participaciones_anticipadas is None:
            raise AttributeError(f'Ya se han repartido todas las participaciones.')
        if len(self.__participaciones_anticipadas) > 0:
            conjunto_nombres_anticipados = set(list(zip(*self.__participaciones_anticipadas))[0])
            for nombre in participantes_anticipados:
                if nombre in conjunto_nombres_anticipados:
                    raise ValueError(f"El participante '{nombre}' ya ha recibido una participación anticipada.")
        self._verificar_nombres(participantes_anticipados)
        if self.reconstruccion - self.longitud_secreto < len(participantes_anticipados) + len(self.__participaciones_anticipadas):
            raise ValueError(f'El numero de participaciones anticipadas ({len(participantes_anticipados) + len(self.__participaciones_anticipadas)}) debe ser menor o igual que el parámetro de privacidad ({self.reconstruccion - self.longitud_secreto}).')

        # Generar las participaciones anticipadas, que son elementos aleatorios del cuerpo
        aleatoriedad = random_array(self.cuerpo.order, len(participantes_anticipados))
        aleatoriedad_b64 = int_to_b64str(aleatoriedad, self.longitud_bytes)
        extras = list(zip(participantes_anticipados, aleatoriedad_b64))
        self.__participaciones_anticipadas.extend(extras)
        return extras

    def codificacion(self, secreto):
        """
        Crea las participaciones de todos los participantes de acuerdo al secreto recibido.
        El formato de las participaciones es: (Identificador, Participación).
        Si se han distribuido participaciones anticipadas, las participaciones serán coherentes con las mismas.
        :param secreto: Secreto que se quiere codificar entre todos los participantes.
        :return: Una lista que contienene las participaciones de cada participante que no ha participado en la distribución anticipada.
        """
        # Condition checks
        if self.__participaciones_anticipadas is None:
            raise AttributeError(f'Ya se han repartido todas las participaciones.')
        if len(secreto) != self.longitud_secreto:
            raise ValueError(f'Se esperaba un secreto de longitud {self.longitud_secreto}, pero se ha recibido uno de longitud {len(secreto)}.')
        secreto_i = bytes_to_int(secreto)
        for i, sec in enumerate(secreto_i, 1):
            if sec >= self.cuerpo.order:
                raise ValueError(f'El secreto proporcionado nº{i}  debe ser menor que el orden del cuerpo de trabajo ({self.cuerpo.order}).')

        alpha = self.cuerpo.Range(0, self.longitud_secreto)
        # Standard procedure
        if len(self.__participaciones_anticipadas) == 0:
            x = np.arange(self.longitud_secreto, len(self.participantes_nombre))
            # Hay que construir un polinomio que interpole al secreto en sus respectivos puntos y que sea de grado r - 1
            lagrange = lagrange_poly(alpha, self.cuerpo(secreto_i))
            polinomio = lagrange + Poly.Roots(alpha, field=self.cuerpo) * random_polynomial(self.cuerpo,
                                                                                            self.reconstruccion - self.longitud_secreto - 1)
        # Advance sharing
        else:
            # Obtener el elemento asociado a cada participante y decodificar su participación
            nombres_anticipados, valores_anticipados_b64 = zip(*self.__participaciones_anticipadas)
            puntos_anticipados = self.cuerpo(list(self.participantes_numero[nombre] for nombre in nombres_anticipados))
            valores_anticipados = self.cuerpo(b64str_to_int(valores_anticipados_b64))
            x = np.setdiff1d(list(self.participantes_numero.values()), puntos_anticipados)
            # Se determina un polinomio de grado r-1 compatible con las participaciones anticipadas
            puntos_lagrange = self.cuerpo(np.concatenate([alpha, puntos_anticipados]))
            valores_lagrange = self.cuerpo(np.concatenate([secreto_i, valores_anticipados]))
            lagrange = lagrange_poly(puntos_lagrange, valores_lagrange)
            if len(puntos_anticipados) < self.reconstruccion - self.longitud_secreto:  # Si el número de participaciones anticipadas es menor que r-l, hay que completar el polinomio con aleatoriedad
                polinomio = lagrange + Poly.Roots(np.concatenate([alpha, puntos_anticipados]), field=self.cuerpo) * random_polynomial(
                    self.cuerpo, self.reconstruccion - self.longitud_secreto - len(puntos_anticipados) - 1)
            else:  # Si no, el único polinomio disponible es el de Lagrange
                polinomio = lagrange

        # Generar el resto de las participaciones
        participaciones_b64 = int_to_b64str(polinomio(x), self.longitud_bytes)
        self.__participaciones_anticipadas = None  # Eliminación de las participaciones anticipadas para mayor seguridad
        return list(zip(self.participantes_nombre[x], participaciones_b64))

    def _decodificacion_alternativa(self, participaciones):
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (Identificador, Participación).
        Esta versión reconstruye primero el polinomio generador y a partir de él, devuelve el secreto.
        :param participaciones: Secuencia con las participaciones de los participantes que desean obtener el secreto.
        :return: El secreto.
        """
        # Condition checks
        if len(participaciones) < self.reconstruccion:
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto')
        nombres, valores_b64 = zip(*participaciones[:self.reconstruccion])
        self._verificar_nombres(nombres)

        # Obtener el elemento asociado a cada participante y decodificar su participación
        puntos = self.cuerpo(list(self.participantes_numero[nombre] for nombre in nombres))
        valores = self.cuerpo(b64str_to_int(valores_b64))
        # Reconstruir el polinomio generador y el secreto como su evaluacion en los elementos alpha_j
        polinomio = lagrange_poly(puntos, valores)
        return int_to_bytes(polinomio(np.arange(self.longitud_secreto)))

    def decodificacion(self, participaciones):
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (Identificador, Participación).
        Esta versión reconstruye el secreto a partir de la fórmula del polinomio interpolador de Lagrange evaluado en 0, ..., l-1.
        :param participaciones: Secuencia con las participaciones de los participantes que desean obtener el secreto.
        :return: El secreto.
        """
        # Condition checks
        r = self.reconstruccion
        if len(participaciones) < r:
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto.')
        nombres, valores_b64 = zip(*participaciones[:r])
        self._verificar_nombres(nombres)

        # Obtener el elemento asociado a cada participante y decodificar su participación
        puntos = self.cuerpo(list(self.participantes_numero[nombre] for nombre in nombres))
        valores = self.cuerpo(b64str_to_int(valores_b64))
        # Calcular el valor del polinomio generador en 0, ..., l-1 sin reconstruirlo
        mascara = ~np.eye(self.reconstruccion, dtype=bool)  # Máscara de los elementos x_h de la fórmula
        coef = self.cuerpo.Zeros((self.longitud_secreto, self.reconstruccion))
        alphas = self.cuerpo.Range(0, self.longitud_secreto)[:, None]
        for i in range(self.reconstruccion):
            numerador = np.prod(alphas - puntos[mascara[i]], axis=1)  # producto del numerador aj - xh
            denominador = np.prod(puntos[i] - puntos[mascara[i]])  # producto del denominador xi - xh
            coef[:, i] = numerador / denominador  # Cálculo de l_i
        return int_to_bytes(np.sum(valores * coef, axis=1))  # Se devuelve la suma y_i * l_i(a_j)

    def _verificar_nombres(self, nombres):
        """
        Verifica que los participantes sean válidos, es decir, que no haya nombres duplicados y todos los nombres estén registrados como participantes.
        :param nombres: La secuencia de nombres que se quiere comprobar
        """
        # Comprobar no elementos duplicados
        if len(nombres) != len(set(nombres)):
            raise ValueError(f'Se han encontrado participantes duplicados.')
        # Comprobar que los participantes existen
        conjunto_nombres = self.participantes_numero
        for nombre in nombres:
            if nombre not in conjunto_nombres:
                raise ValueError(f"El participante '{nombre}' no está registrado.")