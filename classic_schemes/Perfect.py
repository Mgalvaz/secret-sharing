import warnings
import numpy as np
from galois import Poly, lagrange_poly, GF

from utils import random_array, random_polynomial, bytes_to_int, int_to_bytes, int_to_b64str, b64str_to_int

class Shamir:
    r"""
    Shamir's secret sharing scheme over the finite field $\mathbb{F}_{p^m}$.

    Example:
        Creates a (4,6)-threshold Shamir scheme over the field $\mathbb{F}_{3^5}$ for the participants ['a', 'b', 'c', 'd', 'e', 'f'].

        .. ipython:: python

            sh = Shamir(3**5, 4, ['a', 'b', 'c', 'd', 'e', 'f'])
    """
    def __init__(self, order, r, participants):
        r"""
        Creates a Shamir secret sharing scheme over the finite field $\mathbb{F}_{p^m}$.
        :param order: The order of the finite field over which the scheme is constructed.
        :param r: The reconstruction threshold of the scheme, i.e., the minimum number of participants required to reconstruct the secret.
        :param participants: A list containing the unique names of all participants in the scheme.
        """
        # Condition checks
        if order <= len(participants):
            raise ValueError(f'The number of participants ({len(participants)}) must be less than the order of the underlying field ({order}).')
        if len(participants) != len(set(participants)):
            raise ValueError('Duplicate participants were found.')
        if len(participants) < r:
            raise ValueError(f'The reconstruction threshold ({r}) must be less than or equal to the number of participants ({len(participants)}).')
        if r < 2:
            raise ValueError(f'The reconstruction threshold ({r}) must be greater than 1.')

        self.field = GF(order)
        self.reconstruction = r
        self.__advance_shares = []
        self.byte_length = ((order - 1).bit_length() + 7) // 8
        self.participants_name = np.array([None] + participants)
        self.participants_number = {nombre: i for i, nombre in enumerate(participants, 1)}

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
        if self.reconstruction - len(self.__advance_shares) <= len(advance_participants):
            raise ValueError(f'The number of advance shares ({len(advance_participants) + len(self.__advance_shares)}) must be less than or equal to the privacy threshold ({self.reconstruction - 1}).')

        # Generate the advance shares, which are random field elements
        randomness = random_array(self.field.order, len(advance_participants))
        randomness_b64 = int_to_b64str(randomness, self.byte_length)
        extras = list(zip(advance_participants, randomness_b64))
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
        if self.__advance_shares is None:
            raise AttributeError('All shares have already been distributed.')
        secret_int = bytes_to_int(secret)
        if secret_int >= self.field.order:
            raise ValueError(f'The provided secret must be smaller than the order of the underlying field ({self.field.order}).')

        # Standard procedure
        if len(self.__advance_shares) == 0:
            polynomial = Poly(random_array(self.field.order, self.reconstruction - 1) + [secret_int], field=self.field)
            x = np.arange(1, len(self.participants_name))
        # Advance sharing
        else:
            # Obtain the element associated with each participant and decode their advance share
            names, values_b64 = zip(*self.__advance_shares)
            x_advance = self.field(list(self.participants_number[name] for name in names) + [0])
            y_advance = self.field(b64str_to_int(values_b64) + [secret_int])
            x = np.setdiff1d(np.arange(1, len(self.participants_name)), x_advance)

            # Determine a polynomial of degree r-1 consistent with the advance shares
            lagrange = lagrange_poly(x_advance, y_advance)
            if len(x_advance) < self.reconstruction - 1: # If fewer than r-1 advance shares are available, complete the polynomial with randomness
                polynomial = lagrange + Poly.Roots(x_advance, field=self.field) * random_polynomial(self.field, self.reconstruction - len(x_advance) - 2)
            else: # Otherwise, the Lagrange polynomial is the only possible one
                polynomial = lagrange

        # Generate the remaining shares
        shares_b64 = int_to_b64str(polynomial(x), self.byte_length)
        self.__advance_shares = None  # Delete the stored advance shares for further security
        return list(zip(self.participants_name[x], shares_b64))

    def _alternative_reconstruct(self, shares):
        """
        Reconstructs the secret encoded in the provided shares.
        The shares are represented as tuples of the form (name, share).
        This version first reconstructs the generating polynomial and then returns the secret as its constant coefficient. This method is correct but less efficient than ``reconstruct``, which computes the secret directly without creating the Lagrange polynomial.

        :param shares: Sequence containing the shares of the participants who wish to reconstruct the secret.
        :return: The secret.
        """
        warnings.warn('You are using the alternative reconstruction method. The reconstruct method is recommended because it is significantly faster.', UserWarning)

        # Condition checks
        if len(shares) < self.reconstruction:
            raise ValueError('Not enough shares were provided to recover the secret.')
        names, values_b64 = zip(*shares[:self.reconstruction])
        self._validate_names(names)

        # Obtain the element associated with each participant and decode their share
        x = self.field(list(self.participants_number[name] for name in names))
        y = self.field(b64str_to_int(values_b64))
        # Reconstruct the generating polynomial and obtain the secret as its constant coefficient
        polynomial = lagrange_poly(x, y)
        return int_to_bytes(polynomial.coefficients(order="asc")[0])

    def reconstruct(self, shares):
        """
        Reconstructs the secret encoded in the provided shares.
        The shares are represented as tuples of the form (name, share).
        This version reconstructs the secret directly using the Lagrange interpolation formula evaluated at 0.
        :param shares: Sequence containing the shares of the participants who wish to reconstruct the secret.
        :return: The secret.
        """
        # Condition checks
        r = self.reconstruction
        if len(shares) < r:
            raise ValueError('Not enough shares were provided to recover the secret.')
        names, values_b64 = zip(*shares[:r])
        self._validate_names(names)

        # Obtain the element associated with each participant and decode # their share
        x = self.field(list(self.participants_number[name] for name in names))
        y = self.field(b64str_to_int(values_b64))
        # Compute the value of the generating polynomial at 0 without explicitly reconstructing it
        mask = ~np.eye(r, dtype=bool) # Mask for the x_h elements in the formula
        x_matrix = np.broadcast_to(x, (r, r)) # Matrix in which each row is the points array
        x_matrix = x_matrix[mask].reshape((r, r - 1)).T # Applying the mask flattens the matrix, so reshape is required
        numerator = np.prod(x_matrix, axis=0) # Numerator product
        denominator = np.prod(x_matrix - x, axis=0) # Denominator product
        coeffs = numerator / denominator # Compute l_i
        return int_to_bytes(np.sum(y * coeffs)) # Return the sum of y_i * l_i

    def _validate_names(self, names):
        """
        Verifies that the participants are valid, i.e., that there are no duplicate names and that all names correspond to registered participants.
        :param names: Sequence of participant names to validate.
        """
        # Check for duplicate names
        if len(names) != len(set(names)):
            raise ValueError('Duplicate participants were found.')
        # Check that all participants are registered
        names_set = self.participants_number
        for name in names:
            if name not in names_set:
                raise ValueError(f"Participant '{name}' is not registered")


class Additive:
    r"""
    Additive secret sharing scheme over the finite field $\mathbb{F}_{p^m}$.

    Ejemplo:
        Creates a (5,5)-threshold scheme over the field $\mathbb{F}_{5^4}$ for the participants ['a', 'b', 'c', 'd', 'e'].
        .. ipython:: python

            sh = Simplificado(5**4, ['a', 'b', 'c', 'd', 'e'])
    """
    def __init__(self, order, participants):
        r"""
        Creates an additive secret sharing scheme over the finite field $\mathbb{F}_{p^m}$.
        :param order: The order of the finite field over which the scheme is constructed.
        :param participants: A list containing the unique names of all participants in the scheme.
        """
        # Condition checks
        if len(participants) < 2:
            raise ValueError(f'The number of participants ({len(participants)}) must be greater than 1.')
        if len(participants) != len(set(participants)):
            raise ValueError('Duplicate participants were found.')

        self.participants = participants
        self.field = GF(order)
        self.byte_length = ((order - 1).bit_length() + 7) // 8
        self.__advance_shares = []

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
        if len(self.participants) - len(self.__advance_shares) <= len(advance_participants):
            raise ValueError(f'The number of advance shares ({len(advance_participants) + len(self.__advance_shares)}) must be less than or equal to the privacy threshold ({len(self.participants) - 1}).')

        # Generate the advance shares, which are random field elements
        randomness = random_array(self.field.order, len(advance_participants))
        randomness_b64 = int_to_b64str(randomness, self.byte_length)
        extras = list(zip(advance_participants, randomness_b64))
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
        if self.__advance_shares is None:
            raise AttributeError('All shares have already been distributed.')
        secret_int = bytes_to_int(secret)
        if secret_int >= self.field.order:
            raise ValueError(f'The provided secret must be smaller than the order of the underlying field ({self.cuerpo.order}).')

        # Standard procedure
        if len(self.__advance_shares) == 0:
            x = self.participants
            advance_sum = self.field(0)

        # Advance sharing
        else:
            # Obtain the element associated with each participant and decode their advance share
            names, values_b64 = zip(*self.__advance_shares)
            y = self.field(b64str_to_int(values_b64))
            x = np.setdiff1d(self.participants, names).tolist()
            advance_sum = y.sum()

        # Generate the remaining shares
        shares = self.field(random_array(self.field.order, len(x) - 1))
        shares = np.append(shares, self.field(secret_int) - shares.sum() - advance_sum) # The last share is equal to the secret minus the sum of all preceding shares
        shares_b64 = int_to_b64str(shares, self.byte_length)
        self.__advance_shares = None  # Delete the stored advance shares for further security
        return list(zip(x, shares_b64))

    def reconstruct(self, participaciones):
        """
        Reconstruye el secreto codificado en las participaciones proporcionadas.
        El formato de las participaciones es: (nombre, participación).
        :param participaciones: Secuencia con las participaciones de los participantes que desean obtener el secreto.
        :return: El secreto.
        """
        # Condition checks
        if len(participaciones) < len(self.participantes):
            raise ValueError('No se han proporcionado suficientes participaciones para recuperar el secreto.')
        nombres, valores_b64 = zip(*participaciones[:len(self.participantes)])
        self._validate_names(nombres)

        # Obtener las participaciones
        valores = self.cuerpo(b64str_to_int(valores_b64))
        # El secreto es la suma de todas las participaciones
        return int_to_bytes(valores.sum())

    def _validate_names(self, nombres):
        """
        Verifica que los participantes sean válidos, es decir, que no haya nombres duplicados y todos los nombres estén registrados como participantes.
        :param nombres: La secuencia de nombres que se quiere comprobar
        """
        # Comprobar no elementos duplicados
        if len(nombres) != len(set(nombres)):
            raise ValueError(f'Se han encontrado participantes duplicados.')
        # Comprobar que los participantes existen
        conjunto_nombres = set(self.participantes)
        for nombre in nombres:
            if nombre not in conjunto_nombres:
                raise ValueError(f"El participante '{nombre}' no está registrado.")