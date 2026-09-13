import numpy as np
from galois import Poly, lagrange_poly, GF

from utils import random_array, random_polynomial, bytes_to_int, int_to_bytes, int_to_b64str, b64str_to_int


class RampShamir:
    r"""
        Shamir's ramp secret sharing scheme over the finite field $\mathbb{F}_{p^m}$.

        Example:
            Creates a Shamir ramp scheme with parameters r = 4 and l = 3 over the field $\mathbb{F}_{3^5}$ for participants ['a', 'b', 'c', 'd', 'e', 'f'].

            .. ipython:: python

                rsh = RampShamir(3**5, 4, 3, ['a', 'b', 'c', 'd', 'e', 'f'])
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

        # Generate the advance shares, which are random field elements
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
                raise ValueError(f'The provided secret element #{i} must be smaller than the order of the underlying field ({self.field.order}).')

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
            if len(x_advance) < self.reconstruction - self.secret_length:  # If the number of advance shares is less than r-l, the polynomial must be inserted with randomness
                polynomial = lagrange + Poly.Roots(x_advance, field=self.field) * random_polynomial(self.field, self.reconstruction - self.secret_length - len(x_advance) - 1)
            else:  # If not, the only possible polynomial is Lagrange's
                polynomial = lagrange
            # Generate remaining shares
            shares_b64 = int_to_b64str(polynomial_s(x) + self.field(x) ** self.secret_length * polynomial(x), self.byte_length)

        self.__advance_shares = None  # Delete the stored advance shares for further security
        return list(zip(self.participants_name[x], shares_b64))

    def reconstruct(self, shares):
        """
        Reconstructs the secret encoded in the provided shares.
        The shares are represented as tuples of the form (name, share).
        This version reconstructs the secret directly using the Lagrange interpolation formula evaluated at 0.
        :param shares: Sequence containing the shares of the participants who wish to reconstruct the secret.
        :return: The secret.
        """
        # Condition checks
        if len(shares) < self.reconstruction:
            raise ValueError('Not enough shares were provided to recover the secret.')
        names, values_b64 = zip(*shares[:self.reconstruction])
        self._validate_names(names)

        # Obtain the element associated with each participant and decode their share
        x = self.field(list(self.participants_number[name] for name in names))
        y = self.field(b64str_to_int(values_b64))
        # Reconstruct the generating polynomial and obtain the secret as its first l coefficients
        polynomial = lagrange_poly(x, y)
        return int_to_bytes(polynomial.coefficients(order="asc")[:self.secret_length])

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

class McElieceSarwate:
    r"""
        McEliece-Sarwate secret sharing scheme over finite field $\mathbb{F}_{p^m}$.

        Example:
            Creates a McEliece-Sarwate scheme with parameters r = 4 and l = 3 over the field $\mathbb{F}_{3^5}$ for participants ['a', 'b', 'c', 'd', 'e', 'f'].

            .. ipython:: python

                es = McElieceSarwate(3**5, 4, 3, ['a', 'b', 'c', 'd', 'e', 'f'])
        """

    def __init__(self, order, r, l,  participants):
        r"""
        Creates a McEliece-Sarwate secret sharing scheme over the finite field $\mathbb{F}_{p^m}$.
        :param order: The finite field over which the scheme is constructed.
        :param r: The reconstruction threshold of the scheme, i.e., the minimum number of participants required to reconstruct the secret.
        :param l: The length of the secret to be shared.
        :param participants: A list containing the unique identifiers of all participants in the scheme.
        """
        # Condition checks
        if order - l < len(participants):
            raise ValueError(f'The number of participants ({len(participants)}) must be less than or equal to the order of the underlying field minus the secret length ({order - l}).')
        if len(participants) != len(set(participants)):
            raise ValueError('Duplicate participants were found.')
        if  l < 2:
            raise ValueError(f'The secret length ({l}) must be greater than 1.')
        if r <= l:
            raise ValueError(f'The secret length ({l}) must be less than the reconstruction threshold ({r}).')
        if len(participants) < r:
            raise ValueError(f'The reconstruction threshold ({r}) must be less than or equal to the number of participants ({len(participants)}).')

        self.field = GF(order)
        self.reconstruction = r
        self.secret_length = l
        self.__advance_shares = []
        self.byte_length = ((order - 1).bit_length() + 7) // 8
        self.participants_name = np.array([None] * l + participants)
        self.participants_number = {name: i for i, name in enumerate(participants, l)}

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
            raise ValueError(
                f'The number of advance shares ({len(advance_participants) + len(self.__advance_shares)}) must be less than or equal to the privacy threshold ({self.reconstruction - self.secret_length}).')

        # Generate the advance shares, which are random field elements
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
            raise AttributeError('All shares have already been distributed.')
        if len(secret) != self.secret_length:
            raise ValueError(f'A secret of length {self.secret_length} was expected, but a secret of length {len(secret)} was provided.')
        secret_int = bytes_to_int(secret)
        for i, sec in enumerate(secret_int, 1):
            if sec >= self.field.order:
                raise ValueError(f'The provided secret element #{i} must be smaller than the order of the underlying field ({self.field.order}).')

        alpha = self.field.Range(0, self.secret_length) # Numbers where the secret is encoded
        # Standard procedure
        if len(self.__advance_shares) == 0:
            x = np.arange(self.secret_length, len(self.participants_name))
            # Construct a polynomial of degree r-1 that interpolates the secret at its corresponding evaluation points
            lagrange = lagrange_poly(alpha, self.field(secret_int))
            polynomial = lagrange + Poly.Roots(alpha, field=self.field) * random_polynomial(self.field, self.reconstruction - self.secret_length - 1)
        # Advance sharing
        else:
            # Obtain the field element associated with each participant and decode their advance share
            names, values_b64 = zip(*self.__advance_shares)
            x_advance = self.field(list(self.participants_number[name] for name in names))
            y_advance = self.field(b64str_to_int(values_b64))
            x = np.setdiff1d(list(self.participants_number.values()), x_advance)
            # Determine a degree r-1 polynomial compatible with the secret and the advance shares
            x_lagrange = self.field(np.concatenate([alpha, x_advance]))
            y_lagrange = self.field(np.concatenate([secret_int, y_advance]))
            lagrange = lagrange_poly(x_lagrange, y_lagrange)
            if len(x_advance) < self.reconstruction - self.secret_length:  # If fewer than r-l padvance shares are available, complete the polynomial with randomness
                polynomial = lagrange + Poly.Roots(np.concatenate([alpha, x_advance]), field=self.field) * random_polynomial(self.field, self.reconstruction - self.secret_length - len(x_advance) - 1)
            else:  # Otherwise, the Lagrange polynomial is the only possible one
                polynomial = lagrange

        # Generate the remaining shares
        shares_b64 = int_to_b64str(polynomial(x), self.byte_length)
        self.__advance_shares = None  # Delete the stored advance shares for further security
        return list(zip(self.participants_name[x], shares_b64))

    def _alternative_reconstruct(self, shares):
        """
        Reconstructs the secret encoded in the provided shares.
        The shares are represented as tuples of the form (name, share).
        This version first reconstructs the generating polynomial and then returns the secret as its first l coefficients. This method is correct but less efficient than ``reconstruct``, which computes the secret directly without creating the Lagrange polynomial.
        :param shares: Sequence containing the shares of the participants who wish to reconstruct the secret.
        :return: The secret.
        """
        # Condition checks
        if len(shares) < self.reconstruction:
            raise ValueError('Not enough shares were provided to recover the secret.')
        names, values_b64 = zip(*shares[:self.reconstruction])
        self._validate_names(names)

        # Obtain the field element associated with each participant and decode their share
        x = self.field(list(self.participants_number[nombre] for nombre in names))
        y = self.field(b64str_to_int(values_b64))
        # Reconstruct the generating polynomial and evaluate it at the secret evaluation points
        polinomio = lagrange_poly(x, y)
        return int_to_bytes(polinomio(np.arange(self.secret_length)))

    def reconstruct(self, shares):
        """
        Reconstructs the secret encoded in the provided shares.
        The shares are represented as tuples of the form (name, share).
        This implementation reconstructs the secret directly by evaluating the Lagrange interpolation formula at 0, ..., l-1.
        :param shares: Sequence containing the shares of the participants who wish to reconstruct the secret.
        :return: The secret.
        """
        # Condition checks
        r = self.reconstruction
        if len(shares) < r:
            raise ValueError('Not enough shares were provided to recover the secret.')
        names, values_b64 = zip(*shares[:r])
        self._validate_names(names)

        # Obtain the field element associated with each participant and decode their share
        x = self.field(list(self.participants_number[nombre] for nombre in names))
        y = self.field(b64str_to_int(values_b64))
        # Compute the value of the generating polynomial at 0, ..., l-1 without explicitly reconstructing it
        mask = ~np.eye(self.reconstruction, dtype=bool)  # Mask for the x_h elements in the formula
        coeffs = self.field.Zeros((self.secret_length, self.reconstruction)) # Where each l_i(a_j) will be stored (coeffs[j,i] = l_i(a_j))
        alphas = self.field.Range(0, self.secret_length)[:, None]
        for i in range(self.reconstruction):
            numerator = np.prod(alphas - x[mask[i]], axis=1)  # Numerator product of each aj - xh
            denominator = np.prod(x[i] - x[mask[i]])  # Denominator product of each xi - xh
            coeffs[:, i] = numerator / denominator  # Calculate l_i
        return int_to_bytes(np.sum(y * coeffs, axis=1))  # Return the sum of y_i * l_i(a_j)

    def _validate_names(self, names):
        """
        Verifies that the participants are valid, i.e., that there are no duplicate names and that all names correspond to registered participants.
        :param names: Sequence of participant names to validate.
        """
        # Check for duplicate participants
        if len(names) != len(set(names)):
            raise ValueError('Duplicate participants were found.')
        # Check that all participants are registered
        names_set = self.participants_number
        for name in names:
            if name not in names_set:
                raise ValueError(f"Participant '{name}' is not registered.")