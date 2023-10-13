'''
This is an implementation of Hecke algebras
focussing on their relationship to braid groups.
For this reason, a ``BraidHeckeAlgebra`` depends only on
one parameter, namely the number of strands
of the corresponding braid group.

The skein relation used for the
``BraidHeckeAlgebra`` class is
`s_i^2 = v^2 1 + vz s_i`,
where `s_i` is the i-th standard generator
of a braid group and `1` is the trivial braid.
'''

var('z')

class BraidHeckeAlgebra:
    '''
    An algebra useful for computing the homfly polynomial
    of closures of braids.

    EXAMPLES::

        sage: H = BraidHeckeAlgebra(3)
        sage: H
        Braid Hecke algebra on 3 strands
    '''
    def __init__(self, strands, z0=z):
        self.strands = strands
        self.G = Permutations(strands).list()
        self.z0 = z0

    def __call__(self, coeffs, degree=None):
        if len(coeffs) != self.dimension():
            raise Exception("Wrong number of strands")
        return BraidHeckeAlgebraElement(self, coeffs, degree)

    def __repr__(self):
        return "Braid Hecke algebra on " + str(self.strands) + " strands with parameter " + str(self.z0)

    def __eq__(self, other):
        '''
        Checks for equality.

        For two ``BraidHeckeAlgebra``s to be equal,
        they need to have the same number of strands,
        the same index permutations,
        and the same substitutions.
        '''
        return self.strands == other.strands and self.G == other.G and self.z0 == other.z0


    def dimension(self):
        '''
        The rank of ``self`` as a module over the ring of Laurent
        polynomials in one variable ``z`` over ``ZZ``.

        EXAMPLES::

            sage: H = BraidHeckeAlgebra(4)
            sage: H.dimension()
            24
        '''
        return factorial(self.strands)

    def from_braid(self, braid_word):
        '''
        The image under the homomorphism from the braid group
        to the Hecke algebra mapping a generator `s_i`
        to the transposition ``(i, i + 1)``.

        EXAMPLES::

            sage: H = BraidHeckeAlgebra(3)
            sage: H.from_braid([1, 2, 1])
            [0, 0, 0, 0, 0, 1]
            sage: H.from_braid([1, 2, 1, 2])
            [0, 0, 0, 1, 0, z]
        '''
        result = self.one()

        # coeffs[i] is the coefficient
        # of the lexicographically i-th permutation
        for i in braid_word:
            if i > 0:
                result.gen_mul(i)
            else:
                result.gen_mul_inv(-i)
        result.simplify()
        return result

    def one(self):
        '''
        Returns the multiplicative unit of the Hecke algebra.

        Sets the degree to `0`

        EXAMPLES::

            sage: h = BraidHeckeAlgebra(3).one()
            sage: h
            [1, 0, 0, 0, 0, 0]
            sage: h.degree
            0
        '''
        return self([1] + (self.dimension() - 1) * [0], 0)

    def zero(self):
        '''
        Returns the additive unit of the Hecke algebra

        EXAMPLES::

            sage: BraidHeckeAlgebra(3).one()
            [0, 0, 0, 0, 0, 0]
        '''
        return self(self.dimension() * [0])



class BraidHeckeAlgebraElement:
    '''
    An element of a ``BraidHeckeAlgebra``.
    Keeps track of the degree as a braid word
    in order to be able to compute the homfly
    polynomial of its closure.
    The degree can be ``None`` if this is not
    the intention.

    EXAMPLES::

        sage: H = BraidHeckeAlgebra(3)
        sage: h = H.from_braid(2 * [1, -2])
        sage: h
        [z^2, -z, z^3, 1, -z^2, 0]
        sage: var('z')
        z
        sage: H([1, z, z, z^2, z^2, z^3 + z], degree=3)
        [1, z, z, z^2, z^2, z^3 + z]
    '''
    def __init__(self, H, coeffs, degree=None):
        '''
        The default constructor for ``BraidHeckeAlgebraElement``s.

        PARAMETERS::

            - ``H``: A `BraidHeckeAlgebra`
            - ``coeffs``: A list of symbolic expressions in ``z``
            - ``degree``: An integer or ``None``. Default ``None``

        EXAMPLES::

            sage: H([1, z, z, z^2, z^2, z^3 + z], degree=3)
        '''
        if len(coeffs) != factorial(H.strands):
            raise Exception("length of coefficients does not match the number of strands")
        self.H = H
        self.coeffs = coeffs
        self.degree = degree
        self.simplify()
        # will later cache the multiplicative homomorphism if computed
        self.A = None

    def __repr__(self):
        '''
        A string representation.

        EXAMPLES::

            sage: print(BraidHeckeAlgebra(3).one())
            [1, 0, 0, 0, 0, 0]
        '''
        return str(self.coeffs)

    def __add__(self, other):
        '''
        The sum of two ``BraidHeckeAlgebraElement``s.

        Raises an exception if ``other`` has different length.
        Sets the degree of the result to the previous degrees if they agree, and to ``None`` otherwise.

        PARAMETERS::

            - ``other``: A ``BraidHeckeAlgebraElement``

        EXAMPLES::

            sage: H = BraidHeckeAlgebra(2)
            sage: h1 = H([1, 0])
            sage: h1.degree is None
            True
            sage: h2 = H.one()
            sage: h2.degree
            0
            sage: sum = one + h
            sage: sum
            [2, 0]
            sage: sum.degree
            0
        '''
        # adding scalars
        if not isinstance(other, BraidHeckeAlgebraElement):
            other = other * self.H.one()

        if len(self) != len(other):
            raise Exception("Tried to add elements of different lengths")

        result = [self[i] + other[i] for i in range(len(self))]
        if not self.degree is None and not other.degree is None:
            if self.degree == other.degree:
                return self.H(result, self.degree)
        return self.H(result)

    def __neg__(self):
        return self.H([-x for x in self.coeffs], self.degree)

    def __sub__(self, other):
        return self + (-other)

    def copy(self):
        '''
        A deep copy.
        '''
        return self.H(list(self.coeffs), self.degree)

    def __mul__(self, other):
        '''
        The product of two ``BraidHeckeAlgebraElement``s.

        Sets the degree of the product to the sum
        of the degrees of the factors, if they are both defined,
        and to ``None`` otherwise.

        Raises an exception if the parent algebras are distinct.

        PARAMETERS::
            - ``other``: A ``BraidHeckeAlgebraElement``

        EXAMPLES::
            H = BraidHeckeAlgebra(2)
            h = H.one()
            h * h == h
            True
        '''
        n = len(self)
        G_list = self.get_G()
        G = Permutations(self.get_strands())
        result = self.H.zero()


        # compute which summands are needed for the computation
        necessary = set([G_list[i] for i in range(len(self)) if not other[i] == 0])
        # close under subwords
        for g in necessary:
            reduced = g.reduced_word()
            reduced.reverse()
            necessary = necessary.union([G.from_reduced_word(reversed(reduced[:i])) for i in range(len(reduced))])

        # initialize the dictionary of summands
        smnd_dict = { str(G_list[0].reduced_word()) : self.copy() }
        for i in range(n):
            if i == 0:
                result += self.H([other[i] * self[j] for j in range(n)])
                continue

            perm = G_list[i]
            if perm in necessary:
                gens = perm.reduced_word()
                gens.reverse()
                smnd = smnd_dict[str(gens[:perm.length()-1])].copy()
                smnd.gen_mul(gens[perm.length()-1])
                smnd_dict[str(gens)] = smnd
                result += self.H([other[i] * smnd[j] for j in range(n)])
                result.simplify()
        if self.H != other.H:
            raise Exception("Multiplying BraidHeckeAlgebra elements with different parents")

        # set degree of the product
        if not self.degree is None and not other.degree is None:
            result.degree = self.degree + other.degree

        return result

    def __rmul__(self, n):
        '''
        Scalar multiplication.
        '''
        return self.H([n * x for x in self[:]], self.degree)

    def __pow__(self, n):
        '''
        The power of a ``BraidHeckeAlgebraElement`` to an integer.

        Multiplies the degree by the power.

        PARAMETERS::
            - ``n``: An integer

        EXAMPLES::

            sage: H = BraidHeckeAlgebra(2)
            sage: h1 = H.from_braid([1])
            sage: h2 = H.from_braid(3 * [1])
            sage: h1^3 == h2
            True
        '''
        result = self.H.one()
        if n < 0:
            return self.inverse()^-n
        if n == 0:
            return result
        if n == 1:
            return self.copy()

        # compute powers of powers of two
        results = [self, (self*self).simplify()]
        k = 1
        while 2^(k+1) <= n:
            results.append((results[k]*results[k]).simplify())
            k += 1
        n -= 2^k
        result = results[k]

        # multiply with the appropriate ones
        while n > 0:
            k = 0
            while 2^(k+1) <= n:
                k += 1
            result *= results[k]
            result.simplify()
            n -= 2^k
        return result

    def __getitem__(self, index):
        '''
        Access a coordinate using subscripts.

        PARAMETERS::
            - ``index``: An integer or a ``Permutation``

        EXAMPLES::

            sage: h = BraidHeckeAlgebra(3).from_braid([1, 2, 1])
            sage: h[5]
            1
            sage: g = Permutation([3, 2, 1])
            sage: h[g]
            1
        '''
        if isinstance(index, Permutation):
            return self.coeffs[self.get_G().index(index)]
        return self.coeffs[index]

    def __setitem__(self, index, value):
        '''
        Sets a coordinate using subscripts

        PARAMETERS::
            - ``index``: An integer or a ``Permutation``
            - ``value``: A symbolic expression in ``z``
        '''

        self.coeffs[index] = value

    # does not check for degree equality
    def __eq__(self, other):
        '''
        Checks for equality.

        In order for two ``BraidHeckeAlgebraElement``s to be equal,
        they need to have the same parent
        and the same coefficients.
        '''
        if self.H != other.H:
            return False
        n = len(self)
        return all(self.coeffs[i] == other.coeffs[i] for i in range(n))

    def __len__(self):
        '''
        The length of the coefficient list.
        '''
        return len(self.coeffs)

    def mul_homo(self):
        '''
        A homomorphism describing the multiplicative effect.

        EXAMPLES::

            sage: h = BraidHeckeAlgebra(2).from_braid([1, 1])
            sage: h.mul_homo()
            [       1       z]
            [       z z^2 + 1]
        '''
        if not self.A is None:
            return self.A
        n = len(self)
        A = []
        for g in self.get_G():
            result = self.copy()
            # multiply with g from the right
            word = g.reduced_word()
            word.reverse()
            for i in word:
                result.gen_mul(i)
            result.simplify()
            A.append(result.coeffs)
        self.A = Matrix(A).transpose()
        return self.A

    def inverse(self):
        '''
        The multiplicative inverse.

        EXAMPLES::

            sage: h = BraidHeckeAlgebra(3).from_braid(4 * [1, 2])
            sage: h * h.inverse()
            [1, 0, 0, 0, 0, 0]
        '''
        degree = -self.degree
        coeffs = list((self.mul_homo()^-1).transpose()[0])
        return self.H(coeffs, degree)

    def gen_mul(self, i):
        '''
        Multiplies with the generator ``(i, i + 1)`` from the right.

        Modifies the coefficients and the degree of ``self``.

        PARAMETERS::
            - ``i``: An integer
        '''
        # reset the multiplicative homomorphism
        self.A = None
        perm = Permutation((i, i + 1))
        result = self.H(len(self) * [0])
        G = self.get_G()
        for j, c in enumerate(self.coeffs):
            # != sometimes breaks for algebraic numbers (9.6)
            if not self[j] == 0:
                current = G[j]
                product = current * perm
                index = G.index(product)
                if product.length() > current.length():
                    result[index] += self[j]
                else:
                    result[index] += self[j]
                    result[j] += self.get_z0() * self[j]
        self.coeffs = result.coeffs
        if not self.degree is None:
            self.degree += 1
        return self

    def gen_mul_inv(self, i):
        '''
        Multiplies with the
        inverse of the generator ``(i, i + 1)`` from the right.

        Modifies the coefficients and the degree of ``self``.

        PARAMETERS::
            - ``i``: An integer
        '''
        self.A = None
        n = len(self)
        perm = Permutation((i, i + 1)) * self.get_G()[0]
        index = self.get_G().index(perm)
        coeffs = [-self.get_z0()] + (index - 1) * [0] + [1] + (len(self) - index - 1) * [0]
        gen_inv = self.H(coeffs)
        self.coeffs = (self * gen_inv).coeffs
        if not self.degree is None:
            self.degree -= 1
        return self


    def find_invariant(self, power=1):
        """
        Returns all minimal polynomials of
        numbers ``z`` such that ``self^power``
        is a polynomial multiple of 1.

        For such ``z``, multiplication with ``self^power``
        acts as multiplication by a polynomial.

        INPUTS::
            - `power`: an integer. Default `power=1`.

        EXAMPLES::

            sage: H = BraidHeckeAlgebra(2)
            sage: h = H.from_braid([1])
            sage: h.find_invariant(power=6)
            {z^2 + 1, z^2 + 3}
        """
        if not self.get_z0() == z:
            raise Exception("tried to substitute, even though", self.get_z0(), "was already substituted")
        h = self^power
        start = 1
        # start at the first nonzero index
        while h[start] == 0 and start < len(h) - 1:
            start += 1

        # if all components are zero, all substitutions are invariant
        if h[start] == 0:
            return set([0])

        polys = set()
        first = True
        for i in range(start, len(self)):
            # cannot be the case if first = True
            if h[i] == 0:
                continue

            # no way to make a non-zero constant zero
            if h[i] in CC:
                return set()

            # if there are no polynomials left, don't continue computing
            if not first and len(polys) == 0:
                break

            # write h[i] as f(z) * g(v, z) for candidate
            # of maximal degree.
            # then new_polys consists of the factors of f
            new_polys = set()
            candidate = h[i]

            for poly,mul in candidate.factor_list():
                # don't include the constant polynomials or
                # polynomials with negative multiplicities (can only be z)
                variables = poly.variables()
                if len(variables) == 1:
                    new_polys.add(poly)

            if first:
                # initialize potential polynomials.
                # exclude the minimal polynomial of 0
                # because this makes little sense
                # to plug into a Laurent polynomial
                polys = new_polys.difference(set([z]))
                first = False
            else:
                polys = polys.intersection(new_polys)
        return polys


    def simplify(self):
        '''
        Tries to apply ``simplify_full()`` to all
        coefficients.
        '''
        def simplify(c):
            if isinstance(c, Expression):
                if c.is_trivial_zero():
                    return 0
                if self.get_z0().is_rational_expression():
                    return c.rational_expand()
                else:
                    try:
                        return c.simplify_full()
                    except:
                        return c
            else:
                return c
        self.coeffs = [simplify(c) for c in self.coeffs]
        return self

    def homfly_polynomial(self):
        '''
        Computes the Homfly polynomial
        (or various substitutions thereof)
        of the closure.

        EXAMPLES::

            sage: h = BraidHeckeAlgebra(2).from_braid([1, 1])
            sage: h.homfly_polynomial()
            -(v^3 - v*z^2 - v)/z
        '''
        homflys = []
        B = BraidGroup(self.get_strands())
        v = var('v')
        for g in self.get_G():
            length = g.length()
            word = g.reduced_word()
            braid = B(word)
            comps = braid.components_in_closure()
            link = Link(braid)

            # one would like to write
            # homflys.append(link.homfly_polynomial()).
            # but for some reason sage removes non-involved strands.
            # it leaves one component in case of the trivial braid though.

            if g.is_one():
                big_B = BraidGroup(2 * self.get_strands())
                big_word = [2 * i + 1 for i in range(self.get_strands())]
            else:
                difference = comps - link.number_of_components()
                big_B = BraidGroup(self.get_strands() + 2 * difference)
                big_word = word + [self.get_strands() - 1 + 2 * (i+1) for i in range(difference)]

            big_braid = big_B(big_word)
            homfly = Link(big_braid).homfly_polynomial(normalization='vz')
            if homfly in CC:
                homflys.append(homfly)
            else:
                homflys.append(homfly(v, self.get_z0()))
        G = self.get_G()
        homfly_self = [self[i] * v^(self.degree - G[i].length()) for i in range(len(self))]
        result = vector(homfly_self).dot_product(vector(homflys))
        try:
            return result.simplify_full()
        except:
            return result



    def get_strands(self):
        '''
        Returns the number of strands. If the number of strands is ``n``,
        then the dimension is ``factorial(n)``.
        '''
        return self.H.strands

    def get_G(self):
        '''
        Returns the permutation list of the same cardinality
        as the length.
        '''
        return self.H.G

    def get_z0(self):
        '''
        Returns ``z`` if ``z`` was not substituted, and the substitution
        ``z0`` otherwise.
        '''
        return self.H.z0

    def trace(self):
        '''
        The trace as an endomorphism of the Hecke algebra.
        '''
        return self.mul_homo().trace()

    def jordan_form(self, transformation=False):
        '''
        The Jordan normal form as an endomorphism of the Hecke algebra.
        '''
        return self.mul_homo().jordan_form(transformation=transformation)

    def eigenvalues(self):
        '''
        The eigenvalues as an endomorphism of the Hecke algebra.
        '''
        return self.mul_homo().eigenvalues()

    def minimal_polynomial(self):
        '''
        The minimal polynomial.
        '''
        # determine rank
        powers = [self.H.one()]
        for k in range(1, len(self) + 1):
            powers += [self * powers[-1]]
            M = matrix([p[:] for p in powers])
            if M.rank() < len(powers):
                break
        rank = len(powers) - 1

        # make variables
        variables = []
        for k in range(rank + 1):
            variables += [var('a' + str(k))]

        # create system of equations
        equations = [variables[rank] == 1]
        for i in range(len(self)):
            lhs = 0
            for k in range(rank + 1):
                lhs += variables[k] * powers[k][i]
            equations += [lhs == 0]

        solutions = solve(equations, variables)[0]
        coefficients = [sol.rhs() for sol in solutions]

        polynomial = 0
        x = var('x')
        for k in range(rank + 1):
            polynomial += coefficients[k] * x^k
        return polynomial.collect(x)
