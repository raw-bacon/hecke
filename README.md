# About this repository
This repository consists of a basic implementation of the Hecke algebra with a focus on its relation to the HOMFLY polynomial, which is in the file `hecke.sage`, as well as the proof of Proposition 4.6 in my PHD thesis, which asserts the following.

*Proposition 4.6.*
1. $P_{T(3, 3\ell+1)}(v, \sqrt{-3}) = (\ell v^4 + \ell v^2 + (\ell + 1)) v^{6 \ell}$,
2. $P_{T(4, 8\ell+1)}(v, \sqrt{-2}) = (2 \ell v^6 + 2 \ell v^4 + 2 \ell v^2 + (2 \ell + 1)) v^{24 \ell}$,
3. $P_{T(6, 36 \ell + 1)}(v, \sqrt{-1}) = (6 \ell v^{10} + 6 \ell v^8 + 6 \ell v^4 + 6 \ell v^2 + (6 \ell + 1))v^{180\ell}$.

The proof of this result can be run using the command
```sh
sage infinite-orbits-torus.sage
```

# The Hecke algebra
In our context, the _Hecke algebra_ $H_n$ is
the vector space over $\mathbb{Z}[v^{\pm 1}, z^{\pm 1}]$
with basis $S_n$ and multiplication rule
$s_i \cdot s_i = v^2 + vz s_i$, where $s_i$
is the transposition exchanging $i$ and $i+1$.
Hecke algebras can be constructed using
```sage
load("hecke.sage")
H = BraidHeckeAlgebra(3)
```
where `3` may be replaced by any number $n \geq 2$.
This version of the Hecke algebra comes with
two important maps:

1. The natural representation $\omega \colon B_n \to H_n$, which can be accessed using `eta = H.from_braid([1,-2,1,-2,1,-2])`, if `H` was previously constructed as above, and `[1,-2,1,-2,1,-2]` may be replaced by a list representing any braid. Here, a positive entry $i$ stands for a positive crossing of the $i$-th strand with the $(i+1)$-st strand, and a negative entry $-i$ stands for a negative crossing of the same strands.
2. The projection $p \colon H_n \to \mathbb{Z}[v^{\pm 1}, z^{\pm 1}]$ describing the HOMFLY polynomial via $P_{\widehat \beta} = (p \circ \omega)(\beta)$. It can be accessed using `eta.homfly_polynomial()`, whenever `eta` was previously defined using the `from_braid` method. The skein relation used is $v^{-1} P_+ - v P_- = z P_0$.

#

Internally, `BraidHeckeAlgebraElements` such as `eta` above do not store
the `v`-variable, as it can be recovered just from the length of the braid
that defined `eta`,
see Proposition 4.9 in my thesis, which asserts the following.

*Proposition 4.9.* _Let $\eta = \omega(\beta) \in H_n$ be in the image of $\omega$. Then the coefficients of $\eta$ are all of the form $\pi_g(\eta) = \eta_g(z)v^{d(\beta) - \ell(g)}$ for some Laurent polynomial $\eta_g \in \mathbb{Z}[z^{\pm 1}]$. Here, $d \colon B_n \to \mathbb{Z}$ is the homomorphism sending the generator $\sigma_i$ to $1$ for all $i$, and $\ell \colon S_n \to \mathbb{N}$ is the length function._

This length is stored as `eta.degree`.
This can be relevant when adding `BraidHeckeAlgebraElement`s,
and it is not implemented if they have different degrees
(but, of course, they can always be multiplied).
