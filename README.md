# hecke
This repository consists of a basic implementation of the Hecke algebra with a focus on its relation to the HOMFLY polynomial, which is in the file `hecke.sage`, as well as the proof of Proposition 4.6 in my PHD thesis, which asserts the following.

1. $P_{T(3, 3\ell+1)}(v, \sqrt{-3}) = (\ell v^4 + \ell v^2 + (\ell + 1)) v^{6 \ell}$,
2. $P_{T(4, 8\ell+1)}(v, \sqrt{-2}) = (2 \ell v^6 + 2 \ell v^4 + 2 \ell v^2 + (2 \ell + 1)) v^{24 \ell}$,
3. $P_{T(6, 36 \ell + 1)}(v, \sqrt{-1}) = (6 \ell v^{10} + 6 \ell v^8 + 6 \ell v^4 + 6 \ell v^2 + (6 \ell + 1))v^{180\ell}$.

The proof of this result can be run using the command
```sh
sage infinite-orbits-torus.sage
```
