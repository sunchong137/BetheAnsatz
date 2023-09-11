* What is this?
  This program program computes Bethe-Ansatz (BA) solutions to ground
  state of the infinite 1-dimensonal Hubbard model at arbitrary coupling
  strengths and fillings. This is one of the very few (numerically)
  exactly solvable model problems of strong correlation, and as such is
  a valuable reference system for testing new electronic structure
  methods.

  The program is based on the algorithm described in

  [1] Shiba - Magnetic susceptibility at zero temperture for the
      one-dimensional Hubbard model [Phys. Rev. B 6 930 (1972)]
      http://dx.doi.org/10.1103/PhysRevB.6.930

  combined with some generic alternating series accelerator described
  in Cohen, Villegas, Zagier, Experiment. Math. 9 pp. 3-12 (2000)
  (http://projecteuclid.org/euclid.em/1046889587).

  It was developed to compute reference data for the paper

  [2] Knizia, Chan - Density Matrix Embedding: A Simple Alternative
      to Dynamical Mean-Field Theory [Phys Rev Lett 109 186404 (2012)]
      http://dx.doi.org/10.1103/PhysRevLett.109.186404

* Compiling and running
  - The program consists of two parts: One C++ part, and one Python part.

  - The C++ part calculates the BA energy for two given parameters U and Q.
    Q is a quantity in [0,1] which encodes the filling fraction
    (<n>) implicitly. The program should not be used for U < 0.1, as the
    used iterative process may become very slow and/or instable then.

    The program can be compiled with "make", which executes:

      c++ -O3 -DNDEBUG bethe_ansatz.cpp -o bethe_ansatz

    (tested with g++ and clang). A prerequisite for this to work is a
    working C++ boost library (for boost::format) (C++ boost packages should
    be available in packages on all OSes)

  - The Python part employs the C++ part for getting out more general
    information as numerical optimizations and/or derivatives over those U,Q
    parameters. See source code for available functions and to select the
    calculation mode.

    Things which can be computed are:
      - Energy/site (E) and Chemical Potential (mu) at arbitrary coupling
        strengths U and fillings <n>.
        The chemical potential is defined as

            mu(U,n) := d/d[n] E(U,n)

      - 1-particle gap at half-filling (called "HOMO-LUMO"-gap in the
        program, which is technically incorrect but often so called in
        the context of polymers). The 1-particle gap is defined as

            gap(U) := lim(x->0) [1/2](mu(U,1/2+x) - mu(U,1/2-x)).

        I.e., as the difference in the chemical potential just above
        and just below half filling. Note that the gap is zero away
        from half filling (indicating metallicity) and only non-zero
        at half-filling.
      - On-site spin-spin and density-density self-correlation functions:
        <nA(i) nB(i)>, <Sz(i) x Sz(i)>, were nA is alpha-spin density
        and nB is beta-spin density. These quantities are computed as
        derivatives of the energy with respect to U:

            <nA(i) nB(i)> = d/d[U] E(U,n)

        The spin-correlation function is related to this via

            <Sz(i) Sz(i)> = (3./4.)*(<n> - 2*<nA(i) nB(i)>).


* Background on the development
  - Despite multiple extensive searches I was, at the time, unable
    to anywhere find either the Hubbard model reference data itself or
    programs to calculate it. Everone I asked just knew that it could
    be and has been done, but failed to give any pointers to how
    or where---execept for citing the Lieb & Wu solution from 1968.

  - The Lieb & Wu paper from 1968, which introduced the Bethe Ansatz
    solution of the 1D Hubbard model
    (http://dx.doi.org/10.1103/PhysRevLett.21.192.2)
    is well written and a mile stone in electronic structure theory.
    But it is, as written, clearly not applicable to anything except
    half filling, and generalizations of the approach to non-half-filling
    are not trivial. The core of the idea is in there, but developing
    an actual numerically working algorithm from this would be a matter
    of at least weeks, if not months or more.

  - I stumbled across the Shiba paper more or less by accident. Shiba's
    paper offers a simple approach to calculating BA solutions for
    arbitrary fillings by employing a clever integral formulation &
    transformation of ht Lieb & Wu-ansatz. I did not actually implement
    Shiba's core extension---a generalization for arbitrary
    magnetizabilities---but would do so by popular request.

-- 
Gerald Knizia, 2014-06-03
