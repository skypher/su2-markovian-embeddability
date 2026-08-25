# Project: Markovian embeddability of SU(2)-covariant spin channels

## Full paper target

Fix a spin \(j\in \tfrac12\mathbb N_0\) and the irreducible representation
\(U_j\) on \(\mathcal H_j\).  Classify the channels on
\(B(\mathcal H_j)\) of the form \(e^{t\mathcal L}\), where
\(t>0\) and \(\mathcal L\) is an \(SU(2)\)-covariant GKSL generator.

The final paper must contain all of the following.

1. A classification of every \(SU(2)\)-covariant GKSL generator.
2. The exact Wigner-\(6j\) matrix giving its eigenvalues on the
   irreducible tensor-rank sectors.
3. Necessary and sufficient endpoint inequalities for spins \(1\) and
   \(3/2\).  Spin \(2\) is an extension after these two cases.
4. A characterization of invertible differentiable covariant
   CP-divisible paths by their instantaneous logarithmic derivatives.
5. Exact symbolic code and reproducible plots comparing the full channel
   simplex, its positive-spectrum part, and its covariantly embeddable part.

The endpoint problem always means embeddability by a covariant generator;
it does not assert anything about noncovariant logarithms of a covariant
channel.

## Fixed conventions

Set \(d=2j+1\).  Use Hilbert--Schmidt orthonormal irreducible tensor
operators \(T_q^{(k)}\), with

\[
T_0^{(0)}=I/\sqrt d,\qquad
\operatorname{Tr}\!\left((T_q^{(k)})^\dagger T_{q'}^{(k')}\right)
=\delta_{kk'}\delta_{qq'}.
\]

Define

\[
\mathcal D_k(X)=\sum_{q=-k}^k T_q^{(k)}X(T_q^{(k)})^\dagger
-\frac12\left\{\sum_{q=-k}^k(T_q^{(k)})^\dagger T_q^{(k)},X\right\}.
\]

All displayed matrices and inequalities must use this normalization, or
state the conversion explicitly.

For a covariant channel, write \(\eta_\ell\) for the scalar on the
rank-\(\ell\) sector.  The logarithmic endpoint criterion is stated only
when every \(\eta_\ell>0\).  A finite-time exponential is invertible, and
the covariant generator classification makes its sector eigenvalues real.

## Priority order

1. Invariant-generator classification.
2. Exact recoupling/eigenvalue identity with convention audit.
3. Spin-\(1\) and spin-\(3/2\) facet inequalities.
4. Differentiable-path theorem, exact computation, and figures.
5. Spin \(2\) extension and exposition.

Paper proof is the default.  Do not start Lean work unless the user asks
for Lean explicitly.
