# Markovian embeddability of SU(2)-covariant spin channels

This repository contains a paper draft and exact symbolic artifacts for
covariant Markovian embedding on a fixed irreducible spin space.

## Main result

With the Hilbert--Schmidt normalized tensor operators fixed in
`AGENTS.md`, every SU(2)-covariant GKSL generator is uniquely

\[
\mathcal L=\sum_{k=1}^{2j}\gamma_k\mathcal D_k,\qquad \gamma_k\geq0.
\]

On tensor rank \(\ell\), its eigenvalue vector is
\(\lambda=M^{(j)}\gamma\), where

\[
M^{(j)}_{\ell k}
=(2k+1)\left[
(-1)^{2j+k+\ell}
\begin{Bmatrix}j&j&k\\j&j&\ell\end{Bmatrix}
-\frac1{2j+1}\right].
\]

The inverse also has a single-symbol formula:

\[
\bigl((M^{(j)})^{-1}\bigr)_{k\ell}
=(2\ell+1)(-1)^{2j+k+\ell}
\begin{Bmatrix}j&j&k\\j&j&\ell\end{Bmatrix}.
\]

Consequently, a covariant channel with sector eigenvalues
\(\eta_0=1,\eta_1,\ldots,\eta_{2j}\) is embeddable by a covariant
time-homogeneous generator exactly when every \(\eta_\ell>0\) and

\[
(M^{(j)})^{-1}\log\eta\geq0.
\]

The paper gives explicit primitive facets for spins \(1\), \(3/2\),
and \(2\), as well as the instantaneous version for invertible
covariant CP-divisible paths.

## Repository map

- `paper/main.tex`: theorem statements, proofs, low-spin facets,
  geometry, and path criterion.
- `paper/main.pdf`: compiled paper.
- `code/exact_cones.py`: exact 3j/6j reconstruction and facet
  generation.
- `code/generate_plot_data.py`: rational polytope certificates and
  plotting data.
- `paper/generated/geometry_exact.txt`: exact vertices, barycentric
  coordinates, matrices, facets, and volumes.
- `paper/generated/exact_check.txt`: saved direct-reconstruction
  transcript.
- `ROADMAP.md`: theorem-slot acceptance ledger and source map.

## Reproduce

From the repository root:

```sh
make check
make figures
make paper
```

`make check` reconstructs the generator matrices directly from tensor
operators and compares them with the 6j formulas for
\(2j=1,2,3,4\).  All Python entry points use unbuffered output.
