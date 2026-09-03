# Embedding SU(2)-covariant spin channels into SU(2)-covariant quantum Markov semigroups

This repository contains Leslie P. Polzer's manuscript and exact symbolic
artifacts for covariant Markovian embedding on a fixed irreducible spin
space.

**[Read the current manuscript](paper/main.pdf)**

## Scope and main result

The endpoint problem in this project always requires the generator to be
SU(2)-covariant; it makes no claim about noncovariant logarithms of a covariant
channel.  The covariant-generator classification and Wigner-6j rate transform
are translated from the earlier isotropic-spin literature.  The manuscript
then derives the determinant identity, endpoint cone, exact low-spin geometry,
all-spin volume, and fixed-sign CP-divisibility criterion in one normalization.

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
and \(2\).  An arbitrary invertible differentiable covariant path is
CP-divisible exactly when

\[
(M^{(j)})^{-1}\frac{d}{dt}\log|\eta(t)|\geq0,
\]

where the absolute value is componentwise.  This includes fixed
negative-spectrum branches; for paths beginning at the identity the
absolute values may be omitted.

In \(n=2j\) sector-eigenvalue coordinates, the embeddable region has
the exact Euclidean volume

\[
\operatorname{vol}_n=\frac{2j+1}{2^n n!}.
\]

## Repository map

- `paper/main.tex`: theorem statements, proofs, low-spin facets,
  geometry, and path criterion.
- `paper/main.pdf`: compiled paper.
- `code/exact_cones.py`: exact 3j/6j reconstruction, every-component
  tensor residuals, Choi normalization checks, and facet generation.
- `code/generate_plot_data.py`: rational polytope certificates and
  plotting data.
- `paper/generated/geometry_exact.txt`: exact vertices, barycentric
  coordinates, matrices, facets, boundary cycles, determinant terms, and
  volumes.
- `paper/generated/exact_check.txt`: saved direct-reconstruction
  transcript.
- `ROADMAP.md`: theorem-slot acceptance ledger and source map.

The initial public manuscript release is tagged `v1.0.0`; subsequent revisions
are recorded by the repository history.

## Reproduce

From the repository root:

```sh
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r requirements.txt
make check
make figures
make paper
```

Both executable entry points support `-h` and `--help` without generating or
changing artifacts.

`make check` reconstructs the generator matrices directly from tensor
operators, checks every magnetic component and each normalized Choi
projector, and compares the results with the 6j formulas for
\(2j=1,2,3,4\).  All Python entry points use unbuffered output.

The checked environment is Python 3.12.3, SymPy 1.12, mpmath 1.2.1,
and TeX Live 2023 with PGFPlots compatibility level 1.18.  The Python
package versions are pinned in `requirements.txt`.

## Citation and reuse

Machine-readable citation metadata are in `CITATION.cff`.  Repository reuse
terms are stated in `LICENSE`.

## Contact

Leslie P. Polzer<br>
Independent Researcher<br>
[polzer@fastmail.com](mailto:polzer@fastmail.com)
