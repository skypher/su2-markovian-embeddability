# Roadmap

## Authoritative theorem

The target is the full classification stated in `AGENTS.md`: general
covariant generators, their exact \(6j\)-symbol eigenmatrix, endpoint
criteria for spins \(1\) and \(3/2\), the invertible CP-divisible path
criterion, and reproducible exact computations and comparison plots.

## Hard-target theorem slots

| Slot | Acceptance check | Status |
|---|---|---|
| G: invariant generators | The canonical GKS matrix is shown to be a nonnegative scalar on each tensor-rank block, its canonical traceless Hamiltonian vanishes, and the converse is proved. | Accepted: Theorem 3.1; GKS Theorem 2.2 and Lemma 2.3 checked in context. |
| M: exact spectral matrix | The action of every \(\mathcal D_k\) on every tensor-rank sector is derived with the fixed tensor convention and checked against exact tensor matrices. | Accepted: Lemma 4.2 and Theorem 4.4, including the closed inverse and \(\lvert\det M^{(j)}\rvert=2j+1\); direct reconstruction, every magnetic-component residual, the tensor Gram matrix, and normalized Choi projectors pass for \(2j=1,2,3,4\). |
| E: endpoint and path criteria | Endpoint positivity, cone membership, low-spin facets, and the invertible differentiable CP-divisibility equivalence on every fixed-sign spectral branch are proved. | Accepted: Theorem 5.1, Corollaries 5.2 and 5.4--5.7, Theorem 7.3, and Corollary 7.4; the local criterion includes the Choi tangent calculation, the path test uses \(d\log|\eta|/dt\), and spin \(2\) endpoint facets are included. |
| C: channel geometry and artifacts | The Choi-simplex vertices, nonnegative-spectrum intersections, exact outputs, tests, and publication figures agree in spins \(1\) and \(3/2\). | Accepted: Propositions 6.1--6.3; rational boundary cycles, determinant terms, the general volume identity, and both figures regenerate from the checked scripts without floating-point combinatorial decisions. |

Overall hard-target meter: **4/4 accepted theorem slots**.  The basis is
the saved exact checker transcript for \(2j=1,2,3,4\), the rational
geometry certificate, and a 17-page PDF build with resolved references and
no typesetting or bibliography warnings.

## Current manuscript target

Submission integration: the PDF uses `Anonymous` as placeholder author
metadata, and journal formatting is not specified.  No unaccepted
mathematical theorem slot remains in the current mission.

## Source ledger

1. Gorini--Kossakowski--Sudarshan, *J. Math. Phys.* **17** (1976),
   Theorem 2.2 and Lemma 2.3: canonical finite-dimensional generator
   form and uniqueness in a fixed traceless orthonormal basis.
2. Breuer, *J. Phys. A* **38** (2005), Appendix A, Eqs. (A1)--(A6), and
   Eq. (2.14) with Appendix B: normalized spherical tensors and the
   \(6j\) recoupling matrix.
3. Aschieri--Ruba--Solovej, *Commun. Math. Phys.* **405** (2024),
   Article 298,
   Sec. 1, Eqs. (1.34)--(1.37), and Appendix A: the equivariant-channel
   simplex and its irreducible Choi vertices.
4. Holevo, *Rep. Math. Phys.* **32** (1993), abstract and main theorem:
   compact-group covariant dynamical semigroups and their generators.
5. Holevo, *J. Funct. Anal.* **131** (1995), abstract and main structure
   theorem: the general covariant dynamical-semigroup framework specialized
   in Theorem 3.1.
6. Cîrstoiu--Korzekwa--Jennings, *Phys. Rev. X* **10** (2020),
   Theorem 7 and Appendix B: symmetry-adapted Lindblad decompositions and
   irreducible tensor operators.
7. Xu--Jagadish, arXiv:2605.23852 (2026), abstract and main semigroup
   classification: recent Weyl-map Markovianity context only.
8. Rudnicki--Klimov--Muñoz--Leuchs--Sánchez-Soto, *Symmetry* **18**
   (2026), abstract: recent SU(2)-covariant channel structure context.
9. Al Nuwairan, *Internat. J. Math.* **25** (2014), Proposition 4.1:
   EPOSIC channels are exactly the extreme covariant channels.
10. Chang--Kim--Kwak--Lee--Youn, *Rev. Math. Phys.* **34** (2022),
   Secs. 2 and 4: Clebsch--Gordan channel simplex and low-rank
   degradability context.
11. Wolf--Eisert--Cubitt--Cirac, *Phys. Rev. Lett.* **101** (2008):
   finite-dimensional channel embedding context and logarithm-branch
   obstruction.
