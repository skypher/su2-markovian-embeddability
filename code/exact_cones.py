#!/usr/bin/env python3
"""Exact SU(2)-covariant channel and generator eigenvalue matrices."""

from __future__ import annotations

import argparse
from functools import reduce
from math import gcd

from sympy import Matrix, Rational, ilcm, simplify, sqrt, zeros
from sympy.physics.wigner import clebsch_gordan, wigner_3j, wigner_6j


def _spin(twice_j: int) -> Rational:
    if twice_j < 1:
        raise ValueError("twice_j must be a positive integer")
    return Rational(twice_j, 2)


def vertex_matrix(twice_j: int) -> Matrix:
    """Rows ell=1,...,2j; columns k=0,...,2j."""
    j = _spin(twice_j)
    d = twice_j + 1
    return Matrix(
        twice_j,
        twice_j + 1,
        lambda ell0, k: simplify(
            d
            * (-1) ** (twice_j + k + ell0 + 1)
            * wigner_6j(j, j, k, j, j, ell0 + 1)
        ),
    )


def generator_matrix(twice_j: int) -> Matrix:
    """Rows ell=1,...,2j; columns k=1,...,2j."""
    j = _spin(twice_j)
    d = twice_j + 1
    return Matrix(
        twice_j,
        twice_j,
        lambda ell0, k0: simplify(
            (2 * (k0 + 1) + 1)
            * (
                (-1) ** (twice_j + k0 + 1 + ell0 + 1)
                * wigner_6j(j, j, k0 + 1, j, j, ell0 + 1)
                - Rational(1, d)
            )
        ),
    )


def inverse_generator_matrix(twice_j: int) -> Matrix:
    """Closed 6j formula for M^(-1); rows k and columns ell."""
    j = _spin(twice_j)
    return Matrix(
        twice_j,
        twice_j,
        lambda k0, ell0: simplify(
            (2 * (ell0 + 1) + 1)
            * (-1) ** (twice_j + k0 + 1 + ell0 + 1)
            * wigner_6j(j, j, k0 + 1, j, j, ell0 + 1)
        ),
    )


def tensor_operator(twice_j: int, k: int, q: int) -> Matrix:
    """Hilbert--Schmidt normalized tensor from the fixed CG convention."""
    j = _spin(twice_j)
    d = twice_j + 1
    magnetic = [-j + r for r in range(d)]
    tensor = zeros(d)
    for row, m_out in enumerate(magnetic):
        for col, m_in in enumerate(magnetic):
            tensor[row, col] = simplify(
                sqrt(Rational(2 * k + 1, d))
                * clebsch_gordan(j, k, j, m_in, q, m_out)
            )
    return tensor


def tensor_operator_3j(twice_j: int, k: int, q: int) -> Matrix:
    """The paper's tensor convention, Breuer Appendix A, Eq. (A3)."""
    j = _spin(twice_j)
    d = twice_j + 1
    magnetic = [-j + r for r in range(d)]
    return Matrix(
        d,
        d,
        lambda row, col: simplify(
            sqrt(2 * k + 1)
            * (-1) ** (j - magnetic[row])
            * wigner_3j(
                j,
                j,
                k,
                magnetic[row],
                -magnetic[col],
                -q,
            )
        ),
    )


def dual_intertwiner(twice_j: int) -> Matrix:
    """Matrix of beta(<j,m|)=(-1)^(j-m)|j,-m>."""
    j = _spin(twice_j)
    d = twice_j + 1
    magnetic = [-j + r for r in range(d)]
    beta = zeros(d)
    for row, m_out in enumerate(magnetic):
        for col, m_in in enumerate(magnetic):
            if m_out == -m_in:
                beta[row, col] = (-1) ** (j - m_in)
    return beta


def column_vectorize(matrix: Matrix) -> Matrix:
    """Vectorization in the paper's output--dual-input tensor order."""
    rows = matrix.rows
    columns = matrix.cols
    return Matrix(
        rows * columns,
        1,
        lambda index, _column: matrix[index // columns, index % columns],
    )


def direct_generator_matrix(twice_j: int) -> Matrix:
    """Compute M directly from exact tensor matrices, without a 6j symbol."""
    d = twice_j + 1
    result = zeros(twice_j)
    for ell in range(1, twice_j + 1):
        probe = tensor_operator(twice_j, ell, 0)
        for k in range(1, twice_j + 1):
            image = zeros(d)
            for q in range(-k, k + 1):
                tensor = tensor_operator(twice_j, k, q)
                image += tensor * probe * tensor.conjugate().T
            image -= Rational(2 * k + 1, d) * probe
            scalar = simplify((probe.conjugate().T * image).trace())
            residual = (image - scalar * probe).applyfunc(simplify)
            if residual != zeros(d):
                raise AssertionError(
                    f"direct Kraus image is not scalar on the probe for "
                    f"twice_j={twice_j}, ell={ell}, k={k}:\n{residual}"
                )
            result[ell - 1, k - 1] = scalar
    return result


def primitive_facet_rows(twice_j: int) -> list[list[int]]:
    """Primitive integer rows A with embeddability equivalent to A log(eta)>=0."""
    inverse = inverse_generator_matrix(twice_j)
    rows: list[list[int]] = []
    for row in inverse.tolist():
        denominator = reduce(ilcm, [entry.q for entry in row], 1)
        integers = [int(entry * denominator) for entry in row]
        divisor = reduce(gcd, [abs(value) for value in integers if value])
        rows.append([value // divisor for value in integers])
    return rows


def verify(twice_j: int) -> None:
    d = twice_j + 1
    beta = dual_intertwiner(twice_j)
    if beta * beta.conjugate().T != Matrix.eye(d):
        raise AssertionError(f"dual intertwiner is not unitary for twice_j={twice_j}")
    tensors: dict[tuple[int, int], Matrix] = {}
    for k in range(twice_j + 1):
        for q in range(-k, k + 1):
            tensor = tensor_operator(twice_j, k, q)
            tensors[k, q] = tensor
            tensor_3j = tensor_operator_3j(twice_j, k, q)
            if tensor != tensor_3j:
                raise AssertionError(
                    f"CG and 3j tensor conventions disagree for "
                    f"twice_j={twice_j}, k={k}, q={q}"
                )
            dual_tensor = beta * tensor.conjugate() * beta.conjugate().T
            if dual_tensor != (-1) ** k * tensor.conjugate().T:
                raise AssertionError(
                    f"dual tensor phase failed for "
                    f"twice_j={twice_j}, k={k}, q={q}"
                )

    tensor_labels = list(tensors)
    tensor_gram = Matrix(
        d * d,
        d * d,
        lambda row, column: simplify(
            (
                tensors[tensor_labels[row]].conjugate().T
                * tensors[tensor_labels[column]]
            ).trace()
        ),
    )
    if tensor_gram != Matrix.eye(d * d):
        raise AssertionError(
            f"tensor basis is not Hilbert--Schmidt orthonormal for "
            f"twice_j={twice_j}"
        )

    for k in range(twice_j + 1):
        scalar_sum = zeros(d)
        for q in range(-k, k + 1):
            tensor = tensors[k, q]
            scalar_sum += tensor.conjugate().T * tensor
        expected_sum = Rational(2 * k + 1, d) * Matrix.eye(d)
        if (scalar_sum - expected_sum).applyfunc(simplify) != zeros(d):
            raise AssertionError(
                f"scalar Kraus sum failed for twice_j={twice_j}, k={k}"
            )

    formula = generator_matrix(twice_j)
    direct = direct_generator_matrix(twice_j)
    if formula != direct:
        raise AssertionError(
            f"6j and tensor definitions disagree for twice_j={twice_j}:\n"
            f"formula={formula}\ndirect={direct}"
        )

    for ell in range(twice_j + 1):
        for m in range(-ell, ell + 1):
            probe = tensors[ell, m]
            for k in range(1, twice_j + 1):
                image = zeros(d)
                for q in range(-k, k + 1):
                    tensor = tensors[k, q]
                    image += tensor * probe * tensor.conjugate().T
                image -= Rational(2 * k + 1, d) * probe
                scalar = 0 if ell == 0 else formula[ell - 1, k - 1]
                residual = (image - scalar * probe).applyfunc(simplify)
                if residual != zeros(d):
                    raise AssertionError(
                        f"full tensor residual failed for twice_j={twice_j}, "
                        f"ell={ell}, m={m}, k={k}:\n{residual}"
                    )

    for k in range(twice_j + 1):
        projector = zeros(d * d)
        choi_direct = zeros(d * d)
        for q in range(-k, k + 1):
            tensor = tensors[k, q]
            vector = column_vectorize(tensor)
            projector += vector * vector.conjugate().T
        for a in range(d):
            for b in range(d):
                matrix_unit = zeros(d)
                matrix_unit[a, b] = 1
                output = zeros(d)
                for q in range(-k, k + 1):
                    tensor = tensors[k, q]
                    output += tensor * matrix_unit * tensor.conjugate().T
                output *= Rational(d, 2 * k + 1)
                for row in range(d):
                    for column in range(d):
                        choi_direct[row * d + a, column * d + b] = output[
                            row, column
                        ]
        if (projector * projector - projector).applyfunc(simplify) != zeros(
            d * d
        ):
            raise AssertionError(
                f"Choi support is not a projector for twice_j={twice_j}, k={k}"
            )
        choi_projector = Rational(d, 2 * k + 1) * projector
        if (choi_direct - choi_projector).applyfunc(simplify) != zeros(d * d):
            raise AssertionError(
                f"direct and projector Choi matrices disagree for "
                f"twice_j={twice_j}, k={k}"
            )
        if simplify(choi_direct.trace()) != d:
            raise AssertionError(
                f"Choi trace normalization failed for twice_j={twice_j}, k={k}"
            )
        output_partial_trace = Matrix(
            d,
            d,
            lambda a, b: simplify(
                sum(
                    choi_direct[output * d + a, output * d + b]
                    for output in range(d)
                )
            ),
        )
        if output_partial_trace != Matrix.eye(d):
            raise AssertionError(
                f"trace preservation failed for twice_j={twice_j}, k={k}"
            )

    inverse_formula = inverse_generator_matrix(twice_j)
    if formula * inverse_formula != Matrix.eye(twice_j):
        raise AssertionError(
            f"closed inverse formula failed for twice_j={twice_j}"
        )
    if abs(formula.det()) != twice_j + 1:
        raise AssertionError(
            f"determinant magnitude failed for twice_j={twice_j}"
        )
    vertices = vertex_matrix(twice_j)
    for k0 in range(twice_j):
        k = k0 + 1
        expected = Rational(2 * k + 1, twice_j + 1) * (
            vertices[:, k] - Matrix.ones(twice_j, 1)
        )
        if formula[:, k0] != expected:
            raise AssertionError(f"vertex relation failed for k={k}")
    column_decays = [
        simplify(-sum(formula[ell, k] for ell in range(twice_j)))
        for k in range(twice_j)
    ]
    if column_decays != [2 * k for k in range(1, twice_j + 1)]:
        raise AssertionError(
            f"column-decay identity failed for twice_j={twice_j}: "
            f"{column_decays}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "twice_j",
        nargs="+",
        type=int,
        help="one or more values of 2j",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="also reconstruct the matrix from Clebsch--Gordan tensors",
    )
    args = parser.parse_args()

    for twice_j in args.twice_j:
        matrix = generator_matrix(twice_j)
        print(f"2j = {twice_j}", flush=True)
        print(f"channel vertices (columns k=0,...,2j):\n{vertex_matrix(twice_j)}", flush=True)
        print(f"generator matrix M:\n{matrix}", flush=True)
        print(f"det(M) = {matrix.det()}", flush=True)
        print(f"M^(-1):\n{inverse_generator_matrix(twice_j)}", flush=True)
        print(f"primitive facet rows: {primitive_facet_rows(twice_j)}", flush=True)
        if args.verify:
            verify(twice_j)
            print(
                "exact every-m tensor reconstruction and Choi checks: PASS",
                flush=True,
            )


if __name__ == "__main__":
    main()
