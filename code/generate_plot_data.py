#!/usr/bin/env python3
"""Generate exact polytope data and sampled semigroup surfaces for PGFPlots."""

from __future__ import annotations

import itertools
import math
from pathlib import Path

from sympy import Matrix, Rational, simplify

from exact_cones import (
    generator_matrix,
    inverse_generator_matrix,
    primitive_facet_rows,
    vertex_matrix,
)


ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "paper" / "generated"


def _halfspaces(vertices: Matrix, positive_spectrum: bool) -> list[tuple]:
    """Return affine forms a.dot(x)+b >= 0 for a simplex and x_i >= 0."""
    dimension = vertices.rows
    augmented = vertices.col_join(Matrix.ones(1, dimension + 1))
    inverse = augmented.inv()
    halfspaces = []
    for row in inverse.tolist():
        halfspaces.append((tuple(row[:dimension]), row[dimension], "simplex"))
    if positive_spectrum:
        for coordinate in range(dimension):
            normal = [Rational(0) for _ in range(dimension)]
            normal[coordinate] = Rational(1)
            halfspaces.append((tuple(normal), Rational(0), "coordinate"))
    return halfspaces


def _evaluate(halfspace: tuple, point: tuple):
    normal, offset, _ = halfspace
    return simplify(sum(a * x for a, x in zip(normal, point)) + offset)


def _intersection_vertices(vertices: Matrix, positive_spectrum: bool) -> tuple:
    dimension = vertices.rows
    halfspaces = _halfspaces(vertices, positive_spectrum)
    points = set()
    for active in itertools.combinations(range(len(halfspaces)), dimension):
        coefficient = Matrix([halfspaces[index][0] for index in active])
        if coefficient.det() == 0:
            continue
        rhs = Matrix([-halfspaces[index][1] for index in active])
        solution = coefficient.inv() * rhs
        point = tuple(simplify(entry) for entry in solution)
        if all(_evaluate(halfspace, point) >= 0 for halfspace in halfspaces):
            points.add(point)
    return tuple(sorted(points, key=lambda point: tuple(float(x) for x in point)))


def _ordered_face(points: list[tuple], normal: tuple) -> list[tuple]:
    if len(points) <= 2:
        return points
    centroid = [
        sum(float(point[index]) for point in points) / len(points)
        for index in range(3)
    ]
    normal_float = [float(value) for value in normal]
    normal_norm = math.sqrt(sum(value * value for value in normal_float))
    normal_unit = [value / normal_norm for value in normal_float]

    first = [
        float(points[0][index]) - centroid[index]
        for index in range(3)
    ]
    first_norm = math.sqrt(sum(value * value for value in first))
    u = [value / first_norm for value in first]
    v = [
        normal_unit[1] * u[2] - normal_unit[2] * u[1],
        normal_unit[2] * u[0] - normal_unit[0] * u[2],
        normal_unit[0] * u[1] - normal_unit[1] * u[0],
    ]

    def angle(point: tuple) -> float:
        displacement = [
            float(point[index]) - centroid[index]
            for index in range(3)
        ]
        horizontal = sum(displacement[index] * u[index] for index in range(3))
        vertical = sum(displacement[index] * v[index] for index in range(3))
        return math.atan2(vertical, horizontal)

    return sorted(points, key=angle)


def _ordered_polygon(points: tuple) -> list[tuple]:
    centroid = [
        sum(float(point[index]) for point in points) / len(points)
        for index in range(2)
    ]
    return sorted(
        points,
        key=lambda point: math.atan2(
            float(point[1]) - centroid[1],
            float(point[0]) - centroid[0],
        ),
    )


def _polygon_area(points: list[tuple]):
    twice_area = Rational(0)
    for index, point in enumerate(points):
        following = points[(index + 1) % len(points)]
        twice_area += point[0] * following[1] - point[1] * following[0]
    return simplify(abs(twice_area) / 2)


def _polytope_faces(vertices: Matrix, positive_spectrum: bool) -> list[list[tuple]]:
    dimension = vertices.rows
    if dimension != 3:
        raise ValueError("3D face output is defined only for spin 3/2")
    points = _intersection_vertices(vertices, positive_spectrum)
    faces = []
    for halfspace in _halfspaces(vertices, positive_spectrum):
        active = [point for point in points if _evaluate(halfspace, point) == 0]
        if len(active) >= 3:
            faces.append(_ordered_face(active, halfspace[0]))
    return faces


def _polytope_volume(vertices: Matrix, positive_spectrum: bool):
    points = _intersection_vertices(vertices, positive_spectrum)
    faces = _polytope_faces(vertices, positive_spectrum)
    center = Matrix(
        [
            simplify(sum(point[index] for point in points) / len(points))
            for index in range(3)
        ]
    )
    volume = Rational(0)
    for face in faces:
        first = Matrix(face[0])
        for index in range(1, len(face) - 1):
            second = Matrix(face[index])
            third = Matrix(face[index + 1])
            volume += abs(
                Matrix.hstack(first - center, second - center, third - center).det()
            ) / 6
    return simplify(volume)


def _embeddable_volume(matrix: Matrix):
    column_decays = [-sum(matrix[:, column]) for column in range(matrix.cols)]
    volume = abs(matrix.det())
    for decay in column_decays:
        volume /= decay
    return simplify(volume), column_decays


def _number(value) -> str:
    return f"{float(value):.12g}"


def _write_polygon(path: Path, points: list[tuple]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for point in points:
            handle.write(" ".join(_number(value) for value in point) + "\n")


def _write_faces(path: Path, faces: list[list[tuple]], style: str) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for face in faces:
            handle.write(f"\\addplot3[{style}] coordinates {{\n")
            for point in face:
                handle.write(
                    "  (" + ",".join(_number(value) for value in point) + ")\n"
                )
            handle.write("} \\closedcycle;\n")


def _write_markov_surface(
    path: Path,
    matrix: Matrix,
    fixed_rate: int,
    grid_size: int = 33,
    maximum_rate: float = 8.0,
) -> None:
    varying = [index for index in range(3) if index != fixed_rate]
    with path.open("w", encoding="utf-8") as handle:
        for row in range(grid_size):
            rate_a = maximum_rate * row / (grid_size - 1)
            for column in range(grid_size):
                rate_b = maximum_rate * column / (grid_size - 1)
                rates = [0.0, 0.0, 0.0]
                rates[varying[0]] = rate_a
                rates[varying[1]] = rate_b
                logarithms = [
                    sum(float(matrix[ell, k]) * rates[k] for k in range(3))
                    for ell in range(3)
                ]
                handle.write(
                    " ".join(f"{math.exp(value):.12g}" for value in logarithms)
                    + "\n"
                )
            if row + 1 < grid_size:
                handle.write("\n")


def _write_exact_summary(path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for twice_j in (2, 3, 4):
            vertices = vertex_matrix(twice_j)
            handle.write(f"2j={twice_j}\n")
            handle.write(f"simplex_vertices={vertices}\n")
            handle.write(f"generator_matrix={generator_matrix(twice_j)}\n")
            handle.write(
                f"inverse_generator_matrix={inverse_generator_matrix(twice_j)}\n"
            )
            handle.write(f"primitive_facets={primitive_facet_rows(twice_j)}\n")
            volume, column_decays = _embeddable_volume(
                generator_matrix(twice_j)
            )
            handle.write(f"column_decays={column_decays}\n")
            handle.write(f"embeddable_eigenvalue_volume={volume}\n")
            if twice_j in (2, 3):
                points = _intersection_vertices(vertices, True)
                augmented_inverse = vertices.col_join(
                    Matrix.ones(1, twice_j + 1)
                ).inv()
                handle.write(f"positive_spectrum_vertices={points}\n")
                handle.write("positive_vertex_barycentric_coordinates=\n")
                for point in points:
                    barycentric = augmented_inverse * Matrix([*point, 1])
                    handle.write(f"  {point}: {tuple(barycentric)}\n")
                if twice_j == 2:
                    full_points = [tuple(vertices[:, index]) for index in range(3)]
                    handle.write(
                        f"full_simplex_volume={_polygon_area(_ordered_polygon(tuple(full_points)))}\n"
                    )
                    handle.write(
                        f"positive_spectrum_volume={_polygon_area(_ordered_polygon(points))}\n"
                    )
                else:
                    handle.write(
                        f"full_simplex_volume={_polytope_volume(vertices, False)}\n"
                    )
                    handle.write(
                        f"positive_spectrum_volume={_polytope_volume(vertices, True)}\n"
                    )
            if twice_j != 4:
                handle.write("\n")


def main() -> None:
    GENERATED.mkdir(parents=True, exist_ok=True)

    spin_one_vertices = vertex_matrix(2)
    positive_polygon = _ordered_polygon(
        _intersection_vertices(spin_one_vertices, True)
    )
    _write_polygon(GENERATED / "spin1_positive_polygon.dat", positive_polygon)

    markov_polygon = []
    samples = 401
    for index in range(samples):
        eta_one = index / (samples - 1)
        markov_polygon.append((eta_one, eta_one**3))
    for index in reversed(range(samples)):
        eta_one = index / (samples - 1)
        markov_polygon.append((eta_one, eta_one ** Rational(3, 5)))
    _write_polygon(GENERATED / "spin1_markov_polygon.dat", markov_polygon)

    spin_three_vertices = vertex_matrix(3)
    full_faces = _polytope_faces(spin_three_vertices, False)
    positive_faces = _polytope_faces(spin_three_vertices, True)
    _write_faces(
        GENERATED / "spin3_full_faces.tex",
        full_faces,
        "draw=gray!80,fill=gray!35,fill opacity=0.22,line width=0.35pt",
    )
    _write_faces(
        GENERATED / "spin3_positive_faces.tex",
        positive_faces,
        "draw=orange!85!black,fill=orange!55,fill opacity=0.28,line width=0.35pt",
    )

    matrix = generator_matrix(3)
    for fixed_rate in range(3):
        _write_markov_surface(
            GENERATED / f"spin3_markov_face_{fixed_rate + 1}.dat",
            matrix,
            fixed_rate,
        )

    _write_exact_summary(GENERATED / "geometry_exact.txt")
    print(
        "generated spin-1 polygons, spin-3/2 polyhedra, and "
        "three Markovian boundary meshes",
        flush=True,
    )


if __name__ == "__main__":
    main()
