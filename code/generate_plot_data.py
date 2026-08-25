#!/usr/bin/env python3
"""Generate exact polytope data and sampled semigroup surfaces for PGFPlots."""

from __future__ import annotations

import itertools
import math
from functools import cmp_to_key
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
    return tuple(sorted(points))


def _dot(left: tuple, right: tuple):
    return simplify(sum(a * b for a, b in zip(left, right)))


def _subtract(left: tuple, right: tuple) -> tuple:
    return tuple(simplify(a - b) for a, b in zip(left, right))


def _cross(left: tuple, right: tuple) -> tuple:
    return (
        simplify(left[1] * right[2] - left[2] * right[1]),
        simplify(left[2] * right[0] - left[0] * right[2]),
        simplify(left[0] * right[1] - left[1] * right[0]),
    )


def _ordered_by_exact_angle(points: list[tuple], coordinates) -> list[tuple]:
    """Order planar points counterclockwise using exact orientation tests."""

    def upper_half(coordinate: tuple) -> bool:
        x_value, y_value = coordinate
        return y_value > 0 or (y_value == 0 and x_value >= 0)

    def compare(left: tuple, right: tuple) -> int:
        left_coordinate = coordinates(left)
        right_coordinate = coordinates(right)
        left_half = upper_half(left_coordinate)
        right_half = upper_half(right_coordinate)
        if left_half != right_half:
            return -1 if left_half else 1
        orientation = simplify(
            left_coordinate[0] * right_coordinate[1]
            - left_coordinate[1] * right_coordinate[0]
        )
        if orientation != 0:
            return -1 if orientation > 0 else 1
        left_radius = _dot(left_coordinate, left_coordinate)
        right_radius = _dot(right_coordinate, right_coordinate)
        if left_radius == right_radius:
            return 0
        return -1 if left_radius < right_radius else 1

    return sorted(points, key=cmp_to_key(compare))


def _ordered_face(points: list[tuple], normal: tuple) -> list[tuple]:
    if len(points) <= 2:
        return points
    centroid = tuple(
        simplify(sum(point[index] for point in points) / len(points))
        for index in range(3)
    )
    horizontal_axis = _subtract(points[0], centroid)
    vertical_axis = _cross(normal, horizontal_axis)
    if horizontal_axis == (0, 0, 0) or vertical_axis == (0, 0, 0):
        raise AssertionError("degenerate exact face coordinates")

    def coordinates(point: tuple) -> tuple:
        displacement = _subtract(point, centroid)
        return (
            _dot(displacement, horizontal_axis),
            _dot(displacement, vertical_axis),
        )

    return _ordered_by_exact_angle(points, coordinates)


def _ordered_polygon(points: tuple) -> list[tuple]:
    centroid = tuple(
        simplify(sum(point[index] for point in points) / len(points))
        for index in range(2)
    )

    def coordinates(point: tuple) -> tuple:
        return _subtract(point, centroid)

    return _ordered_by_exact_angle(list(points), coordinates)


def _polygon_area_certificate(points: list[tuple]):
    determinant_terms = []
    for index, point in enumerate(points):
        following = points[(index + 1) % len(points)]
        determinant_terms.append(
            simplify(point[0] * following[1] - point[1] * following[0])
        )
    area = simplify(abs(sum(determinant_terms, Rational(0))) / 2)
    return area, tuple(points), tuple(determinant_terms)


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


def _polytope_volume_certificate(vertices: Matrix, positive_spectrum: bool):
    points = _intersection_vertices(vertices, positive_spectrum)
    faces = _polytope_faces(vertices, positive_spectrum)
    center = Matrix(
        [
            simplify(sum(point[index] for point in points) / len(points))
            for index in range(3)
        ]
    )
    volume = Rational(0)
    determinant_terms = []
    for face in faces:
        first = Matrix(face[0])
        for index in range(1, len(face) - 1):
            second = Matrix(face[index])
            third = Matrix(face[index + 1])
            term = abs(
                Matrix.hstack(first - center, second - center, third - center).det()
            ) / 6
            determinant_terms.append(simplify(term))
            volume += term
    return (
        simplify(volume),
        tuple(center),
        tuple(tuple(face) for face in faces),
        tuple(determinant_terms),
    )


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
                    full_area_data = _polygon_area_certificate(
                        _ordered_polygon(tuple(full_points))
                    )
                    positive_area_data = _polygon_area_certificate(
                        _ordered_polygon(points)
                    )
                    handle.write(f"full_boundary_cycle={full_area_data[1]}\n")
                    handle.write(f"full_determinant_terms={full_area_data[2]}\n")
                    handle.write(
                        f"full_simplex_volume={full_area_data[0]}\n"
                    )
                    handle.write(
                        f"positive_spectrum_boundary_cycle={positive_area_data[1]}\n"
                    )
                    handle.write(
                        f"positive_spectrum_determinant_terms={positive_area_data[2]}\n"
                    )
                    handle.write(
                        f"positive_spectrum_volume={positive_area_data[0]}\n"
                    )
                else:
                    full_volume_data = _polytope_volume_certificate(vertices, False)
                    positive_volume_data = _polytope_volume_certificate(vertices, True)
                    handle.write(f"full_face_center={full_volume_data[1]}\n")
                    handle.write(f"full_face_cycles={full_volume_data[2]}\n")
                    handle.write(f"full_determinant_terms={full_volume_data[3]}\n")
                    handle.write(
                        f"full_simplex_volume={full_volume_data[0]}\n"
                    )
                    handle.write(
                        f"positive_spectrum_face_center={positive_volume_data[1]}\n"
                    )
                    handle.write(
                        f"positive_spectrum_face_cycles={positive_volume_data[2]}\n"
                    )
                    handle.write(
                        f"positive_spectrum_determinant_terms={positive_volume_data[3]}\n"
                    )
                    handle.write(
                        f"positive_spectrum_volume={positive_volume_data[0]}\n"
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
