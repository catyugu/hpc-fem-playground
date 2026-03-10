#!/usr/bin/env python3
"""
Convert COMSOL-style mphtxt mesh to MFEM v1.0 .mesh.

This converter focuses on robust topology handling for mixed 3D meshes used
in this repository (tet/hex/pyr with tri/quad boundaries). It keeps COMSOL
volume-element local ordering and only aligns explicit boundary face ordering
to adjacent volume faces when an exact face match exists.
"""

import sys
import os
import math
import itertools
from pathlib import Path
from collections import defaultdict


def parse_mphtxt(path):
    with open(path, "r", encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle]

    index = 0
    nlines = len(lines)
    sdim = None

    vertices = []
    while index < nlines:
        if sdim is None and lines[index].strip().endswith("# sdim"):
            try:
                sdim = int(lines[index].split("#")[0].strip())
            except Exception:
                sdim = None

        if lines[index].strip().endswith("# number of mesh vertices"):
            prefix = lines[index].split("#")[0].strip()
            if prefix:
                nv = int(prefix)
                index += 1
            else:
                nv = int(lines[index - 1].strip())
                index += 1

            while index < nlines and "# Mesh vertex coordinates" not in lines[index]:
                index += 1
            if index >= nlines:
                raise RuntimeError("Vertex coordinate section not found")
            index += 1

            while len(vertices) < nv and index < nlines:
                stripped = lines[index].strip()
                index += 1
                if stripped == "" or stripped.startswith("#"):
                    continue
                tokens = stripped.split()
                if len(tokens) >= 3:
                    x, y, z = map(float, tokens[:3])
                elif len(tokens) == 2:
                    x, y = map(float, tokens)
                    z = 0.0
                else:
                    x = float(tokens[0])
                    while index < nlines and lines[index].strip() == "":
                        index += 1
                    y = float(lines[index].split()[0])
                    index += 1
                    while index < nlines and lines[index].strip() == "":
                        index += 1
                    z = float(lines[index].split()[0])
                    index += 1
                vertices.append((x, y, z))
            break
        index += 1

    if not vertices:
        raise RuntimeError("No vertices parsed from mphtxt")

    blocks = []
    while index < nlines:
        if lines[index].strip().startswith("# Type"):
            index += 1
            while index < nlines and lines[index].strip() == "":
                index += 1
            if index >= nlines:
                break

            parts = lines[index].split("#")[0].split()
            type_name = parts[1].strip() if len(parts) >= 2 else parts[0].strip()
            index += 1

            while index < nlines and lines[index].strip() == "":
                index += 1
            if index < nlines and "# number of vertices per element" in lines[index]:
                nverts = int(lines[index].split("#")[0])
                index += 1
            else:
                nverts = int(lines[index].split()[0])
                index += 1

            while index < nlines and lines[index].strip() == "":
                index += 1
            if index < nlines and "# number of elements" in lines[index]:
                nelem = int(lines[index].split("#")[0])
                index += 1
            else:
                nelem = int(lines[index].split()[0])
                index += 1

            while index < nlines and "# Elements" not in lines[index]:
                index += 1
            if index < nlines:
                index += 1

            elems = []
            while len(elems) < nelem and index < nlines:
                stripped = lines[index].strip()
                index += 1
                if stripped == "" or stripped.startswith("#"):
                    continue
                tokens = stripped.split()
                while len(tokens) < nverts and index < nlines:
                    extra = lines[index].strip()
                    index += 1
                    if extra == "" or extra.startswith("#"):
                        continue
                    tokens.extend(extra.split())
                elems.append([int(token) for token in tokens[:nverts]])

            geom_indices = None
            look = index
            while look < nlines and lines[look].strip() == "":
                look += 1
            if look < nlines and lines[look].strip().endswith("# number of geometric entity indices"):
                count = int(lines[look].split("#")[0])
                look += 1
                while look < nlines and lines[look].strip() == "":
                    look += 1
                geom_indices = []
                while len(geom_indices) < count and look < nlines:
                    stripped = lines[look].strip()
                    look += 1
                    if stripped == "" or stripped.startswith("#"):
                        continue
                    for token in stripped.split():
                        if len(geom_indices) < count:
                            geom_indices.append(int(token))
                index = look

            blocks.append(
                {
                    "type_name": type_name,
                    "nverts": nverts,
                    "elems": elems,
                    "geom_indices": geom_indices,
                }
            )
        else:
            index += 1

    return vertices, blocks, sdim


def mfem_geom_from_name(name, nverts, sdim=None):
    lowered = (name or "").lower()
    if lowered.startswith("vtx") or nverts == 1:
        return 0
    if lowered.startswith("edg") or lowered.startswith("line") or nverts == 2:
        return 1
    if lowered.startswith("tri") or nverts == 3:
        return 2
    if lowered.startswith("quad"):
        return 3
    if lowered.startswith("tet"):
        return 4
    if lowered.startswith("hex"):
        return 5
    if lowered.startswith("prism") or lowered.startswith("wedge"):
        return 6
    if lowered.startswith("pyr"):
        return 7

    if nverts == 4 and sdim == 3:
        return 4
    if nverts == 4 and sdim == 2:
        return 3
    if nverts == 5:
        return 7
    if nverts == 6:
        return 6
    if nverts == 8:
        return 5

    raise RuntimeError(f"Cannot map type '{name}' with {nverts} vertices")


def element_faces(geom, vertices):
    v = vertices
    if geom == 4:
        return [[v[0], v[2], v[1]], [v[0], v[1], v[3]], [v[1], v[2], v[3]], [v[2], v[0], v[3]]]
    if geom == 5:
        return [
            [v[0], v[1], v[2], v[3]],
            [v[4], v[7], v[6], v[5]],
            [v[0], v[3], v[7], v[4]],
            [v[1], v[5], v[6], v[2]],
            [v[0], v[4], v[5], v[1]],
            [v[3], v[2], v[6], v[7]],
        ]
    if geom == 6:
        return [
            [v[0], v[2], v[1]],
            [v[3], v[4], v[5]],
            [v[0], v[1], v[4], v[3]],
            [v[1], v[2], v[5], v[4]],
            [v[2], v[0], v[3], v[5]],
        ]
    if geom == 7:
        return [
            [v[0], v[1], v[2], v[3]],
            [v[0], v[1], v[4]],
            [v[1], v[2], v[4]],
            [v[2], v[3], v[4]],
            [v[3], v[0], v[4]],
        ]
    return []


def vector_sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vector_dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def vector_cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def vector_norm(a):
    return math.sqrt(vector_dot(a, a))


def detect_hex_face_sets(local_points):
    local_indices = list(range(8))
    faces = []
    plane_tolerance = 1e-9

    for combination in itertools.combinations(local_indices, 4):
        c0, c1, c2, c3 = combination
        p0 = local_points[c0]
        p1 = local_points[c1]
        p2 = local_points[c2]
        p3 = local_points[c3]

        normal = vector_cross(vector_sub(p1, p0), vector_sub(p2, p0))
        if vector_norm(normal) < plane_tolerance:
            continue
        if abs(vector_dot(normal, vector_sub(p3, p0))) > plane_tolerance:
            continue

        other_indices = [index for index in local_indices if index not in combination]
        side_values = [vector_dot(normal, vector_sub(local_points[index], p0)) for index in other_indices]
        if all(value >= -plane_tolerance for value in side_values) or all(value <= plane_tolerance for value in side_values):
            face = tuple(sorted(combination))
            if face not in faces:
                faces.append(face)

    return faces


def reorder_hex_from_geometry(vertices, all_nodes):
    local_points = [all_nodes[index] for index in vertices]
    face_sets = detect_hex_face_sets(local_points)
    if len(face_sets) != 6:
        return vertices

    neighbors = {index: set() for index in range(8)}
    for i, j in itertools.combinations(range(8), 2):
        shared_face_count = sum(1 for face in face_sets if i in face and j in face)
        if shared_face_count == 2:
            neighbors[i].add(j)
            neighbors[j].add(i)

    # Use deterministic sorting: coordinate sum as primary, index as secondary
    def sort_key(index):
        pt = local_points[index]
        return (sum(pt), index)

    base = min(range(8), key=sort_key)
    base_neighbors = list(neighbors[base])
    if len(base_neighbors) != 3:
        return vertices

    best_triplet = None
    best_score = -1e300
    for a, b, c in itertools.permutations(base_neighbors, 3):
        va = vector_sub(local_points[a], local_points[base])
        vb = vector_sub(local_points[b], local_points[base])
        vc = vector_sub(local_points[c], local_points[base])
        score = vector_dot(vector_cross(va, vb), vc)
        if score > best_score:
            best_score = score
            best_triplet = (a, b, c)

    if best_triplet is None:
        return vertices

    v1, v3, v4 = best_triplet

    def common_neighbor(i, j, excluded):
        candidates = (neighbors[i].intersection(neighbors[j])) - set(excluded)
        if not candidates:
            return None
        return min(candidates)  # Deterministic selection

    v2 = common_neighbor(v1, v3, [base, v4])
    if v2 is None:
        return vertices
    v5 = common_neighbor(v1, v4, [base, v3, v2])
    if v5 is None:
        return vertices
    v7 = common_neighbor(v3, v4, [base, v1, v2, v5])
    if v7 is None:
        return vertices

    v6_candidates = (neighbors[v2].intersection(neighbors[v5]).intersection(neighbors[v7])) - {base, v1, v2, v3, v4, v5, v7}
    if not v6_candidates:
        return vertices
    v6 = min(v6_candidates)  # Deterministic selection

    local_order = [base, v1, v2, v3, v4, v5, v6, v7]
    if len(set(local_order)) != 8:
        return vertices

    return [vertices[index] for index in local_order]


def is_dihedral_quad(base, test):
    if len(base) != 4 or len(test) != 4:
        return False

    for shift in range(4):
        if all(test[i] == base[(i + shift) % 4] for i in range(4)):
            return True

    reversed_base = [base[0], base[3], base[2], base[1]]
    for shift in range(4):
        if all(test[i] == reversed_base[(i + shift) % 4] for i in range(4)):
            return True
    return False


def reorder_pyramid_from_geometry(vertices, all_nodes):
    if len(vertices) != 5:
        return vertices

    local_points = [all_nodes[index] for index in vertices]
    tolerance = 1e-9

    apex = None
    base = None
    for candidate in range(5):
        base_indices = [index for index in range(5) if index != candidate]
        p0 = local_points[base_indices[0]]
        p1 = local_points[base_indices[1]]
        p2 = local_points[base_indices[2]]
        p3 = local_points[base_indices[3]]
        normal = vector_cross(vector_sub(p1, p0), vector_sub(p2, p0))
        if vector_norm(normal) < tolerance:
            continue
        if abs(vector_dot(normal, vector_sub(p3, p0))) > tolerance:
            continue
        apex = candidate
        base = base_indices
        break

    if apex is None or base is None:
        return vertices

    base_center = [0.0, 0.0, 0.0]
    for index in base:
        point = local_points[index]
        base_center[0] += point[0]
        base_center[1] += point[1]
        base_center[2] += point[2]
    base_center = (base_center[0] / 4.0, base_center[1] / 4.0, base_center[2] / 4.0)

    normal = vector_cross(
        vector_sub(local_points[base[1]], local_points[base[0]]),
        vector_sub(local_points[base[2]], local_points[base[0]]),
    )
    if vector_norm(normal) < tolerance:
        return vertices

    if vector_dot(normal, vector_sub(local_points[apex], local_points[base[0]])) > 0.0:
        normal = (-normal[0], -normal[1], -normal[2])

    axis_u = vector_sub(local_points[base[0]], base_center)
    if vector_norm(axis_u) < tolerance:
        return vertices
    axis_v = vector_cross(normal, axis_u)
    if vector_norm(axis_v) < tolerance:
        return vertices

    # Use deterministic sorting: include original index as tiebreaker
    base_angles = []
    for index in base:
        direction = vector_sub(local_points[index], base_center)
        angle = math.atan2(vector_dot(direction, axis_v), vector_dot(direction, axis_u))
        base_angles.append((angle, index))
    base_angles.sort()

    ordered_base = [index for _, index in base_angles]
    reordered_local = ordered_base + [apex]
    return [vertices[index] for index in reordered_local]


def pyramid_variants(vertices):
    base = vertices[:4]
    apex = vertices[4]
    variants = []

    for shift in range(4):
        rotated = base[shift:] + base[:shift]
        variants.append(rotated + [apex])

    reversed_base = [base[0], base[3], base[2], base[1]]
    for shift in range(4):
        rotated = reversed_base[shift:] + reversed_base[:shift]
        variants.append(rotated + [apex])

    unique = []
    seen = set()
    for variant in variants:
        key = tuple(variant)
        if key in seen:
            continue
        seen.add(key)
        unique.append(variant)
    return unique


def write_mfem_mesh(nodes, blocks, outpath, sdim=None):
    hex_perm = None
    raw_hex_perm = os.environ.get("MPHTXT_HEX_PERM", "").strip()
    if raw_hex_perm:
        tokens = raw_hex_perm.replace(",", " ").split()
        if len(tokens) != 8:
            raise RuntimeError("MPHTXT_HEX_PERM must contain exactly 8 integers")
        hex_perm = [int(token) for token in tokens]
        if sorted(hex_perm) != list(range(8)):
            raise RuntimeError("MPHTXT_HEX_PERM must be a permutation of 0..7")

    if sdim is None:
        sdim = 3

    volume_types = {4, 5, 6, 7} if sdim == 3 else ({2, 3} if sdim == 2 else {1})
    boundary_types = {2, 3} if sdim == 3 else ({1} if sdim == 2 else {0})

    elements = []
    explicit_boundaries = []

    for block in blocks:
        geom = mfem_geom_from_name(block["type_name"], block["nverts"], sdim)
        for local_index, elem in enumerate(block["elems"]):
            attr = 1
            if block.get("geom_indices") is not None and local_index < len(block["geom_indices"]):
                raw_attr = block["geom_indices"][local_index]
                attr = raw_attr + 1 if geom in boundary_types else raw_attr

            vertices = list(elem)
            if geom == 5 and hex_perm is not None:
                vertices = [vertices[idx] for idx in hex_perm]
            elif geom == 5:
                vertices = reorder_hex_from_geometry(vertices, nodes)
            elif geom == 7:
                vertices = reorder_pyramid_from_geometry(vertices, nodes)

            if geom in volume_types:
                elements.append((attr, geom, vertices))
            elif geom in boundary_types:
                explicit_boundaries.append((attr, geom, vertices))

    # Align pyramid base orientation with adjacent non-pyramid quad faces.
    quad_reference = {}
    for _, geom, vertices in elements:
        if geom == 7:
            continue
        for face in element_faces(geom, vertices):
            if len(face) != 4:
                continue
            key = tuple(sorted(face))
            if key not in quad_reference:
                quad_reference[key] = face

    aligned_elements = []
    for attr, geom, vertices in elements:
        if geom != 7:
            aligned_elements.append((attr, geom, vertices))
            continue

        key = tuple(sorted(vertices[:4]))
        reference = quad_reference.get(key)
        if reference is None:
            aligned_elements.append((attr, geom, vertices))
            continue

        selected = vertices
        for candidate in pyramid_variants(vertices):
            candidate_base = candidate[:4]
            if is_dihedral_quad(reference, candidate_base):
                selected = candidate
                break

        aligned_elements.append((attr, geom, selected))

    elements = aligned_elements

    # Build oriented face map and incidence count from volume elements.
    orientation_map = {}
    face_incidence = defaultdict(int)
    for _, geom, vertices in elements:
        for face in element_faces(geom, vertices):
            key = tuple(sorted(face))
            orientation_map[key] = face
            face_incidence[key] += 1

    # Map explicit boundary attributes by face key. Keep first seen tag.
    explicit_attr = {}
    for attr, geom, vertices in explicit_boundaries:
        key = tuple(sorted(vertices))
        if key in explicit_attr:
            continue
        explicit_attr[key] = attr

    # MFEM expects boundary entities to be exterior volume faces only.
    boundaries = []
    for key, oriented in orientation_map.items():
        if face_incidence[key] != 1:
            continue

        if len(oriented) == 3:
            geom = 2
        elif len(oriented) == 4:
            geom = 3
        else:
            continue

        attr = explicit_attr.get(key, 1)
        boundaries.append((attr, geom, oriented))

    with open(outpath, "w", encoding="utf-8") as handle:
        handle.write("MFEM mesh v1.0\n\n")
        handle.write("# Converted from mphtxt by mphtxt_to_mfem_mesh.py\n\n")
        handle.write("dimension\n")
        handle.write(f"{sdim}\n\n")

        handle.write("elements\n")
        handle.write(f"{len(elements)}\n")
        for attr, geom, vertices in elements:
            handle.write(f"{attr} {geom} {' '.join(str(v) for v in vertices)}\n")

        handle.write("\nboundary\n")
        handle.write(f"{len(boundaries)}\n")
        for attr, geom, vertices in boundaries:
            handle.write(f"{attr} {geom} {' '.join(str(v) for v in vertices)}\n")

        handle.write("\nvertices\n")
        handle.write(f"{len(nodes)}\n")
        handle.write("3\n")
        for x, y, z in nodes:
            handle.write(f"{x:.16g} {y:.16g} {z:.16g}\n")


def main(argv):
    if len(argv) < 2:
        print("Usage: mphtxt_to_mfem_mesh.py input.mphtxt [output.mesh]")
        return 1

    input_path = Path(argv[1])
    output_path = Path(argv[2]) if len(argv) > 2 else input_path.with_suffix(".mesh")

    print(f"Parsing {input_path} ...")
    nodes, blocks, sdim = parse_mphtxt(input_path)
    print(
        f"Parsed {len(nodes)} nodes and {sum(len(block['elems']) for block in blocks)} "
        f"elements in {len(blocks)} blocks (sdim={sdim})"
    )
    print(f"Writing MFEM mesh to {output_path} ...")
    write_mfem_mesh(nodes, blocks, output_path, sdim)
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
