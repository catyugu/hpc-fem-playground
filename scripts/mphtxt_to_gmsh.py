#!/usr/bin/env python3
"""
Convert a COMSOL-style mphtxt mesh to Gmsh 2.2 ASCII (.msh).

Usage:
  python3 mphtxt_to_gmsh.py <input.mphtxt> <output.msh>

This parser is intentionally conservative and tailored to the structure in
the provided `testdata/mymesh.mphtxt`: it reads the global vertex list and
then all element "Type" blocks (vtx/edg/tri/tet). For each element we also
read the corresponding geometric-entity index list when present and write
that index+1 as the first element tag in the MSH file (so boundary/region
attributes are preserved and can be read by MFEM).

Known mappings (mphtxt -> gmsh element types):
  vtx -> 15 (point)
  edg -> 1  (2-node line)
  tri -> 2  (3-node triangle)
  quad -> 3 (4-node quad) -- not seen in this file but supported
  tet -> 4  (4-node tetrahedron)

The script tries to be robust to extra blank/comment lines. Node indices in
mphtxt are zero-based; Gmsh uses 1-based indexing, so we add +1 to node ids.

Author: generated
"""
import sys
from pathlib import Path


def parse_mphtxt(path):
    with open(path, 'r') as f:
        lines = [ln.rstrip('\n') for ln in f]

    i = 0
    nlines = len(lines)
    sdim = None

    def next_nonempty():
        nonlocal i
        while i < nlines and lines[i].strip() == '':
            i += 1
        if i >= nlines:
            return None
        return lines[i]

    # Find the first occurrence of "# number of mesh vertices"
    verts = []
    while i < nlines:
        # try to capture spatial dimension if present (e.g. '3 # sdim')
        if sdim is None and lines[i].strip().endswith('# sdim'):
            try:
                sdim = int(lines[i].split('#')[0].strip())
            except Exception:
                sdim = None

        if lines[i].strip().endswith('# number of mesh vertices'):
            # previous token (the number) may be on same line before comment
            parts = lines[i].split('#')[0].strip()
            if parts:
                nv = int(parts)
                i += 1
            else:
                # number on previous line
                nv = int(lines[i-1].strip())
                i += 1
            # skip lowest mesh vertex index line if present
            # search for the marker "# Mesh vertex coordinates"
            while i < nlines and '# Mesh vertex coordinates' not in lines[i]:
                i += 1
            if i >= nlines:
                raise RuntimeError('Vertex coordinate section not found')
            # move to next line after marker
            i += 1
            # read nv coordinate lines (may have blank lines or comments between)
            while len(verts) < nv and i < nlines:
                s = lines[i].strip()
                i += 1
                if s == '' or s.startswith('#'):
                    continue
                toks = s.split()
                # some mphtxt occasionally list 1 or 2 numbers per line; assume each vertex provided on one line
                if len(toks) >= 3:
                    x, y, z = map(float, toks[:3])
                elif len(toks) == 2:
                    x, y = map(float, toks)
                    z = 0.0
                elif len(toks) == 1:
                    # unlikely, try next two tokens from subsequent lines
                    x = float(toks[0])
                    # get two more numbers
                    while i < nlines and lines[i].strip() == '':
                        i += 1
                    y = float(lines[i].split()[0]); i += 1
                    while i < nlines and lines[i].strip() == '':
                        i += 1
                    z = float(lines[i].split()[0]); i += 1
                verts.append((x, y, z))
            break
        i += 1

    if not verts:
        raise RuntimeError('No vertices parsed from mphtxt')

    # Now scan for "# Type" blocks and parse element types
    elements = []  # list of dicts with keys: 'type_name','nverts','elems', 'geom_indices'
    # continue from current i
    while i < nlines:
        line = lines[i]
        if line.strip().startswith('# Type'):
            # advance to the type name line
            i += 1
            # read until we find a line that contains type name (e.g., '3 tri # type name')
            while i < nlines and lines[i].strip() == '':
                i += 1
            if i >= nlines:
                break
            # token format: '<num> <name> # type name'
            parts = lines[i].split('#')[0].split()
            if len(parts) >= 2:
                type_name = parts[1].strip()
            else:
                type_name = parts[0].strip()
            i += 1

            # skip empty lines to find "# number of vertices per element"
            while i < nlines and lines[i].strip() == '':
                i += 1
            # read the number-of-vertices-per-element (may be on a line by itself)
            # The file often has the number alone followed by comment
            if i < nlines and '# number of vertices per element' in lines[i]:
                nverts = int(lines[i].split('#')[0])
                i += 1
            else:
                # try previous line or the line itself
                try:
                    nverts = int(lines[i].split()[0])
                    i += 1
                except Exception:
                    # skip until we find the marker
                    while i < nlines and '# number of vertices per element' not in lines[i]:
                        i += 1
                    nverts = int(lines[i].split('#')[0])
                    i += 1

            # read number of elements
            while i < nlines and lines[i].strip() == '':
                i += 1
            if i < nlines and '# number of elements' in lines[i]:
                ne = int(lines[i].split('#')[0])
                i += 1
            else:
                ne = int(lines[i].split()[0]); i += 1

            # skip to '# Elements'
            while i < nlines and '# Elements' not in lines[i]:
                i += 1
            if i < nlines and '# Elements' in lines[i]:
                i += 1

            elems = []
            while len(elems) < ne and i < nlines:
                s = lines[i].strip(); i += 1
                if s == '' or s.startswith('#'):
                    continue
                toks = s.split()
                # only take first nverts integers on the line
                if len(toks) < nverts:
                    # maybe element definitions span multiple lines; gather tokens
                    needed = nverts - len(toks)
                    more = []
                    while needed > 0 and i < nlines:
                        nextt = lines[i].strip()
                        i += 1
                        if nextt == '' or nextt.startswith('#'):
                            continue
                        subtoks = nextt.split()
                        more.extend(subtoks)
                        needed = nverts - (len(toks) + len(more))
                    toks = toks + more
                node_ids = [int(x) for x in toks[:nverts]]
                elems.append(node_ids)

            # After elements, try to read geometric entity indices block
            geom_indices = None
            # skip until maybe we see '# number of geometric entity indices'
            j = i
            while j < nlines and lines[j].strip() == '':
                j += 1
            if j < nlines and lines[j].strip().endswith('# number of geometric entity indices'):
                # number is on same line before '#'
                count = int(lines[j].split('#')[0])
                j += 1
                # skip header line
                while j < nlines and lines[j].strip() == '':
                    j += 1
                # read count integers
                geom_indices = []
                while len(geom_indices) < count and j < nlines:
                    s = lines[j].strip(); j += 1
                    if s == '' or s.startswith('#'):
                        continue
                    for tok in s.split():
                        if len(geom_indices) < count:
                            geom_indices.append(int(tok))
                i = j

            elements.append({'type_name': type_name,
                             'nverts': nverts,
                             'elems': elems,
                             'geom_indices': geom_indices})
        else:
            i += 1

    return verts, elements, sdim
    


def gmsh_type_from_name(name, nverts, sdim=None):
    """
    Map mphtxt element type name and vertex count to Gmsh 2.2 element type code.

    This function is conservative: it uses the textual name when available and
    falls back to a mapping by number of vertices. For ambiguous 4-node blocks
    we use the spatial dimension (sdim) when provided: in 3D a 4-node element is
    usually a tetrahedron (gmsh type 4), while in 2D it's typically a quad
    (gmsh type 3).
    """
    n = nverts
    name = (name or '').lower()
    # points / vertices
    if name.startswith('vtx') or name == 'point' or n == 1:
        return 15
    # lines
    if name.startswith('edg') or name.startswith('line') or n == 2:
        return 1
    # triangles
    if name.startswith('tri') or name.startswith('triangle') or n == 3:
        return 2
    # quads / 4-node faces
    if name.startswith('quad') or name.startswith('quad4'):
        return 3
    # tetrahedra (4-node volumetric)
    if name.startswith('tet') or name.startswith('tetra'):
        return 4
    # prisms, pyramids, hexahedra (common 3D block sizes)
    if name.startswith('prism') or name.startswith('penta') or n == 6:
        return 6  # 6-node prism
    if name.startswith('pyr') or name.startswith('pyramid') or n == 5:
        return 7  # 5-node pyramid
    if name.startswith('hex') or name.startswith('hexa') or n == 8:
        return 5  # 8-node hexahedron

    # If we still have an ambiguous 4-node element, use sdim if available
    if n == 4 and sdim is not None:
        # If ambiguous 4-node elements occur and we don't have a clear name
        # indicating tetra vs quad, prefer quad (surface) in 3D. This is
        # a safer default because interpreting surface faces as tetrahedra
        # often leads to invalid element connectivity/orientation and crashes
        # when read by MFEM.
        return 3

    # fallback mapping for small element sizes
    fallback = {1: 15, 2: 1, 3: 2, 4: 4, 5: 7, 6: 6, 8: 5}
    return fallback.get(n, 0)


def write_msh(nodes, elements, outpath, sdim=None):
    # flatten all elements into a single list and assign gmsh element ids
    flat = []
    for block_idx, block in enumerate(elements):
        gmsh_type = gmsh_type_from_name(block['type_name'], block['nverts'], sdim)
        geom = block.get('geom_indices')
        for idx, e in enumerate(block['elems']):
            tag = 0
            if geom is not None:
                # NOTE: In the COMSOL-style mphtxt files the domain (volumetric)
                # geometric entity indices are 1-based, but boundary/surface
                # geometric entity indices are 0-based. To preserve correct
                # attributes when writing a Gmsh .msh we therefore only add
                # +1 for boundary-like element types (points/lines/tris/quads).
                # Volumetric element "regions" (tets/prisms/hexa) are left as-is.
                if gmsh_type in (15, 1, 2, 3):  # point, line, triangle, quad
                    tag = geom[idx] + 1
                else:
                    tag = geom[idx]
            # validate discovered gmsh_type
            if gmsh_type == 0:
                raise RuntimeError(f"Unknown/unsupported element block type '{block['type_name']}' with {block['nverts']} vertices (block {block_idx})")
            flat.append({'type': gmsh_type, 'nodes': e, 'tag': tag, 'block_idx': block_idx, 'local_idx': idx})

    # Small helper functions for 3D vector math
    def _cross(a, b):
        return (a[1]*b[2] - a[2]*b[1],
                a[2]*b[0] - a[0]*b[2],
                a[0]*b[1] - a[1]*b[0])

    def _dot(a, b):
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

    # Attempt to fix quad orientations (gmsh type 3) for 4-node quads.
    # First, try to align quad ordering with prism (6-node) parent faces
    # when possible. For a standard 6-node prism with nodes [0..5] the
    # three quad side faces (ordered) are typically:
    #   [0,1,4,3], [1,2,5,4], [2,0,3,5]
    # If a quad's node set matches one of these faces (as a set) we replace
    # the quad ordering by the prism face ordering so it matches MFEM's
    # expected orientation coming from volumetric prism elements. If no
    # matching parent is found, fall back to the candidate-rotation method
    # (choose the cyclic/flip variant that maximizes the dot product of the
    # two triangle normals).
    if flat:
        # Build a map of prism-face node sets -> ordered face node list
        prism_face_map = {}
        for el in flat:
            if el['type'] == 6 and len(el['nodes']) >= 6:
                # assume standard prism node ordering [0..5]
                n = el['nodes']
                faces = [
                    [n[0], n[1], n[4], n[3]],
                    [n[1], n[2], n[5], n[4]],
                    [n[2], n[0], n[3], n[5]],
                ]
                for f in faces:
                    prism_face_map[frozenset(f)] = f

        # iterate quads and try to align with parent prism faces first
        for ii, el in enumerate(flat):
            if el['type'] == 3 and len(el['nodes']) == 4:
                nid = el['nodes']
                key = frozenset(nid)
                if key in prism_face_map:
                    flat[ii] = {**el, 'nodes': prism_face_map[key]}
                    continue

                # fallback: pick best candidate rotation/flip by triangle-normal dot
                best = None
                best_dot = None
                candidates = [
                    [nid[0], nid[1], nid[2], nid[3]],
                    [nid[1], nid[2], nid[3], nid[0]],
                    [nid[2], nid[3], nid[0], nid[1]],
                    [nid[3], nid[0], nid[1], nid[2]],
                    [nid[0], nid[3], nid[2], nid[1]],
                    [nid[3], nid[2], nid[1], nid[0]],
                    [nid[2], nid[1], nid[0], nid[3]],
                    [nid[1], nid[0], nid[3], nid[2]],
                ]
                for cand in candidates:
                    # map candidate indices to coordinates
                    pc = []
                    for ni in cand:
                        p = nodes[ni]
                        if len(p) == 2:
                            pc.append((p[0], p[1], 0.0))
                        else:
                            pc.append((p[0], p[1], p[2]))
                    p0, p1, p2, p3 = pc
                    v01 = (p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2])
                    v02 = (p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2])
                    v03 = (p3[0]-p0[0], p3[1]-p0[1], p3[2]-p0[2])
                    n1 = _cross(v01, v02)
                    n2 = _cross(v02, v03)
                    dot = _dot(n1, n2)
                    if best_dot is None or dot > best_dot:
                        best_dot = dot
                        best = cand

                if best is not None and best != nid:
                    flat[ii] = {**el, 'nodes': best}

    with open(outpath, 'w') as f:
        f.write('$MeshFormat\n')
        f.write('2.2 0 8\n')
        f.write('$EndMeshFormat\n')

        f.write('$Nodes\n')
        f.write(f'{len(nodes)}\n')
        for i, (x, y, z) in enumerate(nodes, start=1):
            f.write(f'{i} {x:.16g} {y:.16g} {z:.16g}\n')
        f.write('$EndNodes\n')

        f.write('$Elements\n')
        f.write(f'{len(flat)}\n')
        for eid, el in enumerate(flat, start=1):
            # validate node indices are within bounds
            for nn in el['nodes']:
                if not (0 <= nn < len(nodes)):
                    raise RuntimeError(f"Element {eid} (block {el.get('block_idx')} local {el.get('local_idx')}) references out-of-range node {nn} (0..{len(nodes)-1})")
            nlist = [str(n+1) for n in el['nodes']]  # convert 0-based->1-based
            # number of tags 1 and tag value (0 if unknown)
            f.write(f"{eid} {el['type']} 1 {el['tag']} {' '.join(nlist)}\n")
        f.write('$EndElements\n')


def main(argv):
    if len(argv) < 2:
        print('Usage: mphtxt_to_gmsh.py input.mphtxt [output.msh]')
        return 1
    inp = Path(argv[1])
    out = Path(argv[2]) if len(argv) > 2 else inp.with_suffix('.msh')

    print(f'Parsing {inp} ...')
    nodes, elements, sdim = parse_mphtxt(inp)
    print(f'Parsed {len(nodes)} nodes and {sum(len(b["elems"]) for b in elements)} elements in {len(elements)} blocks (sdim={sdim})')
    print(f'Writing {out} ...')
    write_msh(nodes, elements, out, sdim)
    print('Done.')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
