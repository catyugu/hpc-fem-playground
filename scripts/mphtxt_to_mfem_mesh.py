
#!/usr/bin/env python3
"""
Convert COMSOL-style mphtxt mesh to MFEM v1.0 .mesh (linear, order 1).

Usage:
  python3 mphtxt_to_mfem_mesh.py <input.mphtxt> [output.mesh]

This supports the subset of the mphtxt structure used in the workspace test
files (e.g. `testdata/testmesh1.mphtxt`). It re-uses the conservative
parsing strategy from `mphtxt_to_gmsh.py` to extract the global node list
and element blocks. For this first-phase converter we write a linear
MFEM v1.0 mesh (no "nodes FiniteElementSpace" section) preserving
element attributes when available.

The output format follows MFEM's mesh reader expectations:
  - elements: "<attr> <geom_type> <v0> <v1> ..."
  - boundary: same as elements (attributes first)
  - vertices: count, then spaceDim, then coordinate triples per vertex

Author: generated
"""
import sys
from pathlib import Path


def parse_mphtxt(path):
	# Copy of the robust parser used in mphtxt_to_gmsh.py tailored to the
	# mphtxt files in this repo. It returns (nodes, elements, sdim)
	with open(path, 'r') as f:
		lines = [ln.rstrip('\n') for ln in f]

	i = 0
	nlines = len(lines)
	sdim = None

	verts = []
	while i < nlines:
		if sdim is None and lines[i].strip().endswith('# sdim'):
			try:
				sdim = int(lines[i].split('#')[0].strip())
			except Exception:
				sdim = None

		if lines[i].strip().endswith('# number of mesh vertices'):
			parts = lines[i].split('#')[0].strip()
			if parts:
				nv = int(parts)
				i += 1
			else:
				nv = int(lines[i-1].strip())
				i += 1

			while i < nlines and '# Mesh vertex coordinates' not in lines[i]:
				i += 1
			if i >= nlines:
				raise RuntimeError('Vertex coordinate section not found')
			i += 1
			while len(verts) < nv and i < nlines:
				s = lines[i].strip(); i += 1
				if s == '' or s.startswith('#'):
					continue
				toks = s.split()
				if len(toks) >= 3:
					x, y, z = map(float, toks[:3])
				elif len(toks) == 2:
					x, y = map(float, toks); z = 0.0
				elif len(toks) == 1:
					x = float(toks[0])
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

	elements = []
	while i < nlines:
		line = lines[i]
		if line.strip().startswith('# Type'):
			i += 1
			while i < nlines and lines[i].strip() == '':
				i += 1
			if i >= nlines:
				break
			parts = lines[i].split('#')[0].split()
			if len(parts) >= 2:
				type_name = parts[1].strip()
			else:
				type_name = parts[0].strip()
			i += 1

			while i < nlines and lines[i].strip() == '':
				i += 1
			if i < nlines and '# number of vertices per element' in lines[i]:
				nverts = int(lines[i].split('#')[0]); i += 1
			else:
				try:
					nverts = int(lines[i].split()[0]); i += 1
				except Exception:
					while i < nlines and '# number of vertices per element' not in lines[i]:
						i += 1
					nverts = int(lines[i].split('#')[0]); i += 1

			while i < nlines and lines[i].strip() == '':
				i += 1
			if i < nlines and '# number of elements' in lines[i]:
				ne = int(lines[i].split('#')[0]); i += 1
			else:
				ne = int(lines[i].split()[0]); i += 1

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
				if len(toks) < nverts:
					needed = nverts - len(toks)
					more = []
					while needed > 0 and i < nlines:
						nextt = lines[i].strip(); i += 1
						if nextt == '' or nextt.startswith('#'):
							continue
						subtoks = nextt.split()
						more.extend(subtoks)
						needed = nverts - (len(toks) + len(more))
					toks = toks + more
				node_ids = [int(x) for x in toks[:nverts]]
				elems.append(node_ids)

			geom_indices = None
			j = i
			while j < nlines and lines[j].strip() == '':
				j += 1
			if j < nlines and lines[j].strip().endswith('# number of geometric entity indices'):
				count = int(lines[j].split('#')[0]); j += 1
				while j < nlines and lines[j].strip() == '':
					j += 1
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


def mfem_geom_from_name(name, nverts, sdim=None):
	n = nverts
	name = (name or '').lower()
	if name.startswith('vtx') or name == 'point' or n == 1:
		return 0
	if name.startswith('edg') or name.startswith('line') or n == 2:
		return 1
	if name.startswith('tri') or name.startswith('triangle') or n == 3:
		return 2
	if name.startswith('quad') or name.startswith('square') or n == 4:
		# caution: 4 can also mean tetra in 3D; prefer quad for surface faces
		return 3
	if name.startswith('tet') or name.startswith('tetra'):
		return 4
	if name.startswith('hex') or name.startswith('hexa') or name.startswith('cube'):
		return 5
	if name.startswith('prism') or name.startswith('wedge') or name.startswith('penta'):
		return 6
	if name.startswith('pyr') or name.startswith('pyramid'):
		return 7
	# fallback: use vertex-count heuristics (prefer volumetric for 3D)
	if n == 4 and sdim == 3:
		return 4
	if n == 4 and sdim == 2:
		return 3
	if n == 8:
		return 5
	if n == 3:
		return 2
	if n == 2:
		return 1
	# unknown -> raise
	raise RuntimeError(f"Cannot map element name='{name}' nverts={n} to MFEM geometry")


def write_mfem_mesh(nodes, blocks, outpath, sdim=None):
	# Partition blocks into volume vs boundary elements based on sdim
	vol_types = set()
	bdr_types = set()
	if sdim == 3:
		# include prism/wedge(6) and pyramid(7) as volumetric 3D types
		vol_types = {4, 5, 6, 7}  # tets, hexes, prisms, pyramids
		bdr_types = {2, 3}  # tri, quad
	elif sdim == 2:
		vol_types = {2, 3}
		bdr_types = {1}
	elif sdim == 1:
		vol_types = {1}
		bdr_types = {0}
	else:
		# try 3D-first heuristic
		vol_types = {4, 5}
		bdr_types = {2, 3}

	elements = []
	boundaries = []
	for block in blocks:
		mfem_geom = mfem_geom_from_name(block['type_name'], block['nverts'], sdim)
		for idx, e in enumerate(block['elems']):
			attr = 1
			if block.get('geom_indices') is not None:
				if idx < len(block['geom_indices']):
					attr = block['geom_indices'][idx] + 1
			# Classify strictly for 3D: only triangles/quads are valid boundary
			if mfem_geom in vol_types:
				elements.append((attr, mfem_geom, e))
			elif mfem_geom in bdr_types:
				boundaries.append((attr, mfem_geom, e))
			else:
				# Ambiguous or lower-dimensional entities (points/edges) in a 3D
				# source mesh are not valid MFEM boundary faces. Skip them but
				# record a warning entry so the user can inspect the input.
				# For 2D/1D meshes we can still accept lower-dimensional boundaries.
				if sdim == 3:
					# skip point/edge blocks for 3D output
					# (these typically come from mphtxt as vertex/edge sets)
					# We intentionally do not write them into the MFEM boundary
					# block because MFEM expects face types there.
					continue
				else:
					# For lower-dim meshes, accept the block as a boundary
					boundaries.append((attr, mfem_geom, e))

	# Write the MFEM v1.0 mesh file (linear nodes only)
	# Derive boundary faces from volumetric elements to ensure consistent
	# face types and vertex ordering expected by MFEM. This avoids
	# mismatches between external boundary blocks and element face
	# orientation which can cause MFEM to abort when building face tables.
	def elem_faces(geom, verts_idx):
		# return list of faces (each a list of global vertex indices) for the
		# given element. The ordering below uses common conventions and is
		# sufficient for building face connectivity.
		v = verts_idx
		if geom == 4:  # tetra (4 verts)
			# faces are triangles
			return [ [v[0], v[2], v[1]], [v[0], v[1], v[3]], [v[1], v[2], v[3]], [v[2], v[0], v[3]] ]
		if geom == 5:  # hexahedron (8 verts)
			# 6 quad faces (each 4 verts)
			return [ [v[0],v[1],v[2],v[3]], [v[4],v[7],v[6],v[5]],
				[v[0],v[3],v[7],v[4]], [v[1],v[5],v[6],v[2]],
				[v[0],v[4],v[5],v[1]], [v[3],v[2],v[6],v[7]] ]
		if geom == 6:  # prism / wedge (6 verts)
			# triangular prism: two triangles (0,1,2) and (3,4,5) and three quads
			return [ [v[0],v[2],v[1]], [v[3],v[4],v[5]],
				[v[0],v[1],v[4],v[3]], [v[1],v[2],v[5],v[4]], [v[2],v[0],v[3],v[5]] ]
		if geom == 7:  # pyramid (5 verts)
			# base quad [0,1,2,3] and triangular faces to apex v[4]
			return [ [v[0],v[1],v[2],v[3]], [v[0],v[1],v[4]], [v[1],v[2],v[4]], [v[2],v[3],v[4]], [v[3],v[0],v[4]] ]
		# fallback: try to split into a single face with same verts
		return [ list(v) ]

	# build a lookup of explicit face attributes from any parsed face blocks
	face_attr_map = {}  # key(sorted tuple) -> attribute
	for block in blocks:
		bgeom = mfem_geom_from_name(block['type_name'], block['nverts'], sdim)
		if bgeom in bdr_types:
			# block represents faces (tri/quad); capture their attributes if present
			for idx, fv in enumerate(block['elems']):
				key = tuple(sorted(fv))
				attr = 1
				if block.get('geom_indices') is not None and idx < len(block['geom_indices']):
					attr = block['geom_indices'][idx] + 1
				face_attr_map[key] = attr

	# build face map from elements
	face_map = {}  # key(sorted tuple) -> list of (face_vertices_list, parent_attr)
	for attr, geom, verts_idx in elements:
		faces = elem_faces(geom, verts_idx)
		for fv in faces:
			key = tuple(sorted(fv))
			face_map.setdefault(key, []).append((fv, attr))

	# boundary faces are those that appear only once among all element faces
	boundaries = []
	for key, entries in face_map.items():
		if len(entries) == 1:
			fv, parent_attr = entries[0]
			# prefer explicit attribute from parsed face blocks when present
			attr = face_attr_map.get(key, parent_attr)
			# determine face geom by number of vertices
			fg = 2 if len(fv) == 3 else (3 if len(fv) == 4 else (1 if len(fv) == 2 else 0))
			boundaries.append((attr, fg, fv))

	with open(outpath, 'w') as f:
		f.write('MFEM mesh v1.0\n\n')
		f.write('# Converted from mphtxt by mphtxt_to_mfem_mesh.py\n\n')
		f.write('dimension\n')
		# infer geometry dimension
		if sdim is None:
			# infer from presence of volume elements
			sdim = 3 if any(el[1] in (4,5) for el in elements) else (2 if any(el[1] in (2,3) for el in elements) else 1)
		f.write(f'{sdim}\n\n')

		f.write('elements\n')
		f.write(f'{len(elements)}\n')
		for attr, geom, verts_idx in elements:
			# MFEM element line: <attr> <geom> <v0> <v1> ...
			line = f"{attr} {geom}"
			for vi in verts_idx:
				line += f" {vi}"
			f.write(line + '\n')

		f.write('\n')
		f.write('boundary\n')
		f.write(f'{len(boundaries)}\n')
		for attr, geom, verts_idx in boundaries:
			line = f"{attr} {geom}"
			for vi in verts_idx:
				line += f" {vi}"
			f.write(line + '\n')

		f.write('\n')
		f.write('vertices\n')
		f.write(f'{len(nodes)}\n')
		# write spaceDim token expected by MFEM's reader
		sd = 3 if any(len(p) == 3 for p in nodes) else 2
		f.write(f'{sd}\n')
		for (x, y, z) in nodes:
			f.write(f'{x:.16g} {y:.16g} {z:.16g}\n')


def main(argv):
	if len(argv) < 2:
		print('Usage: mphtxt_to_mfem_mesh.py input.mphtxt [output.mesh]')
		return 1
	inp = Path(argv[1])
	out = Path(argv[2]) if len(argv) > 2 else inp.with_suffix('.mesh')

	print(f'Parsing {inp} ...')
	nodes, blocks, sdim = parse_mphtxt(inp)
	print(f'Parsed {len(nodes)} nodes and {sum(len(b["elems"]) for b in blocks)} elements in {len(blocks)} blocks (sdim={sdim})')
	print(f'Writing MFEM mesh to {out} ...')
	write_mfem_mesh(nodes, blocks, out, sdim)
	print('Done.')
	return 0


if __name__ == '__main__':
	sys.exit(main(sys.argv))
