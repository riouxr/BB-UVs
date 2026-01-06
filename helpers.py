import bpy
import bmesh
import math
import os
import json
from math import isclose, atan2, pi, asin
from mathutils import Vector
from bpy.app.handlers import persistent

# ---------- General helpers ----------
def compute_selection_bbox(ctx):
    # Use current projector targets (uvproj_target) instead of current selection,
    # so behavior stays stable even if user changes selection later.
    obs = get_targets()
    if not obs:
        return None, None

    pts = []
    for o in obs:
        m = o.matrix_world
        for v in o.bound_box:
            pts.append(m @ Vector(v))

    mn = Vector((min(p.x for p in pts),
                 min(p.y for p in pts),
                 min(p.z for p in pts)))

    mx = Vector((max(p.x for p in pts),
                 max(p.y for p in pts),
                 max(p.z for p in pts)))

    return mn, mx


def apply_axis_rotation(obj, axis):
    obj.rotation_euler = (0.0, 0.0, 0.0)

    if axis == 'Z':
        return
    elif axis == 'X':
        obj.rotation_euler = (-math.pi / 2.0, -math.pi, math.pi / 2.0)
    elif axis == 'Y':
        obj.rotation_euler = (math.pi / 2.0, 0.0, 0.0)
        

def fit_projector(ctx, proj, mode):
    mn, mx = compute_selection_bbox(ctx)
    if not mn:
        return

    dx = max(mx.x - mn.x, 1e-6)
    dy = max(mx.y - mn.y, 1e-6)
    dz = max(mx.z - mn.z, 1e-6)

    center = (mn + mx) * 0.5
    proj.location = center

    S = ctx.scene
    keep = getattr(S, "uvproj_keep_ratio", True)
    axis = getattr(S, "uvproj_axis", "Z")

    proj.scale = (1, 1, 1)
    proj.rotation_euler = (0, 0, 0)

    base = 2.0   # Blender primitive = 2 units across

    if mode == "PLANE":
        if axis == 'Z':
            sx, sy = dx, dy
        elif axis == 'X':
            sx, sy = dy, dz
        else:
            sx, sy = dx, dz

        if keep:
            L = max(sx, sy)
            sx = sy = L

        proj.scale = (sx/base, sy/base, 1.0)

    elif mode == "CYL":
        if axis == 'Z':
            h = dz;  r1, r2 = dx, dy
        elif axis == 'X':
            h = dx;  r1, r2 = dy, dz
        else:
            h = dy;  r1, r2 = dx, dz

        diameter = max(r1, r2)

        if keep:
            L = max(dx, dy, dz)
            diameter = L
            h = L

        proj.scale = (diameter/base, diameter/base, h/base)

    elif mode == "SPHERE":
        if keep:
            L = max(dx, dy, dz)
            proj.scale = (L/base, L/base, L/base)
        else:
            proj.scale = (dx/base, dy/base, dz/base)

    elif mode == "CUBE":
        if keep:
            L = max(dx, dy, dz)
            proj.scale = (L/base, L/base, L/base)
        else:
            proj.scale = (dx/base, dy/base, dz/base)
    if mode == "CUBE" and axis in {"X", "Y"}:
        axis = "Z"
    apply_axis_rotation(proj, axis)
    
def _S(context):
    """Return scene-level settings (direct scene props)."""
    return context.scene

def _object_is_uv_mesh(o):
    return (o and o.type == 'MESH' and o.data and len(o.data.uv_layers) > 0 and len(o.data.polygons) > 0)

def _ensure_image(path):
    if not path:
        return None
    abspath = os.path.abspath(bpy.path.abspath(path))
    if not os.path.exists(abspath):
        return None
    for img in bpy.data.images:
        try:
            if os.path.abspath(bpy.path.abspath(img.filepath)) == abspath:
                return img
        except Exception:
            pass
    try:
        return bpy.data.images.load(abspath, check_existing=True)
    except Exception:
        return None

def _ensure_grid_material(img, mat_name="BB_UVs_GridMat"):
    mat = bpy.data.materials.get(mat_name)
    if mat is None:
        mat = bpy.data.materials.new(mat_name)
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    n_out = nt.nodes.new("ShaderNodeOutputMaterial"); n_out.location = (500, 0)
    n_diff = nt.nodes.new("ShaderNodeBsdfDiffuse");   n_diff.location = (250, 0)
    n_tex = nt.nodes.new("ShaderNodeTexImage");       n_tex.location = (0, 0)
    n_map = nt.nodes.new("ShaderNodeMapping");        n_map.location = (-200, 0)
    n_coord = nt.nodes.new("ShaderNodeTexCoord");     n_coord.location = (-400, 0)
    n_tex.interpolation = 'Linear'
    n_tex.extension = 'REPEAT'
    n_tex.image = img
    nt.links.new(n_coord.outputs["UV"], n_map.inputs["Vector"])
    nt.links.new(n_map.outputs["Vector"], n_tex.inputs["Vector"])
    nt.links.new(n_tex.outputs["Color"], n_diff.inputs["Color"])
    nt.links.new(n_diff.outputs["BSDF"], n_out.inputs["Surface"])
    return mat

def _set_all_faces_to_mat(obj, mat_index):
    if obj.mode == 'EDIT':
        bm = bmesh.from_edit_mesh(obj.data)
        for f in bm.faces:
            f.material_index = mat_index
        bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)
    else:
        for p in obj.data.polygons:
            p.material_index = mat_index
        obj.data.update()

def _apply_per_face_indices(obj, indices):
    if not indices:
        return
    if obj.mode == 'EDIT':
        bm = bmesh.from_edit_mesh(obj.data)
        bm.faces.ensure_lookup_table()
        count = min(len(indices), len(bm.faces))
        for i in range(count):
            bm.faces[i].material_index = int(indices[i])
        bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)
    else:
        count = min(len(indices), len(obj.data.polygons))
        for i in range(count):
            obj.data.polygons[i].material_index = int(indices[i])
        obj.data.update()

def _clear_and_rebuild_slots(obj, mat_names):
    mats = obj.data.materials
    mats.clear()
    for name in mat_names:
        if not name:
            m = bpy.data.materials.new("BB_Placeholder")
            mats.append(m)
            continue
        m = bpy.data.materials.get(name)
        if m is None:
            m = bpy.data.materials.new(name + " (Missing)")
        mats.append(m)

def _backup_material_mapping(obj):
    if not _object_is_uv_mesh(obj):
        return
    if obj.get("bb_uvs_prev_mslots") is not None and obj.get("bb_uvs_prev_face_indices") is not None:
        return
    slot_names = [m.name if m else "" for m in obj.data.materials]
    face_indices = [p.material_index for p in obj.data.polygons]
    obj["bb_uvs_prev_mslots"] = json.dumps(slot_names)
    obj["bb_uvs_prev_face_indices"] = json.dumps(face_indices)

def _restore_material_mapping(obj):
    slots_s = obj.get("bb_uvs_prev_mslots", None)
    faces_s = obj.get("bb_uvs_prev_face_indices", None)
    if slots_s is None or faces_s is None:
        return False
    try:
        slot_names = json.loads(slots_s)
        face_indices = json.loads(faces_s)
    except Exception:
        return False
    _clear_and_rebuild_slots(obj, slot_names)
    max_idx = max(0, len(obj.data.materials) - 1)
    face_indices = [min(int(i), max_idx) for i in face_indices]
    _apply_per_face_indices(obj, face_indices)
    try:
        del obj["bb_uvs_prev_mslots"]
        del obj["bb_uvs_prev_face_indices"]
    except Exception:
        pass
    return True

# ---------- TD helpers ----------

def Calculate_TD_Area_To_List():
    active_obj = bpy.context.active_object
    if not _object_is_uv_mesh(active_obj):
        return []
    w = h = 2048
    aspect_ratio = w / h
    if aspect_ratio < 1:
        aspect_ratio = 1 / aspect_ratio
    largest_side = max(w, h)
    src_me = active_obj.data
    tmp_me = src_me.copy()
    bm = bmesh.new()
    bm.from_mesh(tmp_me)
    uv_layer = bm.loops.layers.uv.active
    bm.faces.ensure_lookup_table()
    result = []
    for f in bm.faces:
        if uv_layer is None or len(f.loops) < 3:
            result.append([0.0, 0.0])
            continue
        loops = [loop[uv_layer].uv.copy() for loop in f.loops]
        area_uv = 0.0
        a = len(loops) - 1
        for b in range(len(loops)):
            area_uv += (loops[a].x + loops[b].x) * (loops[a].y - loops[b].y)
            a = b
        area_uv = abs(0.5 * area_uv)
        area_geo = f.calc_area()
        if area_geo > 0 and area_uv > 0:
            td_base = ((largest_side / math.sqrt(aspect_ratio)) * math.sqrt(area_uv)) / (math.sqrt(area_geo) * 100)
            td_base /= bpy.context.scene.unit_settings.scale_length
        else:
            td_base = 0.0001
        result.append([td_base, area_uv])
    bm.free()
    bpy.data.meshes.remove(tmp_me)
    return result

def _move_uvs_object(o, dx, dy):
    bm = bmesh.new()
    bm.from_mesh(o.data)
    uv_layer = bm.loops.layers.uv.active
    if uv_layer is None:
        bm.free()
        return
    bm.faces.ensure_lookup_table()
    for f in bm.faces:
        for loop in f.loops:
            uv = loop[uv_layer].uv
            uv.x += dx; uv.y += dy
    bm.to_mesh(o.data); o.data.update(); bm.free()

def _move_uvs_edit(o, dx, dy, only_selected=False):
    bm = bmesh.from_edit_mesh(o.data)
    uv_layer = bm.loops.layers.uv.active
    if uv_layer is None:
        return
    bm.faces.ensure_lookup_table()
    if not only_selected:
        for f in bm.faces:
            for loop in f.loops:
                uv = loop[uv_layer].uv
                uv.x += dx; uv.y += dy
        bmesh.update_edit_mesh(o.data, loop_triangles=False, destructive=False)
        return
    tool = bpy.context.tool_settings
    uv_sync = tool.use_uv_select_sync
    if not uv_sync:
        for f in bm.faces:
            for loop in f.loops:
                if not loop[uv_layer].select:
                    continue
                uv = loop[uv_layer].uv
                uv.x += dx; uv.y += dy
        bmesh.update_edit_mesh(o.data, loop_triangles=False, destructive=False)
        return
    v_mode, e_mode, f_mode = tool.mesh_select_mode
    if f_mode:
        for f in bm.faces:
            if not f.select:
                continue
            for loop in f.loops:
                uv = loop[uv_layer].uv
                uv.x += dx; uv.y += dy
    else:
        for f in bm.faces:
            for loop in f.loops:
                if not loop.vert.select:
                    continue
                uv = loop[uv_layer].uv
                uv.x += dx; uv.y += dy
    bmesh.update_edit_mesh(o.data, loop_triangles=False, destructive=False)

# ---------- Normalize helpers ----------

def _poly_uv_area(face, uv_layer):
    loops = [l[uv_layer].uv.copy() for l in face.loops]
    if len(loops) < 3:
        return 0.0
    area = 0.0
    a = len(loops) - 1
    for b in range(len(loops)):
        area += (loops[a].x + loops[b].x) * (loops[a].y - loops[b].y)
        a = b
    return abs(0.5 * area)

def _uv_key(v, tol=1e-6):
    return (round(v.x / tol) * tol, round(v.y / tol) * tol)

def _collect_uv_islands(bm, uv_layer, only_selected):
    bm.faces.ensure_lookup_table()
    faces = [f for f in bm.faces if (not only_selected or f.select)]
    if not faces:
        return []
    uv_to_faces = {}
    for f in faces:
        for loop in f.loops:
            uv = loop[uv_layer].uv
            key = _uv_key(uv)
            uv_to_faces.setdefault(key, set()).add(f)
    face_adj = {f: set() for f in faces}
    for shared_faces in uv_to_faces.values():
        shared_faces = list(shared_faces)
        for i in range(len(shared_faces)):
            fi = shared_faces[i]
            for j in range(i + 1, len(shared_faces)):
                fj = shared_faces[j]
                face_adj[fi].add(fj)
                face_adj[fj].add(fi)
    islands = []
    seen = set()
    for f in faces:
        if f in seen:
            continue
        stack = [f]
        comp = set()
        seen.add(f)
        while stack:
            cur = stack.pop()
            comp.add(cur)
            for nb in face_adj[cur]:
                if nb not in seen:
                    seen.add(nb)
                    stack.append(nb)
        islands.append(comp)
    return islands

def _center_of_island_uv(island_faces, uv_layer):
    sx = sy = cnt = 0.0
    for f in island_faces:
        for l in f.loops:
            uv = l[uv_layer].uv
            sx += uv.x; sy += uv.y; cnt += 1.0
    if cnt == 0:
        return 0.0, 0.0
    return sx / cnt, sy / cnt

def _normalize_islands_in_bmesh(bm, xy_independent=True, only_selected=True, respect_aspect=True):
    uv_layer = bm.loops.layers.uv.active
    if not uv_layer:
        return
    islands = _collect_uv_islands(bm, uv_layer, only_selected)
    if not islands:
        return
    tot_area_uv = 0.0
    tot_area_3d = 0.0
    for faces in islands:
        for f in faces:
            a3 = f.calc_area()
            au = _poly_uv_area(f, uv_layer)
            tot_area_3d += a3; tot_area_uv += au
    if tot_area_uv <= 1e-12 or tot_area_3d <= 1e-12:
        return
    tot_fac = tot_area_3d / tot_area_uv
    for faces in islands:
        area3 = 0.0; areauv = 0.0
        for f in faces:
            area3 += f.calc_area(); areauv += _poly_uv_area(f, uv_layer)
        if areauv <= 1e-12 or area3 <= 1e-12:
            continue
        fac = area3 / areauv
        scale = math.sqrt(fac / tot_fac)
        if isclose(scale, 1.0, abs_tol=1e-6):
            continue
        cx, cy = _center_of_island_uv(faces, uv_layer)
        for f in faces:
            for l in f.loops:
                uv = l[uv_layer].uv
                uv.x = (uv.x - cx) * scale + cx
                uv.y = (uv.y - cy) * scale + cy

def _rotate_uvs_group_in_bmesh(bm, angle_deg=90.0, only_selected=True):
    uv_layer = bm.loops.layers.uv.active
    if not uv_layer: return
    loops_to_rotate = []
    tool = bpy.context.tool_settings
    uv_sync = getattr(tool, "use_uv_select_sync", False)
    if not only_selected:
        for f in bm.faces:
            for l in f.loops: loops_to_rotate.append(l)
    else:
        if not uv_sync:
            for f in bm.faces:
                for l in f.loops:
                    # Blender 4.0+: BMLoopUV.select no longer exists, use vert.select
                    if l.vert.select: loops_to_rotate.append(l)
        else:
            v_mode, e_mode, f_mode = tool.mesh_select_mode
            if f_mode:
                for f in bm.faces:
                    if not f.select: continue
                    for l in f.loops: loops_to_rotate.append(l)
            else:
                for f in bm.faces:
                    for l in f.loops:
                        if l.vert.select: loops_to_rotate.append(l)
    if not loops_to_rotate: return
    sx = sy = 0.0; n = 0
    for l in loops_to_rotate:
        uv = l[uv_layer].uv; sx += uv.x; sy += uv.y; n += 1
    if n == 0: return
    cx = sx / n; cy = sy / n
    ang = angle_deg % 360.0
    if ang in (90.0, -270.0):
        def rot(u, v): return -v, u
    elif ang in (270.0, -90.0):
        def rot(u, v): return v, -u
    elif ang in (180.0, -180.0):
        def rot(u, v): return -u, -v
    else:
        rad = math.radians(ang); c, s = math.cos(rad), math.sin(rad)
        def rot(u, v): return u * c - v * s, u * s + v * c
    for l in loops_to_rotate:
        uv = l[uv_layer].uv; u, v = uv.x - cx, uv.y - cy
        ru, rv = rot(u, v); uv.x = ru + cx; uv.y = rv + cy

def _flip_uvs_group_in_bmesh(bm, axis='H', only_selected=True):
    uv_layer = bm.loops.layers.uv.active
    if not uv_layer: return
    loops_to_flip = []
    tool = bpy.context.tool_settings
    uv_sync = getattr(tool, "use_uv_select_sync", False)
    if not only_selected:
        for f in bm.faces:
            for l in f.loops: loops_to_flip.append(l)
    else:
        if not uv_sync:
            for f in bm.faces:
                for l in f.loops:
                    # Blender 4.0+: BMLoopUV.select no longer exists, use vert.select
                    if l.vert.select: loops_to_flip.append(l)
        else:
            v_mode, e_mode, f_mode = tool.mesh_select_mode
            if f_mode:
                for f in bm.faces:
                    if f.select:
                        for l in f.loops: loops_to_flip.append(l)
            else:
                for f in bm.faces:
                    for l in f.loops:
                        if l.vert.select: loops_to_flip.append(l)
    if not loops_to_flip: return
    sx = sy = 0.0
    for l in loops_to_flip:
        uv = l[uv_layer].uv; sx += uv.x; sy += uv.y
    cx = sx / len(loops_to_flip); cy = sy / len(loops_to_flip)
    for l in loops_to_flip:
        uv = l[uv_layer].uv
        if axis == 'H': uv.x = 2*cx - uv.x
        elif axis == 'V': uv.y = 2*cy - uv.y

# ---------- Projector helpers ----------

def ensure_uv(obj):
    if not obj.data.uv_layers:
        obj.data.uv_layers.new(name="UVMap")


def projector_local(obj, proj, co):
    return (proj.matrix_world.inverted() @ obj.matrix_world) @ co


def get_targets():
    return [o for o in bpy.data.objects if o.get("uvproj_target")]


def _unwrap_loop_sequence(values, threshold=0.5):
    """Make a scalar UV sequence continuous by adding integer offsets.
    This lets UVs cross 0–1 borders cleanly (UDIM wrapping)."""
    if not values:
        return values
    out = [values[0]]
    offset = 0.0
    prev = values[0]
    for u in values[1:]:
        base = u + offset
        # Shift by whole tiles until neighbor difference is small
        while base - prev > threshold:
            offset -= 1.0
            base = u + offset
        while prev - base > threshold:
            offset += 1.0
            base = u + offset
        out.append(base)
        prev = base
    return out


# ---------- Projector projections ----------

def apply_planar(obj, proj):
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    uv = bm.loops.layers.uv.verify()

    for f in bm.faces:
        for l in f.loops:
            p = projector_local(obj, proj, l.vert.co)
            l[uv].uv = (p.x * 0.5 + 0.5, p.y * 0.5 + 0.5)

    bm.to_mesh(obj.data)
    bm.free()
    obj.data.update()


def apply_cyl(obj, proj):
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    uv = bm.loops.layers.uv.verify()

    for f in bm.faces:
        # First compute raw cylindrical coordinates
        us = []
        vs = []
        loops = []
        for l in f.loops:
            p = projector_local(obj, proj, l.vert.co)
            a = atan2(p.y, p.x)
            u = (a / (2 * pi)) + 0.5
            v = p.z * 0.5 + 0.5
            us.append(u)
            vs.append(v)
            loops.append(l)

        # Make U continuous across seam (wrap over UDIM tiles)
        us_cont = _unwrap_loop_sequence(us)

        for l, u, v in zip(loops, us_cont, vs):
            l[uv].uv = (u, v)

    bm.to_mesh(obj.data)
    bm.free()
    obj.data.update()


def apply_sphere(obj, proj):
    """Spherical projection using projector as origin.
    Horizontal (U) wraps so seams can cross UDIM tiles cleanly."""
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    uv = bm.loops.layers.uv.verify()

    for f in bm.faces:
        us = []
        vs = []
        loops = []
        for l in f.loops:
            p = projector_local(obj, proj, l.vert.co)
            if p.length_squared == 0.0:
                u = 0.5
                v = 0.5
            else:
                p_n = p.normalized()
                u = (atan2(p_n.y, p_n.x) / (2 * pi)) + 0.5
                z = max(-1.0, min(1.0, p_n.z))
                # z in [-1,1] → v in [0,1]
                v = 0.5 + (asin(z) / pi)
            us.append(u)
            vs.append(v)
            loops.append(l)

        us_cont = _unwrap_loop_sequence(us)

        for l, u, v in zip(loops, us_cont, vs):
            l[uv].uv = (u, v)

    bm.to_mesh(obj.data)
    bm.free()
    obj.data.update()


def apply_cube(obj, proj):
    import bmesh

    bm = bmesh.new()
    bm.from_mesh(obj.data)
    uv = bm.loops.layers.uv.verify()

    grid_w = 1.0 / 3.0
    grid_h = 1.0 / 2.0
    pad = 0.02

    tiles = {
        '+X': (0, 1),
        '-X': (1, 1),
        '+Y': (2, 1),
        '-Y': (0, 0),
        '+Z': (1, 0),
        '-Z': (2, 0),
    }

    min_side = min(grid_w, grid_h)

    for f in bm.faces:
        n = f.normal
        ax = max(range(3), key=lambda i: abs(n[i]))
        sign = '+' if n[ax] >= 0 else '-'
        face = sign + 'XYZ'[ax]

        if ax == 0:
            coord_u, coord_v = 1, 2  # y, z
        elif ax == 1:
            coord_u, coord_v = 0, 2  # x, z
        else:
            coord_u, coord_v = 0, 1  # x, y

        tx, ty = tiles[face]
        offset_u = (grid_w - min_side) / 2
        offset_v = (grid_h - min_side) / 2

        for l in f.loops:
            p = projector_local(obj, proj, l.vert.co)
            u = p[coord_u] * 0.5 + 0.5
            v = p[coord_v] * 0.5 + 0.5

            # pad inside tile
            u = pad + u * (1 - 2 * pad)
            v = pad + v * (1 - 2 * pad)

            U = tx * grid_w + offset_u + u * min_side
            V = ty * grid_h + offset_v + v * min_side

            l[uv].uv = (U, V)

    bm.to_mesh(obj.data)
    bm.free()
    obj.data.update()

# ---------- Projector HANDLER ----------

@persistent
def projector_update(scene=None, depsgraph=None):
    import bpy

    # When called from depsgraph_update_post Blender passes only `depsgraph`
    if scene is None or not hasattr(scene, "uvproj_running"):
        scene = bpy.context.scene

    if not getattr(scene, "uvproj_running", False):
        return

    proj = scene.uvproj_projector
    if not proj:
        return

    current_matrix = [list(row) for row in proj.matrix_world]
    last_matrix = scene.get("uvproj_last_matrix")

    # Only skip if we actually HAVE a previous value
    if last_matrix is not None and last_matrix == current_matrix:
        return

    scene["uvproj_last_matrix"] = current_matrix

    mode = proj.get("uvproj_mode", "PLANE")

    from .helpers import get_targets, ensure_uv, apply_planar, apply_cyl, apply_cube, apply_sphere

    for ob in get_targets():
        ensure_uv(ob)

        if mode == "PLANE":
            apply_planar(ob, proj)
        elif mode == "CYL":
            apply_cyl(ob, proj)
        elif mode == "CUBE":
            apply_cube(ob, proj)
        elif mode == "SPHERE":
            apply_sphere(ob, proj)
        else:
            apply_planar(ob, proj)
