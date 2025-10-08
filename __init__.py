# SPDX-License-Identifier: GPL-3.0-or-later
# -----------------------------------------------------------------------------
# BB UVs (Minimal + Move UVs — Selected/Highlighted Defaults)
# Derived work based on:
# • Texel Density Checker — https://github.com/mrven/Blender-Texel-Density-Checker
# Original author: Ivan Vostrikov (mrven)
# • Univ (Universal UV Tools) — https://extensions.blender.org/add-ons/univ/
# Portions of the normalize/pack workflow and related ideas are adapted.
#
# This file is part of a Blender add-on that includes code and ideas derived
# from the projects “Texel Density Checker” by Ivan Vostrikov (mrven) and
# “Univ” (Blender Extensions). Both are licensed under the GNU General Public
# License version 3.0 or (at your option) any later version.
#
# Copyright (C) 2025 BB (modifications and adaptations)
# Copyright (C) Ivan Vostrikov (mrven) and contributors — Texel Density Checker
# Copyright (C) The Univ project authors and contributors — Univ extension
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.
#
# Summary of modifications relative to the upstream projects:
# - Minimal “Get/Set Texel Density” implementation tailored for fast workflows.
# - Independent per-object TD calculation without changing selections.
# - Normalize-all then single-pack workflow (Edit/Object mode) with snap-to-first-UDIM.
# - Move/Rotate/Flip UVs with contextual Selected/Highlighted toggle.
# - Compact UI panels for UV Editor and 3D View.
#
# Third‑party notices:
# - Portions of the normalize/pack logic and flow are adapted from the Univ extension;
# identifiers and structure may have been renamed/simplified to fit this add-on.
# - Attribution above satisfies GPL 3.0-or-later requirements; licensing remains GPL‑3.0‑or‑later.

# -----------------------------------------------------------------------------
# BB UVs (Minimal + Move UVs + Normalize+Pack in one click)
# -----------------------------------------------------------------------------

bl_info = {
    "name": "BB UVs (Grid Helper + UV tools)",
    "author": "BB",
    "version": (1, 6, 0),
    "blender": (4, 5, 0),
    "location": "UV Editor & 3D View > Sidebar (N) > BB UVs",
    "description": "Texel Density + Move/Rotate/Flip + Normalize/Pack + Grid preview (per-face revert). Uses Scene properties directly (no PropertyGroup).",
    "category": "UV",
}

import bpy
import bmesh
import math
import os
import json
from math import isclose
from bpy.props import StringProperty, BoolProperty, FloatProperty, EnumProperty
from bpy.types import Panel, Operator

# ------------------------------ Utils ------------------------------

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
    n_out = nt.nodes.new("ShaderNodeOutputMaterial"); n_out.location = (400, 0)
    n_diff = nt.nodes.new("ShaderNodeBsdfDiffuse");   n_diff.location = (150, 0)
    n_tex = nt.nodes.new("ShaderNodeTexImage");       n_tex.location = (-150, 0)
    n_tex.interpolation = 'Linear'
    n_tex.extension = 'REPEAT'
    n_tex.image = img
    nt.links.new(n_tex.outputs["Color"], n_diff.inputs["Color"])
    nt.links.new(n_diff.outputs["BSDF"], n_out.inputs["Surface"])
    return mat

def _set_all_faces_to_mat(obj, mat_index):
    if obj.data.is_editmode:
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
    if obj.data.is_editmode:
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
                    if l[uv_layer].select: loops_to_rotate.append(l)
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
                    if l[uv_layer].select: loops_to_flip.append(l)
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

# ------------------------------ TD / Move / Rotate / Flip Operators ------------------------------

class BB_Texel_Density_Check(Operator):
    bl_idname = "bb_uvs.texel_density_check"
    bl_label = "Get TD"
    bl_options = {'REGISTER', 'UNDO'}
    def execute(self, context):
        S = _S(context)
        active_obj = context.active_object
        if not _object_is_uv_mesh(active_obj):
            self.report({'WARNING'}, "No active mesh with UVs"); return {'CANCELLED'}
        lst = Calculate_TD_Area_To_List()
        if not lst: S.bb_density = 0.0; return {'CANCELLED'}
        total_area = sum(a[1] for a in lst)
        if total_area == 0: S.bb_density = 0.0; return {'CANCELLED'}
        weighted_td = sum(val * (area / total_area) for val, area in lst)
        S.bb_density = float(weighted_td); return {'FINISHED'}

class BB_Texel_Density_Set(Operator):
    bl_idname = "bb_uvs.texel_density_set"
    bl_label = "Set TD"
    bl_options = {'REGISTER', 'UNDO'}
    def execute(self, context):
        S = _S(context); target_td = float(S.bb_density)
        if target_td <= 0: self.report({'WARNING'}, "No valid TD. Run Get TD first."); return {'CANCELLED'}
        start_active = context.view_layer.objects.active
        start_mode = getattr(start_active, "mode", "OBJECT") if start_active else "OBJECT"
        start_selected = context.selected_objects[:]
        if start_mode == 'EDIT':
            targets = [o for o in context.objects_in_mode if _object_is_uv_mesh(o)]
        else:
            targets = [o for o in start_selected if _object_is_uv_mesh(o)]
        if not targets: self.report({'WARNING'}, "No mesh with UVs selected"); return {'CANCELLED'}
        for o in targets:
            prev = context.view_layer.objects.active; context.view_layer.objects.active = o
            lst = Calculate_TD_Area_To_List(); context.view_layer.objects.active = prev
            total_area = sum(a[1] for a in lst)
            if total_area <= 0: continue
            current_td = sum(val * (area / total_area) for val, area in lst) or 0.0001
            scale = target_td / current_td
            if o.data.is_editmode:
                bm = bmesh.from_edit_mesh(o.data)
                uv_layer = bm.loops.layers.uv.active
                if uv_layer:
                    bm.faces.ensure_lookup_table()
                    for f in bm.faces:
                        for loop in f.loops:
                            uv = loop[uv_layer].uv; uv.x *= scale; uv.y *= scale
                    bmesh.update_edit_mesh(o.data, loop_triangles=False, destructive=False)
            else:
                bm = bmesh.new(); bm.from_mesh(o.data)
                uv_layer = bm.loops.layers.uv.active
                if uv_layer:
                    bm.faces.ensure_lookup_table()
                    for f in bm.faces:
                        for loop in f.loops:
                            uv = loop[uv_layer].uv; uv.x *= scale; uv.y *= scale
                    bm.to_mesh(o.data); o.data.update()
                bm.free()
        return {'FINISHED'}

class BB_UVs_SetMoveContext(Operator):
    bl_idname = "bb_uvs.set_move_context"
    bl_label = "Set Move Context"
    bl_options = {'INTERNAL', 'UNDO'}
    set_selected: BoolProperty(name="Set Selected Mode", default=True)
    def execute(self, context):
        _S(context).bb_move_selected_only = bool(self.set_selected); return {'FINISHED'}

class BB_UVs_ToggleMoveHighlight(Operator):
    bl_idname = "bb_uvs.toggle_move_highlight"
    bl_label = "Toggle Highlighted Only"
    bl_options = {'INTERNAL', 'UNDO'}
    def execute(self, context):
        S = _S(context); S.bb_move_selected_only = not S.bb_move_selected_only; return {'FINISHED'}

class BB_UVs_MoveUVs(Operator):
    bl_idname = "bb_uvs.move_uvs"
    bl_label = "Move UVs"
    bl_options = {'REGISTER', 'UNDO'}
    direction: EnumProperty(
        items=[
            ('UP', "Up", ""), ('DOWN', "Down", ""), ('LEFT', "Left", ""), ('RIGHT', "Right", ""),
            ('UP_LEFT', "Up Left", ""), ('UP_RIGHT', "Up Right", ""),
            ('DOWN_LEFT', "Down Left", ""), ('DOWN_RIGHT', "Down Right", ""),
        ], name="Direction", default='UP'
    )
    def execute(self, context):
        S = _S(context); amt = float(S.bb_move_amount)
        if amt == 0.0: return {'CANCELLED'}
        dx = dy = 0.0
        if self.direction == 'UP': dy = amt
        elif self.direction == 'DOWN': dy = -amt
        elif self.direction == 'LEFT': dx = -amt
        elif self.direction == 'RIGHT': dx = amt
        elif self.direction == 'UP_LEFT': dx = -amt; dy = amt
        elif self.direction == 'UP_RIGHT': dx = amt; dy = amt
        elif self.direction == 'DOWN_LEFT': dx = -amt; dy = -amt
        elif self.direction == 'DOWN_RIGHT': dx = amt; dy = -amt
        active = context.view_layer.objects.active
        mode = getattr(active, "mode", "OBJECT") if active else "OBJECT"
        if mode == 'EDIT':
            if S.bb_move_selected_only:
                targets = [o for o in context.objects_in_mode if _object_is_uv_mesh(o)]
                for o in targets: _move_uvs_edit(o, dx, dy, only_selected=True)
            else:
                if not _object_is_uv_mesh(active):
                    self.report({'WARNING'}, "Active object is not a mesh with UVs"); return {'CANCELLED'}
                _move_uvs_edit(active, dx, dy, only_selected=False)
        else:
            selected = context.selected_objects[:]
            if S.bb_move_selected_only:
                if not _object_is_uv_mesh(active):
                    self.report({'WARNING'}, "Active object is not a mesh with UVs"); return {'CANCELLED'}
                targets = [active]
            else:
                targets = [o for o in selected if _object_is_uv_mesh(o)]
            if not targets:
                self.report({'WARNING'}, "No mesh with UVs to move"); return {'CANCELLED'}
            for o in targets: _move_uvs_object(o, dx, dy)
        try: bpy.ops.bb_uvs.texel_density_check()
        except Exception: pass
        return {'FINISHED'}

class BB_UVs_RotateUVs(Operator):
    bl_idname = "bb_uvs.rotate_uvs"
    bl_label = "Rotate UVs"
    bl_options = {'REGISTER', 'UNDO'}
    direction: EnumProperty(
        items=[('CCW', "CCW", "Rotate 90° counter-clockwise"),
               ('CW',  "CW",  "Rotate 90° clockwise")],
        name="Direction", default='CCW'
    )
    def execute(self, context):
        angle = 90.0 if self.direction == 'CCW' else -90.0
        active = context.view_layer.objects.active
        if not _object_is_uv_mesh(active):
            self.report({'WARNING'}, "No active mesh with UVs"); return {'CANCELLED'}
        mode = getattr(active, "mode", "OBJECT")
        if mode == 'EDIT':
            targets = [o for o in context.objects_in_mode if _object_is_uv_mesh(o)]
            for o in targets:
                bm = bmesh.from_edit_mesh(o.data)
                _rotate_uvs_group_in_bmesh(bm, angle_deg=angle, only_selected=True)
                bmesh.update_edit_mesh(o.data, loop_triangles=False, destructive=False)
        else:
            bm = bmesh.new(); bm.from_mesh(active.data)
            _rotate_uvs_group_in_bmesh(bm, angle_deg=angle, only_selected=False)
            bm.to_mesh(active.data); active.data.update(); bm.free()
        try: bpy.ops.bb_uvs.texel_density_check()
        except Exception: pass
        return {'FINISHED'}

class BB_UVs_FlipUVs(Operator):
    bl_idname = "bb_uvs.flip_uvs"
    bl_label = "Flip UVs"
    bl_options = {'REGISTER', 'UNDO'}
    axis: EnumProperty(
        items=[('H', "Horizontal", "Flip left-right"),
               ('V', "Vertical",   "Flip top-bottom")],
        name="Axis", default='H'
    )
    def execute(self, context):
        active = context.view_layer.objects.active
        if not _object_is_uv_mesh(active):
            self.report({'WARNING'}, "No active mesh with UVs"); return {'CANCELLED'}
        mode = getattr(active, "mode", "OBJECT")
        if mode == 'EDIT':
            targets = [o for o in context.objects_in_mode if _object_is_uv_mesh(o)]
            for o in targets:
                bm = bmesh.from_edit_mesh(o.data)
                _flip_uvs_group_in_bmesh(bm, axis=self.axis, only_selected=True)
                bmesh.update_edit_mesh(o.data, loop_triangles=False, destructive=False)
        else:
            bm = bmesh.new(); bm.from_mesh(active.data)
            _flip_uvs_group_in_bmesh(bm, axis=self.axis, only_selected=False)
            bm.to_mesh(active.data); active.data.update(); bm.free()
        return {'FINISHED'}

# ------------------------------ Normalize + Pack ------------------------------

class BB_UVs_NormalizePack(Operator):
    bl_idname = "bb_uvs.normalize_pack"
    bl_label = "Pack Together"
    bl_options = {'REGISTER', 'UNDO'}
    @staticmethod
    def _call_pack(args):
        try: return bpy.ops.uv.pack_islands('EXEC_DEFAULT', **args)
        except TypeError: pass
        args2 = dict(args); args2.pop("shape_method", None)
        try: return bpy.ops.uv.pack_islands('EXEC_DEFAULT', **args2)
        except TypeError: pass
        args3 = dict(args2); args3.pop("rotate", None)
        try: return bpy.ops.uv.pack_islands('EXEC_DEFAULT', **args3)
        except TypeError: pass
        args4 = {"udim_source": args.get("udim_source", "ACTIVE_UDIM"),
                 "margin_method": args.get("margin_method", "FRACTION"),
                 "margin": args.get("margin", 0.00390625)}
        return bpy.ops.uv.pack_islands('EXEC_DEFAULT', **args4)
    @staticmethod
    def _snap_uvs_to_first_udim_all():
        objs = [o for o in bpy.context.objects_in_mode if o.type == 'MESH']
        if not objs: return
        min_u = float('inf'); min_v = float('inf')
        for obj in objs:
            bm = bmesh.from_edit_mesh(obj.data)
            uv_layer = bm.loops.layers.uv.active
            if not uv_layer: continue
            for f in bm.faces:
                if not f.select: continue
                for l in f.loops:
                    uv = l[uv_layer].uv
                    if uv.x < min_u: min_u = uv.x
                    if uv.y < min_v: min_v = uv.y
        if min_u == float('inf'): return
        du = -math.floor(min_u); dv = -math.floor(min_v)
        if du == 0 and dv == 0: return
        for obj in objs:
            bm = bmesh.from_edit_mesh(obj.data)
            uv_layer = bm.loops.layers.uv.active
            if not uv_layer: continue
            for f in bm.faces:
                if not f.select: continue
                for l in f.loops:
                    uv = l[uv_layer].uv
                    uv.x += du; uv.y += dv
            bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)
    @staticmethod
    def _compute_global_tot_fac(only_selected=True):
        total_area_3d = 0.0; total_area_uv = 0.0
        objs = [o for o in bpy.context.objects_in_mode if o.type == 'MESH']
        for obj in objs:
            bm = bmesh.from_edit_mesh(obj.data)
            uv_layer = bm.loops.layers.uv.active
            if not uv_layer: continue
            for f in bm.faces:
                if only_selected and not f.select: continue
                a3 = f.calc_area(); au = _poly_uv_area(f, uv_layer)
                total_area_3d += a3; total_area_uv += au
        if total_area_uv <= 1e-12 or total_area_3d <= 1e-12: return None
        return total_area_3d / total_area_uv
    @staticmethod
    def _normalize_all_edit_meshes_to_fac(tot_fac, only_selected=True):
        objs = [o for o in bpy.context.objects_in_mode if o.type == 'MESH']
        for obj in objs:
            bm = bmesh.from_edit_mesh(obj.data)
            uv_layer = bm.loops.layers.uv.active
            if not uv_layer: continue
            islands = _collect_uv_islands(bm, uv_layer, only_selected=only_selected)
            if not islands: continue
            for faces in islands:
                area3 = 0.0; areauv = 0.0
                for f in faces:
                    area3 += f.calc_area(); areauv += _poly_uv_area(f, uv_layer)
                if areauv <= 1e-12 or area3 <= 1e-12: continue
                fac = area3 / areauv
                scale = math.sqrt(fac / tot_fac)
                if isclose(scale, 1.0, abs_tol=1e-6): continue
                cx, cy = _center_of_island_uv(faces, uv_layer)
                for f in faces:
                    for l in f.loops:
                        uv = l[uv_layer].uv
                        uv.x = (uv.x - cx) * scale + cx
                        uv.y = (uv.y - cy) * scale + cy
            bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)
    def execute(self, context):
        args = {
            "udim_source": "ACTIVE_UDIM", "rotate": False, "scale": True,
            "margin_method": "FRACTION", "margin": 0.00390625, "pin": False,
            "merge_overlap": False, "pin_method": "LOCKED", "shape_method": "CONCAVE",
        }
        active = context.view_layer.objects.active
        if not _object_is_uv_mesh(active):
            self.report({'WARNING'}, "No active mesh with UVs"); return {'CANCELLED'}
        start_mode = getattr(active, "mode", "OBJECT")
        if start_mode == 'EDIT':
            ts = context.tool_settings
            ts.mesh_select_mode = (False, False, True)
            bpy.ops.mesh.select_all(action='SELECT')
            try: bpy.ops.uv.select_all(action='SELECT')
            except Exception: pass
            tot_fac = self._compute_global_tot_fac(only_selected=True)
            if tot_fac is None:
                self.report({'WARNING'}, "Not enough UV/mesh data to normalize"); return {'CANCELLED'}
            self._normalize_all_edit_meshes_to_fac(tot_fac, only_selected=True)
            self._call_pack(args); self._snap_uvs_to_first_udim_all()
            return {'FINISHED'}
        prev_active = active; prev_selection = context.selected_objects[:]
        targets = [o for o in prev_selection if _object_is_uv_mesh(o)]
        if not targets:
            self.report({'WARNING'}, "No selected meshes with UVs to pack"); return {'CANCELLED'}
        for o in bpy.data.objects: o.select_set(False)
        for o in targets: o.select_set(True)
        context.view_layer.objects.active = targets[0]
        bpy.ops.object.mode_set(mode='EDIT')
        ts = context.tool_settings
        ts.mesh_select_mode = (False, False, True)
        bpy.ops.mesh.select_all(action='SELECT')
        try: bpy.ops.uv.select_all(action='SELECT')
        except Exception: pass
        tot_fac = self._compute_global_tot_fac(only_selected=False)
        if tot_fac is None:
            bpy.ops.object.mode_set(mode='OBJECT')
            for o in bpy.data.objects: o.select_set(False)
            for o in prev_selection:
                if o and o.name in bpy.data.objects: o.select_set(True)
            context.view_layer.objects.active = prev_active
            self.report({'WARNING'}, "Not enough UV/mesh data to normalize"); return {'CANCELLED'}
        self._normalize_all_edit_meshes_to_fac(tot_fac, only_selected=False)
        self._call_pack(args); self._snap_uvs_to_first_udim_all()
        bpy.ops.object.mode_set(mode='OBJECT')
        for o in bpy.data.objects: o.select_set(False)
        for o in prev_selection:
            if o and o.name in bpy.data.objects: o.select_set(True)
        context.view_layer.objects.active = prev_active
        return {'FINISHED'}

class BB_UVs_PackIndividually(Operator):
    bl_idname = "bb_uvs.normalize_pack_individually"
    bl_label = "Pack Individually"
    bl_options = {'REGISTER', 'UNDO'}
    def execute(self, context):
        args = {
            "udim_source": "ACTIVE_UDIM", "rotate": False, "scale": True,
            "margin_method": "FRACTION", "margin": 0.00390625, "pin": False,
            "merge_overlap": False, "pin_method": "LOCKED", "shape_method": "CONCAVE",
        }
        active = context.view_layer.objects.active
        if not _object_is_uv_mesh(active):
            self.report({'WARNING'}, "No active mesh with UVs"); return {'CANCELLED'}
        start_mode = getattr(active, "mode", "OBJECT")
        original_selection = context.selected_objects[:]
        original_active = active
        original_edit_objs = list(context.objects_in_mode) if start_mode == 'EDIT' else []
        if start_mode == 'EDIT':
            targets = [o for o in original_edit_objs if _object_is_uv_mesh(o)]
        else:
            targets = [o for o in original_selection if _object_is_uv_mesh(o)]
        if not targets:
            self.report({'WARNING'}, "No selected meshes with UVs to pack"); return {'CANCELLED'}
        def _normalize_single(obj):
            bm = bmesh.from_edit_mesh(obj.data)
            _normalize_islands_in_bmesh(bm, xy_independent=True, only_selected=False, respect_aspect=True)
            bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)
        for obj in targets:
            if getattr(context.view_layer.objects.active, "mode", "OBJECT") == 'EDIT':
                bpy.ops.object.mode_set(mode='OBJECT')
            for o in bpy.data.objects: o.select_set(False)
            obj.select_set(True); context.view_layer.objects.active = obj
            bpy.ops.object.mode_set(mode='EDIT')
            ts = context.tool_settings
            ts.mesh_select_mode = (False, False, True)
            bpy.ops.mesh.select_all(action='SELECT')
            try: bpy.ops.uv.select_all(action='SELECT')
            except Exception: pass
            _normalize_single(obj)
            BB_UVs_NormalizePack._call_pack(args)
            BB_UVs_NormalizePack._snap_uvs_to_first_udim_all()
            bpy.ops.object.mode_set(mode='OBJECT')
        for o in bpy.data.objects: o.select_set(False)
        for o in original_selection:
            if o and o.name in bpy.data.objects: o.select_set(True)
        context.view_layer.objects.active = original_active
        if start_mode == 'EDIT':
            for o in bpy.data.objects: o.select_set(False)
            for o in original_edit_objs:
                if o and o.name in bpy.data.objects: o.select_set(True)
            if original_edit_objs:
                context.view_layer.objects.active = original_edit_objs[0]
                bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}

# ------------------------------ Grid Helper Operators ------------------------------

class BB_UVs_ChooseGridImage(Operator):
    bl_idname = "bb_uvs.choose_grid_image"
    bl_label = "Browse File"
    bl_options = {'REGISTER', 'UNDO'}
    filepath: StringProperty(name="File Path", subtype='FILE_PATH')
    def execute(self, context):
        _S(context).bb_grid_image_path = bpy.path.abspath(self.filepath)
        self.report({'INFO'}, f"Grid image set: {os.path.basename(_S(context).bb_grid_image_path)}")
        return {'FINISHED'}
    def invoke(self, context, event):
        context.window_manager.fileselect_add(self); return {'RUNNING_MODAL'}

class BB_UVs_ApplyGrid(Operator):
    bl_idname = "bb_uvs.apply_grid"
    bl_label = "Apply Grid"
    bl_options = {'REGISTER', 'UNDO'}
    def execute(self, context):
        S = _S(context)
        img = _ensure_image(S.bb_grid_image_path)
        if img is None:
            self.report({'WARNING'}, "Invalid grid image path. Click 'Browse File' first.")
            return {'CANCELLED'}
        grid_mat = _ensure_grid_material(img, mat_name="BB_UVs_GridMat")
        mode = getattr(context.view_layer.objects.active, "mode", "OBJECT")
        if mode == 'EDIT':
            targets = [o for o in context.objects_in_mode if _object_is_uv_mesh(o)]
        else:
            targets = [o for o in context.selected_objects if _object_is_uv_mesh(o)]
        if not targets:
            self.report({'WARNING'}, "No mesh with UVs selected")
            return {'CANCELLED'}
        for obj in targets:
            _backup_material_mapping(obj)
            mats = obj.data.materials
            try:
                grid_idx = list(mats).index(grid_mat)
            except ValueError:
                mats.append(grid_mat); grid_idx = len(mats) - 1
            _set_all_faces_to_mat(obj, grid_idx)
        self.report({'INFO'}, f"Applied grid to {len(targets)} object(s).")
        return {'FINISHED'}

class BB_UVs_RevertMaterials(Operator):
    bl_idname = "bb_uvs.revert_materials"
    bl_label = "Revert"
    bl_options = {'REGISTER', 'UNDO'}
    def execute(self, context):
        mode = getattr(context.view_layer.objects.active, "mode", "OBJECT")
        if mode == 'EDIT':
            targets = [o for o in context.objects_in_mode if _object_is_uv_mesh(o)]
        else:
            targets = [o for o in context.selected_objects if _object_is_uv_mesh(o)]
        if not targets:
            self.report({'WARNING'}, "No mesh with UVs selected"); return {'CANCELLED'}
        restored = 0
        for obj in targets:
            if _restore_material_mapping(obj):
                restored += 1
        if restored == 0:
            self.report({'WARNING'}, "No backups found on selected objects."); return {'CANCELLED'}
        self.report({'INFO'}, f"Reverted {restored} object(s)."); return {'FINISHED'}

# ------------------------------ UI Helper ------------------------------

class BB_UVs_SetUI(Operator):
    bl_idname = "bb_uvs.set_ui"
    bl_label = "Set UI"
    bl_options = {'REGISTER', 'UNDO'}
    def execute(self, context):
        image_space = None
        for area in context.screen.areas:
            if area.type == 'IMAGE_EDITOR':
                image_space = area.spaces.active; break
        if image_space is None:
            self.report({'WARNING'}, "No UV/Image Editor found"); return {'CANCELLED'}
        image_space.uv_editor.edge_display_type = 'BLACK'
        image_space.uv_editor.tile_grid_shape = (10, 10)
        image_space.overlay.show_grid_background = True
        ts = context.tool_settings
        try: ts.mesh_select_mode = (False, True, False)
        except Exception: pass
        try: ts.uv_select_mode = 'EDGE'
        except Exception: pass
        self.report({'INFO'}, "UI set: UV grid on, selection set to Edge."); return {'FINISHED'}

# ------------------------------ UI Panels ------------------------------

def draw_bb_uvs_panel(self, context):
    S = context.scene
    layout = self.layout

    # Grid helper (BEFORE Set UI)
    layout.label(text="Grid Helper")
    row = layout.row(align=True)
    row.operator("bb_uvs.choose_grid_image", text="Browse File", icon='FILEBROWSER')
    if getattr(S, "bb_grid_image_path", ""):
        row = layout.row(align=True)
        row.label(text=os.path.basename(bpy.path.abspath(S.bb_grid_image_path)))
    row = layout.row(align=True)
    row.operator("bb_uvs.apply_grid", text="Apply Grid", icon='IMAGE_DATA')
    row.operator("bb_uvs.revert_materials", text="Revert", icon='LOOP_BACK')

    layout.separator()
    layout.label(text="Set UI")
    layout.operator("bb_uvs.set_ui", text="Set UI")

    layout.separator()
    layout.label(text="Pack")
    row = layout.row(align=True)
    row.operator("bb_uvs.normalize_pack_individually", text="Pack Individually")
    row.operator("bb_uvs.normalize_pack", text="Pack Together")

    layout.separator()
    row = layout.row(align=True)
    row.prop(S, "bb_density", text="Texel Density")
    row = layout.row(align=True)
    row.operator("bb_uvs.texel_density_check", text="Get TD")
    row.operator("bb_uvs.texel_density_set", text="Set TD")

    layout.separator()
    layout.label(text="Move")
    active = context.view_layer.objects.active
    mode = getattr(active, "mode", "OBJECT") if active else "OBJECT"
    row = layout.row(align=True); row.label(text="Context:")
    if mode == 'EDIT':
        left = row.row(align=True)
        op = left.operator("bb_uvs.set_move_context", text="Selected", icon='RESTRICT_SELECT_OFF', depress=S.bb_move_selected_only)
        op.set_selected = True
        right = row.row(align=True)
        op = right.operator("bb_uvs.set_move_context", text="Highlighted", icon='RADIOBUT_ON', depress=not S.bb_move_selected_only)
        op.set_selected = False
    else:
        # Object Mode: "Selected" = all selected objects; "Highlighted" = active only
        left = row.row(align=True)
        op = left.operator("bb_uvs.set_move_context", text="Selected", icon='RESTRICT_SELECT_OFF', depress=not S.bb_move_selected_only)
        op.set_selected = False  # move all selected objects

        right = row.row(align=True)
        op = right.operator("bb_uvs.set_move_context", text="Highlighted", icon='RADIOBUT_ON', depress=S.bb_move_selected_only)
        op.set_selected = True   # move only the active (highlighted) object

    row = layout.row(align=True)
    col = row.column(align=True); col.operator_context = 'EXEC_DEFAULT'
    op = col.operator("bb_uvs.move_uvs", text="↖"); op.direction = 'UP_LEFT'
    op = col.operator("bb_uvs.move_uvs", text="←");  op.direction = 'LEFT'
    op = col.operator("bb_uvs.move_uvs", text="↙");  op.direction = 'DOWN_LEFT'

    col = row.column(align=True)
    op = col.operator("bb_uvs.move_uvs", text="↑");  op.direction = 'UP'
    col.prop(S, "bb_move_amount", text="")
    op = col.operator("bb_uvs.move_uvs", text="↓");  op.direction = 'DOWN'

    col = row.column(align=True)
    op = col.operator("bb_uvs.move_uvs", text="↗");  op.direction = 'UP_RIGHT'
    op = col.operator("bb_uvs.move_uvs", text="→");  op.direction = 'RIGHT'
    op = col.operator("bb_uvs.move_uvs", text="↘");  op.direction = 'DOWN_RIGHT'

    layout.separator()
    layout.label(text="Rotate")
    r = layout.row(); r.scale_y = 1.1
    op = r.operator("bb_uvs.rotate_uvs", text="90° CCW", icon='LOOP_BACK');    op.direction = 'CCW'
    op = r.operator("bb_uvs.rotate_uvs", text="90° CW",  icon='LOOP_FORWARDS'); op.direction = 'CW'

    layout.separator()
    layout.label(text="Flip")
    f = layout.row(); f.scale_y = 1.1
    op = f.operator("bb_uvs.flip_uvs", text="Flip Horizontal", icon='ARROW_LEFTRIGHT'); op.axis = 'H'

class UV_PT_BB_UVs(Panel):
    bl_label = "BB UVs"
    bl_space_type = 'IMAGE_EDITOR'
    bl_region_type = 'UI'
    bl_category = "BB UVs"
    @classmethod
    def poll(cls, context):
        area = context.area
        return area and area.type == 'IMAGE_EDITOR' and getattr(context.space_data, "ui_type", "UV") == "UV"
    draw = draw_bb_uvs_panel

class VIEW3D_PT_BB_UVs(Panel):
    bl_label = "BB UVs"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "BB UVs"
    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return (context.mode in {'EDIT_MESH', 'OBJECT'} and
                obj and obj.type == 'MESH' and obj.data and
                len(obj.data.uv_layers) > 0)
    draw = draw_bb_uvs_panel

# ------------------------------ Registration (no PropertyGroup; add/remove Scene props) ------------------------------

def register():
    # Define Scene-level properties (these do not require a PropertyGroup class)
    bpy.types.Scene.bb_density = FloatProperty(
        name="Texel Density", default=0.0, min=0.0, precision=3, step=1,
        description="Texel density (px/cm)"
    )
    bpy.types.Scene.bb_move_amount = FloatProperty(
        name="Amount", default=1.0, min=0.0, soft_max=10.0,
        description="Distance in UDIM tiles (1.0 = one tile)"
    )
    bpy.types.Scene.bb_move_selected_only = BoolProperty(
        name="Contextual Toggle", default=True,
        description="Edit Mode: only selected UVs. Object Mode: only active mesh."
    )
    bpy.types.Scene.bb_grid_image_path = StringProperty(
        name="Grid Image", default="", subtype='FILE_PATH',
        description="Path to the grid/UV texture"
    )

    # Register operators/panels
    for cls in (
        BB_Texel_Density_Check, BB_Texel_Density_Set,
        BB_UVs_MoveUVs, BB_UVs_RotateUVs, BB_UVs_NormalizePack,
        BB_UVs_SetUI, BB_UVs_FlipUVs, BB_UVs_SetMoveContext, BB_UVs_ToggleMoveHighlight,
        BB_UVs_ChooseGridImage, BB_UVs_ApplyGrid, BB_UVs_RevertMaterials,
        UV_PT_BB_UVs, BB_UVs_PackIndividually, VIEW3D_PT_BB_UVs,
    ):
        try:
            bpy.utils.register_class(cls)
        except Exception:
            try:
                bpy.utils.unregister_class(cls)
            except Exception:
                pass
            bpy.utils.register_class(cls)

def unregister():
    # Remove Scene-level properties if they exist
    for prop in ("bb_density", "bb_move_amount", "bb_move_selected_only", "bb_grid_image_path"):
        if hasattr(bpy.types.Scene, prop):
            try:
                delattr(bpy.types.Scene, prop)
            except Exception:
                pass

    # Unregister operators/panels
    for cls in reversed((
        BB_Texel_Density_Check, BB_Texel_Density_Set,
        BB_UVs_MoveUVs, BB_UVs_RotateUVs, BB_UVs_NormalizePack,
        BB_UVs_SetUI, BB_UVs_FlipUVs, BB_UVs_SetMoveContext, BB_UVs_ToggleMoveHighlight,
        BB_UVs_ChooseGridImage, BB_UVs_ApplyGrid, BB_UVs_RevertMaterials,
        UV_PT_BB_UVs, BB_UVs_PackIndividually, VIEW3D_PT_BB_UVs,
    )):
        try:
            bpy.utils.unregister_class(cls)
        except Exception:
            pass

if __name__ == "__main__":
    register()
