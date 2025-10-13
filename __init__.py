# SPDX-License-Identifier: GPL-3.0-or-later
# -----------------------------------------------------------------------------
# BB UVs (Minimal + Move UVs — Selected/Highlighted Defaults)
# Derived work based on:
# • Texel Density Checker — https://github.com/mrven/Blender-Texel-Density-Checker
# Original author: Ivan Vostrikov (mrven)
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

# New operator for setting grid scale
class BB_UVs_SetGridScale(Operator):
    bl_idname = "bb_uvs.set_grid_scale"
    bl_label = "Set Grid Scale"
    bl_options = {'REGISTER', 'UNDO'}
    scale: EnumProperty(
        items=[('1', "1x", ""), ('2', "2x", ""), ('4', "4x", ""), ('6', "6x", "")],
        name="Scale", default='1'
    )
    def execute(self, context):
        mat = bpy.data.materials.get("BB_UVs_GridMat")
        if not mat or not mat.node_tree:
            self.report({'WARNING'}, "Grid material not found")
            return {'CANCELLED'}
        fac = float(self.scale)
        nt = mat.node_tree
        for n in nt.nodes:
            if n.type == 'MAPPING':
                n.inputs['Scale'].default_value = (fac, fac, 1.0)
                break
        context.view_layer.update()
        return {'FINISHED'}

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
            if o.mode == 'EDIT':
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

class BB_UVs_SetMoveCollection(Operator):
    bl_idname = "bb_uvs.set_move_collection"
    bl_label = "Set Move Collection"
    bl_options = {'INTERNAL', 'UNDO'}

    def execute(self, context):
        S = _S(context)
        S.bb_move_collection = not bool(S.bb_move_collection)
        if S.bb_move_collection:
            S.bb_move_selected_only = False
            active = context.view_layer.objects.active
            if active:
                colls = active.users_collection
                if colls:
                    col_objs = set(o for c in colls for o in c.objects if _object_is_uv_mesh(o))
                    #print(f"[Collection Toggle] Active: {active.name}, Collection objects with UVs: {[o.name for o in col_objs]}")
                    if col_objs:
                        bpy.ops.object.select_all(action='DESELECT')
                        for o in col_objs:
                            o.select_set(True)
                        context.view_layer.objects.active = active
               # else:
                    #print("[Collection Toggle] Active object is not in any collection")
            #else:
               # print("[Collection Toggle] No active object")
        #else:
            #print("[Collection Toggle] Collection mode OFF")
        return {'FINISHED'}

class BB_UVs_SetMoveContext(Operator):
    bl_idname = "bb_uvs.set_move_context"
    bl_label = "Set Move Context"
    bl_options = {'INTERNAL', 'UNDO'}
    set_selected: BoolProperty(name="Set Selected Mode", default=True)
    def execute(self, context):
        S = _S(context)
        S.bb_move_selected_only = bool(self.set_selected)
        S.bb_move_collection = False   # ← add this line here
        return {'FINISHED'}

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
        S = _S(context)
        active = context.active_object

        # Movement vector based on direction
        amt = S.bb_move_amount
        dx = dy = 0.0
        if self.direction == 'UP': dy = amt
        elif self.direction == 'DOWN': dy = -amt
        elif self.direction == 'LEFT': dx = -amt
        elif self.direction == 'RIGHT': dx = amt
        elif self.direction == 'UP_LEFT': dx, dy = -amt, amt
        elif self.direction == 'UP_RIGHT': dx, dy = amt, amt
        elif self.direction == 'DOWN_LEFT': dx, dy = -amt, -amt
        elif self.direction == 'DOWN_RIGHT': dx, dy = amt, -amt

        targets = []

        # 🟢 Match logic from working Move UVs operator
        if S.bb_move_selected_only:
            if _object_is_uv_mesh(active):
                targets = [active]

        elif S.bb_move_collection:
            # Only move UVs of objects in the same collection as the active object
            if active and active.users_collection:
                active_coll = active.users_collection[0]  # same logic as working button
                targets = [
                    o for o in active_coll.objects
                    if _object_is_uv_mesh(o)
                ]

        else:
            # Default: all selected mesh objects
            targets = [
                o for o in context.selected_objects
                if _object_is_uv_mesh(o)
            ]

        # Apply the move
        for o in targets:
            _move_uvs_object(o, dx, dy)

        if targets:
            context.view_layer.update()

        return {'FINISHED'}


class BB_UVs_SelectSimilarTopology(Operator):
    bl_idname = "bb_uvs.select_similar_topology"
    bl_label = "Sim"
    bl_options = {'REGISTER', 'UNDO'}
    include_scale: BoolProperty(default=False)

    def execute(self, context):
        active = context.active_object
        if not _object_is_uv_mesh(active):
            self.report({'WARNING'}, "Active object is not a valid mesh with UVs")
            return {'CANCELLED'}
        ref_obj = active
        ref_verts = len(ref_obj.data.vertices)
        ref_faces = len(ref_obj.data.polygons)
        ref_scale = ref_obj.scale.copy()

        bpy.ops.object.select_all(action='DESELECT')
        ref_obj.select_set(True)
        context.view_layer.objects.active = ref_obj

        for obj in context.view_layer.objects:
            if obj == ref_obj or not _object_is_uv_mesh(obj):
                continue
            if len(obj.data.vertices) != ref_verts or len(obj.data.polygons) != ref_faces:
                continue
            if self.include_scale and not all(isclose(obj.scale[i], ref_scale[i], rel_tol=1e-6) for i in range(3)):
                continue
            obj.select_set(True)
        return {'FINISHED'}

class BB_UVs_SelectCrossUDIMs(Operator):
    bl_idname = "bb_uvs.select_cross_udims"
    bl_label = "Cross UDIMs"
    bl_description = "Select UV faces that span across multiple UDIM tiles"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return context.mode == 'EDIT_MESH' and (obj := context.active_object) and obj.type == 'MESH'

    def execute(self, context):
        total_counter = 0
        sync = context.tool_settings.use_uv_select_sync

        # Deselect based on sync mode
        if sync:
            bpy.ops.mesh.select_all(action='DESELECT')
        else:
            bpy.ops.uv.select_all(action='DESELECT')

        for obj in context.objects_in_mode:
            if not _object_is_uv_mesh(obj):
                continue

            bm = bmesh.from_edit_mesh(obj.data)
            uv_layer = bm.loops.layers.uv.active
            if not uv_layer:
                continue

            local_counter = 0
            for f in bm.faces:
                if sync and f.hide:
                    continue

                # Collect all UDIM tile coordinates (floor of UVs)
                tiles = set()
                for loop in f.loops:
                    uv = loop[uv_layer].uv
                    tile_u = math.floor(uv.x)
                    tile_v = math.floor(uv.y)
                    tiles.add((tile_u, tile_v))
                    if len(tiles) > 1:
                        # Early exit: already crosses UDIMs
                        break

                if len(tiles) > 1:
                    if sync:
                        f.select = True
                    else:
                        for loop in f.loops:
                            loop[uv_layer].select = True
                            loop[uv_layer].select_edge = True
                    local_counter += 1

            total_counter += local_counter
            bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)

        if total_counter == 0:
            self.report({'INFO'}, "No faces crossing UDIM boundaries found")
        else:
            self.report({'INFO'}, f"Selected {total_counter} face(s) crossing UDIMs")

        return {'FINISHED'}
    
class BB_UVs_SelectSimilarTopologyScale(Operator):
    bl_idname = "bb_uvs.select_similar_topology_scale"
    bl_label = "Sim All"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        bpy.ops.bb_uvs.select_similar_topology(include_scale=True)
        return {'FINISHED'}
    
class BB_UVs_MoveUVsInteractive(bpy.types.Operator):
    bl_idname = "bb_uvs.move_uvs_interactive"
    bl_label = "Move UVs (Interactive Toggle)"
    bl_options = {'REGISTER', 'UNDO', 'GRAB_CURSOR', 'BLOCKING'}

    _dragging = False
    _last_mouse = None
    _total_dx = 0.0
    _total_dy = 0.0

    def invoke(self, context, event):
        S = _S(context)
        if context.mode != 'OBJECT':
            self.report({'WARNING'}, "Interactive move only available in Object Mode")
            return {'CANCELLED'}
        self._dragging = False
        self._last_mouse = (event.mouse_x, event.mouse_y)
        self._total_dx = 0.0
        self._total_dy = 0.0
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}

    def modal(self, context, event):
        S = _S(context)
        active = context.view_layer.objects.active

        if event.type == 'MOUSEMOVE':
            if not self._dragging:
                return {'RUNNING_MODAL'}
            delta_x = event.mouse_x - self._last_mouse[0]
            delta_y = event.mouse_y - self._last_mouse[1]
            self._last_mouse = (event.mouse_x, event.mouse_y)
            sensitivity = 0.02  # Adjust as needed (pixels to UDIM fraction)
            dx = delta_x * sensitivity
            dy = delta_y * sensitivity
            self._total_dx += dx
            self._total_dy += dy

            # Apply incremental move
            if S.bb_move_collection and active:
                colls = active.users_collection
                if colls:
                    col_objs = set(o for c in colls for o in c.objects if _object_is_uv_mesh(o))
                    for o in col_objs:
                        _move_uvs_object(o, dx, dy)
                    context.view_layer.update()

            elif S.bb_move_selected_only:
                if _object_is_uv_mesh(active):
                    _move_uvs_object(active, dx, dy)
                    context.view_layer.update()

            else:
                targets = [o for o in context.selected_objects if _object_is_uv_mesh(o)]
                for o in targets:
                    _move_uvs_object(o, dx, dy)
                context.view_layer.update()

        elif event.type == 'LEFTMOUSE':
            if event.value == 'PRESS':
                self._dragging = True
                self._last_mouse = (event.mouse_x, event.mouse_y)
                return {'RUNNING_MODAL'}
            elif event.value == 'RELEASE':
                if self._dragging:
                    self._dragging = False
                return {'FINISHED'}

        elif event.type in {'RIGHTMOUSE', 'ESC'}:
            # Revert total move
            neg_dx = -self._total_dx
            neg_dy = -self._total_dy
            if S.bb_move_collection and active:
                colls = active.users_collection
                if colls:
                    col_objs = set(o for c in colls for o in c.objects if _object_is_uv_mesh(o))
                    for o in col_objs:
                        _move_uvs_object(o, neg_dx, neg_dy)
                    context.view_layer.update()

            elif S.bb_move_selected_only:
                if _object_is_uv_mesh(active):
                    _move_uvs_object(active, neg_dx, neg_dy)
                    context.view_layer.update()

            else:
                targets = [o for o in context.selected_objects if _object_is_uv_mesh(o)]
                for o in targets:
                    _move_uvs_object(o, neg_dx, neg_dy)
                context.view_layer.update()

            self._dragging = False
            return {'CANCELLED'}

        return {'RUNNING_MODAL'}
# ------------------------------ Rotate / Flip / Normalize + Pack ------------------------------

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
        mode = context.mode
        if mode == 'EDIT_MESH':
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
        mode = context.mode
        if mode == 'EDIT_MESH':
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
                 "margin": args.get("margin", 0.00195312)}
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
        S = _S(context)
        args = {
            "udim_source": "ACTIVE_UDIM", "rotate": S.bb_pack_rotate, "scale": not S.bb_keep_td,
            "margin_method": "FRACTION", "margin": 0.00390625, "pin": False,
            "merge_overlap": False, "pin_method": "LOCKED", "shape_method": "CONCAVE",
        }
        active = context.view_layer.objects.active
        if not _object_is_uv_mesh(active):
            self.report({'WARNING'}, "No active mesh with UVs"); return {'CANCELLED'}
        start_mode = context.mode
        if start_mode == 'EDIT_MESH':
            ts = context.tool_settings
            ts.mesh_select_mode = (False, False, True)
            bpy.ops.mesh.select_all(action='SELECT')
            try: bpy.ops.uv.select_all(action='SELECT')
            except Exception: pass
            if not S.bb_keep_td:
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
        if not S.bb_keep_td:
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
        S = _S(context)
        args = {
            "udim_source": "ACTIVE_UDIM", "rotate": S.bb_pack_rotate, "scale": not S.bb_keep_td,
            "margin_method": "FRACTION", "margin": 0.00390625, "pin": False,
            "merge_overlap": False, "pin_method": "LOCKED", "shape_method": "CONCAVE",
        }
        active = context.view_layer.objects.active
        if not _object_is_uv_mesh(active):
            self.report({'WARNING'}, "No active mesh with UVs"); return {'CANCELLED'}
        start_mode = context.mode
        original_selection = context.selected_objects[:]
        original_active = active
        original_edit_objs = list(context.objects_in_mode) if start_mode == 'EDIT_MESH' else []
        if start_mode == 'EDIT_MESH':
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
            if context.mode == 'EDIT_MESH':
                bpy.ops.object.mode_set(mode='OBJECT')
            for o in bpy.data.objects: o.select_set(False)
            obj.select_set(True); context.view_layer.objects.active = obj
            bpy.ops.object.mode_set(mode='EDIT')
            ts = context.tool_settings
            ts.mesh_select_mode = (False, False, True)
            bpy.ops.mesh.select_all(action='SELECT')
            try: bpy.ops.uv.select_all(action='SELECT')
            except Exception: pass
            if not S.bb_keep_td:
                _normalize_single(obj)
            BB_UVs_NormalizePack._call_pack(args)
            BB_UVs_NormalizePack._snap_uvs_to_first_udim_all()
            bpy.ops.object.mode_set(mode='OBJECT')
        for o in bpy.data.objects: o.select_set(False)
        for o in original_selection:
            if o and o.name in bpy.data.objects: o.select_set(True)
        context.view_layer.objects.active = original_active
        if start_mode == 'EDIT_MESH':
            for o in bpy.data.objects: o.select_set(False)
            for o in original_edit_objs:
                if o and o.name in bpy.data.objects: o.select_set(True)
            if original_edit_objs:
                context.view_layer.objects.active = original_edit_objs[0]
                bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}

class BB_UVs_SelectFlipped(Operator):
    bl_idname = "bb_uvs.select_flipped"
    bl_label = "Select Flipped"
    bl_description = "Select flipped UV faces"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return context.mode == 'EDIT_MESH' and (obj := context.active_object) and obj.type == 'MESH'

    def execute(self, context):
        total_counter = 0
        
        # Get sync setting
        sync = context.tool_settings.use_uv_select_sync
        
        # Deselect all first based on sync mode
        if sync:
            bpy.ops.mesh.select_all(action='DESELECT')
        else:
            bpy.ops.uv.select_all(action='DESELECT')
        
        # Process all objects in edit mode
        for obj in context.objects_in_mode:
            if not _object_is_uv_mesh(obj):
                continue
                
            bm = bmesh.from_edit_mesh(obj.data)
            uv_layer = bm.loops.layers.uv.active
            if not uv_layer:
                continue
                
            # Calculate loop triangles for accurate area detection
            loop_triangles = bm.calc_loop_triangles()
            
            local_counter = 0
            for tris in loop_triangles:
                # Check if we should process this face based on sync mode
                if sync:
                    if tris[0].face.hide:
                        continue
                else:
                    # In non-sync mode, we want to check ALL faces, not just selected ones
                    # So we don't filter by selection here
                    pass
                
                # Get UV coordinates for the triangle
                a = tris[0][uv_layer].uv
                b = tris[1][uv_layer].uv
                c = tris[2][uv_layer].uv

                # Use the exact same area calculation as UniV
                area = a.cross(b) + b.cross(c) + c.cross(a)
                
                # If area is negative, the face is flipped
                if area < 0:
                    if sync:
                        # In sync mode, select the mesh face
                        tris[0].face.select = True
                    else:
                        # In non-sync mode, select the UVs
                        for loop in tris[0].face.loops:
                            loop[uv_layer].select = True
                            loop[uv_layer].select_edge = True  # Important: also select edges
                    local_counter += 1
            
            total_counter += local_counter
            bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)
        
        if not total_counter:
            self.report({'INFO'}, 'Flipped faces not found')
        else:
            self.report({'INFO'}, f'Found {total_counter} flipped faces')
            
        return {'FINISHED'}

class BB_UVs_SelectBoundaries(bpy.types.Operator):
    bl_idname = "bb_uvs.select_boundaries"
    bl_label = "Select Boundaries"
    bl_description = "Select internal UV island boundary edges (holes)"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return context.mode == 'EDIT_MESH' and (obj := context.active_object) and obj.type == 'MESH'

    def execute(self, context):
        from collections import defaultdict
        import math

        total_selected = 0
        sync = context.tool_settings.use_uv_select_sync
        if sync:
            self.report({'WARNING'}, "Disable UV Sync for UV boundary selection")
            return {'CANCELLED'}

        bpy.ops.uv.select_all(action='DESELECT')

        for obj in context.objects_in_mode:
            if not _object_is_uv_mesh(obj):
                continue

            bm = bmesh.from_edit_mesh(obj.data)
            uv_layer = bm.loops.layers.uv.active
            if not uv_layer:
                continue

            islands = _collect_uv_islands(bm, uv_layer, only_selected=False)
            for island in islands:
                # --- build oriented-edge map and usage counts (unordered) ---
                counts = defaultdict(int)           # keyed by frozenset({k1,k2})
                oriented = {}                       # keyed by (k1,k2) -> (face, i)
                uv_by_key = {}                      # keyed by uv key -> Vector2 (uv coords)

                for face in island:
                    loops = face.loops
                    n = len(loops)
                    for i in range(n):
                        l1 = loops[i]
                        l2 = loops[(i + 1) % n]
                        k1 = _uv_key(l1[uv_layer].uv)
                        k2 = _uv_key(l2[uv_layer].uv)
                        counts[frozenset((k1, k2))] += 1
                        oriented[(k1, k2)] = (face, i)
                        # store uv coords for both keys (overwrite is fine)
                        uv_by_key[k1] = l1[uv_layer].uv.copy()
                        uv_by_key[k2] = l2[uv_layer].uv.copy()

                # boundary = unordered edges used exactly once
                boundary_oriented = {}
                for (a, b), (face, i) in oriented.items():
                    if counts[frozenset((a, b))] == 1:
                        boundary_oriented[(a, b)] = (face, i)

                if not boundary_oriented:
                    continue

                # adjacency of boundary graph using oriented vertex keys
                adj = defaultdict(list)
                for (a, b) in boundary_oriented.keys():
                    adj[a].append(b)
                    # note: do not add adj[b].append(a) here as we keep orientation edges present
                    # but for navigation we want to know all neighbors at a vertex:
                    # add both so we can examine candidate outgoing edges from any vertex
                    adj[b]  # ensure key exists

                # However we do need neighbor sets for choosing next directed edge:
                # Build neighbor list per vertex from the oriented boundary edges (outgoing)
                outgoing = defaultdict(list)
                for (a, b) in boundary_oriented.keys():
                    outgoing[a].append(b)

                # used directed edges set
                used = set()
                loops = []

                # helper: compute signed angle from v_prev->v_curr to v_curr->v_next in UV space
                def angle_between(prev_k, curr_k, next_k):
                    p = uv_by_key[prev_k]
                    c = uv_by_key[curr_k]
                    npt = uv_by_key[next_k]
                    v1x = c.x - p.x
                    v1y = c.y - p.y
                    v2x = npt.x - c.x
                    v2y = npt.y - c.y
                    # compute atan2 of cross/dot to get signed angle
                    cross = v1x * v2y - v1y * v2x
                    dot = v1x * v2x + v1y * v2y
                    return math.atan2(cross, dot)

                # Walk every directed boundary edge that hasn't been used yet.
                for start_edge in list(boundary_oriented.keys()):
                    if start_edge in used:
                        continue

                    # initialize traversal from this directed edge
                    a, b = start_edge
                    cur_prev = a
                    cur = b
                    path = [start_edge]
                    used.add(start_edge)

                    # continue walking
                    safety = 0
                    while True:
                        safety += 1
                        if safety > 10000:
                            break  # avoid infinite loops in degenerate cases

                        # if we returned to the starting vertex and path forms a closed ring -> done
                        if cur == start_edge[0] and len(path) > 2:
                            # convert path of directed edges to list of (face,i) using boundary_oriented
                            loop_face_idx = [boundary_oriented[e] for e in path]
                            loops.append(loop_face_idx)
                            break

                        # select outgoing candidates from 'cur' that haven't been used (directed)
                        candidates = [v for v in outgoing[cur] if (cur, v) not in used]

                        if not candidates:
                            # no unused outgoing edges: maybe there is an incoming edge continuing the loop
                            # look for any neighbor (either outgoing or incoming) that forms an unused directed edge
                            # incoming neighbors are vertices 'x' for which (x, cur) exists in boundary_oriented
                            incoming = [x for (x, y) in boundary_oriented.keys() if y == cur and (x, cur) not in used]
                            candidates = incoming

                        if not candidates:
                            # dead end: abandon this path
                            break

                        # if multiple candidates pick the one with smallest left-turn angle (prefer continuity)
                        best = None
                        best_angle = None
                        for cand in candidates:
                            try:
                                ang = angle_between(cur_prev, cur, cand)
                            except Exception:
                                ang = 0.0
                            # We prefer angles that keep us moving forward around the loop:
                            # choose the candidate with the smallest absolute angle (continuation)
                            if best is None or abs(ang) < best_angle:
                                best = cand
                                best_angle = abs(ang)

                        next_edge = (cur, best)
                        # if next_edge isn't an oriented boundary edge (maybe it's an incoming variant), try the other direction
                        if next_edge not in boundary_oriented:
                            # try to find an oriented key that exists (best might be incoming)
                            possible = None
                            if (best, cur) in boundary_oriented and (best, cur) not in used:
                                # reverse direction available (we'll follow that directed edge)
                                possible = (best, cur)
                            if possible is None:
                                break
                            next_edge = possible

                        # mark and append
                        used.add(next_edge)
                        path.append(next_edge)
                        cur_prev, cur = next_edge[0], next_edge[1]

                if len(loops) < 2:
                    # no holes in this island
                    continue

                # compute signed areas to distinguish outer vs holes
                def signed_area(coords):
                    area = 0.0
                    n = len(coords)
                    for i in range(n - 1):
                        area += coords[i].x * coords[i+1].y - coords[i+1].x * coords[i].y
                    return area * 0.5

                signed_areas = []
                for loop in loops:
                    coords = []
                    for face, i in loop:
                        coords.append(face.loops[i][uv_layer].uv.copy())
                    # close loop
                    if coords and coords[0] != coords[-1]:
                        coords.append(coords[0])
                    if len(coords) >= 3:
                        signed_areas.append(signed_area(coords))
                    else:
                        signed_areas.append(0.0)

                abs_areas = [abs(a) for a in signed_areas]
                outer_idx = max(range(len(abs_areas)), key=lambda i: abs_areas[i])
                outer_sign = math.copysign(1, signed_areas[outer_idx])

                # select edges that are opposite winding to the outer loop (holes)
                for idx, area in enumerate(signed_areas):
                    if idx == outer_idx:
                        continue
                    if math.copysign(1, area) != outer_sign:
                        for face, i in loops[idx]:
                            face.loops[i][uv_layer].select_edge = True
                            total_selected += 1

            bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)

        if total_selected == 0:
            self.report({'INFO'}, "No internal boundaries selected")
        else:
            self.report({'INFO'}, f"Selected {total_selected} internal boundary edges")
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
        mode = context.mode
        if mode == 'EDIT_MESH':
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
        mode = context.mode
        if mode == 'EDIT_MESH':
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

class BB_UVs_SelectZeroArea(Operator):
    bl_idname = "bb_uvs.select_zero_area"
    bl_label = "Select 0 Area"
    bl_description = "Select UV faces with zero or near-zero surface area in UV space"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return context.mode == 'EDIT_MESH' and (obj := context.active_object) and obj.type == 'MESH'

    def execute(self, context):
        sync = context.tool_settings.use_uv_select_sync
        if sync:
            self.report({'WARNING'}, "Disable UV Sync for UV area selection")
            return {'CANCELLED'}

        bpy.ops.uv.select_all(action='DESELECT')

        total = 0
        for obj in context.objects_in_mode:
            if not _object_is_uv_mesh(obj):
                continue
            bm = bmesh.from_edit_mesh(obj.data)
            uv_layer = bm.loops.layers.uv.active
            if not uv_layer:
                continue

            for face in bm.faces:
                area = _poly_uv_area(face, uv_layer)
                if area <= 1e-8:  # near-zero UV area
                    for loop in face.loops:
                        loop[uv_layer].select = True
                        loop[uv_layer].select_edge = True
                    total += 1

            bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)

        if total == 0:
            self.report({'INFO'}, "No zero-area UV faces found")
        else:
            self.report({'INFO'}, f"Selected {total} zero-area UV face(s)")
        return {'FINISHED'}
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

    if context.space_data.type == 'IMAGE_EDITOR':
        layout.separator()
        layout.label(text="Set UI")
        layout.operator("bb_uvs.set_ui", text="Set UI")
    
    layout.separator()
    if context.mode == 'OBJECT':
        layout.label(text="Select Similar")
        row = layout.row(align=True)
        row.operator("bb_uvs.select_similar_topology", text="Topology")
        row.operator("bb_uvs.select_similar_topology_scale", text="Topology + Scale")
    
    if context.space_data.type == 'VIEW_3D':    
        layout.label(text="Grid Helper")
        row = layout.row(align=True)
        row.operator("bb_uvs.choose_grid_image", text="Browse File", icon='FILEBROWSER')
        if getattr(S, "bb_grid_image_path", ""):
            row = layout.row(align=True)
            row.label(text=os.path.basename(bpy.path.abspath(S.bb_grid_image_path)))
        row = layout.row(align=True)
        row.operator("bb_uvs.apply_grid", text="Apply Grid", icon='IMAGE_DATA')
        row.operator("bb_uvs.revert_materials", text="Revert", icon='LOOP_BACK')
        row = layout.row(align=True)
        op = row.operator("bb_uvs.set_grid_scale", text="1x"); op.scale = '1'
        op = row.operator("bb_uvs.set_grid_scale", text="2x"); op.scale = '2'
        op = row.operator("bb_uvs.set_grid_scale", text="4x"); op.scale = '4'
        op = row.operator("bb_uvs.set_grid_scale", text="6x"); op.scale = '6'

    #layout.separator()
    layout.label(text="Pack")
    row = layout.row(align=True)
    row.prop(S, "bb_keep_td", text="Keep TD")
    row.prop(S, "bb_pack_rotate", text="Rotation")  # ← Add this line
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

    # 🔹 New interactive drag button (visible only in UV editor, disabled in Edit Mode)
    if context.area.type == 'IMAGE_EDITOR':
        row = layout.row(align=True)
        row.enabled = (context.mode == 'OBJECT')
        row.operator("bb_uvs.move_uvs_interactive",
                     text="Move UVs (Drag)",
                     icon='ARROW_LEFTRIGHT')

    # --- Existing move context and arrows ---
    active = getattr(context, "active_object", None)
    mode = context.mode

    layout.separator()
    row = layout.row(align=True)
    row.label(text="Context:")
    if mode == 'EDIT_MESH':
        left = row.row(align=True)
        op = left.operator("bb_uvs.set_move_context", text="Selected", icon='RESTRICT_SELECT_OFF', depress=not S.bb_move_selected_only)
        op.set_selected = False  # Selected mode = move all objects

        right = row.row(align=True)  
        op = right.operator("bb_uvs.set_move_context", text="Highlighted", icon='RADIOBUT_ON', depress=S.bb_move_selected_only)
        op.set_selected = True 
    else:
        # Object Mode: choose one context: Selected / Highlighted / Collection
        row = layout.row(align=True)

        # Selected (move all selected objects)
        op = row.operator("bb_uvs.set_move_context", text="Selected", icon='RESTRICT_SELECT_OFF',
                          depress=(not S.bb_move_selected_only and not S.bb_move_collection))
        op.set_selected = False  # move all selected objects

        # Highlighted (move only the active)
        op = row.operator("bb_uvs.set_move_context", text="Highlighted", icon='RADIOBUT_ON',
                          depress=(S.bb_move_selected_only and not S.bb_move_collection))
        op.set_selected = True   # move only the active (highlighted) object

        # Collection (toggle) — now an operator so it behaves like the other buttons
        op = row.operator("bb_uvs.set_move_collection", text="Collection", icon='OUTLINER_COLLECTION',
                          depress=S.bb_move_collection)


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
    # Only show "Check and fix" in the UV Editor (IMAGE_EDITOR)
    if context.space_data.type == 'IMAGE_EDITOR':
        layout.separator()
        layout.label(text="Check and fix")
        if context.mode == 'EDIT_MESH':
            f = layout.row(); f.scale_y = 1.1
            f = layout.row(); f.scale_y = 1.1
            f.operator("bb_uvs.select_cross_udims", text="", icon='SNAP_GRID')
            f.operator("bb_uvs.select_boundaries", text="", icon='EDGESEL')
            f.operator("uv.select_overlap", text="", icon='SELECT_INTERSECT')
            f.operator("bb_uvs.select_zero_area", text="", icon='CANCEL')
            f = layout.row(); f.scale_y = 1.1
            f.operator("bb_uvs.select_flipped", text="Select Flipped", icon='UV_SYNC_SELECT')
            op = f.operator("bb_uvs.flip_uvs", text="Flip Horizontal", icon='ARROW_LEFTRIGHT')
            op.axis = 'H'

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

import bpy
from bpy.props import FloatProperty, BoolProperty, StringProperty

# -------------------------------------------------------------------------
# REGISTER / UNREGISTER
# -------------------------------------------------------------------------

def register():
    # Define Scene-level properties (these do not require a PropertyGroup class)
    bpy.types.Scene.bb_density = FloatProperty(
        name="Texel Density",
        default=0.0, min=0.0, precision=3, step=1,
        description="Texel density (px/cm)"
    )
    bpy.types.Scene.bb_move_amount = FloatProperty(
        name="Amount",
        default=1.0, min=0.0, soft_max=10.0,
        description="Distance in UDIM tiles (1.0 = one tile)"
    )
    bpy.types.Scene.bb_move_selected_only = BoolProperty(
        name="Contextual Toggle",
        default=True,
        description="Edit Mode: only selected UVs. Object Mode: only active mesh."
    )
    bpy.types.Scene.bb_grid_image_path = StringProperty(
        name="Grid Image",
        default="", subtype='FILE_PATH',
        description="Path to the grid/UV texture"
    )
    bpy.types.Scene.bb_keep_td = BoolProperty(
        name="Keep TD",
        default=False,
        description="Keep current texel density (skip normalization and scaling during pack)"
    )
    bpy.types.Scene.bb_move_collection = BoolProperty(
        name="Collection",
        default=False,
        description="When ON, move all UVs from objects in the active object's collection"
    )
    bpy.types.Scene.bb_pack_rotate = BoolProperty(
        name="Rotation",
        default=True,
        description="Allow rotation during packing for better efficiency"
    )

    # ---------------------------------------------------------------------
    # Register all operators and panels
    # ---------------------------------------------------------------------
    classes = (
        BB_Texel_Density_Check, BB_Texel_Density_Set,
        BB_UVs_MoveUVs, BB_UVs_RotateUVs, BB_UVs_NormalizePack,
        BB_UVs_SetUI, BB_UVs_FlipUVs, BB_UVs_SetMoveContext, BB_UVs_ToggleMoveHighlight,
        BB_UVs_ChooseGridImage, BB_UVs_ApplyGrid, BB_UVs_RevertMaterials,
        BB_UVs_SetGridScale, BB_UVs_MoveUVsInteractive,
        UV_PT_BB_UVs, BB_UVs_PackIndividually, VIEW3D_PT_BB_UVs,
        BB_UVs_SetMoveCollection,
        BB_UVs_SelectSimilarTopology,
        BB_UVs_SelectSimilarTopologyScale,
        BB_UVs_SelectFlipped,
        BB_UVs_SelectCrossUDIMs,
        BB_UVs_SelectBoundaries,
        BB_UVs_SelectZeroArea,
    )

    for cls in classes:
        try:
            bpy.utils.register_class(cls)
        except ValueError:
            # If already registered (e.g. live-reload), unregister first
            bpy.utils.unregister_class(cls)
            bpy.utils.register_class(cls)


def unregister():
    # ---------------------------------------------------------------------
    # Remove Scene-level properties
    # ---------------------------------------------------------------------
    props = (
        "bb_density",
        "bb_move_amount",
        "bb_move_selected_only",
        "bb_grid_image_path",
        "bb_keep_td",
        "bb_move_collection",
        "bb_pack_rotate",
    )
    for prop in props:
        if hasattr(bpy.types.Scene, prop):
            try:
                delattr(bpy.types.Scene, prop)
            except Exception:
                pass

    # ---------------------------------------------------------------------
    # Unregister classes in reverse order
    # ---------------------------------------------------------------------
    classes = (
        BB_Texel_Density_Check, BB_Texel_Density_Set,
        BB_UVs_MoveUVs, BB_UVs_RotateUVs, BB_UVs_NormalizePack,
        BB_UVs_SetUI, BB_UVs_FlipUVs, BB_UVs_SetMoveContext, BB_UVs_ToggleMoveHighlight,
        BB_UVs_ChooseGridImage, BB_UVs_ApplyGrid, BB_UVs_RevertMaterials,
        BB_UVs_SetGridScale, BB_UVs_MoveUVsInteractive,
        UV_PT_BB_UVs, BB_UVs_PackIndividually, VIEW3D_PT_BB_UVs,
        BB_UVs_SetMoveCollection,
        BB_UVs_SelectSimilarTopology,
        BB_UVs_SelectSimilarTopologyScale,
        BB_UVs_SelectFlipped,
        BB_UVs_SelectCrossUDIMs,
        BB_UVs_SelectBoundaries,
        BB_UVs_SelectZeroArea,
    )
    for cls in reversed(classes):
        try:
            bpy.utils.unregister_class(cls)
        except Exception:
            pass


# Allow running directly from Text Editor (optional)
if __name__ == "__main__":
    register()
