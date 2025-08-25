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
    "name": "BB UVs (Minimal + Move/Pack)",
    "author": "BB",
    "version": (1, 4, 1),
    "blender": (4, 5, 0),
    "location": "UV Editor & 3D View > Sidebar (N) > BB UVs",
    "description": "Minimal Get/Set Texel Density + Move UVs + Normalize+Pack (hard-coded)",
    "category": "UV",
}

import bpy
import bmesh
import math
from math import isclose
from bpy.props import (
    PointerProperty, StringProperty, BoolProperty, FloatProperty, EnumProperty
)
from bpy.types import PropertyGroup, Panel, Operator

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def _object_is_uv_mesh(o):
    return (o and o.type == 'MESH' and o.data and len(o.data.uv_layers) > 0 and len(o.data.polygons) > 0)

def Calculate_TD_Area_To_List():
    """Return [texel_density, uv_area] per polygon of ACTIVE object without changing selections."""
    active_obj = bpy.context.active_object
    if not _object_is_uv_mesh(active_obj):
        return []

    # Use a fixed square texture for TD math (px/cm)
    w = h = 2048  # was 1024
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

        # Shoelace area in UV space
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
    """Move UVs of a mesh in object data (not edit mode)."""
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
            uv.x += dx
            uv.y += dy
    bm.to_mesh(o.data)
    o.data.update()
    bm.free()

def _move_uvs_edit(o, dx, dy, only_selected=False):
    """Move UVs of an edit-mode mesh. If only_selected=True, move only selected UV loops."""
    bm = bmesh.from_edit_mesh(o.data)
    uv_layer = bm.loops.layers.uv.active
    if uv_layer is None:
        return
    bm.faces.ensure_lookup_table()

    if not only_selected:
        for f in bm.faces:
            for loop in f.loops:
                uv = loop[uv_layer].uv
                uv.x += dx
                uv.y += dy
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
                uv.x += dx
                uv.y += dy
        bmesh.update_edit_mesh(o.data, loop_triangles=False, destructive=False)
        return

    v_mode, e_mode, f_mode = tool.mesh_select_mode
    if f_mode:
        for f in bm.faces:
            if not f.select:
                continue
            for loop in f.loops:
                uv = loop[uv_layer].uv
                uv.x += dx
                uv.y += dy
    else:
        for f in bm.faces:
            for loop in f.loops:
                if not loop.vert.select:
                    continue
                uv = loop[uv_layer].uv
                uv.x += dx
                uv.y += dy

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
    """Return list of sets of faces; each set is a UV island."""
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
            sx += uv.x
            sy += uv.y
            cnt += 1.0
    if cnt == 0:
        return 0.0, 0.0
    return sx / cnt, sy / cnt

def _normalize_islands_in_bmesh(bm, xy_independent=True, only_selected=True, respect_aspect=True):
    """Normalize UV islands scale based on 3D/UV area ratio."""
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
            tot_area_3d += a3
            tot_area_uv += au

    if tot_area_uv <= 1e-12 or tot_area_3d <= 1e-12:
        return

    tot_fac = tot_area_3d / tot_area_uv

    for faces in islands:
        area3 = 0.0
        areauv = 0.0
        for f in faces:
            area3 += f.calc_area()
            areauv += _poly_uv_area(f, uv_layer)

        if areauv <= 1e-12 or area3 <= 1e-12:
            continue

        fac = area3 / areauv
        if xy_independent:
            scale = math.sqrt(fac / tot_fac)
            su, sv = scale, scale
        else:
            scale = math.sqrt(fac / tot_fac)
            su = sv = scale

        cx, cy = _center_of_island_uv(faces, uv_layer)
        if isclose(su, 1.0, abs_tol=1e-6) and isclose(sv, 1.0, abs_tol=1e-6):
            continue

        for f in faces:
            for l in f.loops:
                uv = l[uv_layer].uv
                uv.x = (uv.x - cx) * su + cx
                uv.y = (uv.y - cy) * sv + cy


def _rotate_uvs_group_in_bmesh(bm, angle_deg=90.0, only_selected=True):
    """Rotate UVs around the centroid of the affected UV coordinates as ONE group."""
    uv_layer = bm.loops.layers.uv.active
    if not uv_layer:
        return

    loops_to_rotate = []
    tool = bpy.context.tool_settings
    uv_sync = getattr(tool, "use_uv_select_sync", False)

    if not only_selected:
        for f in bm.faces:
            for l in f.loops:
                loops_to_rotate.append(l)
    else:
        if not uv_sync:
            for f in bm.faces:
                for l in f.loops:
                    if l[uv_layer].select:
                        loops_to_rotate.append(l)
        else:
            v_mode, e_mode, f_mode = tool.mesh_select_mode
            if f_mode:
                for f in bm.faces:
                    if not f.select:
                        continue
                    for l in f.loops:
                        loops_to_rotate.append(l)
            else:
                for f in bm.faces:
                    for l in f.loops:
                        if l.vert.select:
                            loops_to_rotate.append(l)

    if not loops_to_rotate:
        return

    sx = sy = 0.0
    n = 0
    for l in loops_to_rotate:
        uv = l[uv_layer].uv
        sx += uv.x; sy += uv.y
        n += 1
    if n == 0:
        return
    cx = sx / n
    cy = sy / n

    ang = angle_deg % 360.0
    if ang in (90.0, -270.0):
        def rot(u, v): return -v, u
    elif ang in (270.0, -90.0):
        def rot(u, v): return v, -u
    elif ang in (180.0, -180.0):
        def rot(u, v): return -u, -v
    else:
        rad = math.radians(ang)
        c, s = math.cos(rad), math.sin(rad)
        def rot(u, v): return u * c - v * s, u * s + v * c

    for l in loops_to_rotate:
        uv = l[uv_layer].uv
        u, v = uv.x - cx, uv.y - cy
        ru, rv = rot(u, v)
        uv.x = ru + cx
        uv.y = rv + cy

def _flip_uvs_group_in_bmesh(bm, axis='H', only_selected=True):
    """Flip UVs as one group (keeps layout). axis: 'H' (horizontal) or 'V' (vertical)."""
    uv_layer = bm.loops.layers.uv.active
    if not uv_layer:
        return

    loops_to_flip = []
    tool = bpy.context.tool_settings
    uv_sync = getattr(tool, "use_uv_select_sync", False)

    if not only_selected:
        for f in bm.faces:
            for l in f.loops:
                loops_to_flip.append(l)
    else:
        if not uv_sync:
            for f in bm.faces:
                for l in f.loops:
                    if l[uv_layer].select:
                        loops_to_flip.append(l)
        else:
            v_mode, e_mode, f_mode = tool.mesh_select_mode
            if f_mode:
                for f in bm.faces:
                    if f.select:
                        for l in f.loops:
                            loops_to_flip.append(l)
            else:
                for f in bm.faces:
                    for l in f.loops:
                        if l.vert.select:
                            loops_to_flip.append(l)

    if not loops_to_flip:
        return

    sx = sy = 0.0
    for l in loops_to_flip:
        uv = l[uv_layer].uv
        sx += uv.x; sy += uv.y
    cx = sx / len(loops_to_flip)
    cy = sy / len(loops_to_flip)

    for l in loops_to_flip:
        uv = l[uv_layer].uv
        if axis == 'H':
            uv.x = 2*cx - uv.x
        elif axis == 'V':
            uv.y = 2*cy - uv.y
        
# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

class TDSettings(PropertyGroup):
    density: FloatProperty(
        name="Texel Density",
        description="Texel density (px/cm)",
        default=0.0,
        min=0.0,
        precision=3,
        step=1
    )
    move_amount: FloatProperty(
        name="Amount",
        description="Distance in UDIM tiles (1.0 = one tile)",
        default=1.0,
        min=0.0,
        soft_max=10.0
    )
    move_selected_only: BoolProperty(
        name="Contextual Toggle",
        description="Edit Mode: only selected UVs. Object Mode: only active mesh.",
        default=True
    )

# ---------------------------------------------------------------------------
# Operators — TD
# ---------------------------------------------------------------------------

class BB_Texel_Density_Check(Operator):
    bl_idname = "bb_uvs.texel_density_check"
    bl_label = "Get TD"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        td = context.scene.bb_td
        active_obj = context.active_object
        if not _object_is_uv_mesh(active_obj):
            self.report({'WARNING'}, "No active mesh with UVs")
            return {'CANCELLED'}

        lst = Calculate_TD_Area_To_List()
        if not lst:
            td.density = 0.0
            return {'CANCELLED'}

        total_area = sum(a[1] for a in lst)
        if total_area == 0:
            td.density = 0.0
            return {'CANCELLED'}

        weighted_td = sum(val * (area / total_area) for val, area in lst)
        td.density = float(weighted_td)
        return {'FINISHED'}

class BB_Texel_Density_Set(Operator):
    bl_idname = "bb_uvs.texel_density_set"
    bl_label = "Set TD"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        td = context.scene.bb_td
        target_td = float(td.density)
        if target_td <= 0:
            self.report({'WARNING'}, "No valid TD. Run Get TD first.")
            return {'CANCELLED'}

        start_active = context.view_layer.objects.active
        start_mode = getattr(start_active, "mode", "OBJECT") if start_active else "OBJECT"
        start_selected = context.selected_objects[:]

        if start_mode == 'EDIT':
            targets = [o for o in context.objects_in_mode if _object_is_uv_mesh(o)]
        else:
            # Object Mode: always apply to all selected meshes with UVs
            targets = [o for o in start_selected if _object_is_uv_mesh(o)]

        if not targets:
            self.report({'WARNING'}, "No mesh with UVs selected")
            return {'CANCELLED'}

        for o in targets:
            prev = context.view_layer.objects.active
            context.view_layer.objects.active = o
            lst = Calculate_TD_Area_To_List()
            context.view_layer.objects.active = prev

            total_area = sum(a[1] for a in lst)
            if total_area <= 0:
                continue
            current_td = sum(val * (area / total_area) for val, area in lst)
            if current_td <= 0:
                current_td = 0.0001

            scale = target_td / current_td

            if o.data.is_editmode:
                bm = bmesh.from_edit_mesh(o.data)
                uv_layer = bm.loops.layers.uv.active
                if uv_layer:
                    bm.faces.ensure_lookup_table()
                    for f in bm.faces:
                        for loop in f.loops:
                            uv = loop[uv_layer].uv
                            uv.x *= scale
                            uv.y *= scale
                    bmesh.update_edit_mesh(o.data, loop_triangles=False, destructive=False)
            else:
                bm = bmesh.new()
                bm.from_mesh(o.data)
                uv_layer = bm.loops.layers.uv.active
                if uv_layer:
                    bm.faces.ensure_lookup_table()
                    for f in bm.faces:
                        for loop in f.loops:
                            uv = loop[uv_layer].uv
                            uv.x *= scale
                            uv.y *= scale
                    bm.to_mesh(o.data)
                    o.data.update()
                bm.free()

        return {'FINISHED'}

# ---------------------------------------------------------------------------
# Operators — Move / Rotate / Flip
# ---------------------------------------------------------------------------

class BB_UVs_SetMoveContext(Operator):
    bl_idname = "bb_uvs.set_move_context"
    bl_label = "Set Move Context"
    bl_options = {'INTERNAL', 'UNDO'}

    set_selected: BoolProperty(name="Set Selected Mode", default=True)

    def execute(self, context):
        td = context.scene.bb_td
        td.move_selected_only = bool(self.set_selected)
        return {'FINISHED'}

class BB_UVs_ToggleMoveHighlight(Operator):
    bl_idname = "bb_uvs.toggle_move_highlight"
    bl_label = "Toggle Highlighted Only"
    bl_options = {'INTERNAL', 'UNDO'}

    def execute(self, context):
        td = context.scene.bb_td
        td.move_selected_only = not td.move_selected_only
        return {'FINISHED'}

class BB_UVs_MoveUVs(Operator):
    bl_idname = "bb_uvs.move_uvs"
    bl_label = "Move UVs"
    bl_options = {'REGISTER', 'UNDO'}

    direction: EnumProperty(
        items=[
            ('UP', "Up", ""),
            ('DOWN', "Down", ""),
            ('LEFT', "Left", ""),
            ('RIGHT', "Right", ""),
            ('UP_LEFT', "Up Left", ""),
            ('UP_RIGHT', "Up Right", ""),
            ('DOWN_LEFT', "Down Left", ""),
            ('DOWN_RIGHT', "Down Right", ""),
        ],
        name="Direction",
        default='UP'
    )

    def execute(self, context):
        td = context.scene.bb_td
        amt = float(td.move_amount)
        if amt == 0.0:
            return {'CANCELLED'}

        dx = dy = 0.0
        if self.direction == 'UP':
            dy = amt
        elif self.direction == 'DOWN':
            dy = -amt
        elif self.direction == 'LEFT':
            dx = -amt
        elif self.direction == 'RIGHT':
            dx = amt
        elif self.direction == 'UP_LEFT':
            dx = -amt; dy = amt
        elif self.direction == 'UP_RIGHT':
            dx = amt; dy = amt
        elif self.direction == 'DOWN_LEFT':
            dx = -amt; dy = -amt
        elif self.direction == 'DOWN_RIGHT':
            dx = amt; dy = -amt

        active = context.view_layer.objects.active
        mode = getattr(active, "mode", "OBJECT") if active else "OBJECT"

        if mode == 'EDIT':
            if td.move_selected_only:
                targets = [o for o in context.objects_in_mode if _object_is_uv_mesh(o)]
                for o in targets:
                    _move_uvs_edit(o, dx, dy, only_selected=True)
            else:
                if not _object_is_uv_mesh(active):
                    self.report({'WARNING'}, "Active (highlighted) object is not a mesh with UVs")
                    return {'CANCELLED'}
                _move_uvs_edit(active, dx, dy, only_selected=False)
        else:
            selected = context.selected_objects[:]
            if td.move_selected_only:
                if not _object_is_uv_mesh(active):
                    self.report({'WARNING'}, "Active (highlighted) object is not a mesh with UVs")
                    return {'CANCELLED'}
                targets = [active]
            else:
                targets = [o for o in selected if _object_is_uv_mesh(o)]

            if not targets:
                self.report({'WARNING'}, "No mesh with UVs to move")
                return {'CANCELLED'}

            for o in targets:
                _move_uvs_object(o, dx, dy)

        try:
            bpy.ops.bb_uvs.texel_density_check()
        except Exception:
            pass

        return {'FINISHED'}

class BB_UVs_RotateUVs(Operator):
    """Rotate UVs 90° CCW/CW as a single group (respects selection in Edit Mode)."""
    bl_idname = "bb_uvs.rotate_uvs"
    bl_label = "Rotate UVs"
    bl_options = {'REGISTER', 'UNDO'}

    direction: EnumProperty(
        items=[('CCW', "CCW", "Rotate 90° counter-clockwise"),
               ('CW',  "CW",  "Rotate 90° clockwise")],
        name="Direction",
        default='CCW'
    )

    def execute(self, context):
        angle = 90.0 if self.direction == 'CCW' else -90.0
        active = context.view_layer.objects.active
        if not _object_is_uv_mesh(active):
            self.report({'WARNING'}, "No active mesh with UVs")
            return {'CANCELLED'}

        mode = getattr(active, "mode", "OBJECT")

        if mode == 'EDIT':
            targets = [o for o in context.objects_in_mode if _object_is_uv_mesh(o)]
            for o in targets:
                bm = bmesh.from_edit_mesh(o.data)
                _rotate_uvs_group_in_bmesh(bm, angle_deg=angle, only_selected=True)
                bmesh.update_edit_mesh(o.data, loop_triangles=False, destructive=False)
        else:
            bm = bmesh.new()
            bm.from_mesh(active.data)
            _rotate_uvs_group_in_bmesh(bm, angle_deg=angle, only_selected=False)
            bm.to_mesh(active.data)
            active.data.update()
            bm.free()

        try:
            bpy.ops.bb_uvs.texel_density_check()
        except Exception:
            pass
        return {'FINISHED'}

class BB_UVs_FlipUVs(Operator):
    """Flip UVs as a single group (Horizontal or Vertical)."""
    bl_idname = "bb_uvs.flip_uvs"
    bl_label = "Flip UVs"
    bl_options = {'REGISTER', 'UNDO'}

    axis: EnumProperty(
        items=[('H', "Horizontal", "Flip left-right"),
               ('V', "Vertical",   "Flip top-bottom")],
        name="Axis",
        default='H'
    )

    def execute(self, context):
        active = context.view_layer.objects.active
        if not _object_is_uv_mesh(active):
            self.report({'WARNING'}, "No active mesh with UVs")
            return {'CANCELLED'}

        mode = getattr(active, "mode", "OBJECT")
        if mode == 'EDIT':
            targets = [o for o in context.objects_in_mode if _object_is_uv_mesh(o)]
            for o in targets:
                bm = bmesh.from_edit_mesh(o.data)
                _flip_uvs_group_in_bmesh(bm, axis=self.axis, only_selected=True)
                bmesh.update_edit_mesh(o.data, loop_triangles=False, destructive=False)
        else:
            bm = bmesh.new()
            bm.from_mesh(active.data)
            _flip_uvs_group_in_bmesh(bm, axis=self.axis, only_selected=False)
            bm.to_mesh(active.data)
            active.data.update()
            bm.free()

        return {'FINISHED'}

# ---------------------------------------------------------------------------
# Operators — Normalize + Pack (normalize ALL objects first; then one pack)
# ---------------------------------------------------------------------------

class BB_UVs_NormalizePack(Operator):
    """Normalize islands across *all* edit-mode objects, then pack once (snap to first UDIM)."""
    bl_idname = "bb_uvs.normalize_pack"
    bl_label = "Pack"
    bl_options = {'REGISTER', 'UNDO'}

    @staticmethod
    def _call_pack(args):
        try:
            return bpy.ops.uv.pack_islands('EXEC_DEFAULT', **args)
        except TypeError:
            pass
        args2 = dict(args); args2.pop("shape_method", None)
        try:
            return bpy.ops.uv.pack_islands('EXEC_DEFAULT', **args2)
        except TypeError:
            pass
        args3 = dict(args2); args3.pop("rotate", None)
        try:
            return bpy.ops.uv.pack_islands('EXEC_DEFAULT', **args3)
        except TypeError:
            pass
        args4 = {
            "udim_source": args.get("udim_source", "ACTIVE_UDIM"),
            "margin_method": args.get("margin_method", "FRACTION"),
            "margin": args.get("margin", 0.00390625),
        }
        return bpy.ops.uv.pack_islands('EXEC_DEFAULT', **args4)

    @staticmethod
    def _normalize_all_edit_meshes(only_selected=True):
        for obj in bpy.context.objects_in_mode:
            if obj.type != 'MESH':
                continue
            bm = bmesh.from_edit_mesh(obj.data)
            _normalize_islands_in_bmesh(bm, xy_independent=True, only_selected=only_selected, respect_aspect=True)
            bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)

    @staticmethod
    def _snap_uvs_to_first_udim_all():
        """Translate selected UVs so their minimum UDIM lands at (0,0) across *all* edit-mode meshes."""
        objs = [o for o in bpy.context.objects_in_mode if o.type == 'MESH']
        if not objs:
            return

        # Compute global min tile across selection
        min_u = float('inf')
        min_v = float('inf')
        for obj in objs:
            bm = bmesh.from_edit_mesh(obj.data)
            uv_layer = bm.loops.layers.uv.active
            if not uv_layer:
                continue
            for f in bm.faces:
                if not f.select:
                    continue
                for l in f.loops:
                    uv = l[uv_layer].uv
                    if uv.x < min_u: min_u = uv.x
                    if uv.y < min_v: min_v = uv.y
        if min_u == float('inf'):
            return

        du = -math.floor(min_u)
        dv = -math.floor(min_v)
        if du == 0 and dv == 0:
            return

        for obj in objs:
            bm = bmesh.from_edit_mesh(obj.data)
            uv_layer = bm.loops.layers.uv.active
            if not uv_layer:
                continue
            for f in bm.faces:
                if not f.select:
                    continue
                for l in f.loops:
                    uv = l[uv_layer].uv
                    uv.x += du
                    uv.y += dv
            bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)

    def execute(self, context):
        args = {
            "udim_source": "ACTIVE_UDIM",
            "rotate": False,
            "scale": True,
            "margin_method": "FRACTION",
            "margin": 0.00390625,
            "pin": False,
            "merge_overlap": False,
            "pin_method": "LOCKED",
            "shape_method": "CONCAVE",
        }

        active = context.view_layer.objects.active
        if not _object_is_uv_mesh(active):
            self.report({'WARNING'}, "No active mesh with UVs")
            return {'CANCELLED'}

        start_mode = getattr(active, "mode", "OBJECT")

        if start_mode == 'EDIT':
            ts = context.tool_settings
            ts.mesh_select_mode = (False, False, True)
            bpy.ops.mesh.select_all(action='SELECT')
            try:
                bpy.ops.uv.select_all(action='SELECT')
            except Exception:
                pass

            self._normalize_all_edit_meshes(only_selected=True)
            self._call_pack(args)
            self._snap_uvs_to_first_udim_all()
            return {'FINISHED'}

        # OBJECT mode: normalize ALL selected meshes with UVs, then one pack using multi-edit
        prev_active = active
        prev_selection = context.selected_objects[:]

        targets = [o for o in prev_selection if _object_is_uv_mesh(o)]
        if not targets:
            self.report({'WARNING'}, "No selected meshes with UVs to pack")
            return {'CANCELLED'}

        # Enter multi-object edit on targets
        for o in bpy.data.objects:
            o.select_set(False)
        for o in targets:
            o.select_set(True)
        context.view_layer.objects.active = targets[0]

        bpy.ops.object.mode_set(mode='EDIT')

        ts = context.tool_settings
        ts.mesh_select_mode = (False, False, True)
        bpy.ops.mesh.select_all(action='SELECT')
        try:
            bpy.ops.uv.select_all(action='SELECT')
        except Exception:
            pass

        self._normalize_all_edit_meshes(only_selected=False)
        self._call_pack(args)
        self._snap_uvs_to_first_udim_all()

        bpy.ops.object.mode_set(mode='OBJECT')

        # Restore selection and active
        for o in bpy.data.objects:
            o.select_set(False)
        for o in prev_selection:
            if o and o.name in bpy.data.objects:
                o.select_set(True)
        context.view_layer.objects.active = prev_active

        return {'FINISHED'}

# ---------------------------------------------------------------------------
# Operator — UI prefs helper
# ---------------------------------------------------------------------------

class BB_UVs_SetUI(Operator):
    bl_idname = "bb_uvs.set_ui"
    bl_label = "Set UI"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        image_space = None
        for area in context.screen.areas:
            if area.type == 'IMAGE_EDITOR':
                image_space = area.spaces.active
                break
        if image_space is None:
            self.report({'WARNING'}, "No UV/Image Editor found")
            return {'CANCELLED'}

        image_space.uv_editor.edge_display_type = 'BLACK'
        image_space.uv_editor.tile_grid_shape = (10, 10)
        image_space.overlay.show_grid_background = True

        ts = context.tool_settings
        try:
            ts.mesh_select_mode = (False, True, False)
        except Exception:
            pass

        try:
            ts.uv_select_mode = 'EDGE'
        except Exception:
            pass

        self.report({'INFO'}, "UI set: UV grid on, and selection set to Edge (mesh & UV).")
        return {'FINISHED'}

# ---------------------------------------------------------------------------
# UI Panels
# ---------------------------------------------------------------------------

def draw_bb_uvs_panel(self, context):
    td = context.scene.bb_td
    layout = self.layout

    layout.operator("bb_uvs.set_ui", text="Set UI")

    layout.separator()
    layout.label(text="Pack")
    layout.operator("bb_uvs.normalize_pack", text="Pack", icon='UV')

    layout.separator()
    row = layout.row(align=True)
    row.prop(td, "density", text="Texel Density")
    row = layout.row(align=True)
    row.operator("bb_uvs.texel_density_check", text="Get TD")
    row.operator("bb_uvs.texel_density_set", text="Set TD")

    layout.separator()
    layout.label(text="Move")

    active = context.view_layer.objects.active
    mode = getattr(active, "mode", "OBJECT") if active else "OBJECT"

    row = layout.row(align=True)
    row.label(text="Context:")

    if mode == 'EDIT':
        left = row.row(align=True)
        op = left.operator("bb_uvs.set_move_context",
                           text="Selected",
                           icon='RESTRICT_SELECT_OFF',
                           depress=td.move_selected_only)
        op.set_selected = True

        right = row.row(align=True)
        op = right.operator("bb_uvs.set_move_context",
                            text="Highlighted",
                            icon='RADIOBUT_ON',
                            depress=not td.move_selected_only)
        op.set_selected = False
    else:
        row.operator("bb_uvs.toggle_move_highlight",
                     text="Highlighted Only",
                     icon='RESTRICT_SELECT_OFF',
                     depress=td.move_selected_only)

    row = layout.row(align=True)
    col = row.column(align=True)
    col.operator_context = 'EXEC_DEFAULT'
    op = col.operator("bb_uvs.move_uvs", text="↖"); op.direction = 'UP_LEFT'
    op = col.operator("bb_uvs.move_uvs", text="←");  op.direction = 'LEFT'
    op = col.operator("bb_uvs.move_uvs", text="↙");  op.direction = 'DOWN_LEFT'

    col = row.column(align=True)
    op = col.operator("bb_uvs.move_uvs", text="↑");  op.direction = 'UP'
    col.prop(td, "move_amount", text="")
    op = col.operator("bb_uvs.move_uvs", text="↓");  op.direction = 'DOWN'

    col = row.column(align=True)
    op = col.operator("bb_uvs.move_uvs", text="↗");  op.direction = 'UP_RIGHT'
    op = col.operator("bb_uvs.move_uvs", text="→");  op.direction = 'RIGHT'
    op = col.operator("bb_uvs.move_uvs", text="↘");  op.direction = 'DOWN_RIGHT'

    layout.separator()
    layout.label(text="Rotate")
    r = layout.row()
    r.scale_y = 1.1
    op = r.operator("bb_uvs.rotate_uvs", text="90° CCW", icon='LOOP_BACK');    op.direction = 'CCW'
    op = r.operator("bb_uvs.rotate_uvs", text="90° CW",  icon='LOOP_FORWARDS'); op.direction = 'CW'

    layout.separator()
    layout.label(text="Flip")
    f = layout.row()
    f.scale_y = 1.1
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

# ---------------------------------------------------------------------------
# Registration
# ---------------------------------------------------------------------------

classes = (
    TDSettings,
    BB_Texel_Density_Check,
    BB_Texel_Density_Set,
    BB_UVs_MoveUVs,
    BB_UVs_RotateUVs,
    BB_UVs_NormalizePack,
    BB_UVs_SetUI,
    BB_UVs_FlipUVs,
    BB_UVs_SetMoveContext,
    BB_UVs_ToggleMoveHighlight,
    UV_PT_BB_UVs,
    VIEW3D_PT_BB_UVs,
)

def register():
    for cls in classes:
        bpy.utils.register_class(cls)
    bpy.types.Scene.bb_td = PointerProperty(type=TDSettings)

def unregister():
    del bpy.types.Scene.bb_td
    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)

if __name__ == "__main__":
    register()
