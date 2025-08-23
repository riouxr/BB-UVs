# SPDX-License-Identifier: GPL-3.0-or-later
# -----------------------------------------------------------------------------
# BB UVs (Minimal + Move UVs — Selected/Highlighted Defaults)
# Derived work based on: Texel Density Checker — https://github.com/mrven/Blender-Texel-Density-Checker
# Original author: Ivan Vostrikov (mrven)
# 
# This file is part of a Blender add-on that includes code and ideas derived
# from the project “Texel Density Checker” by Ivan Vostrikov (mrven), which is
# licensed under the GNU General Public License version 3.0 or (at your option)
# any later version.
#
# Copyright (C) 2025 BB (modifications and adaptations)
# Copyright (C) Ivan Vostrikov (mrven) and contributors (original project)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# Summary of modifications from the upstream project:
# - Minimal “Get/Set Texel Density” implementation tailored for fast workflows.
# - Independent per-object TD calculation without changing selections.
# - Added Move UVs operator with contextual “Selected Only / Highlighted Only” toggle.
# - Compact UI panels for UV Editor and 3D View.
# -----------------------------------------------------------------------------

bl_info = {
    "name": "BB UVs (Minimal + Move UVs — Selected/Highlighted Defaults)",
    "author": "BB",
    "version": (1, 1, 4),
    "blender": (4, 0, 0),
    "location": "UV Editor > Sidebar (N) > BB UVs and 3D View > Sidebar (N) > BB UVs",
    "description": "Minimal Get/Set Texel Density + Move UVs with Selected/Highlighted toggle",
    "category": "UV",
}

import bpy
import bmesh
import math
from bpy.props import PointerProperty, StringProperty, BoolProperty, FloatProperty, EnumProperty
from bpy.types import PropertyGroup, Panel, Operator

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def _object_is_uv_mesh(o):
    return (o and o.type == 'MESH' and o.data and len(o.data.uv_layers) > 0 and len(o.data.polygons) > 0)

def Calculate_TD_Area_To_List():
    """Return [texel_density, uv_area] for each polygon of ACTIVE object, without touching selections."""
    active_obj = bpy.context.active_object
    if not _object_is_uv_mesh(active_obj):
        return []

    # Fixed square texture for TD math (px/cm)
    w = h = 1024
    aspect_ratio = w / h
    if aspect_ratio < 1:
        aspect_ratio = 1 / aspect_ratio
    largest_side = max(w, h)

    # Work on a temporary mesh copy to avoid mode flips
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

        # Shoelace in UV space
        loops = [loop[uv_layer].uv.copy() for loop in f.loops]
        area_uv = 0.0
        a = len(loops) - 1
        for b in range(len(loops)):
            area_uv += (loops[a].x + loops[b].x) * (loops[a].y - loops[b].y)
            a = b
        area_uv = abs(0.5 * area_uv)

        area_geo = f.calc_area()  # geometric face area

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
    """
    Move UVs of an edit-mode mesh. If only_selected=True, move only selected UV loops.
    Works with both UV Sync OFF (uses UV loop selection flags) and UV Sync ON
    (uses mesh face/vertex selection based on mesh_select_mode).
    """
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

    # only_selected path:
    tool = bpy.context.tool_settings
    uv_sync = tool.use_uv_select_sync

    if not uv_sync:
        # UV Sync OFF: trust UV loop selection flags
        for f in bm.faces:
            for loop in f.loops:
                if not loop[uv_layer].select:
                    continue
                uv = loop[uv_layer].uv
                uv.x += dx
                uv.y += dy
        bmesh.update_edit_mesh(o.data, loop_triangles=False, destructive=False)
        return

    # UV Sync ON: use mesh selection (vertex/edge/face). mesh_select_mode is a 3-tuple (v, e, f).
    v_mode, e_mode, f_mode = tool.mesh_select_mode
    if f_mode:
        # Face select: move loops of selected faces
        for f in bm.faces:
            if not f.select:
                continue
            for loop in f.loops:
                uv = loop[uv_layer].uv
                uv.x += dx
                uv.y += dy
    else:
        # Vertex/edge select: move loops whose verts are selected
        for f in bm.faces:
            for loop in f.loops:
                if not loop.vert.select:
                    continue
                uv = loop[uv_layer].uv
                uv.x += dx
                uv.y += dy

    bmesh.update_edit_mesh(o.data, loop_triangles=False, destructive=False)

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

class TDSettings(PropertyGroup):
    density: StringProperty(
        name="Current TD",
        description="Texel density (px/cm)",
        default="0"
    )
    move_amount: FloatProperty(
        name="Amount",
        description="Distance in UDIM tiles (1.0 = one tile)",
        default=1.0,
        min=0.0,
        soft_max=10.0
    )
    # One toggle used contextually:
    # - Edit Mode: "Selected Only" (move only selected UVs)
    # - Object Mode: "Highlighted Only" (move only the active mesh)
    move_selected_only: BoolProperty(
        name="Contextual Toggle",
        description="Edit Mode: only selected UVs. Object Mode: only active mesh.",
        default=True  # ON by default as requested
    )

# ---------------------------------------------------------------------------
# Operators (TD)
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
            td.density = "0"
            return {'CANCELLED'}

        total_area = sum(a[1] for a in lst)
        if total_area == 0:
            td.density = "0"
            return {'CANCELLED'}

        weighted_td = sum(val * (area / total_area) for val, area in lst)
        td.density = f"{weighted_td:.3f}"
        return {'FINISHED'}

class BB_Texel_Density_Set(Operator):
    bl_idname = "bb_uvs.texel_density_set"
    bl_label = "Set TD"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        td = context.scene.bb_td
        try:
            target_td = float(td.density)
        except:
            target_td = 0.0
        if target_td <= 0:
            self.report({'WARNING'}, "No valid TD. Run Get TD first.")
            return {'CANCELLED'}

        start_active = context.view_layer.objects.active
        start_mode = getattr(start_active, "mode", "OBJECT") if start_active else "OBJECT"
        start_selected = context.selected_objects[:]

        # Determine targets
        if start_mode == 'EDIT':
            targets = [o for o in context.objects_in_mode if _object_is_uv_mesh(o)]
        else:
            targets = [start_active] if td.move_selected_only else [o for o in start_selected if _object_is_uv_mesh(o)]

        if not targets:
            self.report({'WARNING'}, "No mesh with UVs selected")
            return {'CANCELLED'}

        # Scale UVs of each target to match target TD (around origin)
        for o in targets:
            # Compute current TD per object (without touching selections)
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

        # Refresh TD readout (doesn't change selections)
        try:
            bpy.ops.bb_uvs.texel_density_check()
        except:
            pass
        return {'FINISHED'}

# ---------------------------------------------------------------------------
# Operator (Move UVs) — Selected Only (Edit) / Highlighted Only (Object)
# ---------------------------------------------------------------------------

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
            dx = -amt
            dy = amt
        elif self.direction == 'UP_RIGHT':
            dx = amt
            dy = amt
        elif self.direction == 'DOWN_LEFT':
            dx = -amt
            dy = -amt
        elif self.direction == 'DOWN_RIGHT':
            dx = amt
            dy = -amt

        active = context.view_layer.objects.active
        mode = getattr(active, "mode", "OBJECT") if active else "OBJECT"
        selected = context.selected_objects[:]

        if mode == 'EDIT':
            # Edit Mode: operate on all meshes in multi-edit
            targets = [o for o in context.objects_in_mode if _object_is_uv_mesh(o)]
            only_selected_uvs = td.move_selected_only  # "Selected Only" in Edit Mode
            for o in targets:
                _move_uvs_edit(o, dx, dy, only_selected=only_selected_uvs)
        else:
            # Object Mode: "Highlighted Only" => active mesh only (default ON)
            if td.move_selected_only:
                targets = [active] if _object_is_uv_mesh(active) else []
            else:
                targets = [o for o in selected if _object_is_uv_mesh(o)]
            if not targets:
                self.report({'WARNING'}, "No mesh with UVs selected")
                return {'CANCELLED'}
            for o in targets:
                _move_uvs_object(o, dx, dy)

        # Keep TD readout fresh (no selection changes)
        try:
            bpy.ops.bb_uvs.texel_density_check()
        except:
            pass
        return {'FINISHED'}

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
        image_space.overlay.show_grid_background = True  # Enable to show the grid
        return {'FINISHED'}

# ---------------------------------------------------------------------------
# UI Panels
# ---------------------------------------------------------------------------

def draw_bb_uvs_panel(self, context):
    td = context.scene.bb_td
    layout = self.layout

    layout.operator("bb_uvs.set_ui", text="Set UI")

    # Minimal TD block
    layout.label(text=f"Texel Density: {td.density}")
    layout.operator("bb_uvs.texel_density_check", text="Get TD")
    layout.operator("bb_uvs.texel_density_set", text="Set TD")

    # Move UVs section
    layout.separator()

    layout.label(text="Move")

    # Contextual toggle
    active = context.view_layer.objects.active
    mode = getattr(active, "mode", "OBJECT") if active else "OBJECT"
    layout.prop(td, "move_selected_only",
                text=("Selected Only" if mode == 'EDIT' else "Highlighted Only"))

    row = layout.row(align=True)
    # First column
    col = row.column(align=True)
    col.operator_context = 'EXEC_DEFAULT'
    op = col.operator("bb_uvs.move_uvs", text="↖")
    op.direction = 'UP_LEFT'
    op = col.operator("bb_uvs.move_uvs", text="←")
    op.direction = 'LEFT'
    op = col.operator("bb_uvs.move_uvs", text="↙")
    op.direction = 'DOWN_LEFT'
    # Second column
    col = row.column(align=True)
    op = col.operator("bb_uvs.move_uvs", text="↑")
    op.direction = 'UP'
    col.prop(td, "move_amount", text="")
    op = col.operator("bb_uvs.move_uvs", text="↓")
    op.direction = 'DOWN'
    # Third column
    col = row.column(align=True)
    op = col.operator("bb_uvs.move_uvs", text="↗")
    op.direction = 'UP_RIGHT'
    op = col.operator("bb_uvs.move_uvs", text="→")
    op.direction = 'RIGHT'
    op = col.operator("bb_uvs.move_uvs", text="↘")
    op.direction = 'DOWN_RIGHT'
    
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
    BB_UVs_SetUI,
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