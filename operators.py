import bpy
import bmesh
import math
from math import isclose
from .helpers import *

# ---------- TD / Move / Rotate / Flip Operators ----------

class BB_Texel_Density_Check(bpy.types.Operator):
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

class BB_Texel_Density_Set(bpy.types.Operator):
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

class BB_UVs_SetMoveCollection(bpy.types.Operator):
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
                    if col_objs:
                        bpy.ops.object.select_all(action='DESELECT')
                        for o in col_objs:
                            o.select_set(True)
                        context.view_layer.objects.active = active
        return {'FINISHED'}

class BB_UVs_SetMoveContext(bpy.types.Operator):
    bl_idname = "bb_uvs.set_move_context"
    bl_label = "Set Move Context"
    bl_options = {'INTERNAL', 'UNDO'}
    set_selected: bpy.props.BoolProperty(name="Set Selected Mode", default=True)
    def execute(self, context):
        S = _S(context)
        S.bb_move_selected_only = bool(self.set_selected)
        S.bb_move_collection = False
        return {'FINISHED'}

class BB_UVs_ToggleMoveHighlight(bpy.types.Operator):
    bl_idname = "bb_uvs.toggle_move_highlight"
    bl_label = "Toggle Highlighted Only"
    bl_options = {'INTERNAL', 'UNDO'}
    def execute(self, context):
        S = _S(context); S.bb_move_selected_only = not S.bb_move_selected_only; return {'FINISHED'}

class BB_UVs_MoveUVs(bpy.types.Operator):
    bl_idname = "bb_uvs.move_uvs"
    bl_label = "Move UVs"
    bl_options = {'REGISTER', 'UNDO'}
    direction: bpy.props.EnumProperty(
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

        if context.mode == 'EDIT_MESH':
            S.bb_move_collection = False
            if S.bb_move_selected_only:  # Highlighted: move all UVs of active object
                if _object_is_uv_mesh(active):
                    _move_uvs_edit(active, dx, dy, only_selected=False)
            else:  # Selected: move selected UVs of all objects in edit mode
                targets = [o for o in context.objects_in_mode if _object_is_uv_mesh(o)]
                for o in targets:
                    _move_uvs_edit(o, dx, dy, only_selected=True)
        else:  # Object mode
            targets = []
            if S.bb_move_selected_only:
                if _object_is_uv_mesh(active):
                    targets = [active]
            elif S.bb_move_collection:
                if active and active.users_collection:
                    active_coll = active.users_collection[0]
                    targets = [o for o in active_coll.objects if _object_is_uv_mesh(o)]
            else:
                targets = [o for o in context.selected_objects if _object_is_uv_mesh(o)]
            for o in targets:
                _move_uvs_object(o, dx, dy)

        context.view_layer.update()
        return {'FINISHED'}
    
class BB_UVs_SelectSimilarTopology(bpy.types.Operator):
    bl_idname = "bb_uvs.select_similar_topology"
    bl_label = "Sim"
    bl_options = {'REGISTER', 'UNDO'}
    include_scale: bpy.props.BoolProperty(default=False)

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
    
class BB_UVs_SelectSimilarTopologyScale(bpy.types.Operator):
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

            # --- New zoom-aware sensitivity ---
            region = context.region
            v2d = region.view2d if hasattr(region, "view2d") else None
            if v2d:
                x1, y1 = v2d.region_to_view(0, 0)
                x2, y2 = v2d.region_to_view(100, 0)
                uv_per_pixel = abs(x2 - x1) / 100.0
            else:
                uv_per_pixel = 0.005  # fallback

            dx = delta_x * uv_per_pixel
            dy = delta_y * uv_per_pixel
            # -----------------------------------

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

class BB_UVs_RotateUVs(bpy.types.Operator):
    bl_idname = "bb_uvs.rotate_uvs"
    bl_label = "Rotate UVs"
    bl_options = {'REGISTER', 'UNDO'}
    direction: bpy.props.EnumProperty(
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

class BB_UVs_FlipUVs(bpy.types.Operator):
    bl_idname = "bb_uvs.flip_uvs"
    bl_label = "Flip UVs"
    bl_options = {'REGISTER', 'UNDO'}
    axis: bpy.props.EnumProperty(
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

class BB_UVs_NormalizePack(bpy.types.Operator):
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

class BB_UVs_PackIndividually(bpy.types.Operator):
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

class BB_UVs_SelectFlipped(bpy.types.Operator):
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
                            loop[uv_layer].select_edge = True
                    local_counter += 1
            
            total_counter += local_counter
            bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)
        
        if not total_counter:
            self.report({'INFO'}, 'Flipped faces not found')
        else:
            self.report({'INFO'}, f'Found {total_counter} flipped faces')
            
        return {'FINISHED'}

class BB_UVs_SelectCrossUDIMs(bpy.types.Operator):
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

            # Pre-collect all islands first (much faster)
            islands = _collect_uv_islands(bm, uv_layer, only_selected=False)
            
            local_counter = 0
            processed_islands = set()
            
            for island in islands:
                if id(island) in processed_islands:
                    continue
                    
                # Check if any face in this island crosses UDIMs
                island_crosses_udim = False
                for f in island:
                    tiles = set()
                    for loop in f.loops:
                        uv = loop[uv_layer].uv
                        tile_u = math.floor(uv.x)
                        tile_v = math.floor(uv.y)
                        tiles.add((tile_u, tile_v))
                        if len(tiles) > 1:
                            island_crosses_udim = True
                            break
                    if island_crosses_udim:
                        break

                # If island crosses UDIMs, select all faces in it
                if island_crosses_udim:
                    for f in island:
                        if sync:
                            f.select = True
                        else:
                            for loop in f.loops:
                                loop[uv_layer].select = True
                                loop[uv_layer].select_edge = True
                        local_counter += 1
                    
                    processed_islands.add(id(island))

            total_counter += local_counter
            bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)

        if total_counter == 0:
            self.report({'INFO'}, "No faces crossing UDIM boundaries found")
        else:
            self.report({'INFO'}, f"Selected {total_counter} face(s) crossing UDIMs")

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
            bpy.ops.mesh.select_all(action='SELECT')
            context.tool_settings.use_uv_select_sync = False
            # Ensure UVs are fully selected after disabling sync
            bpy.ops.uv.select_all(action='SELECT')
            
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
                            # Select the edge
                            face.loops[i][uv_layer].select_edge = True
                            
                            # Also select the two UV vertices for better visibility
                            face.loops[i][uv_layer].select = True
                            next_i = (i + 1) % len(face.loops)
                            face.loops[next_i][uv_layer].select = True
                            
                            total_selected += 1

            bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)

        # Force UV select mode to EDGE for better visibility of selected edges
        context.tool_settings.uv_select_mode = 'EDGE'

        if total_selected == 0:
            self.report({'INFO'}, "No internal boundaries selected")
        else:
            self.report({'INFO'}, f"Selected {total_selected} internal boundary edges")
            
            
        return {'FINISHED'}

# ------------------------------ Grid Helper Operators ------------------------------

class BB_UVs_ChooseGridImage(bpy.types.Operator):
    bl_idname = "bb_uvs.choose_grid_image"
    bl_label = "Browse File"
    bl_options = {'REGISTER', 'UNDO'}
    filepath: bpy.props.StringProperty(name="File Path", subtype='FILE_PATH')
    def execute(self, context):
        _S(context).bb_grid_image_path = bpy.path.abspath(self.filepath)
        self.report({'INFO'}, f"Grid image set: {os.path.basename(_S(context).bb_grid_image_path)}")
        return {'FINISHED'}
    def invoke(self, context, event):
        context.window_manager.fileselect_add(self); return {'RUNNING_MODAL'}

class BB_UVs_ApplyGrid(bpy.types.Operator):
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

class BB_UVs_RevertMaterials(bpy.types.Operator):
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

class BB_UVs_SetGridScale(bpy.types.Operator):
    bl_idname = "bb_uvs.set_grid_scale"
    bl_label = "Set Grid Scale"
    bl_options = {'REGISTER', 'UNDO'}
    scale: bpy.props.EnumProperty(
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

class BB_UVs_SelectZeroArea(bpy.types.Operator):
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
            # Auto-disable UV sync and continue instead of showing error
            bpy.ops.mesh.select_all(action='SELECT')

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
            
        # Auto-execute unfold after selection
        if total > 0:
            bpy.ops.uv.unwrap(method='ANGLE_BASED', margin=0.001)
            
        return {'FINISHED'}

class BB_UVs_SetUI(bpy.types.Operator):
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
        
        # Use safe attribute setting for edge display type
        if hasattr(image_space.uv_editor, 'edge_display_type'):
            try:
                # Try available edge display types
                available_types = ['BLACK', 'WHITE', 'DASH', 'SOLID']
                for display_type in available_types:
                    try:
                        image_space.uv_editor.edge_display_type = display_type
                        break
                    except TypeError:
                        continue
            except Exception:
                pass
        
        # Set other UI properties safely
        if hasattr(image_space.uv_editor, 'tile_grid_shape'):
            try:
                image_space.uv_editor.tile_grid_shape = (10, 10)
            except Exception:
                pass
        if hasattr(image_space.overlay, 'show_grid_background'):
            try:
                image_space.overlay.show_grid_background = True
            except Exception:
                pass
            
        ts = context.tool_settings
        try: 
            ts.mesh_select_mode = (False, True, False)
        except Exception: 
            pass
        try: 
            if hasattr(ts, 'uv_select_mode'):
                ts.uv_select_mode = 'EDGE'
        except Exception: 
            pass
            
        self.report({'INFO'}, "UI set: UV grid on, selection set to Edge."); return {'FINISHED'}

# ---------- Projector operators ----------

def _prepare_targets(ctx):
    sel = [o for o in ctx.selected_objects if o.type == 'MESH']
    if not sel:
        return None
    for o in sel:
        o["uvproj_target"] = True
        ensure_uv(o)
    return sel

class UVPROJ_OT_add_planar(bpy.types.Operator):
    bl_idname = "uvproj.add_planar"
    bl_label = "Add Planar Projector"

    def execute(self, ctx):
        if not _prepare_targets(ctx):
            self.report({'ERROR'}, "Select meshes first")
            return {'CANCELLED'}

        bpy.ops.mesh.primitive_plane_add(size=2)
        proj = ctx.active_object
        proj["uvproj_mode"] = "PLANE"
        proj.display_type = 'WIRE'
        proj.show_in_front = True

        ctx.scene.uvproj_projector = proj
        ctx.scene.uvproj_running = True
        projector_update(ctx.scene)
        return {'FINISHED'}

class UVPROJ_OT_add_cyl(bpy.types.Operator):
    bl_idname = "uvproj.add_cyl"
    bl_label = "Add Cylindrical Projector"

    def execute(self, ctx):
        if not _prepare_targets(ctx):
            self.report({'ERROR'}, "Select meshes first")
            return {'CANCELLED'}

        bpy.ops.mesh.primitive_cylinder_add(radius=1, depth=2)
        proj = ctx.active_object
        proj["uvproj_mode"] = "CYL"
        proj.display_type = 'WIRE'
        proj.show_in_front = True

        ctx.scene.uvproj_projector = proj
        ctx.scene.uvproj_running = True
        projector_update(ctx.scene)
        return {'FINISHED'}

class UVPROJ_OT_add_cube(bpy.types.Operator):
    bl_idname = "uvproj.add_cube"
    bl_label = "Add Cube Projector"

    def execute(self, ctx):
        if not _prepare_targets(ctx):
            self.report({'ERROR'}, "Select meshes first")
            return {'CANCELLED'}

        bpy.ops.mesh.primitive_cube_add(size=2)
        proj = ctx.active_object
        proj["uvproj_mode"] = "CUBE"
        proj.display_type = 'WIRE'
        proj.show_in_front = True

        ctx.scene.uvproj_projector = proj
        ctx.scene.uvproj_running = True
        projector_update(ctx.scene)
        return {'FINISHED'}

class UVPROJ_OT_apply(bpy.types.Operator):
    bl_idname = "uvproj.apply"
    bl_label = "Apply & Remove Projector"

    def execute(self, ctx):
        targets = get_targets()
        ctx.scene.uvproj_running = False

        if "uvproj_last_matrix" in ctx.scene:
            del ctx.scene["uvproj_last_matrix"]

        for o in targets:
            del o["uvproj_target"]

        proj = ctx.scene.uvproj_projector
        if proj:
            bpy.data.objects.remove(proj, do_unlink=True)

        ctx.scene.uvproj_projector = None

        if targets:
            active = ctx.view_layer.objects.active
            if active is None or active not in targets:
                ctx.view_layer.objects.active = targets[0]
                targets[0].select_set(True)

        return {'FINISHED'}

# Collect all operator classes for registration
operator_classes = (
    BB_Texel_Density_Check,
    BB_Texel_Density_Set,
    BB_UVs_SetMoveCollection,
    BB_UVs_SetMoveContext,
    BB_UVs_ToggleMoveHighlight,
    BB_UVs_MoveUVs,
    BB_UVs_SelectSimilarTopology,
    BB_UVs_SelectSimilarTopologyScale,
    BB_UVs_MoveUVsInteractive,
    BB_UVs_RotateUVs,
    BB_UVs_FlipUVs,
    BB_UVs_NormalizePack,
    BB_UVs_PackIndividually,
    BB_UVs_SelectFlipped,
    BB_UVs_SelectCrossUDIMs,
    BB_UVs_SelectBoundaries,
    BB_UVs_ChooseGridImage,
    BB_UVs_ApplyGrid,
    BB_UVs_RevertMaterials,
    BB_UVs_SetGridScale,
    BB_UVs_SelectZeroArea,
    BB_UVs_SetUI,
    UVPROJ_OT_add_planar,
    UVPROJ_OT_add_cyl,
    UVPROJ_OT_add_cube,
    UVPROJ_OT_apply,
)