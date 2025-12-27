import bpy
from .helpers import *

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
    f = layout.row(); f.scale_y = 1.1
    f.operator("bb_uvs.select_flipped", text="Select Flipped", icon='UV_SYNC_SELECT')
    op = f.operator("bb_uvs.flip_uvs", text="Flip Horizontal", icon='ARROW_LEFTRIGHT')
    op.axis = 'H'

    # Projector section (integrated after Rotate)
    layout.separator()
    layout.label(text="Projector")
    col = layout.column(align=True)
    col.operator("bb_uvs.projector_add_plane")
    col.operator("bb_uvs.projector_add_cylinder")
    col.operator("bb_uvs.projector_add_cube")
    col.separator()
    col.operator("bb_uvs.projector_apply")

    layout.separator()
    # Only show "Check and fix" in the UV Editor (IMAGE_EDITOR)
    if context.space_data.type == 'IMAGE_EDITOR':
        layout.separator()
        
        if context.mode == 'EDIT_MESH':
            layout.label(text="Check and fix")
            f = layout.row(); f.scale_y = 1.1
            f = layout.row(); f.scale_y = 1.1
            f.operator("bb_uvs.select_cross_udims", text="", icon='SNAP_GRID')
            f.operator("bb_uvs.select_boundaries", text="", icon='EDGESEL')
            f.operator("uv.select_overlap", text="", icon='SELECT_INTERSECT')
            f.operator("bb_uvs.select_zero_area", text="", icon='CANCEL')

class UV_PT_BB_UVs(bpy.types.Panel):
    bl_label = "BB UVs"
    bl_space_type = 'IMAGE_EDITOR'
    bl_region_type = 'UI'
    bl_category = "BB UVs"
    @classmethod
    def poll(cls, context):
        area = context.area
        return area and area.type == 'IMAGE_EDITOR' and getattr(context.space_data, "ui_type", "UV") == "UV"
    def draw(self, context):
        draw_bb_uvs_panel(self, context)

class VIEW3D_PT_BB_UVs(bpy.types.Panel):
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
    def draw(self, context):
        draw_bb_uvs_panel(self, context)

# Collect all panel classes for registration
panel_classes = (
    UV_PT_BB_UVs,
    VIEW3D_PT_BB_UVs,
)