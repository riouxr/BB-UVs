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



bl_info = {
    "name": "BB UVs",
    "author": "Robert + ChatGPT",
    "version": (3, 0, 3),
    "blender": (3, 6, 0),
    "category": "UV",
}

import bpy
from .helpers import *
from .operators import *
from .panels import *


def _uvproj_props_update(self, context):
    """Re-fit and re-project when axis / keep-ratio changes."""
    scene = context.scene
    proj = scene.uvproj_projector
    if not proj or proj.name not in bpy.data.objects:
        return

    mode = proj.get("uvproj_mode", "PLANE")

    # Re-fit projector to its targets' bounding box
    fit_projector(context, proj, mode)

    # Re-apply projection to all targets
    projector_update()


# -------------------------------------------------------------------------
# Register
# -------------------------------------------------------------------------
def register():
    # Scene properties (no PropertyGroup needed)
    bpy.types.Scene.bb_density = bpy.props.FloatProperty(
        name="Texel Density",
        default=0.0, min=0.0, precision=3, step=1,
        description="Texel density (px/cm)"
    )

    bpy.types.Scene.bb_move_amount = bpy.props.FloatProperty(
        name="Amount",
        default=1.0, min=0.0, soft_max=10.0,
        description="Distance in UDIM tiles (1.0 = one tile)"
    )

    bpy.types.Scene.bb_move_selected_only = bpy.props.BoolProperty(
        name="Contextual Toggle",
        default=True,
        description="Edit Mode: only selected UVs. Object Mode: only active mesh."
    )

    bpy.types.Scene.bb_grid_image_path = bpy.props.StringProperty(
        name="Grid Image",
        default="", subtype='FILE_PATH',
        description="Path to the grid/UV texture"
    )

    bpy.types.Scene.bb_keep_td = bpy.props.BoolProperty(
        name="Keep TD",
        default=False,
        description="Keep current texel density (skip normalization and scaling during pack)"
    )

    bpy.types.Scene.bb_move_collection = bpy.props.BoolProperty(
        name="Collection",
        default=False,
        description="When ON, move all UVs from objects in the active object's collection"
    )

    bpy.types.Scene.bb_pack_rotate = bpy.props.BoolProperty(
        name="Rotation",
        default=True,
        description="Allow rotation during packing for better efficiency"
    )

    # Projector-specific scene props
    bpy.types.Scene.uvproj_projector = bpy.props.PointerProperty(type=bpy.types.Object)

    bpy.types.Scene.uvproj_running = bpy.props.BoolProperty(default=False)

    bpy.types.Scene.uvproj_keep_ratio = bpy.props.BoolProperty(
        name="Keep Ratio",
        default=True,
        description="Uniformly fit projector to largest bounding box dimension",
        update=_uvproj_props_update,
    )

    bpy.types.Scene.uvproj_axis = bpy.props.EnumProperty(
        name="Axis",
        description="Projection axis",
        items=[
            ('X', 'X', 'Project along X axis'),
            ('Y', 'Y', 'Project along Y axis'),
            ('Z', 'Z', 'Project along Z axis'),
        ],
        default='Z',
        update=_uvproj_props_update,
    )

    # Register classes from submodules
    for cls in operator_classes + panel_classes:
        bpy.utils.register_class(cls)

    # Register handler ONCE
    if projector_update not in bpy.app.handlers.depsgraph_update_post:
        bpy.app.handlers.depsgraph_update_post.append(projector_update)


# -------------------------------------------------------------------------
# Unregister
# -------------------------------------------------------------------------
def unregister():

    # Remove projector handler once
    if projector_update in bpy.app.handlers.depsgraph_update_post:
        bpy.app.handlers.depsgraph_update_post.remove(projector_update)

    # Unregister classes in reverse
    for cls in reversed(operator_classes + panel_classes):
        bpy.utils.unregister_class(cls)

    # Clean up scene props
    del bpy.types.Scene.bb_density
    del bpy.types.Scene.bb_move_amount
    del bpy.types.Scene.bb_move_selected_only
    del bpy.types.Scene.bb_grid_image_path
    del bpy.types.Scene.bb_keep_td
    del bpy.types.Scene.bb_move_collection
    del bpy.types.Scene.bb_pack_rotate

    del bpy.types.Scene.uvproj_projector
    del bpy.types.Scene.uvproj_running
    del bpy.types.Scene.uvproj_keep_ratio
    del bpy.types.Scene.uvproj_axis


if __name__ == "__main__":
    register()
