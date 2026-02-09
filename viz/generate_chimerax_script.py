import os

def create_cxc_script(pdb_path, output_png):
    """Generates a ChimeraX script to visualize binder-target interface."""
    script = f"""
# Open the structure
open {pdb_path}

# Color the target (Chain A) and binder (Chain B)
color /A white
color /B cyan

# Show the interface
surface /A
style /B stick
graphics shadows true
lighting soft

# Snapshot
save {output_png} width 1920 height 1080
exit
"""
    with open("render.cxc", "w") as f:
        f.write(script)
    print(f"✅ ChimeraX script minted for {pdb_path}")
