import nbformat as nbf

nb = nbf.v4.new_notebook()

# Cell 1: Setup
nb.cells.append(nbf.v4.new_code_cell("""
!git clone https://github.com/prateekm1007/control-atlas.git
%cd control-atlas
!pip install rdkit-pypi networkx
!pip install kaggle
"""))

# Cell 2: RFAA Execution
nb.cells.append(nbf.v4.new_code_cell("""
# Run the worker (Mocked for template)
!chmod +x entries/051_gpu_orchestration/spot_worker.sh
!./entries/051_gpu_orchestration/spot_worker.sh batch_001.json
"""))

# Cell 3: Persistence
nb.cells.append(nbf.v4.new_code_cell("""
# Sync results
!python entries/052_kaggle_adapter/kaggle_sync.py rfaa_results/ username/control-atlas-results
"""))

with open('entries/052_kaggle_adapter/Control_Atlas_RFAA_Runner.ipynb', 'w') as f:
    nbf.write(nb, f)
    
print("[+] Notebook generated.")
