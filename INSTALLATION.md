# Installation Guide

## Quick Start

Clone and install:
- git clone https://github.com/prateekm1007/control-atlas.git
- cd control-atlas
- conda create -n control-atlas python=3.11 -y
- conda activate control-atlas
- pip install numpy scipy biopython openmm pdbfixer

## Verify Installation

python tools/native_audit.py structures/champ005/structure.cif --list-chains

Expected output: A: 935 atoms, B: 185 atoms
