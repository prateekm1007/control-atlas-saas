#!/usr/bin/env python3
"""
Entry 045 — Conjecture System
Generates falsifiable conjectures from the Knowledge Graph.
"""

import json
import networkx as nx
from pathlib import Path
from collections import defaultdict

class ConjectureEngine:
    def __init__(self, kg_path):
        with open(kg_path, 'r') as f:
            data = json.load(f)
        # Handle NetworkX JSON format (node-link)
        self.graph = nx.node_link_graph(data, edges="links")
        
    def generate_conjectures(self):
        conjectures = []
        
        # 1. Near-Miss Conjectures
        # Hypothesis failed exactly ONE constraint
        hypotheses = [n for n, d in self.graph.nodes(data=True) if d.get('type') == 'Hypothesis']
        
        for hypo in hypotheses:
            if self.graph.nodes[hypo].get('result') != 'PROVEN FALSE':
                continue
                
            violations = []
            for u, v, d in self.graph.out_edges(hypo, data=True):
                if d.get('relation') == 'VIOLATES':
                    violations.append(v)
            
            if len(violations) == 1:
                constraint_id = violations[0]
                constraint_desc = self.graph.nodes[constraint_id].get('description', 'Unknown')
                conjectures.append(f"NEAR-MISS: {hypo} failed only on {constraint_desc}. Relaxing this constraint may rescue the hypothesis.")

        # 2. Dominant Failure Modes
        # Which constraints kill the most hypotheses?
        constraint_failures = defaultdict(int)
        for u, v, d in self.graph.edges(data=True):
            if d.get('relation') == 'VIOLATES':
                constraint_failures[v] += 1
        
        sorted_failures = sorted(constraint_failures.items(), key=lambda x: x[1], reverse=True)
        for const, count in sorted_failures:
            desc = self.graph.nodes[const].get('description', const)
            conjectures.append(f"DOMINANT FAILURE: {desc} caused {count} rejections.")

        return conjectures

if __name__ == "__main__":
    # Test with the graph generated in Entry 044
    # Ensure this path matches where Entry 044 saves the file
    kg_path = Path(__file__).resolve().parents[2] / "entries/044_knowledge_graph/test_kg.json"
    
    if kg_path.exists():
        system = ConjectureEngine(kg_path)
        print("Generated Conjectures:")
        for c in system.generate_conjectures():
            print(f"- {c}")
    else:
        print(f"[!] Graph file not found at {kg_path}")
        print("    Run Entry 044 (graph_builder.py) first.")
