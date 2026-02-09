#!/usr/bin/env python3
"""
Entry 044 — Scientific Knowledge Graph
Stores proofs and constraints in a graph structure.
"""

import json
import networkx as nx
from pathlib import Path
from datetime import datetime

class KnowledgeGraph:
    def __init__(self):
        self.graph = nx.DiGraph()
        
    def add_proof(self, target_id, compound_id, proof_data):
        """
        Ingest a Proof object (dict) into the graph.
        """
        # Nodes
        self.graph.add_node(target_id, type="Target")
        self.graph.add_node(compound_id, type="Compound")
        
        # The Hypothesis Node (Intersection)
        hypothesis_id = f"{target_id}_{compound_id}"
        self.graph.add_node(hypothesis_id, type="Hypothesis", 
                            result=proof_data["final_conclusion"])
        
        self.graph.add_edge(target_id, hypothesis_id, relation="PART_OF")
        self.graph.add_edge(compound_id, hypothesis_id, relation="PART_OF")
        
        # Constraint Nodes
        for step in proof_data["steps"]:
            constraint_id = f"Constraint_{step['axiom']}"
            self.graph.add_node(constraint_id, type="Constraint", 
                                description=step['premise'])
            
            # Edge: Hypothesis -> Constraint
            relation = "SATISFIES" if step['conclusion'] == "True" else "VIOLATES"
            self.graph.add_edge(hypothesis_id, constraint_id, 
                                relation=relation, 
                                value=step['evidence'])

    def export_json(self, path):
        # FIX: Explicit edges kwarg to silence FutureWarning
        data = nx.node_link_data(self.graph, edges="links")
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
            
    def stats(self):
        return {
            "nodes": self.graph.number_of_nodes(),
            "edges": self.graph.number_of_edges(),
            "violations": len([u for u, v, d in self.graph.edges(data=True) if d['relation'] == 'VIOLATES'])
        }

if __name__ == "__main__":
    # Test
    kg = KnowledgeGraph()
    
    # Mock Proof Data (from Entry 043)
    proof = {
        "final_conclusion": "PROVEN FALSE",
        "steps": [
            {"axiom": "Existence", "premise": "Vol > 150", "evidence": 120, "conclusion": "False"},
            {"axiom": "Stability", "premise": "pLDDT > 70", "evidence": 85, "conclusion": "True"}
        ]
    }
    
    kg.add_proof("EGFR", "Mol_X", proof)
    print(kg.stats())
    kg.export_json("test_kg.json")
