import os
import subprocess
import sys

class ChaiBackend:
    """
    Chai-1 inference wrapper.
    Prepares FASTA inputs and executes structure prediction.
    """
    
    VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
    
    def __init__(self, output_dir="results/v6_discovery"):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def validate_sequence(self, seq, name="sequence"):
        """Verify valid amino acid characters."""
        seq_upper = seq.upper().replace(" ", "").replace("\n", "")
        invalid = set(seq_upper) - self.VALID_AA
        if invalid:
            raise ValueError(f"Invalid amino acids in {name}: {invalid}")
        if len(seq_upper) < 5:
            raise ValueError(f"Sequence {name} too short: {len(seq_upper)} aa")
        return seq_upper

    def prepare_fasta(self, target_seq, binder_seq, name):
        """Creates multi-entity FASTA for Chai-1."""
        target_clean = self.validate_sequence(target_seq, "target")
        binder_clean = self.validate_sequence(binder_seq, "binder")
        
        fasta_content = f">protein|name=TARGET\n{target_clean}\n>protein|name={name}\n{binder_clean}\n"
        path = os.path.join(self.output_dir, f"{name}.fasta")
        
        with open(path, "w") as f:
            f.write(fasta_content)
        
        print(f"[FASTA] Written: {path}")
        print(f"  Target: {len(target_clean)} aa")
        print(f"  Binder: {len(binder_clean)} aa")
        return path

    def run_inference(self, fasta_path, num_models=1):
        """
        Execute Chai-1 prediction.
        Requires chai-lab to be installed in environment.
        """
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"FASTA not found: {fasta_path}")
        
        output_subdir = os.path.join(
            self.output_dir, 
            os.path.splitext(os.path.basename(fasta_path))[0]
        )
        os.makedirs(output_subdir, exist_ok=True)
        
        cmd = [
            sys.executable, "-m", "chai_lab.chai1", 
            "--input", fasta_path,
            "--output", output_subdir,
            "--num-models", str(num_models)
        ]
        
        print(f"[INFERENCE] Running: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            print(f"[INFERENCE] Complete. Output: {output_subdir}")
            return output_subdir
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Inference failed:\n{e.stderr}")
            raise
        except FileNotFoundError:
            raise RuntimeError("chai_lab not installed. Run: pip install chai-lab")


if __name__ == "__main__":
    # Self-test
    backend = ChaiBackend()
    
    # Test validation
    try:
        backend.validate_sequence("ACDEFGHIKLMNPQRSTVWY", "test_valid")
        print("[TEST] Valid sequence: PASS")
    except ValueError as e:
        print(f"[TEST] Valid sequence: FAIL - {e}")
    
    try:
        backend.validate_sequence("ACDEFXXX", "test_invalid")
        print("[TEST] Invalid sequence: FAIL - should have raised error")
    except ValueError:
        print("[TEST] Invalid sequence: PASS (correctly rejected)")
