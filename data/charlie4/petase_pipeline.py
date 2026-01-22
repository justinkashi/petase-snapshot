"""
PETase Variant Prediction Pipeline - NDCG@10 Optimized
========================================================

Complete zero-shot prediction pipeline optimized for NDCG@10 ranking metric.
Generates multiple submission strategies to maximize top 10% ranking accuracy.

Key Features:
- All 4 phases with 31 features
- Ensemble ranking optimization
- Multiple submission strategies
- Top 10% confidence analysis
- Automatic checkpointing

Usage in Colab:
    !python petase_pipeline_ndcg_optimized.py --input test-2025.csv --cds pet-2025-wildtype-cds.csv --output submission.csv --token YOUR_TOKEN

Author: Charlie
Date: 2025-01
Optimized for: NDCG@10 competition metric
"""

import pandas as pd
import numpy as np
import warnings
from pathlib import Path
from tqdm.auto import tqdm
import json
import pickle
from typing import Dict, List, Tuple, Optional
import argparse
import sys
from scipy.stats import rankdata

warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

class Config:
    """Configuration for the pipeline"""
    
    # ESM3 API Configuration
    ESM3_MODEL = "esm3-open"
    ESM3_URL = "https://forge.evolutionaryscale.ai"
    ESM3_TOKEN = None
    
    # Catalytic triad positions
    CATALYTIC_POSITIONS = {'Ser': 160, 'His': 237, 'Asp': 206}
    
    # Amino acid properties
    AA_PKA = {
        'D': {'pKa': 3.9, 'charge_at_pH5.5': -0.8, 'charge_at_pH9.0': -1.0},
        'E': {'pKa': 4.2, 'charge_at_pH5.5': -0.7, 'charge_at_pH9.0': -1.0},
        'H': {'pKa': 6.0, 'charge_at_pH5.5': +0.3, 'charge_at_pH9.0': 0.0},
        'K': {'pKa': 10.5, 'charge_at_pH5.5': +1.0, 'charge_at_pH9.0': +0.7},
        'R': {'pKa': 12.5, 'charge_at_pH5.5': +1.0, 'charge_at_pH9.0': +1.0},
    }
    
    # Chou-Fasman propensities
    HELIX_PROPENSITY = {
        'A': 1.42, 'E': 1.53, 'L': 1.21, 'M': 1.45, 'Q': 1.17, 'K': 1.16, 'R': 0.98,
        'H': 1.00, 'V': 1.06, 'I': 1.08, 'Y': 0.69, 'C': 0.70, 'W': 1.08, 'F': 1.13,
        'T': 0.83, 'S': 0.77, 'N': 0.67, 'D': 1.01, 'G': 0.57, 'P': 0.57
    }
    
    SHEET_PROPENSITY = {
        'V': 1.70, 'I': 1.60, 'Y': 1.47, 'F': 1.38, 'W': 1.37, 'L': 1.30, 'T': 1.19,
        'C': 1.19, 'M': 1.05, 'A': 0.83, 'G': 0.75, 'R': 0.93, 'K': 0.74, 'S': 0.75,
        'Q': 1.10, 'H': 0.87, 'N': 0.89, 'P': 0.55, 'D': 0.54, 'E': 0.37
    }
    
    # Hydrophobicity scale
    HYDROPHOBICITY = {
        'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
        'G': -0.4, 'T': -0.7, 'W': -0.9, 'S': -0.8, 'Y': -1.3, 'P': -1.6,
        'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
    }
    
    # Amino acid sizes
    AA_SIZE = {
        'G': 1, 'A': 2, 'S': 2, 'C': 2, 'D': 2, 'P': 2, 'N': 2, 'T': 3,
        'E': 3, 'V': 3, 'Q': 3, 'H': 3, 'M': 4, 'I': 4, 'L': 4, 'K': 4,
        'R': 4, 'F': 5, 'Y': 6, 'W': 6
    }
    
    # E. coli codon usage
    ECOLI_CODON_USAGE = {
        'GCG': 0.36, 'GCA': 0.21, 'GCT': 0.16, 'GCC': 0.27,
        'TGC': 0.55, 'TGT': 0.45,
        'GAT': 0.63, 'GAC': 0.37,
        'GAA': 0.68, 'GAG': 0.32,
        'TTT': 0.57, 'TTC': 0.43,
        'GGT': 0.35, 'GGC': 0.40, 'GGA': 0.11, 'GGG': 0.15,
        'CAT': 0.57, 'CAC': 0.43,
        'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11,
        'AAA': 0.74, 'AAG': 0.26,
        'CTG': 0.49, 'TTG': 0.13, 'CTT': 0.10, 'CTC': 0.10, 'CTA': 0.04, 'TTA': 0.13,
        'ATG': 1.00,
        'AAT': 0.45, 'AAC': 0.55,
        'CCG': 0.52, 'CCA': 0.19, 'CCT': 0.16, 'CCC': 0.12,
        'CAA': 0.34, 'CAG': 0.66,
        'CGT': 0.38, 'CGC': 0.40, 'CGA': 0.07, 'CGG': 0.10, 'AGA': 0.04, 'AGG': 0.02,
        'AGC': 0.28, 'TCT': 0.15, 'TCC': 0.15, 'AGT': 0.15, 'TCA': 0.12, 'TCG': 0.15,
        'ACC': 0.44, 'ACT': 0.17, 'ACA': 0.13, 'ACG': 0.27,
        'GTG': 0.37, 'GTT': 0.26, 'GTC': 0.20, 'GTA': 0.17,
        'TGG': 1.00,
        'TAT': 0.57, 'TAC': 0.43,
        'TAA': 0.61, 'TAG': 0.09, 'TGA': 0.30,
    }
    
    RARE_CODONS = {codon for codon, freq in ECOLI_CODON_USAGE.items() if freq < 0.10}
    
    # Base model configurations
    MODELS = {
        'activity_1': {
            'features': ['mutation_llr', 'charge_change_pH5.5', 'introduces_charged_at_pH5.5',
                        'dist_to_active_site', 'structure_risk_score'],
            'weights': [0.30, 0.20, 0.10, 0.20, 0.20],
            'invert': [False, False, False, False, True],
            'output_range': (0, 10)
        },
        'activity_2': {
            'features': ['mutation_llr', 'charge_change_pH9.0', 'introduces_charged_at_pH9.0',
                        'pH_differential', 'dist_to_active_site', 'structure_risk_score'],
            'weights': [0.25, 0.20, 0.10, 0.15, 0.15, 0.15],
            'invert': [False, False, False, False, False, True],
            'output_range': (0, 10)
        },
        'expression': {
            'features': ['pseudo_likelihood', 'cai', 'rare_codon_freq',
                        'rare_codon_cluster_size', 'structure_risk_score'],
            'weights': [0.35, 0.20, 0.20, 0.15, 0.10],
            'invert': [False, False, True, True, True],
            'output_range': (0, 2)
        }
    }


# ============================================================================
# PHASE 1-3: SAME AS BEFORE (Feature Engineering)
# ============================================================================

class Phase1ESM3:
    """ESM3-based evolutionary feature extraction"""
    
    def __init__(self, token: str):
        print("🧬 Initializing ESM3 client...")
        try:
            from esm.models.esm3 import ESM3 
            from esm.sdk.api import ESMProtein, LogitsConfig

            self.model = ESM3.from_pretrained("esm3-open").to("cuda").eval()
            self.ESMProtein = ESMProtein
            self.LogitsConfig = LogitsConfig
            print("✅ ESM3 LOCALLY THIS TIME! initialized successfully")
        except Exception as e:
            print(f"❌ Failed to initialize ESM3 client: {e}")
            raise
    
    def find_parent_wildtype(self, mut_seq: str, wt_sequences: List[str]) -> str:
        min_distance = float('inf')
        parent_wt = None
        for wt in wt_sequences:
            if len(wt) != len(mut_seq):
                continue
            distance = sum(a != b for a, b in zip(wt, mut_seq))
            if distance < min_distance:
                min_distance = distance
                parent_wt = wt
        return parent_wt
    
    def find_mutation_details(self, wt_seq: str, mut_seq: str) -> Tuple[int, str, str]:
        for i, (wt_aa, mut_aa) in enumerate(zip(wt_seq, mut_seq)):
            if wt_aa != mut_aa:
                return i, wt_aa, mut_aa
        return -1, '', ''
    
    def get_position_logits(self, sequence: str) -> np.ndarray:
        """
        Get amino acid logits for each position in a sequence.
        Uses two-step approach: encode then logits.
        """
        try:
            import torch
            
            # Step 1: Encode the protein
            protein = self.ESMProtein(sequence=sequence)
            encoded = self.model.encode(protein)
            output = self.model.logits(encoded, self.LogitsConfig(sequence=True, return_embeddings=True))
                    
            # Extract logits tensor
            logits = output.logits
            if hasattr(logits, 'sequence'):
                logits = logits.sequence
            
            # Convert to numpy array
            # FIX: ESM3 returns BFloat16, which numpy doesn't support
            # Must convert to float32 first
            if isinstance(logits, torch.Tensor):
                logits = logits.detach().to(torch.float32).cpu().numpy()
            else:
                logits = np.array(logits, dtype=np.float32)
            
            # Handle shape: could be (1, L, V) or (L, V)
            if logits.ndim == 3:
                logits = logits[0]  # Remove batch dimension
            
            return logits.astype(np.float32)
            
        except Exception as e:
            print(f"⚠️  ESM3 API error: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    @staticmethod
    def logits_to_log_probs(logits):
        """Convert logits to log probabilities using log-softmax (numerically stable)."""
        max_logits = np.max(logits, axis=-1, keepdims=True)
        shifted_logits = logits - max_logits
        exp_logits = np.exp(shifted_logits)
        sum_exp = np.sum(exp_logits, axis=-1, keepdims=True)
        log_probs = shifted_logits - np.log(sum_exp)
        return log_probs
    
    def compute_mutation_llr(self, wt_seq: str, mut_seq: str, position: int) -> float:
        wt_logits = self.get_position_logits(wt_seq)
        mut_logits = self.get_position_logits(mut_seq)
        
        if wt_logits is None or mut_logits is None:
            return 0.0
        
        # Convert to log probabilities (numerically stable)
        wt_log_probs = self.logits_to_log_probs(wt_logits)
        mut_log_probs = self.logits_to_log_probs(mut_logits)
        
        aa_to_idx = {aa: i for i, aa in enumerate('ACDEFGHIKLMNPQRSTVWY')}
        wt_aa = wt_seq[position]
        mut_aa = mut_seq[position]
        
        wt_idx = aa_to_idx.get(wt_aa, 0)
        mut_idx = aa_to_idx.get(mut_aa, 0)
        
        llr = mut_log_probs[position][mut_idx] - wt_log_probs[position][wt_idx]
        return float(llr)
    
    def compute_pseudo_likelihood(self, sequence: str) -> float:
        logits = self.get_position_logits(sequence)
        if logits is None:
            return -6.0
        
        # Convert to log probabilities
        log_probs = self.logits_to_log_probs(logits)
        
        aa_to_idx = {aa: i for i, aa in enumerate('ACDEFGHIKLMNPQRSTVWY')}
        total_log_prob = 0
        valid_positions = 0
        
        for i, aa in enumerate(sequence):
            if aa not in aa_to_idx:
                continue
            aa_idx = aa_to_idx[aa]
            total_log_prob += log_probs[i][aa_idx]
            valid_positions += 1
        
        return total_log_prob / valid_positions if valid_positions > 0 else -6.0
    
    def process_sequences(self, df: pd.DataFrame, wt_sequences: List[str], 
                         checkpoint_file: str = 'phase1_checkpoint.pkl') -> pd.DataFrame:
        print("\n" + "="*70)
        print("PHASE 1: ESM3 Evolutionary Features")
        print("="*70)
        
        checkpoint_path = Path(checkpoint_file)
        if checkpoint_path.exists():
            print(f"📂 Loading checkpoint from {checkpoint_file}")
            with open(checkpoint_file, 'rb') as f:
                checkpoint = pickle.load(f)
            df = checkpoint['df']
            start_idx = checkpoint['last_index'] + 1
            print(f"✅ Resuming from index {start_idx}")
        else:
            df['is_wildtype'] = False
            df['parent_wt'] = ''
            df['mutation_position'] = -1
            df['wt_aa'] = ''
            df['mut_aa'] = ''
            df['mutation_llr'] = 0.0
            df['pseudo_likelihood'] = -6.0
            start_idx = 0
        
        for idx in tqdm(range(start_idx, len(df)), desc="Processing sequences"):
            seq = df.loc[idx, 'sequence']
            
            if seq in wt_sequences:
                df.loc[idx, 'is_wildtype'] = True
                df.loc[idx, 'parent_wt'] = seq
                df.loc[idx, 'pseudo_likelihood'] = self.compute_pseudo_likelihood(seq)
            else:
                parent_wt = self.find_parent_wildtype(seq, wt_sequences)
                if parent_wt:
                    df.loc[idx, 'parent_wt'] = parent_wt
                    position, wt_aa, mut_aa = self.find_mutation_details(parent_wt, seq)
                    df.loc[idx, 'mutation_position'] = position
                    df.loc[idx, 'wt_aa'] = wt_aa
                    df.loc[idx, 'mut_aa'] = mut_aa
                    
                    if position >= 0:
                        df.loc[idx, 'mutation_llr'] = self.compute_mutation_llr(parent_wt, seq, position)
                    df.loc[idx, 'pseudo_likelihood'] = self.compute_pseudo_likelihood(seq)
            
            if (idx + 1) % 100 == 0:
                checkpoint = {'df': df, 'last_index': idx}
                with open(checkpoint_file, 'wb') as f:
                    pickle.dump(checkpoint, f)
        
        checkpoint = {'df': df, 'last_index': len(df) - 1}
        with open(checkpoint_file, 'wb') as f:
            pickle.dump(checkpoint, f)
        
        print(f"\n✅ Phase 1 complete: {len(df)} sequences processed")
        return df


class Phase2Structure:
    """Structure-based feature extraction"""
    
    @staticmethod
    def estimate_distance(mutation_pos: int) -> float:
        catalytic_positions = list(Config.CATALYTIC_POSITIONS.values())
        seq_distances = [abs(mutation_pos - cat_pos) for cat_pos in catalytic_positions]
        min_seq_dist = min(seq_distances)
        
        if min_seq_dist < 5:
            correction = 0.9
        elif min_seq_dist < 15:
            correction = 1.2
        else:
            correction = 1.5
        
        return min_seq_dist * 3.8 * correction
    
    @staticmethod
    def predict_secondary_structure(sequence: str, position: int) -> str:
        if position < 3 or position >= len(sequence) - 3:
            return 'coil'
        
        window = sequence[position-3:position+4]
        helix_score = sum(Config.HELIX_PROPENSITY.get(aa, 1.0) for aa in window) / len(window)
        sheet_score = sum(Config.SHEET_PROPENSITY.get(aa, 1.0) for aa in window) / len(window)
        
        if helix_score > 1.03 and helix_score > sheet_score:
            return 'helix'
        elif sheet_score > 1.05:
            return 'sheet'
        else:
            return 'coil'
    
    @staticmethod
    def estimate_burial(sequence: str, position: int, aa: str) -> str:
        h = Config.HYDROPHOBICITY.get(aa, 0)
        ss = Phase2Structure.predict_secondary_structure(sequence, position)
        d_term = min(position, len(sequence) - position)
        
        if d_term < 10:
            return 'exposed'
        elif ss in ['helix', 'sheet'] and h > 2.0:
            return 'buried'
        elif h > 0:
            return 'partial'
        else:
            return 'exposed'
    
    @staticmethod
    def compute_structure_risk(row: pd.Series) -> float:
        risk = 0.0
        
        dist = row['dist_to_active_site']
        if dist < 10:
            risk += 0.30
        elif dist < 15:
            risk += 0.20
        elif dist < 30:
            risk += 0.10
        
        burial = row.get('burial_status', 'exposed')
        if burial == 'buried':
            risk += 0.20
        elif burial == 'partial':
            risk += 0.10
        
        ss = row.get('secondary_structure', 'coil')
        if ss in ['helix', 'sheet']:
            risk += 0.15
        
        if not row['is_wildtype']:
            wt_aa = row['wt_aa']
            mut_aa = row['mut_aa']
            wt_h = Config.HYDROPHOBICITY.get(wt_aa, 0)
            mut_h = Config.HYDROPHOBICITY.get(mut_aa, 0)
            delta_h = abs(mut_h - wt_h)
            
            if delta_h > 3.0:
                risk += 0.15
            elif delta_h > 1.5:
                risk += 0.08
            
            wt_size = Config.AA_SIZE.get(wt_aa, 3)
            mut_size = Config.AA_SIZE.get(mut_aa, 3)
            delta_size = abs(mut_size - wt_size)
            
            if delta_size >= 3:
                risk += 0.10
            elif delta_size >= 2:
                risk += 0.05
        
        return min(risk, 1.0)
    
    @staticmethod
    def process_features(df: pd.DataFrame) -> pd.DataFrame:
        print("\n" + "="*70)
        print("PHASE 2: Structure Features")
        print("="*70)
        
        print("📏 Computing distances to active site...")
        df['dist_to_active_site'] = df.apply(
            lambda row: Phase2Structure.estimate_distance(row['mutation_position'])
            if not row['is_wildtype'] else 0,
            axis=1
        )
        
        print("🧬 Predicting secondary structure...")
        df['secondary_structure'] = df.apply(
            lambda row: Phase2Structure.predict_secondary_structure(
                row['sequence'], row['mutation_position']
            ) if not row['is_wildtype'] else 'coil',
            axis=1
        )
        
        print("🏔️  Estimating burial status...")
        df['burial_status'] = df.apply(
            lambda row: Phase2Structure.estimate_burial(
                row['sequence'], row['mutation_position'],
                row['mut_aa'] if not row['is_wildtype'] else ''
            ),
            axis=1
        )
        
        df['near_active_site'] = df['dist_to_active_site'] < 15
        
        print("⚠️  Computing structure risk scores...")
        df['structure_risk_score'] = df.apply(Phase2Structure.compute_structure_risk, axis=1)
        
        print(f"✅ Phase 2 complete: {len(df.columns)} features total")
        return df


class Phase3Enhancement:
    """pH-dependent and codon optimization features"""
    
    @staticmethod
    def compute_charge_change(wt_aa: str, mut_aa: str, pH: float) -> float:
        pH_key = f'charge_at_pH{pH}'
        wt_charge = Config.AA_PKA.get(wt_aa, {}).get(pH_key, 0)
        mut_charge = Config.AA_PKA.get(mut_aa, {}).get(pH_key, 0)
        return mut_charge - wt_charge
    
    @staticmethod
    def compute_cai(cds: str) -> float:
        if len(cds) % 3 != 0:
            return 0.5
        
        codons = [cds[i:i+3] for i in range(0, len(cds), 3)]
        weights = []
        
        for codon in codons:
            codon_upper = codon.upper()
            if codon_upper in Config.ECOLI_CODON_USAGE:
                weights.append(Config.ECOLI_CODON_USAGE[codon_upper])
            else:
                weights.append(0.01)
        
        if not weights:
            return 0.5
        
        from scipy.stats import gmean
        return gmean(weights)
    
    @staticmethod
    def compute_rare_codon_freq(cds: str) -> float:
        if len(cds) % 3 != 0:
            return 0
        
        codons = [cds[i:i+3].upper() for i in range(0, len(cds), 3)]
        rare_count = sum(1 for codon in codons if codon in Config.RARE_CODONS)
        return rare_count / len(codons) if codons else 0
    
    @staticmethod
    def find_rare_codon_clusters(cds: str) -> int:
        if len(cds) % 3 != 0:
            return 0
        
        codons = [cds[i:i+3].upper() for i in range(0, len(cds), 3)]
        max_cluster = 0
        current_cluster = 0
        
        for codon in codons:
            if codon in Config.RARE_CODONS:
                current_cluster += 1
                max_cluster = max(max_cluster, current_cluster)
            else:
                current_cluster = 0
        
        return max_cluster
    
    @staticmethod
    def process_features(df: pd.DataFrame, cds_df: pd.DataFrame) -> pd.DataFrame:
        print("\n" + "="*70)
        print("PHASE 3: pH and Codon Features")
        print("="*70)
        
        print("🧪 Computing pH-dependent charges...")
        df['charge_change_pH5.5'] = df.apply(
            lambda row: Phase3Enhancement.compute_charge_change(
                row['wt_aa'], row['mut_aa'], 5.5
            ) if not row['is_wildtype'] else 0,
            axis=1
        )
        
        df['charge_change_pH9.0'] = df.apply(
            lambda row: Phase3Enhancement.compute_charge_change(
                row['wt_aa'], row['mut_aa'], 9.0
            ) if not row['is_wildtype'] else 0,
            axis=1
        )
        
        df['pH_differential'] = df['charge_change_pH9.0'] - df['charge_change_pH5.5']
        
        df['introduces_charged_at_pH5.5'] = df.apply(
            lambda row: abs(Config.AA_PKA.get(row['mut_aa'], {}).get('charge_at_pH5.5', 0)) > 0.5
            if not row['is_wildtype'] else False,
            axis=1
        ).astype(int)
        
        df['introduces_charged_at_pH9.0'] = df.apply(
            lambda row: abs(Config.AA_PKA.get(row['mut_aa'], {}).get('charge_at_pH9.0', 0)) > 0.5
            if not row['is_wildtype'] else False,
            axis=1
        ).astype(int)
        
        print("🧬 Computing codon usage features...")
        cds_dict = dict(zip(cds_df['sequence'], cds_df['cds']))
        df['cds'] = df['parent_wt'].map(cds_dict)
        
        df['cai'] = df['cds'].apply(
            lambda cds: Phase3Enhancement.compute_cai(cds) if pd.notna(cds) else 0.35
        )
        
        df['rare_codon_freq'] = df['cds'].apply(
            lambda cds: Phase3Enhancement.compute_rare_codon_freq(cds) if pd.notna(cds) else 0
        )
        
        df['gc_content'] = df['cds'].apply(
            lambda cds: (cds.upper().count('G') + cds.upper().count('C')) / len(cds)
            if pd.notna(cds) and len(cds) > 0 else 0.5
        )
        
        df['rare_codon_cluster_size'] = df['cds'].apply(
            lambda cds: Phase3Enhancement.find_rare_codon_clusters(cds) if pd.notna(cds) else 0
        )
        
        print(f"✅ Phase 3 complete: {len(df.columns)} features total")
        df = df.drop('cds', axis=1)
        
        return df


# ============================================================================
# PHASE 4: NDCG-OPTIMIZED PREDICTIONS
# ============================================================================

class Phase4NDCGOptimized:
    """Final predictions optimized for NDCG@10 metric"""
    
    @staticmethod
    def predict_base_model(df: pd.DataFrame, model_config: Dict) -> np.ndarray:
        """Original linear model predictions"""
        features = model_config['features']
        weights = np.array(model_config['weights'])
        invert = model_config['invert']
        output_range = model_config['output_range']
        
        X = df[features].values
        X_min = np.nanmin(X, axis=0)
        X_max = np.nanmax(X, axis=0)
        X_norm = (X - X_min) / (X_max - X_min + 1e-10)
        X_norm = np.nan_to_num(X_norm, nan=0.5)
        
        for i, should_invert in enumerate(invert):
            if should_invert:
                X_norm[:, i] = 1 - X_norm[:, i]
        
        scores = np.dot(X_norm, weights)
        predictions = scores * (output_range[1] - output_range[0]) + output_range[0]
        
        return predictions
    
    @staticmethod
    def compute_confidence(df: pd.DataFrame, target: str) -> np.ndarray:
        """
        Compute prediction confidence for each variant
        High confidence = features strongly agree
        """
        confidence = np.zeros(len(df))
        
        # Strong evolutionary signal
        strong_llr = np.abs(df['mutation_llr']) > 5
        confidence += strong_llr * 0.3
        
        # Low structural risk
        low_risk = df['structure_risk_score'] < 0.3
        confidence += low_risk * 0.2
        
        # Features agreement
        if 'activity' in target:
            # For activity: LLR positive + low risk + moderate charge change
            good_llr = df['mutation_llr'] > 0
            low_charge = np.abs(df['charge_change_pH5.5']) < 1.0
            feature_agree = good_llr & low_risk & low_charge
        else:
            # For expression: high PLL + good CAI + low rare codons
            good_pll = df['pseudo_likelihood'] > -5.5
            good_cai = df['cai'] > 0.33
            feature_agree = good_pll & good_cai & low_risk
        
        confidence += feature_agree * 0.5
        
        return confidence
    
    @staticmethod
    def create_ranking_strategies(df: pd.DataFrame, target: str) -> Dict[str, np.ndarray]:
        """
        Create multiple ranking strategies for ensemble
        """
        rankings = {}
        
        # Strategy 1: Base model
        base_pred = df[f'{target}_base']
        rankings['base_model'] = rankdata(-base_pred)  # Negative for descending
        
        # Strategy 2: Key feature only (most reliable)
        if 'activity' in target:
            key_feature = df['mutation_llr']
        else:
            key_feature = df['pseudo_likelihood']
        rankings['key_feature'] = rankdata(-key_feature)
        
        # Strategy 3: Confidence-weighted
        confidence = Phase4NDCGOptimized.compute_confidence(df, target)
        weighted_pred = base_pred * (1 + confidence)
        rankings['confidence_weighted'] = rankdata(-weighted_pred)
        
        # Strategy 4: Multi-feature median rank
        if 'activity' in target:
            features = ['mutation_llr', 'structure_risk_score', 'dist_to_active_site']
            individual_ranks = [
                rankdata(-df['mutation_llr']),
                rankdata(df['structure_risk_score']),  # Lower is better
                rankdata(df['dist_to_active_site'])    # Farther is better
            ]
        else:
            features = ['pseudo_likelihood', 'cai', 'rare_codon_freq']
            individual_ranks = [
                rankdata(-df['pseudo_likelihood']),
                rankdata(-df['cai']),
                rankdata(df['rare_codon_freq'])  # Lower is better
            ]
        
        median_rank = np.median(individual_ranks, axis=0)
        rankings['median_rank'] = median_rank
        
        # Strategy 5: Biology-constrained
        bio_score = base_pred.copy()
        
        # Penalize risky mutations
        high_risk = df['structure_risk_score'] > 0.6
        bio_score[high_risk] *= 0.5
        
        # Penalize active site mutations with negative LLR
        active_site_bad = (df['dist_to_active_site'] < 15) & (df['mutation_llr'] < 0)
        bio_score[active_site_bad] *= 0.3
        
        rankings['bio_constrained'] = rankdata(-bio_score)
        
        return rankings
    
    @staticmethod
    def ensemble_rankings(rankings: Dict[str, np.ndarray], 
                         strategy: str = 'balanced') -> np.ndarray:
        """
        Ensemble multiple rankings with different weighting strategies
        """
        if strategy == 'balanced':
            # Equal weight to all strategies
            weights = {
                'base_model': 0.30,
                'key_feature': 0.25,
                'confidence_weighted': 0.20,
                'median_rank': 0.15,
                'bio_constrained': 0.10
            }
        elif strategy == 'conservative':
            # Trust key feature more
            weights = {
                'base_model': 0.20,
                'key_feature': 0.40,
                'confidence_weighted': 0.15,
                'median_rank': 0.15,
                'bio_constrained': 0.10
            }
        elif strategy == 'aggressive':
            # Trust model more
            weights = {
                'base_model': 0.50,
                'key_feature': 0.15,
                'confidence_weighted': 0.20,
                'median_rank': 0.10,
                'bio_constrained': 0.05
            }
        else:  # 'diverse'
            # More even distribution
            weights = {
                'base_model': 0.20,
                'key_feature': 0.20,
                'confidence_weighted': 0.20,
                'median_rank': 0.20,
                'bio_constrained': 0.20
            }
        
        final_rank = sum(weights[name] * rank for name, rank in rankings.items())
        return final_rank
    
    @staticmethod
    def rank_to_prediction(ranks: np.ndarray, output_range: Tuple[float, float]) -> np.ndarray:
        """
        Convert ranks back to predictions in biological range
        Uses smooth percentile-based mapping
        """
        percentiles = rankdata(ranks) / len(ranks)
        predictions = (1 - percentiles) * (output_range[1] - output_range[0]) + output_range[0]
        return predictions
    
    @staticmethod
    def analyze_top_variants(df: pd.DataFrame, target: str, top_k: int = 500):
        """
        Analyze characteristics of top 10% variants
        """
        top_variants = df.nlargest(top_k, target)
        
        print(f"\n📊 Top {top_k} variants for {target}:")
        print(f"   Mean LLR: {top_variants['mutation_llr'].mean():.2f}")
        print(f"   Mean structure risk: {top_variants['structure_risk_score'].mean():.3f}")
        print(f"   Mean distance: {top_variants['dist_to_active_site'].mean():.1f} Å")
        
        if 'activity' in target:
            print(f"   Near active site (<20Å): {(top_variants['dist_to_active_site'] < 20).sum()} ({100*(top_variants['dist_to_active_site'] < 20).sum()/top_k:.1f}%)")
        else:
            print(f"   Mean CAI: {top_variants['cai'].mean():.3f}")
            print(f"   Mean rare codon freq: {top_variants['rare_codon_freq'].mean():.3f}")
        
        high_confidence = (np.abs(top_variants['mutation_llr']) > 5).sum()
        print(f"   High confidence (|LLR|>5): {high_confidence} ({100*high_confidence/top_k:.1f}%)")
    
    @staticmethod
    def generate_predictions(df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
        """
        Generate multiple submission files with different strategies
        """
        print("\n" + "="*70)
        print("PHASE 4: NDCG-Optimized Predictions")
        print("="*70)
        
        submissions = {}
        
        # First, generate base predictions
        print("\n🎯 Generating base model predictions...")
        for target_name, model_config in Config.MODELS.items():
            df[f'{target_name}_base'] = Phase4NDCGOptimized.predict_base_model(df, model_config)
        
        # Strategy variations
        strategies = ['balanced', 'conservative', 'aggressive', 'diverse']
        
        for strategy in strategies:
            print(f"\n🔄 Creating {strategy} ensemble...")
            
            submission = df[['sequence']].copy()
            
            for target_name, model_config in Config.MODELS.items():
                # Create rankings
                rankings = Phase4NDCGOptimized.create_ranking_strategies(df, target_name)
                
                # Ensemble
                final_rank = Phase4NDCGOptimized.ensemble_rankings(rankings, strategy)
                
                # Convert to predictions
                predictions = Phase4NDCGOptimized.rank_to_prediction(
                    final_rank, 
                    model_config['output_range']
                )
                
                submission[target_name] = predictions
            
            submissions[strategy] = submission
            
            # Analyze top 10%
            top_k = int(len(df) * 0.1)
            for target_name in Config.MODELS.keys():
                Phase4NDCGOptimized.analyze_top_variants(submission, target_name, top_k)
        
        # Also create a "base_only" submission for comparison
        print(f"\n📋 Creating base model submission (no ranking optimization)...")
        base_submission = df[['sequence']].copy()
        for target_name in Config.MODELS.keys():
            base_submission[target_name] = df[f'{target_name}_base']
        submissions['base_only'] = base_submission
        
        print(f"\n✅ Generated {len(submissions)} submission variants")
        return submissions


# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description='PETase Prediction Pipeline - NDCG Optimized')
    parser.add_argument('--input', type=str, required=True, help='Input CSV file')
    parser.add_argument('--cds', type=str, required=True, help='CDS file')
    parser.add_argument('--output', type=str, default='submission', help='Output prefix (will create multiple files)')
    parser.add_argument('--token', type=str, default=None, help='ESM3 API token')
    parser.add_argument('--skip-phase1', action='store_true', help='Skip Phase 1')
    parser.add_argument('--strategy', type=str, default='all', 
                       choices=['all', 'balanced', 'conservative', 'aggressive', 'diverse', 'base_only'],
                       help='Which strategy to use (default: all)')
    
    args = parser.parse_args()
    
    print("\n" + "="*70)
    print(" "*10 + "PETase Prediction Pipeline - NDCG@10 Optimized")
    print("="*70)
    print(f"\n📁 Input: {args.input}")
    print(f"📁 CDS: {args.cds}")
    print(f"📁 Output prefix: {args.output}")
    print(f"🎯 Strategy: {args.strategy}")
    
    if args.token is None:
        print("\n🔑 ESM3 API token required")
        args.token = input("Enter your ESM3 token: ").strip()
    
    print("\n📂 Loading data...")
    test_df = pd.read_csv(args.input)
    cds_df = pd.read_csv(args.cds)
    print(f"✅ Loaded {len(test_df)} sequences")
    print(f"✅ Loaded {len(cds_df)} wildtype CDS")
    
    # Auto-detect column names
    print(f"\n🔍 Detecting column names...")
    print(f"   Test CSV columns: {test_df.columns.tolist()}")
    print(f"   CDS CSV columns: {cds_df.columns.tolist()}")
    
    # Find sequence column in test data
    test_seq_col = None
    for col in ['sequence', 'seq', 'protein_sequence', 'aa_sequence', 'amino_acid_sequence']:
        if col in test_df.columns:
            test_seq_col = col
            break
    
    if test_seq_col is None:
        raise ValueError(f"❌ Cannot find sequence column in test file. Available columns: {test_df.columns.tolist()}")
    
    if test_seq_col != 'sequence':
        print(f"   Renaming test '{test_seq_col}' → 'sequence'")
        test_df = test_df.rename(columns={test_seq_col: 'sequence'})
    
    # Find sequence column in CDS data
    cds_seq_col = None
    for col in ['sequence', 'seq', 'protein_sequence', 'aa_sequence', 'amino_acid_sequence']:
        if col in cds_df.columns:
            cds_seq_col = col
            break
    
    if cds_seq_col is None:
        raise ValueError(f"❌ Cannot find sequence column in CDS file. Available columns: {cds_df.columns.tolist()}")
    
    if cds_seq_col != 'sequence':
        print(f"   Renaming CDS '{cds_seq_col}' → 'sequence'")
        cds_df = cds_df.rename(columns={cds_seq_col: 'sequence'})
    
    # Find CDS column
    cds_col = None
    for col in ['cds', 'coding_sequence', 'nucleotide_sequence', 'dna_sequence', 'nt_sequence']:
        if col in cds_df.columns:
            cds_col = col
            break
    
    if cds_col is None:
        raise ValueError(f"❌ Cannot find CDS column in CDS file. Available columns: {cds_df.columns.tolist()}")
    
    if cds_col != 'cds':
        print(f"   Renaming CDS '{cds_col}' → 'cds'")
        cds_df = cds_df.rename(columns={cds_col: 'cds'})
    
    print(f"✅ Column names detected and standardized")
    
    wt_sequences = cds_df['sequence'].tolist()
    
    # Phase 1-3: Feature Engineering
    if not args.skip_phase1:
        phase1 = Phase1ESM3(args.token)
        test_df = phase1.process_sequences(test_df, wt_sequences)
        test_df.to_csv('phase1_complete.csv', index=False)
        print(f"💾 Saved phase1_complete.csv")
    else:
        print("\n⏭️  Skipping Phase 1")
        test_df = pd.read_csv('phase1_complete.csv')
    
    test_df = Phase2Structure.process_features(test_df)
    test_df.to_csv('phase2_complete.csv', index=False)
    print(f"💾 Saved phase2_complete.csv")
    
    test_df = Phase3Enhancement.process_features(test_df, cds_df)
    test_df.to_csv('phase3_complete.csv', index=False)
    print(f"💾 Saved phase3_complete.csv")
    
    # Phase 4: NDCG-Optimized Predictions
    submissions = Phase4NDCGOptimized.generate_predictions(test_df)
    
    # Save submissions
    print(f"\n💾 Saving submission files...")
    
    if args.strategy == 'all':
        for strategy_name, submission in submissions.items():
            filename = f"{args.output}_{strategy_name}.csv"
            submission.to_csv(filename, index=False)
            print(f"   ✅ {filename}")
    else:
        submission = submissions[args.strategy]
        filename = f"{args.output}.csv"
        submission.to_csv(filename, index=False)
        print(f"   ✅ {filename}")
    
    print(f"\n{'='*70}")
    print(f"✅ PIPELINE COMPLETE")
    print(f"{'='*70}")
    
    if args.strategy == 'all':
        print(f"\n📊 Generated {len(submissions)} submission variants:")
        print(f"   1. {args.output}_balanced.csv      (Recommended: balanced ensemble)")
        print(f"   2. {args.output}_conservative.csv  (Trust key features more)")
        print(f"   3. {args.output}_aggressive.csv    (Trust model more)")
        print(f"   4. {args.output}_diverse.csv       (Even weight distribution)")
        print(f"   5. {args.output}_base_only.csv     (Original model, no optimization)")
        print(f"\n💡 Recommendation: Submit 'balanced' first, then 'conservative' if needed")
    else:
        print(f"\n📊 Submission file: {args.output}.csv")
    
    print(f"\n🎉 Ready for competition submission!")


if __name__ == "__main__":
    main()
