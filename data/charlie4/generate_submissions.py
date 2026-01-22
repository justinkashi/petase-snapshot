"""
Phase 4 Standalone: Generate Submissions from Phase 3 Data
===========================================================

This script reads phase3_complete.csv and generates all 5 submission strategies.
No need to re-run Phase 1-3!

Usage:
    python generate_submissions.py

Output:
    - submission_balanced.csv       (Recommended for 1st submission)
    - submission_conservative.csv   (Recommended for 2nd submission)
    - submission_aggressive.csv
    - submission_diverse.csv
    - submission_base_only.csv
"""

import pandas as pd
import numpy as np
from scipy.stats import rankdata
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("Phase 4: NDCG-Optimized Predictions (Standalone)")
print("="*70)

# ============================================================================
# Configuration
# ============================================================================

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
# Phase 4: Prediction Functions
# ============================================================================

def predict_base_model(df, model_config):
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

def compute_confidence(df, target):
    """Compute prediction confidence"""
    confidence = np.zeros(len(df))
    
    # Strong evolutionary signal
    strong_llr = np.abs(df['mutation_llr']) > 5
    confidence += strong_llr * 0.3
    
    # Low structural risk
    low_risk = df['structure_risk_score'] < 0.3
    confidence += low_risk * 0.2
    
    # Features agreement
    if 'activity' in target:
        good_llr = df['mutation_llr'] > 0
        low_charge = np.abs(df['charge_change_pH5.5']) < 1.0
        feature_agree = good_llr & low_risk & low_charge
    else:
        good_pll = df['pseudo_likelihood'] > -5.5
        good_cai = df['cai'] > 0.33
        feature_agree = good_pll & good_cai & low_risk
    
    confidence += feature_agree * 0.5
    
    return confidence

def create_ranking_strategies(df, target):
    """Create multiple ranking strategies"""
    rankings = {}
    
    # Strategy 1: Base model
    base_pred = df[f'{target}_base']
    rankings['base_model'] = rankdata(-base_pred)
    
    # Strategy 2: Key feature only
    if 'activity' in target:
        key_feature = df['mutation_llr']
    else:
        key_feature = df['pseudo_likelihood']
    rankings['key_feature'] = rankdata(-key_feature)
    
    # Strategy 3: Confidence-weighted
    confidence = compute_confidence(df, target)
    weighted_pred = base_pred * (1 + confidence)
    rankings['confidence_weighted'] = rankdata(-weighted_pred)
    
    # Strategy 4: Median rank
    if 'activity' in target:
        individual_ranks = [
            rankdata(-df['mutation_llr']),
            rankdata(df['structure_risk_score']),
            rankdata(df['dist_to_active_site'])
        ]
    else:
        individual_ranks = [
            rankdata(-df['pseudo_likelihood']),
            rankdata(-df['cai']),
            rankdata(df['rare_codon_freq'])
        ]
    
    median_rank = np.median(individual_ranks, axis=0)
    rankings['median_rank'] = median_rank
    
    # Strategy 5: Biology-constrained
    bio_score = base_pred.copy()
    high_risk = df['structure_risk_score'] > 0.6
    bio_score[high_risk] *= 0.5
    active_site_bad = (df['dist_to_active_site'] < 15) & (df['mutation_llr'] < 0)
    bio_score[active_site_bad] *= 0.3
    rankings['bio_constrained'] = rankdata(-bio_score)
    
    return rankings

def ensemble_rankings(rankings, strategy='balanced'):
    """Ensemble rankings with different strategies"""
    if strategy == 'balanced':
        weights = {
            'base_model': 0.30,
            'key_feature': 0.25,
            'confidence_weighted': 0.20,
            'median_rank': 0.15,
            'bio_constrained': 0.10
        }
    elif strategy == 'conservative':
        weights = {
            'base_model': 0.20,
            'key_feature': 0.40,
            'confidence_weighted': 0.15,
            'median_rank': 0.15,
            'bio_constrained': 0.10
        }
    elif strategy == 'aggressive':
        weights = {
            'base_model': 0.50,
            'key_feature': 0.15,
            'confidence_weighted': 0.20,
            'median_rank': 0.10,
            'bio_constrained': 0.05
        }
    else:  # 'diverse'
        weights = {
            'base_model': 0.20,
            'key_feature': 0.20,
            'confidence_weighted': 0.20,
            'median_rank': 0.20,
            'bio_constrained': 0.20
        }
    
    final_rank = sum(weights[name] * rank for name, rank in rankings.items())
    return final_rank

def rank_to_prediction(ranks, output_range):
    """Convert ranks to predictions"""
    percentiles = rankdata(ranks) / len(ranks)
    predictions = (1 - percentiles) * (output_range[1] - output_range[0]) + output_range[0]
    return predictions

def analyze_top_variants(submission_df, features_df, target, top_k=500):
    """Analyze top 10% variants"""
    # Get indices of top variants
    top_indices = submission_df.nlargest(top_k, target).index
    # Get corresponding rows from features_df
    top_variants = features_df.loc[top_indices]
    
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

# ============================================================================
# Main Processing
# ============================================================================

print("\n📂 Loading phase3_complete.csv...")
df = pd.read_csv('phase3_complete.csv')
print(f"✅ Loaded {len(df)} sequences with {len(df.columns)} features")

# Generate base predictions
print("\n🎯 Generating base model predictions...")
for target_name, model_config in MODELS.items():
    df[f'{target_name}_base'] = predict_base_model(df, model_config)
    print(f"   ✅ {target_name}: mean={df[f'{target_name}_base'].mean():.2f}, range=[{df[f'{target_name}_base'].min():.2f}, {df[f'{target_name}_base'].max():.2f}]")

# Generate submissions for each strategy
strategies = {
    'balanced': 'Balanced ensemble (Recommended for 1st submission)',
    'conservative': 'Trust LLR more (Recommended for 2nd submission)',
    'aggressive': 'Trust model more',
    'diverse': 'Equal weights',
}

submissions = {}

for strategy, description in strategies.items():
    print(f"\n🔄 Creating {strategy} ensemble ({description})...")
    
    submission = df[['sequence']].copy()
    
    for target_name, model_config in MODELS.items():
        # Create rankings
        rankings = create_ranking_strategies(df, target_name)
        
        # Ensemble
        final_rank = ensemble_rankings(rankings, strategy)
        
        # Convert to predictions
        predictions = rank_to_prediction(final_rank, model_config['output_range'])
        
        submission[target_name] = predictions
    
    submissions[strategy] = submission
    
    # Analyze top 10%
    top_k = int(len(df) * 0.1)
    for target_name in MODELS.keys():
        analyze_top_variants(submission, df, target_name, top_k)

# Base only submission
print(f"\n📋 Creating base model submission (no ranking optimization)...")
base_submission = df[['sequence']].copy()
for target_name in MODELS.keys():
    base_submission[target_name] = df[f'{target_name}_base']
submissions['base_only'] = base_submission

print(f"\n✅ Generated {len(submissions)} submission variants")

# ============================================================================
# Save Submissions
# ============================================================================

print(f"\n💾 Saving submission files...")

for strategy_name, submission in submissions.items():
    filename = f"submission_{strategy_name}.csv"
    submission.to_csv(filename, index=False)
    print(f"   ✅ {filename} ({len(submission)} sequences)")

print("\n" + "="*70)
print("✅ PHASE 4 COMPLETE!")
print("="*70)

print(f"\n📊 Generated submission files:")
print(f"   1. submission_balanced.csv       (Recommended for 1st submission)")
print(f"   2. submission_conservative.csv   (Recommended for 2nd submission)")
print(f"   3. submission_aggressive.csv")
print(f"   4. submission_diverse.csv")
print(f"   5. submission_base_only.csv      (Control group)")

print(f"\n💡 Recommendation: Submit 'balanced' first, then 'conservative' if needed")
print(f"\n🎉 Ready for competition submission!")
