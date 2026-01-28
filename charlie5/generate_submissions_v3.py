"""
Phase 4 V3: Target-Specific Diverse Submissions
================================================

Key Improvements:

- Use different ranking strategies for the three targets (optimize Activity1, Activity2, and Expression separately).
- Adopt truly discriminative features (use pseudo-likelihood instead of CAI).
- Introduce physicochemical change features (e.g., hydrophobicity and volume).
- Optimize specifically for pH-dependent differences between Activity 1 and Activity 2.

Official Metrics:
- Activity 1: pH 5.5, 30°C, Citrate buffer
- Activity 2: pH 9.0, 30°C, Glycine buffer  
- Expression: E. coli BL21(DE3), 18°C, IPTG-induced, His-tag purification

Usage:
    python generate_submissions_v3.py

"""

import pandas as pd
import numpy as np
from scipy.stats import rankdata
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("Phase 4 V3: Target-Specific Diverse Submissions")
print("="*70)

# ============================================================================
# Amino Acid Properties (for computing mutation effects)
# ============================================================================

# Kyte-Doolittle hydrophobicity scale
HYDROPHOBICITY = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8,
    'G': -0.4, 'T': -0.7, 'W': -0.9, 'S': -0.8, 'Y': -1.3, 'P': -1.6,
    'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
}

# Amino acid molecular volume (Å³)
AA_VOLUME = {
    'G': 60, 'A': 88, 'S': 89, 'C': 108, 'D': 111, 'P': 112, 'N': 114,
    'T': 116, 'E': 138, 'V': 140, 'Q': 143, 'H': 153, 'M': 162, 'I': 166,
    'L': 166, 'K': 168, 'R': 173, 'F': 189, 'Y': 193, 'W': 227
}

# Amino acid flexibility (B-factor proxy)
AA_FLEXIBILITY = {
    'G': 1.0, 'S': 0.9, 'D': 0.8, 'N': 0.8, 'P': 0.7, 'K': 0.7, 'E': 0.6,
    'Q': 0.6, 'R': 0.5, 'T': 0.5, 'A': 0.4, 'H': 0.4, 'C': 0.3, 'M': 0.3,
    'V': 0.2, 'I': 0.2, 'L': 0.2, 'F': 0.1, 'Y': 0.1, 'W': 0.1
}

# ============================================================================
# Load Data
# ============================================================================

print("\n📂 Loading phase3_complete.csv...")
df = pd.read_csv('phase3_complete.csv')
print(f"✅ Loaded {len(df)} sequences")

# ============================================================================
# Compute Additional Mutation Features
# ============================================================================

print("\n🔬 Computing additional mutation features...")

def get_delta_hydrophobicity(row):
    """Compute change in hydrophobicity upon mutation"""
    if row['is_wildtype'] or pd.isna(row['wt_aa']) or pd.isna(row['mut_aa']):
        return 0.0
    wt_h = HYDROPHOBICITY.get(row['wt_aa'], 0)
    mut_h = HYDROPHOBICITY.get(row['mut_aa'], 0)
    return mut_h - wt_h

def get_delta_volume(row):
    """Compute change in side chain volume upon mutation"""
    if row['is_wildtype'] or pd.isna(row['wt_aa']) or pd.isna(row['mut_aa']):
        return 0.0
    wt_v = AA_VOLUME.get(row['wt_aa'], 140)
    mut_v = AA_VOLUME.get(row['mut_aa'], 140)
    return mut_v - wt_v

def get_delta_flexibility(row):
    """Compute change in local flexibility upon mutation"""
    if row['is_wildtype'] or pd.isna(row['wt_aa']) or pd.isna(row['mut_aa']):
        return 0.0
    wt_f = AA_FLEXIBILITY.get(row['wt_aa'], 0.5)
    mut_f = AA_FLEXIBILITY.get(row['mut_aa'], 0.5)
    return mut_f - wt_f

df['delta_hydrophobicity'] = df.apply(get_delta_hydrophobicity, axis=1)
df['delta_volume'] = df.apply(get_delta_volume, axis=1)
df['delta_flexibility'] = df.apply(get_delta_flexibility, axis=1)

# Absolute changes (for conservativeness)
df['abs_delta_hydro'] = df['delta_hydrophobicity'].abs()
df['abs_delta_volume'] = df['delta_volume'].abs()

print(f"   ✅ delta_hydrophobicity: {df['delta_hydrophobicity'].nunique()} unique values")
print(f"   ✅ delta_volume: {df['delta_volume'].nunique()} unique values")
print(f"   ✅ delta_flexibility: {df['delta_flexibility'].nunique()} unique values")

# ============================================================================
# Define Feature Dimensions with Good Discriminability
# ============================================================================

print("\n🎯 Creating high-discriminability ranking dimensions...")

# Dimension 1: Evolutionary Signal (ESM3 LLR)
# Higher LLR = mutation is evolutionarily acceptable = GOOD
evo_score = df['mutation_llr'].values
print(f"   ✅ Evolutionary (mutation_llr): {len(np.unique(evo_score))} unique values")

# Dimension 2: Protein Foldability (pseudo_likelihood)
# Higher PLL = sequence is more "protein-like" = GOOD for expression
fold_score = df['pseudo_likelihood'].values
print(f"   ✅ Foldability (pseudo_likelihood): {len(np.unique(fold_score))} unique values")

# Dimension 3: Structural Safety
# Lower risk + farther from active site = GOOD for activity
# Normalize distance to [0,1] range
dist_normalized = df['dist_to_active_site'] / df['dist_to_active_site'].max()
struct_score = -df['structure_risk_score'] + dist_normalized * 0.5
print(f"   ✅ Structural safety: {len(np.unique(struct_score))} unique values")

# Dimension 4: Physicochemical Conservation  
# Smaller changes in hydrophobicity and volume = more conservative = GOOD
# Use negative so that smaller absolute changes get higher scores
physico_score = -(df['abs_delta_hydro'] / 9.0 + df['abs_delta_volume'] / 167.0)  # Normalized
print(f"   ✅ Physicochemical conservation: {len(np.unique(physico_score))} unique values")

# Dimension 5: Active Site Proximity (for activity only)
# Being near active site is GOOD for activity (if mutation is favorable)
# This creates a risk-reward tradeoff
proximity_score = -df['dist_to_active_site']  # Closer = higher score
print(f"   ✅ Active site proximity: {len(np.unique(proximity_score))} unique values")

# Dimension 6: pH-specific charge effects
# For pH 5.5: H is partially charged (+0.3)
# For pH 9.0: K is partially uncharged, H is neutral
charge_score_pH55 = -df['charge_change_pH5.5'].abs()  # Less change = better
charge_score_pH90 = -df['charge_change_pH9.0'].abs()  # Less change = better
print(f"   ✅ Charge stability pH5.5: {len(np.unique(charge_score_pH55))} unique values")
print(f"   ✅ Charge stability pH9.0: {len(np.unique(charge_score_pH90))} unique values")

# ============================================================================
# Verify Independence
# ============================================================================

print("\n📊 Verifying dimension independence (Spearman correlations):")
dimensions = {
    'evolutionary': evo_score,
    'foldability': fold_score, 
    'structural': struct_score,
    'physicochemical': physico_score,
}

for i, (n1, s1) in enumerate(dimensions.items()):
    for n2, s2 in list(dimensions.items())[i+1:]:
        corr = pd.Series(s1).corr(pd.Series(s2), method='spearman')
        status = "✅" if abs(corr) < 0.5 else "⚠️"
        print(f"   {status} {n1} vs {n2}: {corr:.3f}")

# ============================================================================
# Target-Specific Ranking Functions
# ============================================================================

def rank_for_activity1(df, strategy='balanced'):
    """
    Activity 1: pH 5.5, 30°C, Citrate buffer
    - Enzyme needs to be stable and active at mildly acidic pH
    - Histidine residues are partially protonated
    """
    
    # Core: evolutionary signal is most important
    evo_rank = rankdata(-evo_score)
    struct_rank = rankdata(-struct_score)
    physico_rank = rankdata(-physico_score)
    charge_rank = rankdata(-charge_score_pH55)
    
    if strategy == 'evolutionary':
        # Trust ESM3 evolutionary signal most
        weights = {'evo': 0.55, 'struct': 0.25, 'physico': 0.15, 'charge': 0.05}
    elif strategy == 'structural':
        # Emphasize structural safety
        weights = {'evo': 0.30, 'struct': 0.45, 'physico': 0.15, 'charge': 0.10}
    elif strategy == 'conservative':
        # Prefer conservative mutations
        weights = {'evo': 0.30, 'struct': 0.25, 'physico': 0.35, 'charge': 0.10}
    elif strategy == 'balanced':
        weights = {'evo': 0.40, 'struct': 0.30, 'physico': 0.20, 'charge': 0.10}
    else:  # aggressive - trust model, accept more risk
        weights = {'evo': 0.50, 'struct': 0.20, 'physico': 0.20, 'charge': 0.10}
    
    combined = (weights['evo'] * evo_rank + 
                weights['struct'] * struct_rank + 
                weights['physico'] * physico_rank +
                weights['charge'] * charge_rank)
    
    return combined

def rank_for_activity2(df, strategy='balanced'):
    """
    Activity 2: pH 9.0, 30°C, Glycine buffer
    - Enzyme needs to be stable at alkaline pH
    - Lysine residues start to deprotonate
    - Histidine is neutral
    - Alkaline conditions can destabilize some proteins
    """
    
    evo_rank = rankdata(-evo_score)
    struct_rank = rankdata(-struct_score)
    physico_rank = rankdata(-physico_score)
    charge_rank = rankdata(-charge_score_pH90)
    fold_rank = rankdata(-fold_score)  # Foldability more important at harsh pH
    
    if strategy == 'evolutionary':
        weights = {'evo': 0.50, 'struct': 0.20, 'physico': 0.10, 'charge': 0.10, 'fold': 0.10}
    elif strategy == 'structural':
        weights = {'evo': 0.25, 'struct': 0.40, 'physico': 0.10, 'charge': 0.10, 'fold': 0.15}
    elif strategy == 'conservative':
        weights = {'evo': 0.25, 'struct': 0.20, 'physico': 0.30, 'charge': 0.10, 'fold': 0.15}
    elif strategy == 'balanced':
        weights = {'evo': 0.35, 'struct': 0.25, 'physico': 0.15, 'charge': 0.10, 'fold': 0.15}
    else:  # aggressive
        weights = {'evo': 0.45, 'struct': 0.20, 'physico': 0.15, 'charge': 0.10, 'fold': 0.10}
    
    combined = (weights['evo'] * evo_rank + 
                weights['struct'] * struct_rank + 
                weights['physico'] * physico_rank +
                weights['charge'] * charge_rank +
                weights['fold'] * fold_rank)
    
    return combined

def rank_for_expression(df, strategy='balanced'):
    """
    Expression: E. coli BL21(DE3), 18°C, IPTG-induced, His-tag purification
    - Low temperature (18°C) favors proper folding
    - Key factors: protein solubility, folding, no aggregation
    - pseudo_likelihood is excellent predictor of expressibility
    """
    
    fold_rank = rankdata(-fold_score)  # Most important for expression!
    evo_rank = rankdata(-evo_score)
    struct_rank = rankdata(-struct_score)
    physico_rank = rankdata(-physico_score)
    
    # CAI still matters for expression, even though it's per-parent
    cai_rank = rankdata(-df['cai'].values)
    
    if strategy == 'evolutionary':
        weights = {'fold': 0.35, 'evo': 0.35, 'struct': 0.15, 'physico': 0.10, 'cai': 0.05}
    elif strategy == 'structural':
        weights = {'fold': 0.30, 'evo': 0.20, 'struct': 0.30, 'physico': 0.10, 'cai': 0.10}
    elif strategy == 'conservative':
        weights = {'fold': 0.30, 'evo': 0.20, 'struct': 0.15, 'physico': 0.25, 'cai': 0.10}
    elif strategy == 'balanced':
        weights = {'fold': 0.35, 'evo': 0.25, 'struct': 0.20, 'physico': 0.10, 'cai': 0.10}
    else:  # foldability-focused
        weights = {'fold': 0.50, 'evo': 0.20, 'struct': 0.15, 'physico': 0.10, 'cai': 0.05}
    
    combined = (weights['fold'] * fold_rank + 
                weights['evo'] * evo_rank + 
                weights['struct'] * struct_rank + 
                weights['physico'] * physico_rank +
                weights['cai'] * cai_rank)
    
    return combined

# ============================================================================
# Utility Functions
# ============================================================================

def rank_to_prediction(ranks, output_range):
    """Convert ranks to predictions (lower rank = higher prediction)"""
    percentiles = rankdata(ranks) / len(ranks)
    predictions = (1 - percentiles) * (output_range[1] - output_range[0]) + output_range[0]
    return predictions

def add_noise_to_break_ties(ranks, noise_level=0.001):
    """Add tiny noise to break ties while preserving overall ranking"""
    noise = np.random.RandomState(42).randn(len(ranks)) * noise_level
    return ranks + noise

# ============================================================================
# Generate Submissions
# ============================================================================

OUTPUT_RANGES = {
    'activity_1': (0, 10),
    'activity_2': (0, 10),
    'expression': (0, 2)
}

STRATEGIES = {
    'evolutionary': 'Trust ESM3 evolutionary signal',
    'structural': 'Emphasize structural safety', 
    'conservative': 'Prefer conservative mutations',
    'balanced': 'Balanced approach',
    'aggressive': 'Accept more risk for potential gain'
}

submissions = {}

print("\n📝 Generating target-specific submissions...")

for strategy_name, description in STRATEGIES.items():
    print(f"\n   🔄 {strategy_name}: {description}")
    
    submission = df[['sequence']].copy()
    
    # Activity 1 (pH 5.5)
    ranks_a1 = rank_for_activity1(df, strategy_name)
    ranks_a1 = add_noise_to_break_ties(ranks_a1)
    submission['activity_1'] = rank_to_prediction(ranks_a1, OUTPUT_RANGES['activity_1'])
    
    # Activity 2 (pH 9.0) - DIFFERENT ranking!
    ranks_a2 = rank_for_activity2(df, strategy_name)
    ranks_a2 = add_noise_to_break_ties(ranks_a2)
    submission['activity_2'] = rank_to_prediction(ranks_a2, OUTPUT_RANGES['activity_2'])
    
    # Expression - DIFFERENT ranking!
    ranks_expr = rank_for_expression(df, strategy_name)
    ranks_expr = add_noise_to_break_ties(ranks_expr)
    submission['expression'] = rank_to_prediction(ranks_expr, OUTPUT_RANGES['expression'])
    
    submissions[strategy_name] = submission
    
    # Show top 5
    print(f"      Activity_1 top 5: {submission['activity_1'].nlargest(5).values.round(2)}")
    print(f"      Activity_2 top 5: {submission['activity_2'].nlargest(5).values.round(2)}")
    print(f"      Expression top 5: {submission['expression'].nlargest(5).values.round(2)}")

# ============================================================================
# Verify Diversity
# ============================================================================

print("\n" + "="*70)
print("📊 Diversity Verification")
print("="*70)

print("\n1. Activity_1 Ranking correlation (Spearman):")
strategy_names = list(STRATEGIES.keys())
for i, n1 in enumerate(strategy_names):
    for n2 in strategy_names[i+1:]:
        corr = submissions[n1]['activity_1'].corr(submissions[n2]['activity_1'], method='spearman')
        print(f"   {n1:15} vs {n2:15}: {corr:.3f}")

print("\n2. Activity_1 vs Activity_2 correlation (same strategy):")
for name in strategy_names:
    corr = submissions[name]['activity_1'].corr(submissions[name]['activity_2'], method='spearman')
    print(f"   {name:15}: activity_1 vs activity_2 = {corr:.3f}")

print("\n3. Top 100 Overlap rate (Activity_1):")
for i, n1 in enumerate(strategy_names):
    for n2 in strategy_names[i+1:]:
        top100_1 = set(submissions[n1].nlargest(100, 'activity_1').index)
        top100_2 = set(submissions[n2].nlargest(100, 'activity_1').index)
        overlap = len(top100_1 & top100_2)
        print(f"   {n1:15} vs {n2:15}: {overlap}%")

print("\n4. Activity_1 vs Activity_2 Top 100 Overlap rate (same strategy):")
for name in strategy_names:
    top100_a1 = set(submissions[name].nlargest(100, 'activity_1').index)
    top100_a2 = set(submissions[name].nlargest(100, 'activity_2').index)
    overlap = len(top100_a1 & top100_a2)
    print(f"   {name:15}: {overlap}%")

# ============================================================================
# Analyze Top Variants
# ============================================================================

print("\n" + "="*70)
print("🔬 Top Variant Analysis (balanced strategy)")
print("="*70)

balanced = submissions['balanced']
for target in ['activity_1', 'activity_2', 'expression']:
    top_idx = balanced.nlargest(100, target).index
    top_variants = df.loc[top_idx]
    
    print(f"\n   {target} Top 100:")
    print(f"      Mean LLR: {top_variants['mutation_llr'].mean():.2f}")
    print(f"      Mean PLL: {top_variants['pseudo_likelihood'].mean():.2f}")
    print(f"      Mean struct_risk: {top_variants['structure_risk_score'].mean():.3f}")
    print(f"      Mean dist_to_active: {top_variants['dist_to_active_site'].mean():.1f} Å")
    print(f"      Wildtype count: {top_variants['is_wildtype'].sum()}")

# ============================================================================
# Save Submissions
# ============================================================================

print("\n💾 Saving submission files...")

for strategy_name, submission in submissions.items():
    filename = f"submission_v3_{strategy_name}.csv"
    submission.to_csv(filename, index=False)
    print(f"   ✅ {filename}")

# Also save the best two as primary submissions
print("\n📋 Creating primary submission files...")
submissions['evolutionary'].to_csv('submission_primary_1.csv', index=False)
submissions['balanced'].to_csv('submission_primary_2.csv', index=False)
print("   ✅ submission_primary_1.csv (evolutionary strategy)")
print("   ✅ submission_primary_2.csv (balanced strategy)")

print("\n" + "="*70)
print("✅ V3 COMPLETE!")
print("="*70)

print("""
📋 Generated files:
   1. submission_v3_evolutionary.csv  - Trust ESM3 signal
   2. submission_v3_structural.csv    - Emphasize structural safety
   3. submission_v3_conservative.csv  - Prefer conservative mutations
   4. submission_v3_balanced.csv      - Balanced approach
   5. submission_v3_aggressive.csv    - Accept more risk
   
   Primary submissions:
   - submission_primary_1.csv (evolutionary)
   - submission_primary_2.csv (balanced)

💡 Key improvements in V3:
   1. Target-specific ranking (Activity1 ≠ Activity2 ≠ Expression)
   2. High-discriminability features (pseudo_likelihood, delta_hydrophobicity)
   3. pH-specific charge effects (pH 5.5 vs 9.0)
   4. Physicochemical conservation metrics
   
🎯 Recommended submission order:
   1st: submission_primary_1.csv (evolutionary - trust ESM3)
   2nd: submission_primary_2.csv (balanced - hedge bets)
""")
