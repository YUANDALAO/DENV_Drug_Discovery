#!/usr/bin/env python3
"""
LibInvent Vocabulary 诊断脚本
快速检查prior和agent模型的vocabulary状态
"""

import torch
from pathlib import Path
import sys

# 配置路径
REINVENT_HOME = Path("/home/august/aidd/REINVENT4")
PRIOR_PATH = REINVENT_HOME / "priors" / "libinvent.prior"
AGENT_PATH = REINVENT_HOME / "priors" / "denv_libinvent_model_v2.model"
SCAFFOLD_PATH = REINVENT_HOME / "data" / "pyrrolidine_dual_aryl.smi"

def load_model_vocabulary(model_path):
    """加载模型并提取vocabulary"""
    try:
        model = torch.load(model_path, map_location='cpu')
        return {
            'scaffold_vocab': model.get('scaffold_vocabulary'),
            'decorator_vocab': model.get('decorator_vocabulary'),
            'model': model
        }
    except Exception as e:
        print(f"❌ 无法加载模型 {model_path}: {e}")
        return None

def analyze_vocabulary(vocab, name):
    """分析vocabulary内容"""
    print(f"\n{'='*60}")
    print(f"📊 {name} Vocabulary Analysis")
    print('='*60)
    
    if vocab is None:
        print("⚠️  Vocabulary not found!")
        return False
    
    tokens = vocab.tokens()
    print(f"Total tokens: {len(tokens)}")
    print(f"Tokens: {tokens}")
    
    # 关键tokens检查
    critical_tokens = ['[*]', '*', '(', ')', 'C', 'N', 'O', 'S']
    print(f"\n🔍 Critical Tokens Check:")
    for token in critical_tokens:
        status = '✅' if token in tokens else '❌'
        print(f"  {status} {token}")
    
    # 特别关注[*]
    has_attachment = '[*]' in tokens
    if has_attachment:
        print(f"\n✅ Attachment point token [*] is present")
    else:
        print(f"\n❌ CRITICAL: Attachment point token [*] is MISSING!")
        print(f"   This will cause the error you're seeing!")
    
    return has_attachment

def check_scaffold_file(scaffold_path):
    """检查scaffold文件格式"""
    print(f"\n{'='*60}")
    print(f"📄 Scaffold File Analysis")
    print('='*60)
    
    if not scaffold_path.exists():
        print(f"❌ Scaffold file not found: {scaffold_path}")
        return False
    
    print(f"File: {scaffold_path}")
    
    with open(scaffold_path, 'r') as f:
        lines = f.readlines()
    
    print(f"Total scaffolds: {len(lines)}")
    print(f"\nFirst 5 scaffolds:")
    
    has_bracket_star = False
    for i, line in enumerate(lines[:5], 1):
        smi = line.strip()
        print(f"  {i}. {smi}")
        if '[*]' in smi:
            print(f"     ⚠️  Contains [*] - should use * instead")
            has_bracket_star = True
    
    if has_bracket_star:
        print(f"\n❌ PROBLEM: Scaffold file contains [*] tokens")
        print(f"   LibInvent expects * not [*] in scaffold SMILES")
        print(f"   This might cause issues during generation")
    else:
        print(f"\n✅ Scaffold format looks correct (using * not [*])")
    
    return not has_bracket_star

def main():
    print("="*60)
    print("🔬 REINVENT4 LibInvent Vocabulary Diagnostic Tool")
    print("="*60)
    
    # Check files exist
    print("\n📁 File Existence Check:")
    files_ok = True
    for path, name in [
        (PRIOR_PATH, "Prior"),
        (AGENT_PATH, "Agent"),
        (SCAFFOLD_PATH, "Scaffold"),
    ]:
        if path.exists():
            print(f"  ✅ {name}: {path}")
        else:
            print(f"  ❌ {name}: {path} NOT FOUND")
            files_ok = False
    
    if not files_ok:
        print("\n❌ Some required files are missing!")
        print("   Please check the paths in this script.")
        sys.exit(1)
    
    # Load and analyze vocabularies
    print("\n" + "="*60)
    print("Loading models...")
    print("="*60)
    
    prior_data = load_model_vocabulary(PRIOR_PATH)
    agent_data = load_model_vocabulary(AGENT_PATH)
    
    if not prior_data or not agent_data:
        print("\n❌ Failed to load models!")
        sys.exit(1)
    
    # Analyze Prior
    prior_scaffold_ok = analyze_vocabulary(
        prior_data['scaffold_vocab'], 
        "PRIOR - Scaffold"
    )
    prior_decorator_ok = analyze_vocabulary(
        prior_data['decorator_vocab'],
        "PRIOR - Decorator"
    )
    
    # Analyze Agent
    agent_scaffold_ok = analyze_vocabulary(
        agent_data['scaffold_vocab'],
        "AGENT - Scaffold"
    )
    agent_decorator_ok = analyze_vocabulary(
        agent_data['decorator_vocab'],
        "AGENT - Decorator"
    )
    
    # Check scaffold file
    scaffold_ok = check_scaffold_file(SCAFFOLD_PATH)
    
    # Summary and recommendations
    print(f"\n{'='*60}")
    print(f"📋 DIAGNOSIS SUMMARY")
    print('='*60)
    
    issues = []
    
    if not agent_decorator_ok:
        issues.append({
            'severity': 'CRITICAL',
            'issue': 'Agent decorator vocabulary missing [*] token',
            'solution': 'Retrain agent model or use libinvent.prior as agent'
        })
    
    if not agent_scaffold_ok:
        issues.append({
            'severity': 'CRITICAL',
            'issue': 'Agent scaffold vocabulary missing [*] token',
            'solution': 'Retrain agent model or use libinvent.prior as agent'
        })
    
    if not scaffold_ok:
        issues.append({
            'severity': 'WARNING',
            'issue': 'Scaffold file contains [*] instead of *',
            'solution': 'Replace [*] with * in scaffold SMILES'
        })
    
    if issues:
        print(f"\n🚨 ISSUES FOUND: {len(issues)}")
        for i, issue in enumerate(issues, 1):
            print(f"\n{i}. [{issue['severity']}] {issue['issue']}")
            print(f"   → Solution: {issue['solution']}")
        
        print(f"\n{'='*60}")
        print(f"💡 RECOMMENDED ACTION")
        print('='*60)
        
        if any(i['severity'] == 'CRITICAL' for i in issues):
            print("""
IMMEDIATE FIX (Quick workaround):
---------------------------------
Modify your TOML configuration to use the same prior for both:

[parameters]
prior_file = "priors/libinvent.prior"
agent_file = "priors/libinvent.prior"  # Changed from denv_libinvent_model_v2.model

This will allow you to continue running, but you'll lose the benefit
of transfer learning.

PROPER FIX (Best practice):
---------------------------
1. Prepare proper training data with correct SMILES format
2. Retrain transfer learning model ensuring vocabulary consistency
3. Verify new model has [*] token before using in RL

Commands:
# Check what tokens are needed in your data
grep -o '\[.*\]' data/chembl_denv.smi | sort -u

# Retrain with proper vocabulary
reinvent -l configs/transfer_learning.toml

# Verify new model
python diagnostic_script.py
""")
    else:
        print("\n✅ No critical issues found!")
        print("   Your vocabularies appear to be compatible.")
    
    print(f"\n{'='*60}")
    print("Diagnostic complete!")
    print('='*60)

if __name__ == "__main__":
    main()
