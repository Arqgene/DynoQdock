import os
import re
from typing import Dict, List, Tuple, Optional

def verify_pdb_structure(pdb_file: str) -> Dict:
    """
    Verify PDB file structure and extract statistics.
    
    Args:
        pdb_file: Path to PDB file
    
    Returns:
        Dictionary with verification results and statistics
    """
    if not os.path.exists(pdb_file):
        return {
            'valid': False,
            'error': 'File not found',
            'statistics': {}
        }
    
    try:
        file_size = os.path.getsize(pdb_file)
        
        if file_size == 0:
            return {
                'valid': False,
                'error': 'File is empty',
                'statistics': {'file_size_kb': 0}
            }
        
        atom_count = 0
        chains = set()
        residues = set()
        has_coordinates = False
        warnings = []
        
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atom_count += 1
                    
                    if len(line) >= 54:
                        try:
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            has_coordinates = True
                        except ValueError:
                            pass
                    
                    if len(line) >= 22:
                        chain = line[21:22].strip()
                        if chain:
                            chains.add(chain)
                    
                    if len(line) >= 26:
                        res_num = line[22:26].strip()
                        res_name = line[17:20].strip()
                        if res_num and res_name:
                            residues.add((res_name, res_num))
        
        if atom_count == 0:
            return {
                'valid': False,
                'error': 'No atoms found in PDB file',
                'statistics': {'file_size_kb': round(file_size / 1024, 2)}
            }
        
        if not has_coordinates:
            warnings.append('No valid 3D coordinates found')
        
        if atom_count < 100:
            warnings.append(f'Very few atoms ({atom_count}) - structure may be incomplete')
        
        if len(chains) == 0:
            warnings.append('No chain identifiers found')
        
        return {
            'valid': True,
            'error': None,
            'statistics': {
                'file_size_kb': round(file_size / 1024, 2),
                'atom_count': atom_count,
                'chain_count': len(chains),
                'chains': sorted(list(chains)),
                'residue_count': len(residues),
                'has_coordinates': has_coordinates
            },
            'warnings': warnings
        }
    
    except Exception as e:
        return {
            'valid': False,
            'error': f'Error reading PDB file: {str(e)}',
            'statistics': {}
        }

def verify_pdbqt_structure(pdbqt_file: str, is_protein: bool = True) -> Dict:
    """
    Verify PDBQT file structure and extract statistics.
    
    Args:
        pdbqt_file: Path to PDBQT file
        is_protein: True for protein, False for ligand
    
    Returns:
        Dictionary with verification results and statistics
    """
    if not os.path.exists(pdbqt_file):
        return {
            'valid': False,
            'error': 'File not found',
            'statistics': {}
        }
    
    try:
        file_size = os.path.getsize(pdbqt_file)
        
        if file_size == 0:
            return {
                'valid': False,
                'error': 'File is empty',
                'statistics': {'file_size_kb': 0}
            }
        
        atom_count = 0
        has_charges = False
        has_atom_types = False
        has_coordinates = False
        root_found = False
        warnings = []
        
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atom_count += 1
                    
                    if len(line) >= 54:
                        try:
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            has_coordinates = True
                        except ValueError:
                            pass
                    
                    if len(line) >= 70:
                        try:
                            charge = line[70:76].strip()
                            if charge:
                                has_charges = True
                        except:
                            pass
                    
                    if len(line) >= 79:
                        atom_type = line[77:79].strip()
                        if atom_type:
                            has_atom_types = True
                
                if line.startswith('ROOT'):
                    root_found = True
        
        if atom_count == 0:
            return {
                'valid': False,
                'error': 'No atoms found in PDBQT file',
                'statistics': {'file_size_kb': round(file_size / 1024, 2)}
            }
        
        if not has_coordinates:
            warnings.append('No valid 3D coordinates found')
        
        if not has_charges:
            warnings.append('No partial charges found - may affect docking accuracy')
        
        if not has_atom_types:
            warnings.append('No atom types found in PDBQT format')
        
        if is_protein:
            if atom_count < 100:
                warnings.append(f'Very few atoms ({atom_count}) for a protein')
            if atom_count > 50000:
                warnings.append(f'Very large protein ({atom_count} atoms) - docking may be slow')
        else:
            if not root_found:
                warnings.append('No ROOT found - ligand may not be properly formatted')
            if atom_count < 5:
                warnings.append(f'Very small ligand ({atom_count} atoms)')
            if atom_count > 150:
                warnings.append(f'Very large ligand ({atom_count} atoms) - consider fragmentation')
        
        return {
            'valid': True,
            'error': None,
            'statistics': {
                'file_size_kb': round(file_size / 1024, 2),
                'atom_count': atom_count,
                'has_coordinates': has_coordinates,
                'has_partial_charges': has_charges,
                'has_atom_types': has_atom_types,
                'has_root': root_found
            },
            'warnings': warnings
        }
    
    except Exception as e:
        return {
            'valid': False,
            'error': f'Error reading PDBQT file: {str(e)}',
            'statistics': {}
        }

def verify_protein_preparation(pdb_file: Optional[str], pdbqt_file: str) -> Dict:
    """
    Comprehensive verification of protein preparation pipeline.
    
    Args:
        pdb_file: Path to intermediate PDB file (cleaned/hydrogenated), can be None
        pdbqt_file: Path to final PDBQT file
    
    Returns:
        Dictionary with comprehensive verification results
    """
    results = {
        'pdb': None,
        'pdbqt': None,
        'overall_valid': False,
        'summary': ''
    }
    
    if pdb_file and os.path.exists(pdb_file):
        results['pdb'] = verify_pdb_structure(pdb_file)
    
    results['pdbqt'] = verify_pdbqt_structure(pdbqt_file, is_protein=True)
    
    if results['pdbqt']['valid']:
        results['overall_valid'] = True
        stats = results['pdbqt']['statistics']
        warnings = results['pdbqt']['warnings']
        
        summary_parts = [
            f"✅ Protein prepared successfully",
            f"📊 {stats['atom_count']} atoms",
            f"💾 {stats['file_size_kb']} KB"
        ]
        
        if stats.get('has_partial_charges'):
            summary_parts.append("⚡ Charges assigned")
        
        if warnings:
            summary_parts.append(f"⚠️ {len(warnings)} warning(s)")
        
        results['summary'] = " | ".join(summary_parts)
    else:
        results['summary'] = f"❌ Verification failed: {results['pdbqt']['error']}"
    
    return results

def verify_ligand_preparation(pdbqt_file: str, original_file: Optional[str] = None) -> Dict:
    """
    Comprehensive verification of ligand preparation pipeline.
    
    Args:
        pdbqt_file: Path to final PDBQT file
        original_file: Path to original file (SDF, MOL, etc.), can be None
    
    Returns:
        Dictionary with comprehensive verification results
    """
    results = {
        'original': None,
        'pdbqt': None,
        'overall_valid': False,
        'summary': ''
    }
    
    results['pdbqt'] = verify_pdbqt_structure(pdbqt_file, is_protein=False)
    
    if results['pdbqt']['valid']:
        results['overall_valid'] = True
        stats = results['pdbqt']['statistics']
        warnings = results['pdbqt']['warnings']
        
        summary_parts = [
            f"✅ Ligand prepared successfully",
            f"📊 {stats['atom_count']} atoms",
            f"💾 {stats['file_size_kb']} KB"
        ]
        
        if stats.get('has_partial_charges'):
            summary_parts.append("⚡ Charges assigned")
        
        if stats.get('has_root'):
            summary_parts.append("🎯 Rotatable bonds defined")
        
        if warnings:
            summary_parts.append(f"⚠️ {len(warnings)} warning(s)")
        
        results['summary'] = " | ".join(summary_parts)
    else:
        results['summary'] = f"❌ Verification failed: {results['pdbqt']['error']}"
    
    return results

def estimate_molecular_weight(pdbqt_file: str) -> Optional[float]:
    """
    Estimate molecular weight from PDBQT file based on atom types.
    
    Args:
        pdbqt_file: Path to PDBQT file
    
    Returns:
        Estimated molecular weight in Daltons, or None if unable to estimate
    """
    atomic_weights = {
        'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
        'F': 18.998, 'P': 30.974, 'S': 32.065, 'Cl': 35.453,
        'Br': 79.904, 'I': 126.904
    }
    
    try:
        total_weight = 0.0
        
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    if len(line) >= 79:
                        atom_type = line[77:79].strip()
                        
                        element = atom_type[0].upper()
                        
                        if element in atomic_weights:
                            total_weight += atomic_weights[element]
        
        return round(total_weight, 2) if total_weight > 0 else None
    
    except Exception as e:
        return None
