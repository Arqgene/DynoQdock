import os
import requests
import subprocess
import urllib.parse

def fetch_smiles_from_pubchem(compound_name):
    """
    Fetch canonical SMILES from PubChem using compound name.
    
    Args:
        compound_name: Name of the compound (e.g., 'aspirin', 'caffeine')
    
    Returns:
        tuple: (smiles_string, compound_cid, error_message)
    """
    encoded_name = urllib.parse.quote(compound_name)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/property/CanonicalSMILES/TXT"
    
    try:
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            smiles = response.text.strip()
            
            cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/cids/TXT"
            cid_response = requests.get(cid_url, timeout=10)
            cid = cid_response.text.strip() if cid_response.status_code == 200 else "unknown"
            
            return smiles, cid, None
        elif response.status_code == 404:
            return None, None, f"Compound '{compound_name}' not found in chemical database"
        else:
            return None, None, f"Chemical database error: HTTP {response.status_code}"
    except Exception as e:
        return None, None, f"Failed to fetch compound data: {str(e)}"

def smiles_to_3d_sdf(smiles, output_sdf):
    """
    Convert SMILES to 3D SDF structure using OpenBabel.
    
    Args:
        smiles: SMILES string
        output_sdf: Path to output SDF file
    
    Returns:
        tuple: (success, error_message)
    """
    try:
        cmd = ['obabel', f'-:{smiles}', '-O', output_sdf, '--gen3d', '-h']
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        if result.returncode != 0:
            error_msg = result.stderr if result.stderr else result.stdout
            return False, f"3D structure generation failed: {error_msg}"
        
        if not os.path.exists(output_sdf) or os.path.getsize(output_sdf) == 0:
            return False, "3D structure generation produced empty file"
        
        return True, None
    
    except subprocess.TimeoutExpired:
        return False, "3D structure generation timed out"
    except Exception as e:
        return False, f"Error generating 3D structure: {str(e)}"

def convert_to_pdbqt(input_file, output_pdbqt, is_protein=False):
    """
    Convert molecular format to PDBQT using OpenBabel.
    
    Args:
        input_file: Path to input file (SDF, MOL, PDB, etc.)
        output_pdbqt: Path to output PDBQT file
        is_protein: Whether the molecule is a protein (default: False)
    
    Returns:
        tuple: (success, error_message)
    """
    try:
        if is_protein:
            cmd = ['obabel', input_file, '-O', output_pdbqt, '-xr']
        else:
            cmd = ['obabel', input_file, '-O', output_pdbqt, '-p', '7.4']
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        if result.returncode != 0:
            error_msg = result.stderr if result.stderr else result.stdout
            return False, f"PDBQT conversion failed: {error_msg}"
        
        if not os.path.exists(output_pdbqt) or os.path.getsize(output_pdbqt) == 0:
            return False, "PDBQT conversion produced empty file"
        
        return True, None
    
    except subprocess.TimeoutExpired:
        return False, "PDBQT conversion timed out"
    except Exception as e:
        return False, f"Error converting to PDBQT: {str(e)}"

def prepare_ligand_from_name(compound_name, output_pdbqt):
    """
    Complete ligand preparation pipeline from compound name:
    1. Fetch SMILES from PubChem
    2. Generate 3D structure (SDF)
    3. Convert to PDBQT
    
    Args:
        compound_name: Name of the compound
        output_pdbqt: Path to output PDBQT file
    
    Returns:
        tuple: (success, error_message, smiles, cid, sdf_path)
    """
    try:
        smiles, cid, error = fetch_smiles_from_pubchem(compound_name)
        if error:
            return False, error, None, None, None
        
        base_dir = os.path.dirname(output_pdbqt)
        sdf_path = os.path.join(base_dir, f'ligand_{cid}.sdf')
        
        success, error = smiles_to_3d_sdf(smiles, sdf_path)
        if not success:
            return False, error, smiles, cid, None
        
        success, error = convert_to_pdbqt(sdf_path, output_pdbqt, is_protein=False)
        if not success:
            return False, error, smiles, cid, sdf_path
        
        return True, None, smiles, cid, sdf_path
    
    except Exception as e:
        return False, f"Ligand preparation error: {str(e)}", None, None, None

def prepare_ligand_from_file(input_file, output_pdbqt):
    """
    Prepare ligand from uploaded file by converting to PDBQT.
    
    Args:
        input_file: Path to input file (SDF, MOL, MOL2, PDB, PDBQT)
        output_pdbqt: Path to output PDBQT file
    
    Returns:
        tuple: (success, error_message)
    """
    try:
        if input_file.endswith('.pdbqt'):
            import shutil
            shutil.copy(input_file, output_pdbqt)
            return True, None
        else:
            return convert_to_pdbqt(input_file, output_pdbqt, is_protein=False)
    except Exception as e:
        return False, f"Ligand file preparation error: {str(e)}"
