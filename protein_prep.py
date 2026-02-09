import os
import requests
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO, Select
import subprocess


def detect_file_format(file_path):
    """
    Detect if a file is PDB or mmCIF format.
    
    Returns:
        str: 'pdb', 'mmcif', or 'unknown'
    """
    try:
        with open(file_path, 'r') as f:
            first_lines = [f.readline() for _ in range(5)]
            content = ''.join(first_lines).lower()

            if content.startswith(
                    'data_') or '_entry.id' in content or 'loop_' in content:
                return 'mmcif'
            elif content.startswith('header') or content.startswith(
                    'atom') or content.startswith(
                        'hetatm') or 'remark' in content:
                return 'pdb'
            else:
                for line in first_lines:
                    if line.startswith('ATOM') or line.startswith(
                            'HETATM') or line.startswith('HEADER'):
                        return 'pdb'
                return 'unknown'
    except:
        return 'unknown'


class ProteinOnlySelect(Select):
    """Select only protein atoms (remove water, ligands, etc.)"""

    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0


class ChainASelect(Select):
    """Select only Chain A"""

    def accept_chain(self, chain):
        return chain.id == 'A'


def fetch_uniprot_fasta(uniprot_id):
    """
    Fetch FASTA sequence from UniProt.
    
    Args:
        uniprot_id: UniProt accession ID (e.g., 'P04637')
    
    Returns:
        tuple: (fasta_sequence, error_message)
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"

    try:
        response = requests.get(url, timeout=10)

        if response.status_code == 200:
            return response.text, None
        elif response.status_code == 404:
            return None, f"Molecular sequence ID '{uniprot_id}' not found"
        else:
            return None, f"Sequence database error: HTTP {response.status_code}"
    except Exception as e:
        return None, f"Failed to fetch molecular sequence: {str(e)}"


def search_uniprot_by_name(protein_name, require_alphafold=False):
    """
    Search protein database for a protein by name and get the first result's ID.
    Prioritizes reviewed entries and human proteins.
    """
    queries = [
        f"(protein_name:{protein_name}) AND (reviewed:true) AND (organism_id:9606)",
        f"(protein_name:{protein_name}) AND (reviewed:true)",
        f"(protein_name:{protein_name}) AND (organism_id:9606)",
        f"protein_name:{protein_name}"
    ]

    for query in queries:
        url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&format=json&size=5"

        try:
            response = requests.get(url, timeout=10)

            if response.status_code == 200:
                data = response.json()
                results = data.get('results', [])

                if results:
                    if require_alphafold:
                        for entry in results:
                            uniprot_id = entry.get('primaryAccession')
                            protein_full_name = entry.get(
                                'proteinDescription',
                                {}).get('recommendedName',
                                        {}).get('fullName',
                                                {}).get('value', protein_name)

                            af_check_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
                            af_response = requests.head(af_check_url,
                                                        timeout=5)

                            if af_response.status_code == 200:
                                return uniprot_id, protein_full_name, None

                    entry = results[0]
                    uniprot_id = entry.get('primaryAccession')
                    protein_full_name = entry.get(
                        'proteinDescription',
                        {}).get('recommendedName',
                                {}).get('fullName',
                                        {}).get('value', protein_name)
                    return uniprot_id, protein_full_name, None
        except Exception:
            continue

    return None, None, f"No results found for protein name '{protein_name}'"


def predict_structure_esmfold(fasta_sequence, output_path):
    """
    Predict protein structure using high-speed sequence-to-structure model.
    """
    sequence_lines = fasta_sequence.strip().split('\n')
    sequence = ''.join(
        [line for line in sequence_lines if not line.startswith('>')])

    if len(sequence) > 400:
        return False, f"Sequence too long ({len(sequence)} residues). Maximum 400 amino acids supported."

    if len(sequence) < 10:
        return False, "Sequence too short. Minimum 10 amino acids required."

    api_url = "https://api.esmatlas.com/foldSequence/v1/pdb/"

    try:
        # ESMFold API expects raw sequence string without newlines or headers
        clean_sequence = "".join(sequence.split())
        response = requests.post(api_url, data=clean_sequence, timeout=120)

        if response.status_code == 200:
            pdb_content = response.text

            if len(pdb_content) < 100:
                return False, "Database returned invalid structure data"

            with open(output_path, 'w') as f:
                f.write(pdb_content)

            return True, None
        elif response.status_code == 400:
            return False, "Invalid sequence format"
        elif response.status_code == 503:
            return False, "Structural prediction service temporarily unavailable. Please try again."
        else:
            return False, f"Structural database error: HTTP {response.status_code}"
    except requests.Timeout:
        return False, "Structure prediction timed out. Try a shorter sequence."
    except Exception as e:
        return False, f"Failed to predict structure: {str(e)}"


def fetch_alphafold_structure(uniprot_id, output_path):
    """
    Download predicted structure from structural database.
    """
    pdb_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"

    try:
        response = requests.get(pdb_url, timeout=30)

        if response.status_code == 200:
            with open(output_path, 'wb') as f:
                f.write(response.content)
            return True, None
        elif response.status_code == 404:
            return False, f"No pre-computed structure available for {uniprot_id}"
        else:
            return False, f"Structural database error: HTTP {response.status_code}"
    except Exception as e:
        return False, f"Failed to download structure: {str(e)}"


def clean_protein_structure_text_based(input_pdb,
                                       output_pdb,
                                       keep_chain='A',
                                       remove_water=True,
                                       remove_hetero=True):
    """
    Text-based PDB cleaner as fallback when BioPython fails.
    Handles malformed PDB files that BioPython cannot parse.
    
    Args:
        input_pdb: Path to input PDB file
        output_pdb: Path to output cleaned PDB file
        keep_chain: Chain to keep (default: 'A', set to None to keep all chains)
        remove_water: Remove water molecules (default: True)
        remove_hetero: Remove all heteroatoms including ligands (default: True)
    
    Returns:
        tuple: (success, error_message)
    """
    standard_residues = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
        'SEC', 'PYL', 'ASX', 'GLX', 'UNK'
    }
    water_residues = {'HOH', 'WAT', 'H2O', 'TIP', 'TIP3', 'SOL'}

    try:
        cleaned_lines = []
        atom_count = 0

        with open(input_pdb, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    if len(line) < 26:
                        continue

                    chain_id = line[21:22].strip() if len(line) > 21 else ''
                    res_name = line[17:20].strip() if len(line) > 19 else ''

                    if keep_chain and chain_id and chain_id != keep_chain:
                        continue

                    if remove_water and res_name in water_residues:
                        continue

                    if remove_hetero:
                        if line.startswith('HETATM'):
                            continue
                        if res_name not in standard_residues:
                            continue

                    cleaned_lines.append(line)
                    atom_count += 1

                elif line.startswith('TER') or line.startswith('END'):
                    cleaned_lines.append(line)
                elif line.startswith('MODEL') or line.startswith('ENDMDL'):
                    cleaned_lines.append(line)
                elif line.startswith('HEADER') or line.startswith(
                        'TITLE') or line.startswith('COMPND'):
                    cleaned_lines.append(line)

        if atom_count == 0:
            return False, "No valid protein atoms found after cleaning"

        with open(output_pdb, 'w') as f:
            f.writelines(cleaned_lines)
            if not any(line.startswith('END') for line in cleaned_lines):
                f.write('END\n')

        return True, None

    except Exception as e:
        return False, f"Text-based cleaning failed: {str(e)}"


def clean_protein_structure(input_pdb,
                            output_pdb,
                            keep_chain='A',
                            remove_water=True,
                            remove_hetero=True):
    """
    Clean protein structure by removing water molecules, heteroatoms, and unwanted chains.
    Supports both PDB and mmCIF formats.
    Uses BioPython first, falls back to text-based parser for problematic files.
    
    Args:
        input_pdb: Path to input PDB file (can be PDB or mmCIF format)
        output_pdb: Path to output cleaned PDB file
        keep_chain: Chain to keep (default: 'A', set to None to keep all chains)
        remove_water: Remove water molecules (default: True)
        remove_hetero: Remove all heteroatoms including ligands (default: True)
    
    Returns:
        tuple: (success, error_message)
    """
    try:
        file_format = detect_file_format(input_pdb)

        if file_format == 'mmcif':
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(PERMISSIVE=True, QUIET=True)

        structure = parser.get_structure("protein", input_pdb)

        io = PDBIO()
        io.set_structure(structure)

        if keep_chain and remove_hetero:

            class ChainAndProteinSelect(Select):

                def accept_chain(self, chain):
                    return chain.id == keep_chain

                def accept_residue(self, residue):
                    return 1 if residue.id[0] == " " else 0

            io.save(output_pdb, ChainAndProteinSelect())
        elif keep_chain:

            class ChainSelect(Select):

                def accept_chain(self, chain):
                    return chain.id == keep_chain

                def accept_residue(self, residue):
                    if remove_water:
                        return 1 if residue.id[0] != "W" else 0
                    return 1

            io.save(output_pdb, ChainSelect())
        elif remove_hetero:
            io.save(output_pdb, ProteinOnlySelect())
        else:

            class NoWaterSelect(Select):

                def accept_residue(self, residue):
                    return 1 if residue.id[0] != "W" else 0

            io.save(output_pdb, NoWaterSelect() if remove_water else Select())

        if not os.path.exists(output_pdb) or os.path.getsize(output_pdb) == 0:
            return False, "Cleaned PDB file is empty or not created"

        return True, None

    except Exception as e:
        success, error = clean_protein_structure_text_based(
            input_pdb, output_pdb, keep_chain, remove_water, remove_hetero)
        if success:
            return True, None
        return False, f"BioPython failed ({str(e)}), text-based fallback also failed: {error}"


def add_hydrogens_openbabel(input_pdb, output_pdb, ph=7.0):
    """
    Add hydrogens to protein structure using OpenBabel.
    
    Args:
        input_pdb: Path to input PDB file
        output_pdb: Path to output PDB file with hydrogens
        ph: pH value for protonation (default: 7.0)
    
    Returns:
        tuple: (success, error_message)
    """
    try:
        cmd = ['obabel', input_pdb, '-O', output_pdb, '-p', str(ph)]
        result = subprocess.run(cmd,
                                capture_output=True,
                                text=True,
                                timeout=60)

        if result.returncode != 0:
            error_msg = result.stderr if result.stderr else result.stdout
            return False, f"OpenBabel hydrogen addition failed: {error_msg}"

        if not os.path.exists(output_pdb) or os.path.getsize(output_pdb) == 0:
            return False, "Hydrogen addition produced empty file"

        return True, None

    except subprocess.TimeoutExpired:
        return False, "Hydrogen addition timed out"
    except Exception as e:
        return False, f"Error adding hydrogens: {str(e)}"


def prepare_protein(input_pdb, output_pdbqt, keep_chain='A', add_h=True):
    """
    Complete protein preparation pipeline:
    1. Clean structure (remove water, heteroatoms, unwanted chains)
    2. Add hydrogens (optional)
    3. Convert to PDBQT format
    
    Args:
        input_pdb: Path to input PDB file
        output_pdbqt: Path to output PDBQT file
        keep_chain: Chain to keep (default: 'A')
        add_h: Add hydrogens before conversion (default: True)
    
    Returns:
        tuple: (success, error_message, intermediate_pdb_path)
    """
    try:
        base_dir = os.path.dirname(output_pdbqt)
        cleaned_pdb = os.path.join(base_dir, 'cleaned_protein.pdb')

        success, error = clean_protein_structure(input_pdb,
                                                 cleaned_pdb,
                                                 keep_chain=keep_chain)
        if not success:
            return False, error, None

        if add_h:
            hydrogenated_pdb = os.path.join(base_dir, 'protein_with_h.pdb')
            success, error = add_hydrogens_openbabel(cleaned_pdb,
                                                     hydrogenated_pdb)
            if not success:
                print(
                    f"Warning: Hydrogen addition failed: {error}. Proceeding without hydrogens."
                )
                final_pdb = cleaned_pdb
            else:
                final_pdb = hydrogenated_pdb
        else:
            final_pdb = cleaned_pdb

        cmd = ['obabel', final_pdb, '-O', output_pdbqt, '-xr']
        result = subprocess.run(cmd,
                                capture_output=True,
                                text=True,
                                timeout=60)

        if result.returncode != 0:
            error_msg = result.stderr if result.stderr else result.stdout
            return False, f"PDBQT conversion failed: {error_msg}", final_pdb

        if not os.path.exists(output_pdbqt) or os.path.getsize(
                output_pdbqt) == 0:
            return False, "PDBQT conversion produced empty file", final_pdb

        return True, None, final_pdb

    except Exception as e:
        return False, f"Protein preparation error: {str(e)}", None
