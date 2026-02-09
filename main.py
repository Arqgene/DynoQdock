import os
import subprocess
import re
import platform
from flask import Flask, request, jsonify, send_from_directory, redirect, url_for, session
from flask_cors import CORS
from werkzeug.utils import secure_filename
from werkzeug.security import generate_password_hash, check_password_hash
import shutil
import json
import psycopg2
from psycopg2.extras import RealDictCursor
import protein_prep
import ligand_prep
import verify_structures

app = Flask(__name__, static_folder='static')
CORS(app)
app.secret_key = os.environ.get('SESSION_SECRET', 'arqgene-docking-secret-2026')
app.config['UPLOAD_FOLDER'] = 'data'
app.config['POSES_FOLDER'] = 'data/poses'
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024
app.config['DATABASE_URL'] = os.environ.get('DATABASE_URL')

ALLOWED_EXTENSIONS = {'pdb', 'pdbqt', 'sdf', 'mol', 'mol2'}

def get_db_connection():
    conn = psycopg2.connect(app.config['DATABASE_URL'])
    return conn

def get_smina_command():
    """Get the correct Smina executable"""
    path = os.path.join(os.path.dirname(__file__), 'smina.exe')
    if os.path.exists(path):
        return path
    return 'smina'

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def convert_to_pdbqt(input_file, output_file, is_protein=False):
    """Convert any molecular format to PDBQT using OpenBabel"""
    try:
        if is_protein:
            cmd = ['obabel', input_file, '-O', output_file, '-xr']
        else:
            cmd = ['obabel', input_file, '-O', output_file, '-p', '7.4']
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            error_msg = result.stderr if result.stderr else result.stdout
            print(f"OpenBabel conversion failed for {input_file}: {error_msg}")
            raise Exception(f"Conversion failed: {error_msg}")
        
        if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
            print(f"OpenBabel produced empty or missing output file: {output_file}")
            raise Exception("Conversion produced empty file - please check input format")
        
        return True
    except subprocess.SubprocessError as e:
        print(f"Subprocess error during conversion: {e}")
        raise Exception(f"Conversion tool error: {str(e)}")
    except Exception as e:
        print(f"Conversion error: {e}")
        raise

def convert_pdbqt_to_pdb(input_file, output_file):
    """Convert PDBQT to PDB for visualization"""
    try:
        cmd = ['obabel', input_file, '-O', output_file]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            error_msg = result.stderr if result.stderr else result.stdout
            print(f"OpenBabel PDB conversion failed: {error_msg}")
            return False
        
        if not os.path.exists(output_file):
            print(f"OpenBabel produced missing output file: {output_file}")
            return False
        
        return True
    except Exception as e:
        print(f"PDB conversion error: {e}")
        return False

def parse_vina_results(output_file):
    affinities = []
    try:
        with open(output_file, 'r') as f:
            for line in f:
                if 'minimizedAffinity' in line:
                    val = line.split()[-1]
                    affinities.append(float(val))
                elif 'VINA RESULT:' in line:
                    val = line.split()[3]
                    affinities.append(float(val))
        return affinities[:9]
    except Exception as e:
        print(f"Error parsing results: {e}")
        return []

def split_poses(multi_pose_file, output_dir):
    """Split multi-pose PDBQT file into individual pose files"""
    poses = []
    current_pose = []
    pose_num = 1
    
    try:
        with open(multi_pose_file, 'r') as f:
            for line in f:
                if line.startswith('MODEL'):
                    current_pose = [line]
                elif line.startswith('ENDMDL'):
                    current_pose.append(line)
                    pose_file = os.path.join(output_dir, f'pose_{pose_num}.pdbqt')
                    with open(pose_file, 'w') as pf:
                        pf.writelines(current_pose)
                    poses.append(pose_file)
                    pose_num += 1
                    current_pose = []
                elif current_pose:
                    current_pose.append(line)
    except Exception as e:
        print(f"Error splitting poses: {e}")
    
    return poses

def combine_protein_ligand(protein_file, ligand_file, output_file):
    """Combine protein and ligand PDBQT files into a single file for visualization"""
    try:
        with open(output_file, 'w') as out:
            with open(protein_file, 'r') as prot:
                for line in prot:
                    if not line.startswith('END'):
                        out.write(line)
            
            out.write('TER\n')
            
            with open(ligand_file, 'r') as lig:
                for line in lig:
                    if not line.startswith(('MODEL', 'ENDMDL', 'END')):
                        out.write(line)
            
            out.write('END\n')
        
        return True
    except Exception as e:
        print(f"Error combining protein and ligand: {e}")
        return False

@app.before_request
def check_auth():
    # List of endpoints that don't require authentication
    public_endpoints = ['login', 'login_page', 'signup', 'logout', 'static', 'serve_pose', 'serve_data']
    
    # Check if the current endpoint is public or if user is logged in
    if request.endpoint and request.endpoint not in public_endpoints:
        if 'user_id' not in session:
            # List of all API/Action endpoints that should return 401
            api_endpoints = [
                '/api/', '/prepare_protein', '/prepare_ligand', '/dock', 
                '/get_results', '/upload_batch', '/get_fasta', '/predict_structure',
                '/dock_batch'
            ]
            is_api = any(request.path.startswith(p) for p in api_endpoints) or request.path in api_endpoints
            
            if is_api:
                return jsonify({'error': 'Authentication required. Please log in.'}), 401
            return redirect(url_for('login_page'))

@app.route('/')
def index():
    return send_from_directory('static', 'index.html')

@app.route('/login')
def login_page():
    return send_from_directory('static', 'login.html')

@app.route('/api/auth/signup', methods=['POST'])
def signup():
    data = request.get_json()
    email = data.get('email')
    password = data.get('password')
    name = data.get('name')
    institution = data.get('institution')
    
    if not email or not password:
        return jsonify({'error': 'Missing email or password'}), 400
    
    hashed_password = generate_password_hash(password)
    
    try:
        conn = get_db_connection()
        cur = conn.cursor()
        cur.execute('INSERT INTO users (email, password, name, institution) VALUES (%s, %s, %s, %s) RETURNING id', 
                   (email, hashed_password, name, institution))
        user_id = cur.fetchone()[0]
        conn.commit()
        cur.close()
        conn.close()
        session['user_id'] = user_id
        return jsonify({'success': True, 'message': 'Account created successfully'})
    except psycopg2.IntegrityError:
        return jsonify({'error': 'Email already exists'}), 400
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/auth/login', methods=['POST'])
def login():
    data = request.get_json()
    email = data.get('email')
    password = data.get('password')
    
    try:
        conn = get_db_connection()
        cur = conn.cursor(cursor_factory=RealDictCursor)
        cur.execute('SELECT * FROM users WHERE email = %s', (email,))
        user = cur.fetchone()
        cur.close()
        conn.close()
        
        if user and check_password_hash(user['password'], password):
            session['user_id'] = user['id']
            return jsonify({'success': True, 'message': 'Logged in successfully'})
        else:
            return jsonify({'error': 'Invalid email or password'}), 401
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/auth/logout')
def logout():
    session.pop('user_id', None)
    return redirect(url_for('login_page'))

@app.route('/upload_batch', methods=['POST'])
def upload_batch():
    protein_ids = request.form.get('protein_ids', '').split(',')
    protein_names = request.form.get('protein_names', '').split(',')
    ligand_names = request.form.get('ligand_names', '').split(',')
    
    os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
    os.makedirs(app.config['POSES_FOLDER'], exist_ok=True)
    
    protein_paths = []
    ligand_paths = []
    
    if 'proteins' in request.files:
        for protein in request.files.getlist('proteins'):
            if protein and allowed_file(protein.filename):
                filename = secure_filename(protein.filename)
                path = os.path.join(app.config['UPLOAD_FOLDER'], f"batch_prot_{filename}")
                protein.save(path)
                
                pdbqt_path = path + ".pdbqt"
                success, error, cleaned_pdb = protein_prep.prepare_protein(path, pdbqt_path)
                if success:
                    verify_structures.verify_protein_preparation(cleaned_pdb, pdbqt_path)
                    protein_paths.append(pdbqt_path)
                
    if 'ligands' in request.files:
        for ligand in request.files.getlist('ligands'):
            if ligand and allowed_file(ligand.filename):
                filename = secure_filename(ligand.filename)
                path = os.path.join(app.config['UPLOAD_FOLDER'], f"batch_lig_{filename}")
                ligand.save(path)
                
                pdbqt_path = path + ".pdbqt"
                if ligand_prep.prepare_ligand_from_file(path, pdbqt_path)[0]:
                    ligand_paths.append(pdbqt_path)

    for pid in [p.strip() for p in protein_ids if p.strip()]:
        raw_pdb = os.path.join(app.config['UPLOAD_FOLDER'], f"batch_raw_{pid}.pdb")
        pdbqt = os.path.join(app.config['UPLOAD_FOLDER'], f"batch_prot_{pid}.pdbqt")
        if not protein_prep.fetch_alphafold_structure(pid, raw_pdb)[0]:
            fasta, _ = protein_prep.fetch_uniprot_fasta(pid)
            if fasta:
                protein_prep.predict_structure_esmfold(fasta, raw_pdb)
        
        if os.path.exists(raw_pdb):
            success, _, cleaned = protein_prep.prepare_protein(raw_pdb, pdbqt)
            if success:
                verify_structures.verify_protein_preparation(cleaned, pdbqt)
                protein_paths.append(pdbqt)
            
    for pname in [p.strip() for p in protein_names if p.strip()]:
        pid, _, _ = protein_prep.search_uniprot_by_name(pname, require_alphafold=False)
        if pid:
            raw_pdb = os.path.join(app.config['UPLOAD_FOLDER'], f"batch_raw_{pid}.pdb")
            pdbqt = os.path.join(app.config['UPLOAD_FOLDER'], f"batch_prot_{pid}.pdbqt")
            if not protein_prep.fetch_alphafold_structure(pid, raw_pdb)[0]:
                fasta, _ = protein_prep.fetch_uniprot_fasta(pid)
                if fasta:
                    protein_prep.predict_structure_esmfold(fasta, raw_pdb)
            
            if os.path.exists(raw_pdb):
                success, _, cleaned = protein_prep.prepare_protein(raw_pdb, pdbqt)
                if success:
                    verify_structures.verify_protein_preparation(cleaned, pdbqt)
                    protein_paths.append(pdbqt)
                
    for lname in [l.strip() for l in ligand_names if l.strip()]:
        pdbqt = os.path.join(app.config['UPLOAD_FOLDER'], f"batch_lig_{secure_filename(lname)}.pdbqt")
        if ligand_prep.prepare_ligand_from_name(lname, pdbqt)[0]:
            verify_structures.verify_ligand_preparation(pdbqt)
            ligand_paths.append(pdbqt)
            
    return jsonify({
        'message': f'Prepared {len(protein_paths)} proteins and {len(ligand_paths)} ligands',
        'proteins': [os.path.basename(p) for p in protein_paths],
        'ligands': [os.path.basename(p) for p in ligand_paths]
    })

@app.route('/dock_batch', methods=['POST'])
def dock_batch():
    if 'user_id' not in session:
        return jsonify({'error': 'Authentication required'}), 401
    
    data = request.get_json()
    proteins = data.get('proteins', [])
    ligands = data.get('ligands', [])
    
    if not proteins or not ligands:
        return jsonify({'error': 'No proteins or ligands specified for batch docking'}), 400
    
    results = []
    smina_cmd = get_smina_command()
    
    for prot_file in proteins:
        prot_path = os.path.join(app.config['UPLOAD_FOLDER'], prot_file)
        if not os.path.exists(prot_path): continue
            
        for lig_file in ligands:
            lig_path = os.path.join(app.config['UPLOAD_FOLDER'], lig_file)
            if not os.path.exists(lig_path): continue
                
            prot_name = prot_file.replace('batch_prot_', '').replace('.pdbqt', '')
            lig_name = lig_file.replace('batch_lig_', '').replace('.pdbqt', '')
            
            output_file = os.path.join(app.config['UPLOAD_FOLDER'], f'batch_{prot_name}_{lig_name}_out.pdbqt')
            
            # Blind docking for batch (default 30A box at center 0,0,0)
            cmd = [
                smina_cmd, '--receptor', prot_path, '--ligand', lig_path,
                '--num_modes', '9', '--exhaustiveness', '1', 
                '--center_x', '0', '--center_y', '0', '--center_z', '0',
                '--size_x', '30', '--size_y', '30', '--size_z', '30',
                '--out', output_file, '--verbosity', '0'
            ]
            
            try:
                print(f"Running docking for {prot_name} and {lig_name}...")
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
                
                if result.returncode != 0:
                    print(f"Smina failed for {prot_name}/{lig_name}: {result.stderr}")
                    continue

                affinities = parse_vina_results(output_file)
                if affinities:
                    # Verification step for batch (compliance)
                    verify_structures.verify_ligand_preparation(lig_path)
                    
                    complex_pdb = f'batch_{prot_name}_{lig_name}_complex.pdb'
                    complex_pdb_path = os.path.join(app.config['UPLOAD_FOLDER'], complex_pdb)
                    
                    combined_pdbqt = output_file + ".complex.pdbqt"
                    if combine_protein_ligand(prot_path, output_file, combined_pdbqt):
                        if convert_pdbqt_to_pdb(combined_pdbqt, complex_pdb_path):
                            results.append({
                                'protein': prot_name,
                                'ligand': lig_name,
                                'best_affinity': affinities[0],
                                'affinities': affinities,
                                'complex_file': complex_pdb
                            })
                        else:
                            print(f"Failed to convert complex to PDB for {prot_name}/{lig_name}")
                    else:
                        print(f"Failed to combine protein/ligand for {prot_name}/{lig_name}")
                else:
                    print(f"No affinities found for {prot_name}/{lig_name}")
            except subprocess.TimeoutExpired:
                print(f"Docking timed out for {prot_name}/{lig_name}")
            except Exception as e:
                print(f"Batch docking error for {prot_name}/{lig_name}: {e}")
                
    return jsonify({'results': results})

@app.route('/get_fasta', methods=['POST'])
def get_fasta():
    if 'user_id' not in session:
        return jsonify({'error': 'Authentication required'}), 401
    
    uniprot_id = request.form.get('uniprot_id', '').strip()
    protein_name = request.form.get('protein_name', '').strip()
    
    if not uniprot_id and protein_name:
        uniprot_id, _, error = protein_prep.search_uniprot_by_name(protein_name)
        if error:
            return jsonify({'error': error}), 404
            
    if not uniprot_id:
        return jsonify({'error': 'Protein ID or Name required'}), 400
        
    fasta, error = protein_prep.fetch_uniprot_fasta(uniprot_id)
    if error:
        return jsonify({'error': error}), 404
        
    return jsonify({
        'success': True,
        'fasta': fasta,
        'uniprot_id': uniprot_id,
        'protein_name': protein_name
    })

@app.route('/predict_structure', methods=['POST'])
def predict_structure():
    if 'user_id' not in session:
        return jsonify({'error': 'Authentication required'}), 401
    
    fasta = request.form.get('fasta', '').strip()
    uniprot_id = request.form.get('uniprot_id', '').strip()
    protein_name = request.form.get('protein_name', '').strip()
    
    os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
    output_pdb = os.path.join(app.config['UPLOAD_FOLDER'], 'protein.pdb')
    protein_pdbqt = os.path.join(app.config['UPLOAD_FOLDER'], 'protein.pdbqt')
    
    if not fasta:
        if not uniprot_id and protein_name:
            uniprot_id, _, _ = protein_prep.search_uniprot_by_name(protein_name)
        
        if uniprot_id:
            fasta, error = protein_prep.fetch_uniprot_fasta(uniprot_id)
            if error:
                return jsonify({'error': error}), 404
        else:
            return jsonify({'error': 'Sequence or Protein ID required for prediction'}), 400
            
    success, error = protein_prep.predict_structure_esmfold(fasta, output_pdb)
    if not success:
        return jsonify({'error': f'Prediction failed: {error}'}), 500
        
    # After prediction, prepare it for docking
    success, error, cleaned_pdb = protein_prep.prepare_protein(output_pdb, protein_pdbqt)
    if success:
        verification = verify_structures.verify_protein_preparation(cleaned_pdb, protein_pdbqt)
        return jsonify({
            'success': True, 
            'message': 'Structure predicted and prepared successfully',
            'verification': verification
        })
    else:
        return jsonify({'error': f'Preparation of predicted structure failed: {error}'}), 500

@app.route('/prepare_protein', methods=['POST'])
def prepare_protein():
    os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
    os.makedirs(app.config['POSES_FOLDER'], exist_ok=True)
    protein_pdbqt = os.path.join(app.config['UPLOAD_FOLDER'], 'protein.pdbqt')
    
    if 'file' in request.files and request.files['file'].filename:
        protein_file = request.files['file']
        if not protein_file.filename or not allowed_file(protein_file.filename):
            return jsonify({'error': 'Invalid file format'}), 400
        filename = secure_filename(protein_file.filename)
        input_pdb = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        protein_file.save(input_pdb)
        success, error, cleaned_pdb = protein_prep.prepare_protein(input_pdb, protein_pdbqt)
        if success:
            verification = verify_structures.verify_protein_preparation(cleaned_pdb, protein_pdbqt)
            return jsonify({'success': True, 'message': 'Protein structure cleaned and prepared successfully', 'verification': verification})
        else:
            return jsonify({'error': f'Protein preparation failed: {error}'}), 500
    
    elif 'uniprot_id' in request.form and request.form['uniprot_id']:
        uniprot_id = request.form['uniprot_id'].strip()
        raw_pdb = os.path.join(app.config['UPLOAD_FOLDER'], f'raw_{uniprot_id}.pdb')
        success, error = protein_prep.fetch_alphafold_structure(uniprot_id, raw_pdb)
        if not success:
            fasta, fasta_error = protein_prep.fetch_uniprot_fasta(uniprot_id)
            if fasta:
                success, error = protein_prep.predict_structure_esmfold(fasta, raw_pdb)
            else:
                return jsonify({'error': f'Structure retrieval failed: {fasta_error}'}), 404
        if not success:
            return jsonify({'error': f'Structure retrieval failed: {error}'}), 404
        success, error, cleaned_pdb = protein_prep.prepare_protein(raw_pdb, protein_pdbqt)
        if success:
            verification = verify_structures.verify_protein_preparation(cleaned_pdb, protein_pdbqt)
            return jsonify({'success': True, 'message': 'Structure retrieved and prepared successfully', 'uniprot_id': uniprot_id, 'verification': verification})
        else:
            return jsonify({'error': f'Preparation failed: {error}'}), 500
    
    elif 'protein_name' in request.form and request.form['protein_name']:
        protein_name = request.form['protein_name'].strip()
        uniprot_id, full_name, error = protein_prep.search_uniprot_by_name(protein_name, require_alphafold=False)
        if error:
            return jsonify({'error': f'Search failed: {error}'}), 404
        raw_pdb = os.path.join(app.config['UPLOAD_FOLDER'], f'raw_{uniprot_id}.pdb')
        success, error = protein_prep.fetch_alphafold_structure(uniprot_id, raw_pdb)
        if not success:
            fasta, fasta_error = protein_prep.fetch_uniprot_fasta(uniprot_id)
            if fasta:
                success, error = protein_prep.predict_structure_esmfold(fasta, raw_pdb)
            else:
                return jsonify({'error': f'Structure retrieval failed: {fasta_error}'}), 404
        if not success:
            return jsonify({'error': f'Structure retrieval failed: {error}'}), 404
        success, error, cleaned_pdb = protein_prep.prepare_protein(raw_pdb, protein_pdbqt)
        if success:
            verification = verify_structures.verify_protein_preparation(cleaned_pdb, protein_pdbqt)
            return jsonify({'success': True, 'message': f'Structure for "{full_name}" retrieved and prepared successfully', 'uniprot_id': uniprot_id, 'verification': verification})
        else:
            return jsonify({'error': f'Protein preparation failed: {error}'}), 500
    else:
        return jsonify({'error': 'Please provide either a file, ID, or name'}), 400

@app.route('/prepare_ligand', methods=['POST'])
def prepare_ligand():
    os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
    ligand_pdbqt = os.path.join(app.config['UPLOAD_FOLDER'], 'ligand.pdbqt')
    
    if 'file' in request.files and request.files['file'].filename:
        ligand_file = request.files['file']
        if not ligand_file.filename or not allowed_file(ligand_file.filename):
            return jsonify({'error': 'Invalid file format'}), 400
        filename = secure_filename(ligand_file.filename)
        input_file = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        ligand_file.save(input_file)
        success, error = ligand_prep.prepare_ligand_from_file(input_file, ligand_pdbqt)
        if success:
            verification = verify_structures.verify_ligand_preparation(ligand_pdbqt, input_file)
            mol_weight = verify_structures.estimate_molecular_weight(ligand_pdbqt)
            return jsonify({'success': True, 'message': 'Ligand file prepared successfully', 'verification': verification, 'molecular_weight': mol_weight})
        else:
            return jsonify({'error': f'Ligand preparation failed: {error}'}), 500
    
    elif 'compound_name' in request.form and request.form['compound_name']:
        compound_name = request.form['compound_name'].strip()
        success, error, smiles, cid, sdf_path = ligand_prep.prepare_ligand_from_name(compound_name, ligand_pdbqt)
        if success:
            verification = verify_structures.verify_ligand_preparation(ligand_pdbqt, sdf_path)
            mol_weight = verify_structures.estimate_molecular_weight(ligand_pdbqt)
            return jsonify({'success': True, 'message': f'Ligand "{compound_name}" generated successfully', 'verification': verification, 'molecular_weight': mol_weight})
        else:
            return jsonify({'error': f'Ligand generation failed: {error}'}), 500
    else:
        return jsonify({'error': 'Please provide either a file or compound name'}), 400

@app.route('/dock', methods=['POST'])
def run_docking():
    if 'user_id' not in session:
        return jsonify({'error': 'Authentication required'}), 401
    protein_pdbqt = os.path.join(app.config['UPLOAD_FOLDER'], 'protein.pdbqt')
    ligand_pdbqt = os.path.join(app.config['UPLOAD_FOLDER'], 'ligand.pdbqt')
    if not os.path.exists(protein_pdbqt) or not os.path.exists(ligand_pdbqt):
        return jsonify({'error': 'Please upload files first'}), 400
    for f in os.listdir(app.config['POSES_FOLDER']):
        os.remove(os.path.join(app.config['POSES_FOLDER'], f))
    output_file = os.path.join(app.config['UPLOAD_FOLDER'], 'all_poses.pdbqt')
    smina_cmd = get_smina_command()
    cmd = [smina_cmd, '--receptor', protein_pdbqt, '--ligand', ligand_pdbqt, '--num_modes', '9', '--exhaustiveness', '8', '--out', output_file, '--verbosity', '1']
    data = request.get_json() if request.is_json else {}
    grid_mode = data.get('grid_mode', 'manual')
    if grid_mode == 'manual':
        cmd.extend(['--center_x', str(data.get('center_x', 0)), '--center_y', str(data.get('center_y', 0)), '--center_z', str(data.get('center_z', 0)), '--size_x', str(data.get('size_x', 25)), '--size_y', str(data.get('size_y', 25)), '--size_z', str(data.get('size_z', 25))])
    else:
        cmd.extend(['--center_x', '0', '--center_y', '0', '--center_z', '0', '--size_x', '30', '--size_y', '30', '--size_z', '30'])
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            return jsonify({'error': f'Docking failed: {result.stderr}'}), 500
        affinities = parse_vina_results(output_file)
        if not affinities:
            return jsonify({'error': 'No docking results found'}), 500
        pose_files = split_poses(output_file, app.config['POSES_FOLDER'])
        results = []
        for i, (pose_file, affinity) in enumerate(zip(pose_files, affinities), 1):
            complex_pdbqt = os.path.join(app.config['POSES_FOLDER'], f'complex_{i}.pdbqt')
            complex_pdb = os.path.join(app.config['POSES_FOLDER'], f'complex_{i}.pdb')
            if combine_protein_ligand(protein_pdbqt, pose_file, complex_pdbqt):
                if convert_pdbqt_to_pdb(complex_pdbqt, complex_pdb):
                    results.append({'pose': i, 'affinity': affinity, 'path': f'data/poses/complex_{i}.pdb'})
        with open(os.path.join(app.config['UPLOAD_FOLDER'], 'results.json'), 'w') as f:
            json.dump(results, f)
        return jsonify({'results': results})
    except subprocess.TimeoutExpired:
        return jsonify({'error': 'Docking timeout'}), 500
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/results', methods=['GET'])
def get_results():
    results_file = os.path.join(app.config['UPLOAD_FOLDER'], 'results.json')
    if not os.path.exists(results_file):
        return jsonify({'error': 'No results available'}), 404
    with open(results_file, 'r') as f:
        results = json.load(f)
    return jsonify({'results': results})

@app.route('/data/poses/<filename>')
def serve_pose(filename):
    return send_from_directory(app.config['POSES_FOLDER'], filename)

@app.route('/data/<filename>')
def serve_data(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename)

if __name__ == '__main__':
    os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
    os.makedirs(app.config['POSES_FOLDER'], exist_ok=True)
    app.run(host='0.0.0.0', port=5000, debug=False)
