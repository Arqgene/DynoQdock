# 🧬 Protein-Ligand Docking & 3D Visualization App

A web application for molecular docking using **Smina** and interactive 3D visualization with **MolStar**. Upload protein and ligand files, run molecular docking, and explore the top 9 binding poses with their binding affinities.

## 🚀 Features

- **Automated Structure Retrieval**: Fetch proteins from AlphaFold Database and ligands from PubChem
- **AI-Powered Prediction**: ESMFold integration for instant protein structure prediction
- **Structure Verification**: Automated quality checks with detailed statistics (atom count, charges, molecular weight)
- **Manual Grid Box Control**: Targeted docking with custom active site coordinates
- **Multi-format Support**: Upload proteins and ligands in PDB, PDBQT, SDF, MOL, or MOL2 formats
- **Automatic Format Conversion**: Files are automatically converted to PDBQT format using OpenBabel
- **Smina Docking Engine**: Generate top 9 binding poses with binding affinity scores (blind or targeted)
- **Interactive 3D Visualization**: View molecular structures using the MolStar viewer
- **Real-time Feedback**: Live status logging with verification warnings and docking progress
- **Clean Interface**: Professional white background design with step-by-step workflow

## 🛠️ Tech Stack

### Backend
- **Python 3.11** - Core runtime
- **Flask** - Web framework
- **Smina** - Molecular docking engine (static binary)
- **OpenBabel** - Chemical format conversion

### Frontend
- **MolStar** (via CDN) - 3D molecular visualization
- **Vanilla JavaScript** - Interactive UI
- **CSS** - Clean, responsive styling

## 📁 Project Structure

```
.
├── main.py              # Flask backend with docking logic
├── protein_prep.py      # Protein preparation & AlphaFold/ESMFold integration
├── ligand_prep.py       # Ligand preparation & PubChem integration
├── verify_structures.py # Structure verification module (NEW)
├── static/
│   ├── index.html       # Main UI with grid box controls
│   ├── style.css        # Styling
│   └── viewer.js        # MolStar integration & verification display
├── data/
│   ├── poses/           # Generated docking poses
│   └── ...              # Uploaded files
├── smina.static         # Smina binary
├── WINDOWS_SETUP.md     # Detailed Windows setup guide
└── README.md            # This file
```

## 🚀 Quick Start

### Running on Replit

1. Click the **Run** button
2. Flask server starts automatically on port 5000
3. Access the web interface in the preview window

### Local Installation

```bash
# Install dependencies
pip install flask
sudo apt install openbabel

# Download Smina
curl -L https://sourceforge.net/projects/smina/files/smina.static/download -o smina.static
chmod +x smina.static

# Run the app
python main.py
```

## 📖 Usage

### Step 1: Upload Files
1. Click **"Choose File"** for protein (PDB or PDBQT format)
2. Click **"Choose File"** for ligand (PDB, PDBQT, SDF, MOL, or MOL2)
3. Click **"Upload Files"**

### Step 2: Run Docking
1. Click **"Run Smina Docking"**
2. Wait 1-2 minutes for docking to complete
3. Results will appear automatically

### Step 3: Explore Results
- View the **top 9 poses** in the results table
- See **binding affinity** for each pose (lower = stronger binding)
- Click **"View 3D"** to visualize any pose in the molecular viewer
- Pose #1 loads automatically (best binding affinity)

## ⚙️ Docking Parameters

The app uses optimized Smina parameters:
- **Number of modes**: 9 (top 9 poses)
- **Exhaustiveness**: 8 (balance between speed and accuracy)
- **Search space**: Auto-calculated from ligand position (±8 Å)
- **pH**: 7.4 (physiological pH for protonation)

## 📊 Understanding Results

### Binding Affinity (kcal/mol)
- **Negative values** = Favorable binding (more negative = stronger)
- **Typical range**: -15 to 0 kcal/mol
- **Strong binding**: < -9 kcal/mol
- **Moderate binding**: -6 to -9 kcal/mol
- **Weak binding**: > -6 kcal/mol

### Pose Ranking
Poses are ranked by binding affinity (Pose #1 = strongest predicted binding)

## 🔧 Technical Details

### File Conversion Pipeline
1. **Upload** → Any supported format
2. **Convert** → PDBQT format (if needed) using OpenBabel
3. **Dock** → Smina generates 9 poses in PDBQT
4. **Visualize** → Convert poses to PDB for MolStar

### API Endpoints
- `POST /upload` - Upload and convert protein/ligand files
- `POST /dock` - Run Smina docking
- `GET /results` - Retrieve docking results
- `GET /data/poses/<filename>` - Serve pose files for visualization

## 🧪 Example Files

You can test the app with sample protein-ligand pairs:
- **HIV-1 Protease** + **Indinavir** (drug molecule)
- **SARS-CoV-2 Main Protease** + **Inhibitor**
- Any protein from **PDB database** + small molecule ligand

Download structures from:
- [Protein Data Bank (PDB)](https://www.rcsb.org/)
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/) (for ligands)

## ⚠️ Limitations

- Maximum file size: 50 MB
- Docking timeout: 5 minutes
- Designed for small molecules (ligands < 50 heavy atoms work best)
- Development server (not for production use)

## 🔬 Scientific References

- **Smina**: Fork of AutoDock Vina with enhanced scoring
- **OpenBabel**: Chemical toolbox for format conversion
- **MolStar**: Modern molecular visualization toolkit

## 📝 License

This project uses open-source tools:
- Smina (Apache License)
- OpenBabel (GPL)
- MolStar (MIT License)

---

Built with ❤️ for molecular biology and drug discovery research
