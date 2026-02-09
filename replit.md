# Automated Protein-Ligand Docking Platform

## Overview

This project is an advanced web-based molecular docking application designed for researchers and students in computational chemistry and drug discovery. It automates the entire protein-ligand docking workflow, from structure retrieval and preparation to docking execution and 3D visualization. Users can either upload structure files or simply input protein and compound names, with the system automatically fetching data from external databases (AlphaFold DB, UniProt, PubChem), preparing molecules, performing docking with Smina, and visualizing results. The platform aims to provide a user-friendly, automated interface, eliminating the need for command-line expertise or manual structure preparation.

**Key Capabilities:**
- **Automated Structure Handling**: Fetches AlphaFold protein structures and generates 3D ligand structures from PubChem SMILES.
- **Automated Structure Verification**: Performs quality checks on prepared protein and ligand structures.
- **Flexible Docking Modes**: Supports both blind docking and targeted docking with manual grid box control for active site-specific analysis.
- **Smart Docking Workflow**: Integrates preparation, verification, docking, and visualization into a seamless pipeline.
- **Interactive Results**: Provides live status logging, binding affinity analysis, and an interactive 3D pose viewer.
- **ESMFold Integration**: Offers rapid protein structure prediction via ESMFold for proteins not yet in AlphaFold DB.

## User Preferences

Preferred communication style: Simple, everyday language.

## System Architecture

### Application Framework
- **Web Framework**: Flask serves as a lightweight Python web server running on port 5000.
- **Frontend**: A static HTML/CSS/JavaScript single-page application, without complex build processes.
- **Data Storage**: All transient data, including uploads and docking results, are stored in a local filesystem under the `data/` directory.

### Molecular Processing Pipeline
The system orchestrates a multi-step molecular processing pipeline:
1.  **Structure Retrieval/Upload**: Users can upload PDB/PDBQT/SDF/MOL/MOL2 files or provide protein/ligand names for automated fetching from UniProt, AlphaFold DB, or PubChem. ESMFold is used for on-demand protein structure prediction if AlphaFold structures are unavailable.
2.  **Format Conversion & Preparation**: OpenBabel is used to convert all input files to PDBQT format, clean protein structures (e.g., remove water), add hydrogens, and generate 3D ligand structures.
3.  **Structure Verification**: A dedicated module (`verify_structures.py`) performs quality checks, including atom counts, 3D coordinate validation, partial charge assignment, and molecular weight estimation for ligands. Warnings are provided for potential issues.
4.  **Docking Execution**: The Smina static binary performs molecular docking, generating 9 binding poses. Users can specify manual grid box coordinates for targeted docking or use the default blind docking.
5.  **Result Processing**: Binding affinity scores are extracted from Smina output logs, and PDBQT poses are converted to PDB format for visualization. Protein-ligand complex files are generated for each pose.

### Frontend Architecture
- **Technology**: Uses vanilla JavaScript for direct DOM manipulation and `fetch` API for backend communication.
- **Visualization**: Integrates the MolStar 3D molecular viewer via CDN for interactive display of protein-ligand complexes.
- **User Interface**: Features a progressive disclosure UI, revealing sections (upload, docking, results, viewer) sequentially. The design is responsive, optimized for desktop, with a clean white background.

### System Design Choices
- **Cross-Platform Compatibility**: Automatic OS detection for Smina and OpenBabel binary execution.
- **Error Handling**: Comprehensive error detection for file conversion, docking failures, and invalid user inputs.
- **Targeted Docking**: Implemented manual grid box configuration to allow users to define active site coordinates, improving accuracy and speed for known binding sites.

## External Dependencies

### Binary Tools (System-Level)
-   **Smina**: Static binary molecular docking engine (`smina.static`) for generating binding poses and affinity scores.
-   **OpenBabel**: Chemical toolkit for molecular file format conversion (e.g., PDB to PDBQT), adding hydrogens, and generating 3D structures from SMILES.

### Python Dependencies
-   **Flask**: Lightweight WSGI web framework for handling API endpoints and serving static files.
-   **Subprocess Module**: Used to execute Smina and OpenBabel binaries and capture their output.
-   **requests**: For making HTTP requests to external APIs (UniProt, PubChem).
-   **biopython**: For PDB file parsing and manipulation.
-   **Flask-CORS**: Enables Cross-Origin Resource Sharing for the MolStar viewer.

### Frontend Libraries (CDN)
-   **MolStar**: 3D molecular visualization toolkit (version 3.42.0) loaded from CDN, providing interactive rendering of molecular structures in WebGL.

### Automated Data Sources
-   **AlphaFold Database**: For downloading predicted protein structures via UniProt accession IDs.
-   **UniProt REST API**: For fetching protein FASTA sequences and information by protein name or ID.
-   **PubChem REST API**: For retrieving canonical SMILES strings of compounds by name.
-   **ESMFold API**: For on-demand, rapid protein structure prediction based on sequence.