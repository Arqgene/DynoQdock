# 🧬 ArqGene - Local Installation & Troubleshooting Guide

This guide provides a comprehensive, step-by-step process to set up and run the ArqGene Molecular Docking Platform on your local Windows machine.

---

## 🚀 Phase 1: Prerequisites & Preparation
Before starting, ensure your system meets the minimum requirements.

### 1. System Requirements
*   **OS**: Windows 10 or 11 (64-bit).
*   **Storage**: At least 500MB of free space.
*   **Internet**: Required for the initial setup and fetching molecular structures.

### 2. Extracting the Project
1.  Download the project ZIP file.
2.  Right-click the file and select **Extract All...**.
3.  Choose a simple folder path like `C:\ArqGene` to avoid permission issues.

---

## 🛠️ Phase 2: Automated Installation
We have provided a "one-click" script that handles the technical setup for you.

### 1. Run the Setup Script
1.  Open the extracted folder.
2.  Find the file named **`setup_and_run.bat`**.
3.  Double-click it to start.

### 2. What the Script Does Automatically:
*   **Checks for Python 3.11**: Downloads and installs it if missing.
*   **Checks for OpenBabel**: Installs the molecular conversion engine automatically.
*   **Downloads Smina**: Fetches the docking engine required for calculations.
*   **Sets up Environment**: Creates a virtual "sandbox" for the app's libraries.
*   **Launches the App**: Opens your browser to `http://localhost:5000`.

---

## 🧪 Phase 3: Running Your First Docking
Once the browser opens, follow these steps to verify everything is working.

1.  **Login**: Use your credentials (or create a new local account).
2.  **Prepare Protein**: Enter a UniProt ID (e.g., `P00734`) and click **Prepare Protein**.
3.  **Prepare Ligand**: Enter a compound name (e.g., `warfarin`) and click **Prepare Ligand**.
4.  **Run Docking**: Click the **Run Smina Docking** button.
5.  **Visualize**: Click **View Pose** to explore the 3D results.

---

## 🆘 Phase 4: Troubleshooting & Common Errors
If something goes wrong, check these common scenarios.

### ❌ Error: "Python is not recognized"
*   **Cause**: Python was installed but the "Add to PATH" option was missed.
*   **Fix**: The script tries to fix this automatically. If it still fails, restart your computer to refresh system variables.

### ❌ Error: "The system cannot find the path specified"
*   **Cause**: You might be running the script directly inside the ZIP file without extracting it first.
*   **Fix**: Make sure you have **extracted** the ZIP to a folder on your Desktop or C: drive.

### ❌ Error: "Access Denied" or Admin Prompts
*   **Cause**: Windows needs permission to install OpenBabel or Python.
*   **Fix**: If a window pops up asking for permission, click **Yes**. You can also try right-clicking `setup_and_run.bat` and selecting **Run as Administrator**.

### ❌ Error: Browser doesn't open
*   **Cause**: The server took longer than 15 seconds to start.
*   **Fix**: Check the command window. If it says "Server is running," manually type `http://localhost:5000` in your browser.

### ❌ Error: "Smina failed" in the web app
*   **Cause**: Antivirus software might be blocking the `smina.exe` file.
*   **Fix**: Add the project folder to your Antivirus "Exclusions" list.

---

## 📝 Final Checklist for Success
- ✅ **Extract** the files (don't run from ZIP).
- ✅ **Internet** is connected.
- ✅ **Double-click** `setup_and_run.bat`.
- ✅ **Wait** for the "Setup Complete" message.
