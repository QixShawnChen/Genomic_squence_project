from flask import Flask, render_template, request, redirect, url_for, send_file
import requests
import os
import utility  
import numpy as np

import uuid

app = Flask(__name__)





@app.route("/")
def home():
    return render_template("home.html")


@app.route("/enter-sequence", methods=["GET", "POST"])
def enter_sequence():
    results = {}
    if request.method == "POST":
        dna_sequence = request.form["sequence"].upper()
        selected_functions = request.form.getlist("functions")

        if "transcription" in selected_functions:
            results["Transcription"] = utility.transcribe_dna_to_rna(dna_sequence)

        if "translation" in selected_functions:
            results["Translation"] = utility.translate_orf_with_ambiguity(dna_sequence)

        if "find_orfs" in selected_functions:
            orfs = utility.find_orfs(dna_sequence)
            results["ORFs"] = orfs if orfs else "No ORFs found"

        if "gc_content" in selected_functions:
            gc_content = utility.calculate_gc_content(dna_sequence)
            results["GC Content"] = f"{gc_content:.2f}%"

        if "cpg_islands" in selected_functions:
            cpg_islands = utility.find_cpg_islands(dna_sequence)
            results["CpG Islands"] = cpg_islands if cpg_islands else "No CpG islands found"

        if "codon_usage" in selected_functions:
            codon_usage = utility.codon_usage(dna_sequence)
            results["Codon Usage"] = dict(codon_usage)

        if "restriction_sites" in selected_functions:
            restriction_sites = utility.find_restriction_sites(dna_sequence)
            results["Restriction Sites"] = restriction_sites

    return render_template("enter_sequence.html", results=results)


@app.route("/upload-file", methods=["GET", "POST"])
def upload_file():
    results = []
    output_file_path = None  # Initialize output file path
    download_link = None  # Initialize download link
    
    if request.method == "POST":
        if "fasta_file" in request.files and request.files["fasta_file"].filename:
            file = request.files["fasta_file"]
            sequences = utility.parse_fasta(file)  

            selected_functions = request.form.getlist("functions")

            # Process each sequence and generate results
            for seq_name, dna_sequence in sequences.items():
                sequence_results = {"sequence_name": seq_name, "original_sequence": dna_sequence, "analysis": {}}

                if "transcription" in selected_functions:
                    sequence_results["analysis"]["Transcription"] = utility.transcribe_dna_to_rna(dna_sequence)

                if "translation" in selected_functions:
                    sequence_results["analysis"]["Translation"] = utility.translate_orf_with_ambiguity(dna_sequence)

                if "find_orfs" in selected_functions:
                    orfs = utility.find_orfs(dna_sequence)
                    sequence_results["analysis"]["ORFs"] = orfs if orfs else "No ORFs found"

                if "gc_content" in selected_functions:
                    gc_content = utility.calculate_gc_content(dna_sequence)
                    sequence_results["analysis"]["GC Content"] = f"{gc_content:.2f}%"

                if "cpg_islands" in selected_functions:
                    cpg_islands = utility.find_cpg_islands(dna_sequence)
                    sequence_results["analysis"]["CpG Islands"] = cpg_islands if cpg_islands else "No CpG islands found"

                if "codon_usage" in selected_functions:
                    codon_usage = utility.codon_usage(dna_sequence)
                    sequence_results["analysis"]["Codon Usage"] = dict(codon_usage)

                if "restriction_sites" in selected_functions:
                    restriction_sites = utility.find_restriction_sites(dna_sequence)
                    sequence_results["analysis"]["Restriction Sites"] = restriction_sites

                results.append(sequence_results)

            # Create a single text file for all results
            if results:
                # Create output directory
                output_dir = "output"
                os.makedirs(output_dir, exist_ok=True)

                # Generate a unique filename for the results file
                output_file_path = os.path.join(output_dir, f"results_{uuid.uuid4().hex}.txt")

                # Write results to the file
                with open(output_file_path, "w") as output_file:
                    for result in results:
                        output_file.write(f"Results for {result['sequence_name']}:\n")
                        output_file.write(f"Original Sequence: {result['original_sequence']}\n")
                        for key, value in result["analysis"].items():
                            output_file.write(f"{key}: {value}\n")
                        output_file.write("\n")

                # Generate a relative URL for the download
                download_link = f"/download/{os.path.basename(output_file_path)}"

        else:
            # Handle case where no file was uploaded
            results = {"error": "No file was uploaded or file was invalid."}

    # Always return a response, even if no file was uploaded
    return render_template("upload_file.html", results=results, download_link=download_link)





@app.route("/search-uniprot", methods=["GET", "POST"])
def search_uniprot():
    uniprot_results = {}
    pdb_url = None
    png_url = None
    uniprot_accession = ""

    if request.method == "POST":
        uniprot_accession = request.form["uniprot_accession"].strip()
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_accession}.json"
        response = requests.get(url)

        if response.status_code == 200:
            uniprot_data = response.json()
            uniprot_results = {
                "ID": uniprot_data.get("id"),
                "Accession": uniprot_data.get("accession"),
                "Protein Name": uniprot_data.get("protein", {}).get("recommendedName", {}).get("fullName", {}).get("value"),
                "Organism": uniprot_data.get("organism", {}).get("scientificName"),
                "Sequence": uniprot_data.get("sequence", {}).get("value"),
            }

            # Fetch AlphaFold data for the PDB and PNG
            alphafold_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_accession}"
            alphafold_response = requests.get(alphafold_url)

            if alphafold_response.status_code == 200:
                alphafold_data = alphafold_response.json()
                if alphafold_data:
                    pdb_url = alphafold_data[0].get("pdbUrl")
                    png_url = alphafold_data[0].get("paeImageUrl")
        else:
            uniprot_results = {"Error": f"Accession {uniprot_accession} not found."}

    return render_template(
        "alphafold_res.html",
        uniprot_accession = uniprot_accession,
        uniprot_results=uniprot_results,
        pdb_url=pdb_url,
        png_url=png_url,
    )


@app.route("/jaspar_id", methods=["GET", "POST"])
def jaspar():
    results = None
    if request.method == "POST":
        taxonomy_id = request.form["taxonomy_id"]
        try:
            matrices = utility.fetch_matrices_for_species(taxonomy_id)
            return render_template("matrix_list.html", taxonomy_id=taxonomy_id, matrices=matrices)
        except Exception as e:
            results = {"error": str(e)}

    return render_template("jaspar_id.html", results=results)



@app.route("/jaspar/scan", methods=["GET", "POST"])
def jaspar_scan():
    taxonomy_id = request.args.get("taxonomy_id")
    matrices = request.args.getlist("matrices")
    results = None
    pfm = None
    pwm = None
    download_link = None

    if request.method == "POST":
        dna_sequence = request.form.get("sequence", "").strip().upper()
        matrix_id = request.form["matrix_id"]
        file_upload = request.files.get("file_upload")

        file_results = []
        try:
            # Process uploaded file if provided
            if file_upload and file_upload.filename:
                # Parse uploaded sequences
                file_sequences = utility.parse_fasta(file_upload)
                for seq_name, seq in file_sequences.items():
                    try:
                        # Fetch PFM and PWM
                        pfm = utility.fetch_pfm_from_jaspar(matrix_id)
                        pwm = utility.convert_pfm_to_pwm(pfm)

                        # Scan the sequence for binding sites
                        binding_sites = utility.scan_sequence_with_pwm(seq, pwm)

                        file_results.append(f"Results for {seq_name}:\n")
                        file_results.append(f"Original Sequence: {seq}:\n")
                        file_results.append(f"Binding Sites: {binding_sites}\n\n")
                    except ValueError as e:
                        # Handle specific case where sequence is too short
                        file_results.append(f"Results for {seq_name}:\n")
                        file_results.append(f"Error: {str(e)}\n\n")
                        print(f"Skipping sequence {seq_name} due to error: {str(e)}")
                    except Exception as e:
                        # Handle other unexpected errors
                        file_results.append(f"Results for {seq_name}:\n")
                        file_results.append(f"An unexpected error occurred: {str(e)}\n\n")
                        print(f"Skipping sequence {seq_name} due to unexpected error: {str(e)}")

                # Write results to a file
                output_dir = "output"
                os.makedirs(output_dir, exist_ok=True)
                output_file = os.path.join(output_dir, f"results_{uuid.uuid4().hex}.txt")
                with open(output_file, "w") as f:
                    f.writelines(file_results)

                download_link = f"/download/{os.path.basename(output_file)}"

            # Process entered sequence if provided
            if dna_sequence:
                pfm = utility.fetch_pfm_from_jaspar(matrix_id)
                pwm = utility.convert_pfm_to_pwm(pfm)
                binding_sites = utility.scan_sequence_with_pwm(dna_sequence, pwm)

                if not binding_sites:
                    results = {"error": f"No binding sites found for matrix ID {matrix_id}."}
                else:
                    results = {"matrix_id": matrix_id, "binding_sites": binding_sites}
        except Exception as e:
            print(f"Error occurred: {str(e)}")
            results = {"error": f"An error occurred: {str(e)}"}

    return render_template(
        "jaspar_scan.html",
        taxonomy_id=taxonomy_id,
        matrices=matrices,
        results=results,
        pfm=pfm,
        pwm=pwm,
        download_link=download_link,
    )


@app.route("/download/<filename>")
def download_file(filename):
    output_dir = "output"
    file_path = os.path.join(output_dir, filename)
    try:
        return send_file(file_path, as_attachment=True)
    except FileNotFoundError:
        return "File not found", 404



if __name__ == "__main__":
    app.run(debug=True)



