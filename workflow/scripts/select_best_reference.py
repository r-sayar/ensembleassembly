import json
import os # For checking file existence, good practice

# This part is for standalone testing.
# When run by Snakemake, json_file_path would come from snakemake.input.multiqc_report
# and output_file_path from snakemake.output.best_reference
json_file_path = snakemake.input.multiqc_report + "/multiqc_data.json"#"results/qc/multiqc_all/denovo/multiqc_data/multiqc_data.json"
output_file_path = snakemake.output.best_reference#"results/best_reference.txt"

# --- Robust file loading (recommended) ---
if not os.path.exists(json_file_path):
    error_msg = f"CRITICAL ERROR: MultiQC data file not found at {json_file_path}"
    print(error_msg)
    with open(output_file_path, "w") as f_err: # Write error to output for Snakemake
        f_err.write(error_msg + "\n")
    raise FileNotFoundError(error_msg)

try:
    with open(json_file_path) as f:
        multiqc_data = json.load(f)
except json.JSONDecodeError as e:
    error_msg = f"CRITICAL ERROR: Could not decode JSON from {json_file_path}. Error: {e}"
    print(error_msg)
    with open(output_file_path, "w") as f_err:
        f_err.write(error_msg + "\n")
    raise
# --- End of robust file loading ---

references = []

# Access the actual raw data sections
saved_raw_data = multiqc_data.get("report_saved_raw_data", {})
busco_data_dict = saved_raw_data.get("multiqc_busco", {})
quast_data_dict = saved_raw_data.get("multiqc_quast", {})

if not busco_data_dict:
    print("Warning: 'multiqc_busco' data not found in 'report_saved_raw_data'.")
if not quast_data_dict:
    print("Warning: 'multiqc_quast' data not found in 'report_saved_raw_data'.")

# Iterate through BUSCO results
# busco_sample_name will be like "denovo_sample1", "sample1", etc.
# busco_metrics will be the dictionary of BUSCO scores for that sample
for busco_sample_name_original, busco_metrics in busco_data_dict.items():
    if not isinstance(busco_metrics, dict):
        print(f"Warning: Skipping malformed BUSCO entry for {busco_sample_name_original}")
        continue

    # Remove "denovo_" prefix if it exists
    cleaned_sample_name = busco_sample_name_original
    if busco_sample_name_original.startswith("denovo_"):
        cleaned_sample_name = busco_sample_name_original[len("denovo_"):]
        print(f"Info: Cleaned sample name from '{busco_sample_name_original}' to '{cleaned_sample_name}'")


    complete_buscos = busco_metrics.get("complete", 0.0)
    total_buscos = busco_metrics.get("total", 0.0)
    busco_score_percentage = 0.0
    if total_buscos > 0:
        busco_score_percentage = (float(complete_buscos) / float(total_buscos)) * 100.0
    else:
        print(f"Warning: Total BUSCOs is 0 for sample {cleaned_sample_name} (original: {busco_sample_name_original}), completeness set to 0%.")

    # Find corresponding QUAST data for the same original sample name
    # QUAST data keys might still have the "denovo_" prefix from MultiQC's processing
    quast_sample_metrics = quast_data_dict.get(busco_sample_name_original) # Use original name for lookup
    n50 = 0.0 # Default N50 if QUAST data or N50 key is missing

    if quast_sample_metrics and isinstance(quast_sample_metrics, dict):
        n50_value = quast_sample_metrics.get("N50")
        if n50_value is not None:
            try:
                n50 = float(n50_value)
            except ValueError:
                print(f"Warning: Could not convert N50 value '{n50_value}' to float for sample {cleaned_sample_name} (original: {busco_sample_name_original}). Using 0.0.")
        else:
            print(f"Warning: N50 key not found for sample {cleaned_sample_name} (original: {busco_sample_name_original}) in QUAST data. Using 0.0.")
    else:
        print(f"Info: No QUAST data found for BUSCO sample {cleaned_sample_name} (original: {busco_sample_name_original}). N50 will be 0.0.")
        
    references.append((cleaned_sample_name, busco_score_percentage, n50)) # Use cleaned_sample_name
    print(f"Processed: {cleaned_sample_name}, BUSCO %: {busco_score_percentage:.2f}, N50: {n50}")


if not references:
    error_msg = "CRITICAL ERROR: No reference samples could be processed from MultiQC data."
    print(error_msg)
    with open(output_file_path, "w") as f_err:
        f_err.write(error_msg + "\n")
    raise ValueError(error_msg)

# Sort by BUSCO score (descending), then N50 (descending)
best_reference_tuple = max(references, key=lambda x: (x[1], x[2]))

print(f"Selected best reference: {best_reference_tuple}")

with open(output_file_path, "w") as f:
    f.write(f"Best Reference: {best_reference_tuple[0]}\n")
    f.write(f"BUSCO Score (%): {best_reference_tuple[1]:.2f}\n") # Formatted percentage
    f.write(f"N50: {best_reference_tuple[2]}\n")

print(f"Best reference information written to {output_file_path}")