import os
import sys
import pandas as pd
import re
import shutil

def rename_files(csv_file, input_directory, output_directory):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv_file)

    # Iterate through each file in the input directory
    for filename in os.listdir(input_directory):
        # Extract the accession code and read direction from the filename using regex
        accession_match = re.match(r'^(SRR\d+)_(F|R)_paired\.fastq\.gz', filename)
        if accession_match:
            accession_code = accession_match.group(1)
            read_direction = accession_match.group(2)
            
            # Find the matching row in the DataFrame based on the accession code
            row = df[df['Run'] == accession_code]

            if not row.empty:
                country = row['geo_loc_name_country'].values[0]
                collection_date = row['Collection_Date'].values[0]
                
                # Replace spaces with underscores
                country = country.replace(' ', '_')
                
                # Create the new file name (Run_geo_loc_name_country_Collection_Date_readDirection.fastq.gz)
                new_file_name = f"{accession_code}_{country}_{collection_date}_{read_direction}.fastq.gz"

                # Get the current and new file paths
                current_file_path = os.path.join(input_directory, filename)
                new_file_path = os.path.join(output_directory, new_file_name)

                # Copy the renamed file to the output directory
                shutil.copy(current_file_path, new_file_path)
                print(f"Copied {current_file_path} to {new_file_path}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python rename_files.py <csv_file> <input_directory> <output_directory>")
    else:
        csv_file = sys.argv[1]
        input_directory = sys.argv[2]
        output_directory = sys.argv[3]
        rename_files(csv_file, input_directory, output_directory)
