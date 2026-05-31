import os

# Set this to your top-level directory containing the sample subdirectories
top_dir = "/projects/rpci/songyao/pnfioric/arc_project/airr_calculation/fastq_files/"

output_file = "/projects/rpci/songyao/pnfioric/arc_project/fastq_list.txt"

with open(output_file, "w") as out_f:
    # Loop over each sample subdirectory
    for sample_dir in sorted(os.listdir(top_dir)):
        sample_path = os.path.join(top_dir, sample_dir)

        if os.path.isdir(sample_path):

            # Find all R1 FASTQs
            r1_files = sorted([
                f for f in os.listdir(sample_path)
                if f.endswith("_R1_001.fastq.gz")
            ])

            for r1 in r1_files:

                # Infer matching R2 FASTQ
                r2 = r1.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")

                # Optional: skip if R2 does not exist
                if not os.path.exists(os.path.join(sample_path, r2)):
                    print(f"Missing R2 for: {r1}")
                    continue

                # Use the directory name as the sample ID
                sample_id = sample_dir

                # Write output
                out_f.write(
                    f"{sample_dir}/{r1}\t"
                    f"{sample_dir}/{r2}\t"
                    f"{sample_id}\n"
                )
