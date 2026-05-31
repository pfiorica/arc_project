import os
import time
import subprocess

os.chdir("/projects/rpci/songyao/pnfioric/arc_project/airr_calculation")
# Path to the folder containing job scripts
jobs_folder = "/projects/rpci/songyao/pnfioric/arc_project/airr_calculation/jobs_mixcr"

# Get list of job scripts
job_scripts = [os.path.join(jobs_folder, f) for f in os.listdir(jobs_folder) if f.endswith(".sh")]

# Submit jobs in batches of 100 per hour
batch_size = 800
for i in range(0, len(job_scripts), batch_size):
    batch = job_scripts[i:i+batch_size]
    print(f"Submitting batch {i//batch_size + 1}: {len(batch)} jobs")
    
    for job in batch:
        subprocess.run(["sbatch", job])
    
    print("Waiting for 800s before submitting the next batch...")
    time.sleep(800)  # Wait for 1 hour
