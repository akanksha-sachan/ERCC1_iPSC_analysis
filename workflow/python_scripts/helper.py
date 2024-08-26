import os
import sys
import numpy as np
import pandas as pd

import multiprocessing as mp

############ helper functions for processing RNA-seq data ############

def run_fastqc(data_path, to_extract_fasta_file=None, raw_or_trimmed="raw", user_name="Akanksha",
               sbatch_template_file="/ix/cigcore/utils/code/fastqc_template.sbatch",
               do_remove_original_fastq=False, do_submit=True, email="idm17@pitt.edu",
               suffix="fastq.gz", char_num_after_basename="Auto", max_job_num=48, 
               do_remove_bam=False, do_write=True, transfer_prefix=None, minlen=None, 
               hours=20):
    """
    Run FastQC on the raw or trimmed fastq files.
    """
    
    def validate_path(path):
        if not os.path.isabs(path):
            raise ValueError("dataPath must be absolute, i.e., start with '/'!")
    
    def get_basenames(files, suffix, char_num_after_basename):
        if char_num_after_basename == "Auto":
            if suffix == "fastq.gz":
                if "_R1" in files[0] or "_R2" in files[0]:
                    print("R1/R2 exist in file names!")
                else:
                    raise ValueError("File names don't include R1/R2!")
                
                if len(re.findall(r'_R', files[0])) > 1:
                    raise ValueError("File names should not have more than one R1/R2 in their names!")
                
                char_num_after_basename = len(re.split(r'_R', files[0])[-1]) + len('_R')
            elif suffix == "bam":
                char_num_after_basename = 4
            else:
                raise ValueError("charNumAfterBasename could not be determined automatically!")
        
        basenames = list(set([re.sub(f'.{{{char_num_after_basename}}}$', '', f) for f in files]))
        return basenames, char_num_after_basename

    # Validate the data_path
    validate_path(data_path)
    
    # Set paths and create necessary directories
    qc_out_path = os.path.join(data_path, "qc", raw_or_trimmed)
    if raw_or_trimmed == "trimmed":
        data_path = os.path.join(data_path, "trimmed")
    
    split_data_path = data_path.split("/")
    proj_path = "/".join(split_data_path[:5])
    code_path = os.path.join(proj_path, "code", user_name)
    log_path = os.path.join(code_path, "log", "fastqc")
    os.makedirs(log_path, exist_ok=True)
    
    # List files and determine basenames
    files = [f for f in os.listdir(data_path) if f.endswith(suffix)]
    basenames, char_num_after_basename = get_basenames(files, suffix, char_num_after_basename)
    
    print("I found the following samples:")
    print(basenames)
    n_samples = len(basenames)
    print(f"Number of samples: {n_samples}")
    
    # Read and modify the sbatch template
    with open(sbatch_template_file, 'r') as file:
        template = file.read()
    
    fastqc_job = f"fastqc -o {qc_out_path} {data_path}/$FASTQNAMES[${{SLURM_ARRAY_TASK_ID}}]_{suffix}"
    
    script = template
    script = script.replace("LOGPATHOLDER", log_path)
    script = script.replace("hoursPlaceHolder", str(hours))
    script = script.replace("arrayLengthHolder", str(n_samples - 1))
    script = script.replace("emailHolder", email)
    script = script.replace("dataPathPlaceHolder", data_path)
    script = script.replace("qcFolderNameHolder", qc_out_path)
    script = script.replace("jobCommandHolder1", fastqc_job)
    
    sbatch_file_name = f"{split_data_path[-1]}_fastqc_{raw_or_trimmed}.sbatch"
    sbatch_file_path = os.path.join(code_path, sbatch_file_name)
    
    if do_write:
        with open(sbatch_file_path, 'w') as file:
            file.write(script)
        print(f"I wrote: {sbatch_file_path}")
        command_s = f"sbatch {sbatch_file_path}\n"
        print(command_s)
        
        if do_submit:
            sleep_seconds = 3
            while len(subprocess.getoutput("squeue -u $USER").splitlines()) - 1 >= max_job_num:
                print(f"Still too many jobs ({max_job_num}), sleeping for {sleep_seconds / 60} minutes ....")
                time.sleep(sleep_seconds)
                sleep_seconds *= 1.2
            
            subprocess.run(command_s, shell=True, executable='/bin/bash')
            print(f"Submitted at: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    if transfer_prefix:
        to_path = re.sub(r'^.*/proj/', '~/OneDrive\ -\ University\ of\ Pittsburgh/proj', data_path)
        transfer_command = f"scp {transfer_prefix} {os.path.join(data_path, 'qc')} {to_path}"
        print("Copy QC reports: Run the following from the local computer AFTER the cleaning jobs are done.")
        print(transfer_command)
    
    # Save the cleaning process information
    cleaning = {"basenames": basenames, "commands": command_s}
    with open(os.path.join(log_path, "cleaning.json"), 'w') as file:
        import json
        json.dump(cleaning, file)
    
    return cleaning

    pass