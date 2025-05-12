#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import datetime
import glob
import re
from multiprocessing import Pool
import time
import pickle
from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from bs4 import BeautifulSoup
import argparse
import shutil
import re
import csv
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()
parser = argparse.ArgumentParser(
    description="Pan-Genomic-Guided crRNA Design: A sophisticated tool for designing highly specific crRNAs targeting user-defined bacterial genera."
)
required = parser.add_argument_group('Required arguments')
required.add_argument(
    '-g', '--genus',
    required=True,
    type=str,
    metavar='STR',
    help='Target bacterial genus for analysis (required).'
)
required.add_argument(
    '-offg', '--offgenomes',
    required=True,
    type=str,
    metavar='PATH',
    help='Path to the genome file for off-target effect evaluation (required).'
)
required.add_argument(
    '-firedriver', '--firefox_driver_path',
    required=True,
    type=str,
    metavar='PATH',
    help='Path to the Firefox WebDriver executable (required).'
)
required.add_argument(
    '-ffry', '--flashfry_path',
    required=True,
    type=str,
    metavar='PATH',
    help='Path to the FlashFry executable (required).'
)
optional = parser.add_argument_group('Optional arguments')
optional.add_argument(
    '-i', '--input_path',
    required=False,
    type=str,
    metavar='PATH',
    help='Path to the directory containing the genomes for analysis. If provided, the genomes will not be downloaded. This argument can be used independently.'
)
optional.add_argument(
    '-gff', '--gff_file',
    required=False,
    type=str,
    metavar='PATH',
    help='Path to Prokka-generated GFF files (requires --input_path). If provided, genomes will not be annotated.'
)
optional.add_argument(
        '-coregene', '--core_gene_folder',
        required=False,
        type=str,
        metavar='PATH',
        help='Path to the directory containing pre-computed core genes. If provided, the pipeline will skip genome download, annotation, and pan-genome analysis, starting directly from crRNA design (requires --input_path).'
    )
optional.add_argument(
    '-m', '--memory',
    required=False,
    type=str,
    default='12g',
    metavar='STR',
    help='Memory allocation for FlashFry (default: 12g).'
)
optional.add_argument(
    '-maxOff', '--maximumOffTargets',
    required=False,
    type=int,
    metavar='INT',
    default=12000,
    help='Maximum number of off-target sites to consider (default: 12000).'
)
optional.add_argument(
    '-PAM', '--PAM_sequence',
    required=False,
    default='TTTN',
    type=str,
    metavar='SEQ', 
    help='Protospacer Adjacent Motif (PAM) sequence (default: TTTN).'
)
optional.add_argument(
    '-t', '--threads',
    required=False,
    default=8,
    type=int,
    metavar='INT',
    help='Number of threads for parallel processing (default: 8).'
)
args = parser.parse_args()
genus = args.genus
input_path = args.input_path
gff_file = args.gff_file
core_gene_folder = args.core_gene_folder
offgenomes = args.offgenomes
firefox_driver_path = args.firefox_driver_path
threads = args.threads
PAM = args.PAM_sequence
flashfry_path = args.flashfry_path
memory = args.memory
maximumOffTargets = args.maximumOffTargets
max_processes = os.cpu_count() 
def download_genomes(genus_name):
    start_time = time.time()
    print(f"[INFO] Downloading of genomes for {genus_name} from NCBI GenBank database...")
    output_dir = f"{genus_name.replace(' ', '_')}_all_genbank"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    download_cmd = f"ncbi-genome-download -s genbank -F fasta -l all -r 8 -p 4 --flat-output -P bacteria -g \"{genus_name}\" -o \"{output_dir}\""
    subprocess.run(download_cmd, shell=True, check=True)
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"[INFO] Downloading genomes for {genus_name} took {execution_time} seconds.") 
    return os.path.abspath(output_dir)
def extract_gz_files(genus_name, input_dir):
    print(f"[INFO] Checking for *.gz files for {genus_name}...")
    os.chdir(input_dir) 
    gz_files = [file for file in os.listdir(input_dir) if file.endswith('.gz')]
    extracted_file_list = []
    num_genomes = 0
    if not gz_files:
        print(f"[INFO] No *.gz files found for {genus_name}. Proceeding with available files.")
        extracted_file_list = [os.path.abspath(file) for file in os.listdir(input_dir) if file.endswith('.fna')]
        num_genomes = len(extracted_file_list)
        print(f"[INFO] Total number of genomes extracted for {genus_name}: {num_genomes}")
    else:
        start_time = time.time()
        print(f"[INFO] Decompressing genome files for {genus_name}...")
        decompress_cmd = "gunzip -dv *.gz"
        subprocess.run(decompress_cmd, shell=True, check=True)
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"[INFO] Decompressing genome files for {genus_name} took {execution_time} seconds.")
        extracted_file_list = [os.path.abspath(file) for file in os.listdir(input_dir) if file.endswith('.fna')] 
        num_genomes = len(extracted_file_list)
        print(f"[INFO] Total number of genomes extracted for {genus_name}: {num_genomes}")
    return extracted_file_list, num_genomes
def annotate_genome(fasta_path):
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"File {fasta_path} not found")
    filename = os.path.splitext(os.path.basename(fasta_path))[0] 
    output_dir = filename + "_prokka"
    prokka_cmd = f"prokka --prefix {filename} --cpus {threads} --norrna --notrna --outdir {output_dir} {fasta_path}" 
    subprocess.run(prokka_cmd, shell=True, check=True)
def annotate_genomes_sequential(extracted_file_list):
    for fasta_path in extracted_file_list:
        annotate_genome(fasta_path)
def get_gff_files(input_dir, gff_dir):    
    gff_files = glob.glob(os.path.join(input_dir, '**', '*.gff'), recursive=True)
    gff_files = [os.path.abspath(file) for file in gff_files]
    gff_dir = gff_dir.replace(' ', '_')  
    if not os.path.exists(gff_dir):
        os.makedirs(gff_dir)
    for file in gff_files:
        shutil.copy(file, gff_dir)
    return os.path.abspath(gff_dir)
def pan_genome_analysis(gff_dir, num_threads):
    start_time = time.time()
    print("[INFO] Starting pan-genome analysis...")
    if not gff_dir:
        raise ValueError("No input files provided for pan-genome analysis")
    output_dir = "panaroo_output"
    os.chdir(gff_dir)
    panaroo_cmd = f"panaroo -i *.gff -o {output_dir} --clean-mode strict -t {num_threads} --remove-invalid-genes" 
    subprocess.run(panaroo_cmd, shell=True, check=True)
    combined_DNA_CDS = "combined_DNA_CDS.fasta"
    combined_DNA_CDS_path = os.path.join(output_dir, combined_DNA_CDS)
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"[INFO] Pan-genome analysis completed. Total time consumed: {execution_time} seconds.")
    return os.path.abspath(output_dir), os.path.abspath(combined_DNA_CDS_path)
def extract_core_genes(panaroo_output_dir, total_genomes):
    os.chdir(panaroo_output_dir)
    core_gene_folder = "core_gene_folder" 
    if not os.path.exists(core_gene_folder):
        os.makedirs(core_gene_folder)
    file_path_rtab = "gene_presence_absence.Rtab"
    file_path_csv = "gene_presence_absence.csv"
    file_path_simple_csv = "gene_data.csv" 
    ro.globalenv['panaroo_output_dir'] = panaroo_output_dir
    ro.globalenv['file_path_rtab'] = file_path_rtab
    ro.globalenv['file_path_csv'] = file_path_csv
    ro.globalenv['file_path_simple_csv'] = file_path_simple_csv
    ro.globalenv['total_genomes'] = total_genomes
    ro.globalenv['core_gene_folder'] = core_gene_folder
    ro.r('''
    x <- read.delim(file_path_rtab)
    n_coregenes <- dim(as.data.frame(x$Gene[rowSums(x[, -1]) >= total_genomes*0.99])) [1] 
    n_softcoregenes <- dim(as.data.frame(x$Gene[rowSums(x[, -1]) >= total_genomes*0.95 & rowSums(x[, -1])<total_genomes*0.99]))[1] 
    n_shellgenes <- dim(as.data.frame(x$Gene[rowSums(x[, -1]) >= total_genomes*0.15 & rowSums(x[, -1])<total_genomes*0.95]))[1] 
    n_cloudgenes <- dim(as.data.frame(x$Gene[rowSums(x[, -1]) >= 0 & rowSums(x[, -1])<total_genomes*0.15]))[1] 
    target <- as.data.frame(x$Gene[rowSums(x[, -1]) == total_genomes]) 
    y <- read.csv(file_path_csv)
    target_gene <- subset(y,y$Gene%in%target$`x$Gene[rowSums(x[, -1]) == total_genomes]`)
    rows_with_semicolon <- apply(target_gene[, 4:ncol(target_gene)], 1, function(row) any(grepl(";", row)))
    target_gene_cleaned <- target_gene[!rows_with_semicolon, ]
    library(tidyr)
    df_long <- gather(target_gene_cleaned, key = "gff_file", value = "annotation_id", 4:ncol(target_gene_cleaned))
    gene_data_simpe <- read.csv(file_path_simple_csv)
    gene_data_simpe_01 <- gene_data_simpe[,c(3,4)]
    df_long$annotation_id <- gsub("_pseudo", "", df_long$annotation_id)
    z <- merge(df_long,gene_data_simpe_01,by.x = "annotation_id",by.y = "annotation_id",all.x=FALSE)
    genes <- unique(z$Gene)
    for(category in genes) {
    subset_data <- z[z$Gene == category, 6]
    filename <- file.path(core_gene_folder, paste("output_", category, ".csv", sep = ""))
    write.table(subset_data, file = filename, row.names = F,col.names = F,sep = "\t",quote = F)
    }
    ''')
    summary_statistics = {}
    with open("summary_statistics.txt", 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) == 3:
                category = parts[0].strip().lower().replace(' ', '_')
                count = int(parts[2].strip())
                summary_statistics[category] = count
    n_coregenes = ro.globalenv['n_coregenes'][0]
    n_softcoregenes = ro.globalenv['n_softcoregenes'][0]
    n_shellgenes = ro.globalenv['n_shellgenes'][0]
    n_cloudgenes = ro.globalenv['n_cloudgenes'][0]
    matches = {
        'core_genes': n_coregenes == summary_statistics.get('core_genes', 0),
        'soft_core_genes': n_softcoregenes == summary_statistics.get('soft_core_genes', 0),
        'shell_genes': n_shellgenes == summary_statistics.get('shell_genes', 0),
        'cloud_genes': n_cloudgenes == summary_statistics.get('cloud_genes', 0)
    }
    for category, match in matches.items():
        if match:
            print(f"{category.replace('_', ' ').title()} match the summary statistics.") 
        else:
            print(f"{category.replace('_', ' ').title()} do not match the summary statistics.")
    return os.path.abspath(core_gene_folder)
def extract_fasta_from_csv(panaroo_output_dir, core_gene_folder, combined_DNA_CDS):
    os.chdir(panaroo_output_dir)
    csv_files = [os.path.join(core_gene_folder, f) for f in os.listdir(core_gene_folder) if f.endswith('.csv')]
    for file in csv_files:
        filename = os.path.splitext(os.path.basename(file))[0]  
        output_fa_file = os.path.join(core_gene_folder, f'{filename}.fa')
        grep_cmd = f"seqkit grep -n -f {file} {combined_DNA_CDS} -j {threads} -o {output_fa_file}"
        subprocess.run(grep_cmd, shell=True, check=True)
def design_crRNA(core_gene_folder):
    os.chdir(core_gene_folder)
    ro.r('''
        library(CaSilico)
    ''')
    ro.globalenv['num_threads'] = threads
    ro.globalenv['core_gene_folder'] = core_gene_folder
    r_code = '''
        run_CaSilico <- function(file_path) {
            file_name <- tools::file_path_sans_ext(basename(file_path))
            message(paste("Running file:", file_name))
            setwd(core_gene_folder) 
            CaSilico(ResultsFolder = file_name,
                     TargetFasta = file_path,
                     CrisprTypes = "casV_A",
                     ConservationMethod = 1,
                     ConservationThreshold = 1,
                     OffTarget = F,
                     OffAsk = F,
                     Threads = num_threads)
            message(paste("The file", file_name, "has finished running"))
        }
    '''
    ro.r(r_code) 
    run_CaSilico = ro.globalenv['run_CaSilico'] 
    input_list = glob.glob(os.path.join(core_gene_folder, '*.fa'))
    for input_file in input_list:
        run_CaSilico(input_file) 
def web_crawler(folder_path, firefox_driver_path, output_folder):
    firefox_options = Options()
    firefox_options.add_argument("--headless") 
    driver = webdriver.Firefox(executable_path=firefox_driver_path, options=firefox_options)
    gene_name = folder_path.split('output_')[-1] 
    html_file_path = os.path.join(folder_path, 'Spacers_Information.Html')
    driver.get(f"file://{html_file_path}") 
    all_data = []
    headers = []
    prev_rows = None 
    while True: 
        html = driver.page_source
        soup = BeautifulSoup(html, 'html.parser')
        table = soup.find('table')
        if not headers:
            header_rows = table.find_all('tr')[:2] 
            top_row = [cell.text.strip() for cell in header_rows[0].find_all('th')] 
            bottom_row = [cell.text.strip() for cell in header_rows[1].find_all('th')]
            for header in top_row:
                if header == "Off_Target": 
                    for item in bottom_row:
                        headers.append(f"{header}_{item}") 
                else:
                    headers.append(header)
        rows = []
        for row in table.find_all('tr')[1:]: 
            cols = row.find_all('td') 
            if cols:
                row_data = [ele.text.strip() for ele in cols]
                rows.append(row_data)
        if not rows or rows == prev_rows:
            print(f"Gene {gene_name} processing completed.")
            break
        else:
            all_data.extend(rows)
            prev_rows = rows 
        next_page_button = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.LINK_TEXT, "Next"))
        )
        driver.execute_script("arguments[0].scrollIntoView();", next_page_button) 
        next_page_button.click() 
        WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.TAG_NAME, "body"))
        )
    df = pd.DataFrame(all_data, columns=headers)
    df['core gene'] = gene_name
    csv_file_path = os.path.join(output_folder, f'{gene_name}_paged_table_data.csv') 
    df.to_csv(csv_file_path, index=False)
    driver.quit()
def web_crawler_in_parallel(core_gene_folder, driver_path):
    os.chdir(core_gene_folder)
    output_folder = "crawler_data"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder) 
    pattern = re.compile(r'^casV_A-ConservationMethod_1-threshold_100-result_output_.*$')
    folder_paths = [os.path.join(core_gene_folder, name) for name in os.listdir(core_gene_folder) if os.path.isdir(os.path.join(core_gene_folder, name))]
    filtered_folders = [folder_path for folder_path in folder_paths if pattern.match(os.path.basename(folder_path))]
    max_processes = os.cpu_count()
    with Pool(processes=max_processes) as pool:
        results = [pool.apply_async(web_crawler, args=(folder_path, driver_path, output_folder)) for folder_path in filtered_folders]
        for result in results:
            result.get()
def summarize_crawler_results(core_gene_folder):
    start_time = time.time()
    print("[INFO] Initiating the summary analysis of crRNA...")
    os.chdir(os.path.join(core_gene_folder, "crawler_data"))
    current_dir = os.getcwd()
    ro.globalenv['current_dir'] = current_dir
    r_code = '''
    setwd(current_dir)
    file_list <- list.files(path = current_dir) 
    df <- read.csv(file_list[1]) 
    for (file in file_list[-1]) { 
    temp_df <- read.csv(file)
    df <- rbind(df, temp_df)
    }
    y <- df[, -c(1, ncol(df)-1)]
    num_rows <- nrow(y)
    gRNA_values <- paste0("gRNA", sprintf("%04d", 1:num_rows))  
    y[, 2] <- gRNA_values
    write.csv(y, "all predicted gRNA.csv", row.names = FALSE)
    z <- y
    z$GC_content <- as.numeric(gsub("%", "", z$GC_content))
    filtered_z <- z[z$GC_content >= 40 & z$GC_content <= 65,]
    filtered_z <- filtered_z[filtered_z$Rank == 1,]
    filtered_z <- filtered_z[!grepl("Unnatural", filtered_z$Self_Complementarity_Free_Energy), ]
    filtered_z <- filtered_z[filtered_z$Conservation.Efficiency == 1,]
    write.csv(filtered_z, "filtered gRNA.csv", row.names = FALSE)
    format_fasta <- function(df) {
    fasta_lines <- paste0(">", df$Guide_RNA_ID, "\n", df$Guide_RNA, "\n")
    return(fasta_lines)
    }
    fasta_lines <- format_fasta(filtered_z)
    fasta_file_path <- "initial filtered guide RNA.fasta"
    writeLines(fasta_lines, fasta_file_path)
    fasta_file_absolute_path <- normalizePath(fasta_file_path)
    '''
    ro.r(r_code) 
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"[INFO] Summary analysis of crRNA took {execution_time} seconds.")
    initial_filtered_gRNA = str(ro.globalenv['fasta_file_absolute_path'][0])  
    return initial_filtered_gRNA
def add_pam_sequence(gRNA):
    return PAM + gRNA
def generate_fai(fasta_file):
    fai_file = fasta_file + ".fai"
    if not os.path.exists(fai_file):
        start_time = time.time()
        cmd = f"samtools faidx {fasta_file}" 
        subprocess.run(cmd, shell=True, check=True)
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"[INFO] Building index for {fasta_file} took {execution_time} seconds.")
    seqkit_fai_file = fasta_file + ".seqkit.fai" 
    shutil.copy(fai_file, seqkit_fai_file) 
def find_gRNA_targets(gRNA, genome_file):
    cmd = f"seqkit locate -d -i -m 0 -j 1 -p {gRNA} {genome_file}" 
    result = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    return len(result.strip().split('\n')) > 1 
def process_genome(args):
    gRNA, genome_file = args 
    pam_gRNA = add_pam_sequence(gRNA)
    return find_gRNA_targets(pam_gRNA, genome_file)
def target_coverage_analysis(gRNA_file, extracted_genome_list):
    gRNA_records = SeqIO.parse(gRNA_file, "fasta")
    gRNA_list = [(record.id, str(record.seq)) for record in gRNA_records]
    total_genomes = len(extracted_genome_list)
    with Pool(processes=max_processes) as pool:
        pool.map(generate_fai, extracted_genome_list)
    results = []
    full_coverage_gRNAs = []  
    for gRNA_id, gRNA_seq in gRNA_list:
        with Pool(processes=max_processes) as pool:
            matches = pool.map(process_genome, [(gRNA_seq, genome_file) for genome_file in extracted_genome_list])
        successful_targets = sum(matches) 
        coverage = successful_targets / total_genomes
        result = (gRNA_id, gRNA_seq, successful_targets, total_genomes, coverage)
        results.append(result)
        if f"{coverage*100:.2f}%" == "100.00%":
            full_coverage_gRNAs.append(result)
    full_coverage_dict = {gRNA_id: (gRNA_seq, successful_targets, total_genomes, coverage) 
                          for gRNA_id, gRNA_seq, successful_targets, total_genomes, coverage in full_coverage_gRNAs}
    return full_coverage_gRNAs, full_coverage_dict 
def generate_pam_added_fasta(full_coverage_gRNAs, specified_PAM): 
    output_file = f"full_coverage_gRNAs_with_PAM_{specified_PAM}.fasta"
    records = []
    for gRNA_info in full_coverage_gRNAs:
        gRNA_id, gRNA_seq = gRNA_info[0], gRNA_info[1]
        pam_gRNA_seq = (specified_PAM + gRNA_seq).upper()
        record = SeqRecord(
            Seq(pam_gRNA_seq),
            id=gRNA_id,
            description=f"Full coverage gRNA with PAM {specified_PAM}"
        )
        records.append(record)
    SeqIO.write(records, output_file, "fasta")
    return os.path.abspath(output_file)
def convert_to_uppercase(input_file, output_file):
    with open(output_file, 'w') as f_out:
        for record in SeqIO.parse(input_file, "fasta"):
            record.seq = record.seq.upper()
            SeqIO.write(record, f_out, "fasta")
    return os.path.abspath(output_file)
def build_index(reference_path):
    reference_dir = os.path.dirname(reference_path)
    tmp_dir = os.path.join(reference_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    database_name = os.path.splitext(os.path.basename(reference_path))[0] + "_db"
    database_path = os.path.join(reference_dir, database_name)
    cmd = f"java -Xmx{memory} -jar {flashfry_path} index --tmpLocation {tmp_dir} --database {database_path} --reference {reference_path} --enzyme cpf1"
    subprocess.run(cmd, shell=True, check=True)
    shutil.rmtree(tmp_dir)
    return os.path.abspath(database_path)
def run_offtarget_discovery(database_name, grna_PAM_file, output_file):
    cmd = f"java -Xmx{memory} -jar {flashfry_path} discover --database {database_name} --fasta {grna_PAM_file} --output {output_file} --maximumOffTargets {maximumOffTargets}"
    subprocess.run(cmd, shell=True, check=True)
    return os.path.abspath(output_file)
def find_mismatches(target, offtarget):
    mismatches = []
    for i in range(len(target)):
        if target[i] != offtarget[i]:
            mismatches.append(i + 1)  
    return mismatches
def filter_off_targets(input_file_path, output_file_path, offtarget_count_file):
    offtarget_count = {} 
    with open(input_file_path, 'r') as infile, open(output_file_path, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t') 
        writer.writerow(['gRNA ID', 'Target Sequence with PAM', 'Off-target Sequence with PAM', 'Mismatch Positions'])
        next(infile) 
        for line in infile:
            columns = line.strip().split('\t')
            if len(columns) < 9:
                continue
            gRNA = columns[0]
            target = columns[3]
            if columns[8].strip():
                offtargets = columns[8].split(',')
            else:
                offtargets = []
            for seq in offtargets:
                if seq: 
                    match = re.search(r'_(\d+)_\d+$', seq) 
                    if match:
                        offtarget_count_value = int(match.group(1))  
                    cleaned_seq = re.sub(r'_\d+_\d+$', '', seq)
                    mismatches = find_mismatches(target, cleaned_seq)
                    filtered_mismatches = [m for m in mismatches if m != 4]
                    if not filtered_mismatches:  
                        offtarget_count[gRNA] = offtarget_count.get(gRNA, 0) + offtarget_count_value  
                    else:
                        if filtered_mismatches and min(filtered_mismatches) > 10:  
                            offtarget_count[gRNA] = offtarget_count.get(gRNA, 0) + offtarget_count_value  
                    writer.writerow([gRNA, target, cleaned_seq, ",".join(map(str, mismatches))])
    with open(offtarget_count_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(['gRNA ID', 'Potential Off-target Sites'])
        for gRNA, count in offtarget_count.items():
            writer.writerow([gRNA, count])
    with open(offtarget_count_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        existing_gRNAs = set(row[0] for i, row in enumerate(reader) if i > 0)  
    with open(offtarget_count_file, 'a', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        with open(input_file_path, 'r') as infile:
            reader = csv.reader(infile, delimiter='\t')
            next(reader) 
            for line in reader:
                gRNA = line[0]
                if gRNA not in existing_gRNAs:
                    if line[-4] == 'OVERFLOW':
                        offtarget_count = 'uncertain'
                    else:
                        offtarget_count = 0
                    writer.writerow([gRNA, offtarget_count])
    return os.path.abspath(offtarget_count_file)
def merge_off_target_and_coverage(offtarget_count_file, full_coverage_dict):
    with open(offtarget_count_file, 'r') as infile, open('off_target_and_target_coverage.tsv', 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        next(reader)
        new_header = ['gRNA ID', 'gRNA Sequence', 'Number of Targeted Genomes', 'Total Genomes Analyzed', 'Targeting Coverage (%)', 'Potential Off-target Sites']
        writer.writerow(new_header)
        for row in reader:
            gRNA_id = row[0]
            match = re.match(r'^(gRNA\d+)', gRNA_id)
            if match:
                basic_gRNA_id = match.group(1)  
                matched = False
                for key in full_coverage_dict.keys():
                    if basic_gRNA_id == key:
                        matched = True
                        gRNA_seq, successful_targets, total_genomes, coverage = full_coverage_dict[key]
                        new_row = [basic_gRNA_id, gRNA_seq, successful_targets, total_genomes, f"{coverage * 100:.2f}%", row[1]]
                        writer.writerow(new_row)
                        break 
                if not matched: 
                    new_row = [basic_gRNA_id, 'N/A', 'N/A', 'N/A', 'N/A', row[1]]
                    writer.writerow(new_row)
            else:
                new_row = ['N/A', 'N/A', 'N/A', 'N/A', 'N/A', row[1]]
                writer.writerow(new_row)
def print_completion_message():
    current_time = datetime.datetime.now()
    print("\n" + "="*80)
    print("[INFO] Execution completed successfully!")
    print(f"[INFO] Completion time: {current_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("\n[CREDITS] This software was developed by:")
    print("         Yong Sheng, Hengyu Wang, Xuhu Mao, Qian Li,")
    print("         Jingmin Yan, Chenglong Rao and Ling Deng")
    print("         from Army Medical University")
    print("\n[VERSION] Software version: 1.0.0")
    print("[CONTACT] For support, please contact: yongsheng@tmmu.edu.cn")
    print("="*80 + "\n")
def serialize_and_store_data(obj, name): 
    folder_name = "pickle_data"
    os.makedirs(folder_name, exist_ok=True) 
    file_name = f"{folder_name}/{name}.pkl"
    try:
        with open(file_name, "wb") as f: 
            pickle.dump(obj, f) 
        print(f"[INFO] {name} was stored successfully in the folder {folder_name}.")
    except OSError:
        print(f"[ERROR] Failed to store {name} in the folder {folder_name}")
def load_serialized_data(file_name):
    file_path = os.path.join("pickle_data", f"{file_name}.pkl")
    try:
        with open(file_path, "rb") as f: 
            obj = pickle.load(f) 
        return obj
    except FileNotFoundError:
        print(f"[ERROR] File {file_name} not found")
        return None
    except pickle.UnpicklingError:
        print(f"[ERROR] Failed to load data from {file_name}")
        return None
def main():
    print(f"[INFO] Start time: {datetime.datetime.now()}")
    start_time = time.time()
    global core_gene_folder
    if core_gene_folder is None:
        global input_path
        global gff_file
        if gff_file and os.path.exists(gff_file):
            print(f"[INFO] Using existing GFF files: {gff_file}")
            print("[INFO] Skipping genome annotation as GFF files are provided.")
            gff_dir = gff_file
            gff_files = [f for f in os.listdir(gff_dir) if f.lower().endswith('.gff')]
            num_genomes = len(gff_files)
        else:
            if input_path and os.path.exists(input_path):
                print(f"[INFO] Using existing genomes from: {input_path}")
            else:
                input_path = download_genomes(genus)
            extracted_file_list, num_genomes = extract_gz_files(genus, input_path)
            print("[INFO] Starting genome annotation...")
            annotate_genomes_sequential(extracted_file_list)
            gff_dir = get_gff_files(input_path, f"{genus}_all_gff")
            end_time = time.time()
            execution_time = end_time - start_time
            print(f"[INFO] Genome annotation took {execution_time} seconds.")
        panaroo_output_dir, combined_DNA_CDS = pan_genome_analysis(gff_dir, 16)
        core_gene_folder = extract_core_genes(panaroo_output_dir, num_genomes)
        extract_fasta_from_csv(panaroo_output_dir, core_gene_folder, combined_DNA_CDS)
    fasta_extensions = ('.fa', '.fasta', '.fna') 
    fasta_files = [f for f in os.listdir(core_gene_folder) if f.lower().endswith(fasta_extensions)]
    if not fasta_files:
        raise ValueError(f"No fasta files found in the core gene folder: {core_gene_folder}")
    first_fasta_file = os.path.join(core_gene_folder, fasta_files[0])
    num_genomes = sum(1 for _ in SeqIO.parse(first_fasta_file, "fasta")) 
    start_time = time.time()
    print("[INFO] Starting crRNA design...")
    design_crRNA(core_gene_folder)
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"crRNA design completed. Total time consumed: {execution_time} seconds.")
    start_time = time.time()
    print("[INFO] Starting web scraping...")
    web_crawler_in_parallel(core_gene_folder, firefox_driver_path)
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Web scraping completed. Total time consumed: {execution_time} seconds.")
    initial_filtered_gRNA = summarize_crawler_results(core_gene_folder) 
    extracted_genome_list = [os.path.abspath(os.path.join(input_path, file)) for file in os.listdir(input_path) if file.endswith('.fna')]
    start_time = time.time()
    print("[INFO] Initiating gRNA target coverage analysis...")
    full_coverage_gRNAs, full_coverage_dict = target_coverage_analysis(initial_filtered_gRNA, extracted_genome_list)
    os.chdir(core_gene_folder)
    serialize_and_store_data(full_coverage_gRNAs, "full_coverage_gRNAs")
    serialize_and_store_data(full_coverage_dict, "full_coverage_dict")
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"[INFO] Target coverage analysis completed successfully. Total time consumed: {execution_time} seconds.")
    start_time = time.time()
    print("[INFO] Initiating off-target effect evaluation for filtered gRNAs...")
    print(f"[INFO] Evaluating {len(full_coverage_gRNAs)} gRNA sequences against the genome...")
    full_coverage_gRNAs_with_specific_PAM = generate_pam_added_fasta(full_coverage_gRNAs, "TTTA") 
    offgenomes_upper = convert_to_uppercase(offgenomes, "offgenomes_upper.fasta") 
    offgenomes_db = build_index(offgenomes_upper) 
    serialize_and_store_data(offgenomes_db, "offgenomes_db")
    os.makedirs(os.path.join(core_gene_folder, "gRNA_targeting_analysis"), exist_ok=True)
    os.chdir(os.path.join(core_gene_folder, "gRNA_targeting_analysis"))
    output_off_target_analysis = run_offtarget_discovery(offgenomes_db, full_coverage_gRNAs_with_specific_PAM, "output_off_target_analysis.tsv") 
    offtarget_count_file = filter_off_targets(output_off_target_analysis, "off_target_info.txt", "off_target_count.txt") 
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"[INFO] Off-target effect evaluation completed successfully. Total time consumed: {execution_time} seconds.")
    start_time = time.time()
    print("[INFO] Initiating merger of target coverage and off-target information...")
    merge_off_target_and_coverage(offtarget_count_file, full_coverage_dict)
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"[INFO] Merger of target coverage and off-target information completed successfully. Total time consumed: {execution_time} seconds.")
if __name__ == "__main__":
    main()
    print_completion_message()