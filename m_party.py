#!/usr/bin/env python
# run tool main script without indicating python

"""
M-PARTY - Mining Protein dAtasets for Target EnzYmes

by José Freitas

Dec 2026
"""

import sys
import shutil
# sys.path.insert(0, f'{"/".join(sys.path[0].split("/")[:-1])}/share')
sys.path.append(f'{sys.path[0]}/workflow/scripts')
sys.path.append(f'{sys.path[0]}/workflow/pathing_utils')
# sys.path.append(f'{sys.path[0]}/M-PARTY')
import os
from pathlib import Path
import time
import yaml
import pandas as pd
from tqdm import tqdm 
import snakemake
import itertools
import threading
import logging

from workflow.pathing_utils.cli_args import get_parser, process_arguments
from workflow.scripts.hmmsearch_run import run_hmmsearch
from workflow.scripts.hmm_process import *
from workflow.scripts.hmm_vali import concat_final_model, file_generator, exec_testing, hmm_filtration, remove_fp_models, make_paths_dic, delete_inter_files
import workflow.scripts.UPIMAPI_parser as UPIMAPI_parser
from workflow.scripts.seq_download import get_fasta_sequences
from workflow.scripts.CDHIT_seq_download import fasta_retriever_from_cdhit
import workflow.scripts.CDHIT_parser as CDHIT_parser
from workflow.scripts.mparty_util import build_upi_query_db, threshold2clusters, get_tsv_files, save_as_tsv, concat_code_hmm, compress_fasta, return_fasta_content, check_id, check_db_existance
import workflow.scripts.BLAST_parser as BLAST_parser
import workflow.scripts.DIAMOND_parser as DIAMOND_parser
from workflow.scripts.command_run import run_tcoffee, run_hmmbuild, run_hmmemit, concat_fasta, run_sra_download, download_sra_robust
from workflow.scripts.InterPro_retriever import get_IP_sequences
from workflow.scripts.KEGG_retriever import get_kegg_genes
from workflow.scripts.KMA_parser import run_KMA, kma_parser, get_hit_sequences
from config.process_arguments import get_arguments, write_yaml_json, resolve_config
import workflow.scripts.output_scripts.table_report_utils as table_report_utils
import workflow.scripts.output_scripts.text_report_utils as text_report_utils
from workflow.scripts.FASTA_processing import parse_fasta, clean_sequence_ids
from workflow.pathing_utils.fixed_paths import PathManager, declare_fixed_paths
from workflow.pathing_utils.path_generator import dir_generator_from_list, check_results_directory, file_generator


def setup_logging(verbose: bool = False, log_file: Path | None = None, db_name: str | None= None):
    level = logging.DEBUG if verbose else logging.INFO

    handlers = [
        logging.StreamHandler(sys.stdout)
    ]

    if log_file:
        Path(PathManager.log_path).mkdir(parents = True, exist_ok = True)
        if db_name:
            out_path = PathManager.log_path / f'{log_file}_{db_name}'
        else:
            out_path = PathManager.log_path / log_file
        handlers.append(logging.FileHandler(out_path))

    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
        handlers=handlers
    )

logger = logging.getLogger(__name__)


def table_report(
        args: dict,
        dataframe: pd.DataFrame, 
        path: str, 
        type_format: str, 
        db_name: str
    ):
    """Saves a table in a user specified format, with the processed and filtered information from the 
    hmmsearch execution with the HMMs against the query sequences.

    Args:
        args (dict): list of arguments from argparse
        dataframe (pd.DataFrame): Dataframe with only the relevant information from hmmsearch execution 
        for all hmm from all threshold ranges.
        path (str): output path.
        type_format (str): Specify the output format.
        db_name(str): Name of the databases name.

    Raises:
        TypeError: Raises TypeError error if user gives an unsupported output format.
    """
    summary_dic = table_report_utils.create_summary_dict(dataframe=dataframe)
    
    if args.expansion:
        indexes = dataframe.index.values.tolist()
        for i in range(len(indexes)):
            summary_dic["models"][i] = f'{indexes[i]}_{summary_dic["models"][i]}'
    
    df = pd.DataFrame.from_dict(summary_dic)
    list_ids_permodel = {}
    if not args.expansion:
            mother_seqs = f'{sys.path[0]}/resources/Data/FASTA/{db_name}/CDHIT/clusters/'
            for model in tqdm(list(set(summary_dic["models"])), desc = "Tracebacking model's sequences", unit = "model"):
                for file in file_generator(mother_seqs):
                    if file.split(".")[0] == model:
                        if model not in list_ids_permodel:
                            list_ids_permodel[model] = parse_fasta(os.path.join(mother_seqs, file))
                            break
    else:
        mother_seqs = f'{sys.path[0]}/resources/Data/FASTA/{db_name}/CDHIT/'
        for val in summary_dic["models"]:
            thresh = val.split("_")[0]
            model = val.split("_")[-1]
            for folder in os.listdir(mother_seqs):
                if os.path.isdir(os.path.join(mother_seqs, folder)) and folder == thresh:
                    for file in os.listdir(os.path.join(mother_seqs, folder)):
                        if file.endswith(".fasta") and model == file.split(".")[0]:
                            key = thresh + "_" + model
                            if key not in list_ids_permodel:
                                list_ids_permodel[key] = parse_fasta(mother_seqs + folder + "/" + file)
                            # else:
                            #     list_ids_permodel[key].append(parse_fasta(mother_seqs + folder + "/" + file, meta_gen = True if args.input_type == "metagenome" else False))
    
    table_name = "report_table." + type_format
    table_report_utils.check_output(type=type_format, outdir=path, table_name=table_name, dataframe=df, ids_per_model=list_ids_permodel)


def text_report(
        dataframe: pd.DataFrame, 
        path: str, 
        bit_threshold: float, 
        eval_threshold: float, 
        vali: bool = False, 
        kma: bool = False
    ):
    """Write the final report as .txt file, with a summary of the results from the annotation 
    performed with hmmsearch. Starts by calculating the number of in-built HMM profiles, and gives an insight of the 
    filtration thresholds.

    Args:
        dataframe (pd.DataFrame): Dataframe with only the relevant information from hmmsearch execution 
        for all hmm from all threshold ranges.
        path (str): output path.
    """
    # number of initial HMM profiles
    number_init_hmms, number_validated_hmms = 0, 0
    for dir in os.listdir(PathManager.hmm_database_path):
        if os.path.isdir(os.path.join(PathManager.hmm_database_path, dir)):
            for _ in os.listdir(os.path.join(PathManager.hmm_database_path, dir)):
                number_init_hmms += 1

    if vali:
        for dir in os.listdir(PathManager.validated_hmm_dir):
            if os.path.isdir(os.path.join(PathManager.validated_hmm_dir, dir)):
                for _ in os.listdir(os.path.join(PathManager.validated_hmm_dir, dir)):
                    number_validated_hmms += 1

    # get the IDs from all hits after quality check
    query_names = get_match_ids(dataframe, to_list = True, only_relevant = True)

    # get number of hits given for each sequence
    number_hits_perseq = get_number_hits_perseq(query_names)

    # get the unique sequences
    unique_seqs = get_unique_hits(query_names)
    inputed_seqs = config["seqids"]
    
    variables = text_report_utils.write_var_file()
    text_report_utils.write_text_report(config, path, args, variables)


def get_number_hits_perseq(hit_ids_list: list) -> dict:
    """Given a list of sequences IDs from the hits against the hmm models from hmmsearch, counts the number of each ID.

    Args:
        hit_ids_list (list): List of sequence IDs.

    Returns:
        dict: Dictionary containing each ID as key and the respective number of occurrences as value.
    """
    counter = {}
    for i in hit_ids_list:
        counter[i] = counter.get(i, 0) + 1
    return counter


def get_unique_hits(hit_ids_list: list) -> list:
    """Given a list of sequence IDs from the hits against the hmm models from hmmsearch, return a new list with only the unique elements.

    Args:
        hit_IDs_list (list): List of sequence IDs.

    Returns:
        list: List with only a single occurrence of each ID.
    """
    unique_ids_list = []
    for x in hit_ids_list:
        if x not in unique_ids_list:
            unique_ids_list.append(x)
    return unique_ids_list


def get_aligned_seqs(
        config, 
        hit_ids_list: list, 
        path: str, 
        inputed_seqs: str, 
        kma_alignfile: str = None
    ):
    """Writes an ouput Fasta file with the sequences from the input files that had a hit in hmmsearch 
    annotation against the hmm models.

    Args:
        hit_ids_list (list): list of IDs that hit.
        path (str): ouput path.
        inputed_seqs (str): name of the initial input file.
    """
    # returns a list the sequences that hit against the models (only one entry)
    unique_ids = get_unique_hits(hit_ids_list)

    if config.get("seqids") == "too_big":
        check_id(inputed_seqs, path, unique_ids)
    
    else:
        with open(path + "aligned.fasta", "w") as wf:
            if config.get("input_type") == "metagenome":
                input_ids = parse_fasta(kma_alignfile, remove_excess_id = False, kma_res = True)
                inp_seqs = kma_alignfile
            else:
                input_ids = parse_fasta(inputed_seqs, remove_excess_id = False)
                inp_seqs = inputed_seqs
                
            with open(inp_seqs, "r") as rf:
                lines = rf.readlines()
                for x in unique_ids:
                    if x in input_ids:
                        iterador = iter(lines)
                        linha = next(iterador)
                        while linha is not None:
                            if x not in linha:
                                linha = next(iterador, None)
                                continue
                            elif x in linha:
                                wf.write(linha)
                                linha = next(iterador, None)
                                while linha is not None and not linha.startswith(">"):
                                    wf.write(linha)
                                    linha = next(iterador, None)
                            elif x not in linha and linha.startswith(">"):
                                break
                            linha = next(iterador, None)
                    else:
                        continue
            rf.close()
        wf.close()


def generate_output_files(
        dataframe: pd.DataFrame, 
        hit_ids_list: list, 
        inputed_seqs: str,
        config: dict,
        bit_threshold: float = None, 
        eval_threshold: float = None, 
        kma: bool = False,
        kma_alignfile: str = None
    ):
    """Function that initializes the output files creation simultaneously, for now, only two files are generated:
    report and aligned sequences.
    Path will always be the output folder defined by the user when running tool in CLI, so no pat argument is required.

    Args:
        args (dict): list of argumets from argparse
        dataframe (pd.DataFrame): Dataframe with only the relevant information from hmmsearch execution.
        hit_ids_list (list): list of Uniprot IDs that hit.
        inputed_seqs (str): name of the initial input file.
    """
    out_folder = config.get("output") + "/"
    if kma:
        get_aligned_seqs(config, hit_ids_list, out_folder, inputed_seqs, kma_alignfile = kma_alignfile)
        dataframe.to_excel(f'{out_folder}report_table.xlsx', sheet_name = "Table_Report", index = 0)
    else:
        table_report(config, dataframe, out_folder, config.get("output_type"), config.get("hmm_database_name"))
        if config.get("report_text"):
            if config.get("hmm_validation"):
                text_report(dataframe, out_folder, bit_threshold, eval_threshold, vali = True)
            else:
                text_report(dataframe, out_folder, bit_threshold, eval_threshold)
        get_aligned_seqs(config, hit_ids_list, out_folder, inputed_seqs)


def fetch_sra(config):
    """Function to fetch SRA files from Sequence Reads Archive

    Args:
        config (file): The parsed config file object

    Raises:
        ValueError: If input type is metagenome at the same time as the interpro flag is given with an ID
    """
    print("Starting download of SRA files...\n")
    time.sleep(1)
    
    if check_db_existance(config):
        
        if config.get("split_files"):
            print("Samples will be splitted into forward and reverse reads files...\n")

        if not config.get("use_cache"):
            print('Downloading directly to file. If taking to much time, turn on "--use_cache"...\n')

            for accession in tqdm(config.get("fetch_sra")):
                run_sra_download(
                    accession, 
                    str(PathManager.sra_fastq_path),
                    config.get("split_files"),
                    config.get("verbose")
                )
        else:
            for accession in tqdm(config.get("fetch_sra")):
                download_sra_robust(
                    accession, 
                    str(PathManager.sra_fastq_path),
                    config.get("split_files"),
                    config.get("verbose")
                )
            
    return True


def database_construction(config):
    """Pipeline for the database construction workflow

    Args:
        config (file): The parsed config file object

    Raises:
        ValueError: If input type is metagenome at the same time as the interpro flag is given with an ID
    """
    print("HMM database construction workflow from user input started...\n")
    time.sleep(1)
    
    if check_db_existance(config):

        time.sleep(2)

        if config.get("expansion"):
            expand_base_sequences(config=config)

        else:
            # make necessary directories
            dir_generator_from_list(
                [
                    PathManager.tcoffee_path, 
                    PathManager.cdhit_path / "clusters", 
                    PathManager.hmm_database_path, 
                    PathManager.sra_fastq_path
                ]
            )
            
            if config.get("KEGG_ID"):
                # if given ID is Kegg Orthology
                if config.get("KEGG_ID")[0].startswith("K"):
                    kegg_sequences_path = get_kegg_genes(
                        PathManager.fasta_type_dir / Path(config.get("KEGG_ID")[0]).with_suffix(".fasta"), 
                        type_seq = "nuc" if config.get("input_type_db_const") == "nucleic" else "AA",
                        ko = config.get("KEGG_ID"), 
                        verbose = config.get("verbose")
                    )

                # If given ID is an E.C. number
                else:
                    kegg_sequences_path = get_kegg_genes(
                        PathManager.fasta_type_dir / Path(config.get("KEGG_ID")[0]).with_suffix(".fasta"), 
                        type_seq = "nuc" if config.get("input_type_db_const") == "nucleic" else "AA",
                        ec_number = config.get("KEGG_ID"), 
                        verbose = config.get("verbose")
                    )

                # Only build HMMs if input is protein or nucleic
                if config.get("input_type") != "metagenome":
                    build_hmms_from_seqs(config, kegg_sequences_path, "KEGG")

            if config.get("InterPro_ID"):
                # for interpro is only possible to run for aminoacids and so for HMM and not KMA and raw metagenomes
                if config.get("input_type") == "metagenome":
                    raise ValueError("Metagenomic samples cannot be annalyzed with proteins as database")

                # if given ID is a InterProt ID
                elif config.get("InterPro_ID")[0].startswith("IPR") and len(config.get("InterPro_ID")) == 1:
                    inp_seqs_path = get_IP_sequences(
                        PathManager.fasta_type_dir / Path(config.get("InterPro_ID")[0]).with_suffix(".fasta"), 
                        interpro_ID = config.get("InterPro_ID"), 
                        reviewed = config.get("curated"), 
                        verbose = config.get("verbose")
                    )

                # if given ID is a list of proteins from InterProt
                elif config.get("InterPro_ID")[0].startswith("A"):
                    inp_seqs_path = get_IP_sequences(
                        PathManager.fasta_type_dir / Path(config.get("InterPro_ID")).with_suffix(".fasta"), 
                        protein = config.get("InterPro_ID"), 
                        verbose = config.get("verbose")
                    )

                # Start HMM construction
                build_hmms_from_seqs(config, inp_seqs_path, "InP", ident_perc=0.8)

            # if a FASTA file with interest proteins/nucleiotides is given
            if config.get("input_file_db_const"):
                # Will not build HMMs if input is a metagenome
                
                if config.get("input_type") == "metagenome":
                    # instead, copy the file to the same output FASTA dir
                    shutil.copyfile(config.get("input_file_db_const"), PathManager.fasta_type_dir)

                else:
                    # Start HMM construction
                    shutil.copyfile(config.get("input_file_db_const"), PathManager.fasta_type_dir / config.get("input_file_db_const").split("/")[-1].split(".")[0])
                    build_hmms_from_seqs(
                        config,
                        sequences_path=config.get("input_file_db_const"),
                        from_database=config.get("input_file_db_const").split("/")[-1].split(".")[0]
                    )

            # remove files wrongly going to the root dir
            files = [f for f in os.listdir('.') if os.path.isfile(f)]
            for file in files:
                if file.endswith(".dnd"):
                    delete_inter_files(file)

    if config.get("hmm_validation"):
        validate_hmm(config=config)


def expand_base_sequences(config):
    Path("resources/Data/FASTA/DataBases").mkdir(parents = True, exist_ok = True)
    Path(f'resources/Data/Tables/{config.get("hmm_database_name")}').mkdir(parents = True, exist_ok = True)
    query_db = build_upi_query_db("resources/Data/FASTA/DataBases", config = config, verbose = config["verbose"])

    if config["alignment_method"] == "diamond":
        ### FASTA to DMND
        diamond_file = DIAMOND_parser.build_diamond_DB(query_db, "resources/Data/FASTA/", verbose = config["verbose"])  # ver a cena do overwrite para estes passos
        Path(f'resources/Alignments/{config.get("hmm_database_name")}/BLAST/diamond_output/').mkdir(parents = True, exist_ok = True)
        aligned_tsv = DIAMOND_parser.run_DIAMOND(config.get("input_file_db_const"), f'resources/Alignments/{config.get("hmm_database_name")}/{config["alignment_method"].upper()}/diamond_output/out.tsv', diamond_file, config.get("threads"))
        handle = DIAMOND_parser(aligned_tsv)
        dic_enzymes = DIAMOND_parser.DIAMOND_iter_per_sim(handle)
        if config["verbose"]:
            print(f'Saving IDs from the ranges of {config["thresholds"]} percentages of similarity.\n')
        save_as_tsv(dic_enzymes, f'resources/Data/Tables/{config.get("hmm_database_name")}/DIAMOND_results_per_sim.tsv')

    elif config["alignment_method"] == "upimapi":
        # aligned_TSV = run_UPIMAPI(query_DB, f'resources/Alignments/{args.hmm_db_name}/{config["alignment_method"].upper()}/upimapi_results', args.input_seqs_db_const, args.threads)
        aligned_tsv = f'resources/Alignments/{config.get("hmm_database_name")}/{config["alignment_method"].upper()}/upimapi_results/UPIMAPI_results.tsv'
        handle = UPIMAPI_parser.UPIMAPI_parser(aligned_tsv)
        dic_enzymes = UPIMAPI_parser.UPIMAPI_iter_per_sim(handle)
        if config["verbose"]:
            print(f'Saving IDs for the minimum cutoff values of {config["thresholds"]} percentages of similarity.\n')
        save_as_tsv(dic_enzymes, f'resources/Data/Tables/{config.get("hmm_database_name")}/UPIMAPI_results_per_sim.tsv')

    elif config["alignment_method"] == "blast":
        # blastdb_file = build_blast_DB(query_DB, "resources/Data/FASTA/DataBases/BLAST", args.input_type_db_const, verbose = config["verbose"])
        Path(f'resources/Alignments/{config.get("hmm_database_name")}/BLAST/BLAST_results').mkdir(parents = True, exist_ok = True)
        # run_BLAST(args.input_seqs_db_const, f'resources/Alignments/{args.hmm_db_name}/BLAST/BLAST_results/test.tsv', blastdb_file, 8)
        aligned_tsv = f'resources/Alignments/{config.get("hmm_database_name")}/BLAST/BLAST_results/test.tsv'
        handle = BLAST_parser(aligned_tsv)
        dic_enzymes = BLAST_parser.BLAST_iter_per_sim(handle)
        if config["verbose"]:
            print(f'Saving IDs from the ranges of {config["thresholds"]} percentages of similarity.\n')
        Path(f'resources/Data/Tables/{config.get("hmm_database_name")}/').mkdir(parents = True, exist_ok = True)
        save_as_tsv(dic_enzymes, f'resources/Data/Tables/{config.get("hmm_database_name")}/BLAST_results_per_sim.tsv')

    else:
        raise ValueError("--align_method flag only ranges from 'diamond', 'upimapi' or 'blast'. Chose one from the list.")

    Path(f'resources/Data/FASTA/{config.get("hmm_database_name")}/{config["alignment_method"].upper()}/').mkdir(parents = True, exist_ok = True)
    Path(f'resources/Data/Tables/{config.get("hmm_database_name")}/CDHIT_clusters/').mkdir(parents = True, exist_ok = True)
    for thresh in config["thresholds"]:
        if config["verbose"]:
            print(f'Retrieving sequences from {thresh} range\n')
        try:
            get_fasta_sequences(f'resources/Data/Tables/{config.get("hmm_database_name")}/{config["alignment_method"].upper()}_results_per_sim.tsv', f'resources/Data/FASTA/{config.get("hmm_database_name")}/{config["alignment_method"].upper()}/{thresh}.fasta')
        except Exception as exc:
            print(exc)
            raise FileNotFoundError(f'resources/Data/Tables/{config["alignment_method"].upper()} not found.')
        ### run CDHIT
        if config["verbose"]:
            print(f'CDHIT run for {thresh} range\n')
            Path(f'resources/Data/FASTA/{config.get("hmm_database_name")}/CDHIT/{thresh}/').mkdir(parents = True, exist_ok = True)
        try:
            CDHIT_parser.run_CDHIT(f'resources/Data/FASTA/{config.get("hmm_database_name")}/{config["alignment_method"].upper()}/{thresh}.fasta', f'resources/Data/FASTA/{config.get("hmm_database_name")}/CDHIT/cd-hit_after_{config["alignment_method"]}_{thresh}.fasta', 8)
            handle = CDHIT_parser.cdhit_parser(f'resources/Data/FASTA/{config.get("hmm_database_name")}/CDHIT/cd-hit_after_{config["alignment_method"]}_{thresh}.fasta.clstr')
            handle2 = counter(handle, tsv_ready = True, remove_duplicates = True)
            save_as_tsv(handle2, f'resources/Data/Tables/{config.get("hmm_database_name")}/CDHIT_clusters/cdhit_clusters_{thresh}_after{config["alignment_method"]}.tsv')

            if config["verbose"]:
                print("Retrieving sequences divided by clusters from CDHIT\n")
            fasta_retriever_from_cdhit(f'resources/Data/Tables/{config.get("hmm_database_name")}/CDHIT_clusters/cdhit_clusters_{thresh}_after{config["alignment_method"]}.tsv', 
                                        f'resources/Data/FASTA/{config.get("hmm_database_name")}/CDHIT/{thresh}')
        except Exception as exc:
            print(exc)
            if config["verbose"]:
                print(f'[WARNING] Minimum cutoff of {thresh} of similarity not detected.\n')
            time.sleep(2)
            continue

    ### add cluster per threshol to config
    os.remove("config/config.yaml")

    if config_format == "yaml":
        files = get_tsv_files(config)
        threshandclust = threshold2clusters(files)
        print(threshandclust)
        for thresh, cluster in threshandclust.items():
            for c in range(len(cluster)):
                cluster[c] = str(cluster[c])
            config[thresh] = cluster
        newthresh = []
        for thresh in config["thresholds"]:
            if thresh not in threshandclust:
                continue
            else:
                newthresh.append(thresh)
        config["thresholds"] = newthresh

    with open("config/config.yaml", "w") as dump_file:
        yaml.dump(config, dump_file)
        dump_file.close()

    snakemake.main(
        f'-s {config.get("snakefile")} --printshellcmds --cores {config["threads"]} --configfile config/{config.get("config_file")}'
        f'{" --unlock" if config.get("unlock") else ""}')

    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    for file in files:
        if file.endswith(".dnd"):
            delete_inter_files(file)

    print("HMM database created!")
    time.sleep(2)


def build_hmms_from_seqs(
        config: dict,
        sequences_path: str,
        from_database: str,
        ident_perc: float = 0.7,
    ):
    """Sequence that condenses the code responsible for building the hmm from the sequences received from the user.

    Args:
        args (dict): list of arguments from argparse
        sequences_path (str): a path or filename for the sequences to be used for hmm building. Can be from KEGG, InterPro or a custom sequence file
        from_database (str): The database the sequences come from
        ident_perc (float, optional): The identity threshold to be used for CDHIT clustering. Defaults to 0.7
    """
    if from_database == "KEGG" or from_database == "InP":
        filename = sequences_path.stem
    else:
        filename = sequences_path.split("/")[-1].split(".")[0]
    
    # generate cluster file
    CDHIT_parser.run_CDHIT(
        sequences_path, 
        PathManager.cdhit_path / Path(filename).with_suffix(".fasta"), 
        config.get("threads"), 
        identperc=ident_perc,
        type_seq = "AA" if config.get("input_type_db_const") == "protein" else "NUC", 
    )

    # get cluster info -> {cluster_X: [list of IDs]}
    seqs = CDHIT_parser.cdhit_parser(
                PathManager.cdhit_path / Path(filename).with_suffix(".fasta.clstr"), 
                kegg = True if from_database == "KEGG" else False,
                ip = True if from_database == "InP" else False
            )

    # get FASTA IDs from the sequences file that will be used to build the HMMs
    input_ids = parse_fasta(
        sequences_path, 
        kegg = True if from_database == "KEGG" else False, 
        config = config
    )
    
    # generate FASTAS with each cluster sequence
    CDHIT_parser.get_clustered_sequences(
        seqs, 
        PathManager.cdhit_path / "clusters", 
        sequences_path, 
        input_ids, 
        from_database
    )

    # Run multiple sequence alignement with T-COFFEE
    for file in os.listdir(PathManager.cdhit_path / "clusters"):
        try:
            run_tcoffee(
                PathManager.cdhit_path / "clusters" / file,
                PathManager.tcoffee_path / Path(file.split(".")[0]).with_suffix(".clustal_aln"),
                type_seq = "PROTEIN" if config.get("input_type_db_const") == "protein" else "DNA",
                verbose = config.get("verbose")
            )

        except Exception as exc:
            logger.warning("T-COFFEE failed for file %s: %s", file, exc)
            if config.get("verbose"):
                logger.exception("Full traceback")

            continue

    # Build HMMs with msa files
    for msa in os.listdir(PathManager.tcoffee_path):
        run_hmmbuild(
            PathManager.tcoffee_path / msa, 
            PathManager.hmm_database_path / Path(msa.split(".")[0]).with_suffix(".hmm"), 
            # args.verbose,
            True,
            PathManager.hmm_database_path / Path(msa.split(".")[0]).with_suffix(".txt"), 
        )

        # Get consensus sequence for KMA
        if config.get("consensus"):
            run_hmmemit(
                PathManager.hmm_database_path / Path(msa.split(".")[0]).with_suffix(".hmm"), 
                PathManager.consensus_path / Path(msa.split(".")[0]).with_suffix(".fasta")
            )

            # get consensus sequence
            concat_fasta(
                PathManager.consensus_path, 
                PathManager.consensus_path / "consensus"
            )

    # concat hmm models to a single file
    concat_code_hmm(config.get("hmm_database_name"), from_database + "_model")


def validate_hmm(config):
    """Executes the pipeline for the HMMs validation

    Args:
        config (str): The config file
    """
    print("Starting HMM validation procedures...")
    time.sleep(2)

    pathing = make_paths_dic(args.hmm_db_name)
    exec_testing(thresholds = config["thresholds"], path_dictionary = pathing, database = args.negative_db)
    to_remove = hmm_filtration(pathing)
    remove_fp_models(to_remove, pathing)
    concat_final_model(pathing)
    time.sleep(2)
    print("M-PARTY has concluded model validation! Will now switch to the newlly created models (in the validated_HMM folder\n")


def annotation(config):

    print("Annotation workflow started...\n")
    time.sleep(2)

    if config.get("hmm_validation"):

        if not os.path.exists(PathManager.validated_hmm_dir):
            validate_hmm(config=config)
        else:
            print("Validated HMM already up, proceding to annotation...\n")
            time.sleep(1)

    # if a metagenome is given, runs KMA 
    if config.get("input_type") == "metagenome":
        # paths are hardcoded for convinience
        dir_generator_from_list(
            [
                PathManager.tables_path / 'kma_hits', 
                PathManager.databases_path / 'kma_db' / 'KEGG_cons'
            ]
        )

        paired_workflow, second_input = False, None
        if config.get("input_number") == 2:
            paired_workflow = True
            second_input = config.get("input")[1]
        
        # if consensus is requested, go get the sequence file
        if config.get("consensus"):
            print(config.get("input"))
            kma_out = run_KMA(
                PathManager.consensus_path / Path("consensus").with_suffix(".fasta"), 
                PathManager.databases_path / 'kma_db',
                config.get("input"), 
                PathManager.tables_path / 'kma_hits' / config.get("input")[0].split(".")[0], 
                threads = config.get("threads"),
                paired_end=paired_workflow,
                second_input=second_input
            )

        # otherwise, run for the sequence file in the requested hmm_db_name
        else:
            for file in os.listdir(PathManager.fasta_type_dir):
                if os.path.isfile(os.path.join(PathManager.fasta_type_dir, file)):
                    kma_out = run_KMA(
                        PathManager.fasta_type_dir / file, 
                        PathManager.databases_path / 'kma_db',
                        config.get("input"), 
                        PathManager.tables_path / 'kma_hits' / Path(config.get("input").split("/")[-1].split(".")[0]), 
                        threads = config.get("threads"),
                        paired_end=paired_workflow,
                        second_input=second_input
                    )
                    
        df = kma_parser(kma_out.with_suffix(".res"), )
        hit_seqs = get_hit_sequences(df, to_list = True)
        generate_output_files(
            df, 
            hit_seqs, 
            kma_out, 
            config, 
            kma = True, 
            kma_alignfile = kma_out.with_suffix(".fsa")
        )

    # if input file is not a metagenome
    else:
        if config.get("hmm_validation"):
            for hmm_file in file_generator(PathManager.validated_hmm_dir, full_path = True):
                hmmserach_out_file = Path(f'search_{config["input_file"].split("/")[-1].split(".")[0]}_{hmm_file.split("/")[-1].split(".")[0]}').with_suffix("." + args.hmms_output_type)
                run_hmmsearch(
                    config.get("input"), 
                    hmm_file,
                    PathManager.hmmsearch_results_path / hmmserach_out_file,
                    verbose = config.get("verbose"), 
                    eval = 0.00001,
                    out_type = config.get("hmmsearch_out_type")
                )
        else:
        # if models have been concatenated
            if config.get("concat_models"):
                print(PathManager.hmm_database_path)
                for hmm_file in file_generator(PathManager.hmm_database_path / "concat_model", full_path = True):
                    hmmserach_out_file = Path(f'search_{config["input_file"].split("/")[-1].split(".")[0]}_{hmm_file.split("/")[-1].split(".")[0]}').with_suffix("." + args.hmms_output_type)
                    if os.path.exists(PathManager.hmmsearch_results_path / hmmserach_out_file):
                        os.remove(PathManager.hmmsearch_results_path / hmmserach_out_file)
                        run_hmmsearch(
                            config.get("input"), 
                            hmm_file,
                            PathManager.hmmsearch_results_path / hmmserach_out_file,
                            verbose = config.get("verbose"), 
                            eval = 0.00001,
                            out_type = config.get("hmmsearch_out_type")
                        )
                    else:
                        PathManager.hmmsearch_results_path.mkdir(parents = True, exist_ok = True)
                        run_hmmsearch(
                            config.get("input"), 
                            hmm_file,
                            PathManager.hmmsearch_results_path / hmmserach_out_file,
                            verbose = config.get("verbose"), 
                            eval = 0.00001,
                            out_type = config.get("hmmsearch_out_type")
                        )
            else:
                p = os.listdir(PathManager.hmm_database_path)
                for thresh in p:
                    path = os.path.join(PathManager.hmm_database_path, thresh)
                    Path(path).mkdir(parents = True, exist_ok = True)
                    hmmserach_out_file = Path(f'search_{config["input_file"].split("/")[-1].split(".")[0]}_{hmm_file.split("/")[-1].split(".")[0]}').with_suffix("." + args.hmms_output_type)
                    for hmm_file in file_generator(path, full_path = True):
                        run_hmmsearch(
                            config.get("input"), 
                            hmm_file, 
                            path / hmmserach_out_file, 
                            verbose = config.get("verbose"), 
                            eval = 0.00001,
                            out_type = config.get("hmmsearch_out_type")
                        )
                    BLAST_parser.concat_hmmsearch_results(path, PathManager.hmmsearch_results_path)

        if config.get("expansion"):
            lista_dataframes = dict.fromkeys(config["thresholds"])
            for file in file_generator(PathManager.hmmsearch_results_path):
                thresh = file.split("_")[-1].split(".")[0]
                lista_dataframes[thresh] = read_hmmsearch_table(PathManager.hmmsearch_results_path + file)

            final_df = concat_df_byrow(df_dict = lista_dataframes)
            rel_df = relevant_info_df(final_df)
            quality_df, bs_thresh, eval_thresh = quality_check(rel_df, give_params = True)
            hited_seqs = get_match_ids(quality_df, to_list = True, only_relevant = True)
            

        else:
            for file in file_generator(PathManager.hmmsearch_results_path):
                if config.get("input").split("/")[-1].split(".")[0] in file:
                    dataframe = read_hmmsearch_table(PathManager.hmmsearch_results_path / file)
            rel_df = relevant_info_df(dataframe)
            quality_df, bs_thresh, eval_thresh = quality_check(rel_df, give_params = True)
            hited_seqs = get_match_ids(quality_df, to_list = True, only_relevant = True)
        
        # outout files are always generated
        generate_output_files(quality_df, hited_seqs, config.get("input"), config, bs_thresh, eval_thresh)


def clean(args: dict):
    paths = ["Data/FASTA/", "Alignments/", "Data/HMMs/"]
    for path in paths:
        path = f'resources/{path}{args.hmm_db_name}/'
        if Path(path).exists():
            shutil.rmtree(path)
            if args.verbose:
                logger.info(f'Deleted path {path}')


def main_pipeline(args):

    ### Resolve config ###
    config = resolve_config(args)
    # optionally save it for reproducibility
    if args.display_config or not args.config_file:
        write_yaml_json("yaml", args)

    ### Clean database and exit ###
    if args.clean:
        clean(args)
        return

    done = False
    if args.verbose and args.input is not None:
        def animate():
            for c in itertools.cycle(['|', '/', '-', '\\']):
                if done:
                    break
                sys.stdout.write('\rParsing input sequences IDs: ' + c)
                sys.stdout.flush()
                time.sleep(0.1)
            time.sleep(0.5)
            sys.stdout.flush()
            sys.stdout.write('\rParsing input sequences IDs: Done!\n')

        t = threading.Thread(target=animate)
        t.start()

    done = True
    time.sleep(1)

    st = time.time()

    ### VALIDATION ###
    # first only runs for if user flags --validation alone without input sequences, will validate the models inside database only
    if args.hmm_validation and args.workflow != "database_construction" and args.workflow != "both" and args.input == None:

        validate_hmm(config=config)

    ### ANNOTATION ###
    # runs if input sequences are given
    if args.workflow == "annotation" and args.input is not None:
        annotation(config)
    
    elif args.workflow == "fetch":
        fetch_sra(config)

    ### DATABASE CONSTRUCTION ###
    elif args.workflow == "database_construction":
        database_construction(config=config)

    ### DATABASE CONSTRUCTION + ANNOTATION ###
    elif args.workflow == "both":
        database_construction(config=config)
        annotation(config)



    et = time.time()
    elapsed_time = et - st
    miliseconds_time = elapsed_time * 1000
    minutes_time = (elapsed_time / 60)
    print(f'Execution time: {miliseconds_time:.4f} milliseconds and {minutes_time:.2f} minutes!')
    if args.workflow == "database_construction":
        if args.consensus:
            print("Consensus sequences generated!")
        else:
            print("HMMs generated!")
            print("Next step should be annotation with the just created HMM database.\nIf you need further guidance, refer to M-PARTY documentation in GitHub.")
    else:
        print(f'M-PARTY has stoped running! Results are displayed in the {args.output} folder :)')
    print("Thank you for using M-PARTY! ")


def main():
    # get CLI arguments
    parser = get_parser()
    args = parser.parse_args()

    # check arguments
    process_arguments(args)
    
    # initialize paths
    declare_fixed_paths(args)

    # setup logger
    setup_logging(verbose=args.verbose, log_file=args.log_file, db_name=args.hmm_db_name)

    # start pipeline
    main_pipeline(args)

if __name__ == "__main__":
    main()