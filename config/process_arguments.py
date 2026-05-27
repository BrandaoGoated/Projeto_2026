import yaml
import json
from workflow.pathing_utils.path_generator import check_results_directory
import logging

logger = logging.getLogger(__name__)

def get_arguments(args: dict, sequences: list) -> dict:
    """Converts the arguments given by the CLI to a dictionary"""
    arguments = {
        "seqids": sequences,
        "input": args.input,
        "input_number": 0 if args.input == None else len(args.input),
        "database": args.database,
        "input_file": None if sequences == [] else args.input.split("/")[-1],
        "input_file_db_const": args.input_seqs_db_const,
        "input_type_db_const": args.input_type_db_const,
        "consensus": args.consensus,
        "KEGG_ID": args.kegg,
        "InterPro_ID": args.interpro,
        "hmm_database_name": args.hmm_db_name,
        "alignment_method": args.align_method.lower(),
        "msa_aligner": args.aligner,
        "input_type": args.input_type,
        "metagenomic": True if args.input_type == "metagenome" else False,
        "hmm_validation": args.hmm_validation,
        "expansion": args.expansion,
        "concat_models": args.concat_hmm_models,
        "output_directory": args.output,
        "out_table_format": args.output_type,
        "hmmsearch_out_type": args.hmms_output_type,
        "threads": args.threads,
        "workflow": args.workflow,
        "thresholds": [*range(60, 91, 5)] if args.expansion else False,
        "verbose": args.verbose,
        "overwrite": args.overwrite,
        "split_files": args.split_files,
        "use_cache": args.use_cache,
        "fetch_sra": args.sra,
        "curated": args.curated,
        "snakefile": args.snakefile,
        "unlock": args.unlock,
        "config_file": args.config_file
    }
    return arguments


def check_input_arguments_for_proceding(config: dict, kma_res: bool) -> bool:
    """Checks wether to continue the execution of the parent function by the given arguments. If to continue, return False

    Args:
        config (dict): list of arguments in dict from config
        kma_res (bool): comes from KMA execution

    Returns:
        bool: False if parent has to continue 
    """
    if config.get("hmm_validation") == True and config.get("workflow") == "annotation" and config.get("input") == None:
        if config.get("verbose"):
            print("No input file detected. Proceding to validation")
        return False
    
    # elif config.get("workflow") == "database_construction" and config.get("input") == None and config.get("kegg") == None and config.get("interpro") == None and config.get("input_seqs_db_const") == None:
    #     if config.get("verbose"):
    #         print("No input file detected. Proceding to model construction")
    #     return False
    
    elif config.get("input_type") == "metagenome" and kma_res == False:
        return False
    
    else: return True


def read_config(filename: str) -> tuple:

    config_type = filename.split(".")[-1]
    if config_type == "yaml":
        with open(filename) as stream:
            try:
                config_file = yaml.safe_load(stream)
                stream.close()
            except yaml.YAMLError as exc:
                logger.exception(exc)
    elif config_type == "json":
        with open(filename) as stream:
            try:
                config_file == json.load(stream)
                stream.close()
            except json.decoder.JSONDecodeError as exc:
                logger.exception(exc)
    else:
        quit("Config file must be in .yaml or .json format! Get an example config file from ./config folder.")
    return config_file, config_type


def build_config_from_args(args) -> dict:
    """Converts argparse Namespace directly to the config dict, no file I/O."""
    # seq_ids = _resolve_seq_ids(args)
    seq_ids = []
    return get_arguments(args, seq_ids)


def resolve_config(args):
    """Entry point for config file: returns a config dict from either a provided file or CLI args."""
    if args.config_file is not None:
        config, _ = read_config(args.config_file)
    else:
        config = build_config_from_args(args)
    return config

    
def write_yaml_json(config_type: str, args_dict: dict):
    config_filename = "config"
    check_results_directory(args_dict.output)
    if config_type == "yaml":
        with open(f'{args_dict.output}/{config_filename}.yaml', "w") as file:
            yaml.dump(args_dict, file)
            file.close()
    else:
        with open(f'{args_dict.output}/{config_filename}.json', "w") as file:
            document = json.dumps(args_dict)
            file.write(document)
            file.close()