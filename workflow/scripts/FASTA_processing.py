import re
import logging
from time import time
from config.process_arguments import check_input_arguments_for_proceding

logger = logging.getLogger(__name__)

def check_fasta() -> bool:
    return True


def clean_sequence_ids(
        line: str,
        remove_excess_id: bool, 
        ip: bool, 
        kegg: bool, 
        kma_res: bool
    ) -> str:
    """Function that receives a string and cleans it based on predefined patterns

    Args:
        line (str): line string
        remove_excess_id (bool): Decide wether to remove the excess part of UniProt IDs
        ip (bool): Set to True if sequences from filename were retrieved from InterPro, which has a specific nomenclature
        for the FASTA entries
        kegg (bool): Set to True if sequences from filenames were retrieved from KEGG, which has a specific nomenclature
        for the FASTA entries.
        kma_res (bool): Set to True if this function is set to run for the processing of KMA results

    Returns:
        str: Cleaned string
    """
    if kegg:
        return re.search(r">(\S+)", line).group(1)
    elif ip:
        return re.search(r">([^|]+)\|", line).group(1)
    elif kma_res:
        return line.replace(">", "").strip()
    else:
        if not remove_excess_id:
            return line.split(" ")[0][1:]
        else:
            try:
                return re.search(r"\|(.*)\|", line).group(1)
            except Exception:
                identi = line.split(" ")[0]
                return identi.replace(">", "")


def parse_fasta(
        filename: str, 
        remove_excess_id: bool = True, 
        ip: bool = False, 
        kegg: bool = False, 
        kma_res: bool = False, 
        config: dict = None,
    ) -> list:
    """Given a FASTA file, returns the IDs from all sequences in that file. 
    If file not present, program will be quited and TypeError message raised.

    Args:
        filename (str): Name of FASTA file.
        remove_excess_id (bool, optional): Decide wether to remove the excess part of UniProt IDs. Defaults to True.
        ip (bool, optional): Set to True if sequences from filename were retrieved from InterPro, which has a specific nomenclature
        for the FASTA entries.
        kegg (bool, optional): Set to True if sequences from filenames were retrieved from KEGG, which has a specific nomenclature
        for the FASTA entries.
        kma_res (bool, optional): Set to True if this function is set to run for the processing of KMA results. Defaults to False.
        verbose (bool, optional): Set to True to print aditional messages of wath is happening. Defaults to False.

    Returns:
        list: A list containing IDs from all sequences
    """
    uniq_ids = []
    if check_input_arguments_for_proceding(config, kma_res=kma_res) == False:
        return uniq_ids
    else:
        try:
            with open(filename, "r") as handlefile:
                try:
                    for line in handlefile:
                        if line.startswith(">"):
                            uniq_ids.append(clean_sequence_ids(line, remove_excess_id, ip, kegg, kma_res))
                    if config.get("verbose"):
                        print(f'Input file {filename} detected and sequence IDs retrieved\n')
                        time.sleep(1)
                except Exception as exc:
                    logger.warning(exc)
                    quit("File must be in FASTA format.")
        except TypeError:
            raise TypeError("Missing input file! Make sure -i option is filled")
        return uniq_ids

