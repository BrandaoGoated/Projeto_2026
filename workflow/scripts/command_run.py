import logging
from pathlib import Path
import subprocess
from subprocess import run
import sys

logger = logging.getLogger(__name__)

def run_command(bash_command: list[str], output='', mode='w', verbose=True):
    display = ' '.join(bash_command) + (f' > {output}' if output else '')
    logger.debug(display)

    try:
        if output == '':
            run(bash_command, stdout=sys.stdout if verbose else None, check=True)
        else:
            with open(output, mode) as output_file:
                run(bash_command, stdout=output_file, check=True)  # check=True was missing here!
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with exit code {e.returncode}: {display}")
        sys.exit(1)
    except FileNotFoundError as e:
        logger.error(f"Executable not found: {e}")
        sys.exit(1)


def docker_run_tcoffee(volume, input_file, output_type, output_name):
    run_command(f'docker`run`--rm`-v`{volume}`pegi3s/tcoffee:latest`t_coffee`/data/{input_file}`-run_name`/data/{output_name}`-output`{output_type}', sep="`")

def docker_run_hmmbuild(volume, input_file, output_file):
    run_command(f'docker`run`--rm`-v`{volume}`biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1`hmmbuild`{output_file}`{input_file}', sep="`")

def docker_run_hmmsearch(volume, hmm_file, db, output_file):
    run_command(f'docker`run`--rm`-v`{volume}`biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1`hmmsearch`{hmm_file}`{db}`>`{output_file}', sep="`")

def run_tcoffee(input: str, output: str, type_seq: str = "PROTEIN", verbose: bool = False):
    command = ["t_coffee", str(input), "-output", "clustalw_aln", "-outfile", str(output), "-type", type_seq, "-quiet", "stderr" if verbose else "", "-n_core", "4"]
    # run_command(f't_coffee`{input}`-output`clustalw_aln`-outfile`{output}`-type`{type_seq}`-quiet`{"stderr" if verbose else ""}`-n_core`4', sep = "`")
    run_command(command)

def run_hmmbuild(input: str, output: str, verbose: bool = False, stdout_path: str = None):
    if not verbose and stdout_path == "":
        raise ValueError("Tem que adicionar um caminho para o ficheiro de output de hmmbuild")
    message = ""
    if not verbose:
        message = f'`-o`{stdout_path}'
    
    command = ["hmmbuild", str(output), str(input)]
    # run_command(f'hmmbuild`{output}`{input}{message}', sep = "`")
    run_command(command)

def run_hmmemit(input: str, output: str):
    command = ["hmmemit", "-o", str(output), str(input)]
    # run_command(f'hmmemit`-o`{output}`{input}', sep = "`")
    run_command(command)

def concat_hmm(input_path: str, output_path: str):
    command = ["cat", f"{str(input_path)}*.hmm", ">", f"{str(output_path)}.hmm"]
    # run_command(f'cat`{input_path}*.hmm`>`{output_path}.hmm', sep = "`")
    run_command(command)

def run_sra_download(accession: str, output_dir: str, split_files: bool, verbose: bool):
    command = ["fasterq-dump", accession, "--outdir", output_dir, "--threads", "4"]
    if split_files:
        command.append("--split-files")
    if verbose:
        command.append("--progress")

    run_command(command)

def download_sra_robust(accession: str, output_dir: str, split_files: bool, verbose: bool):
    lock_file = Path(f"/app/{accession}/{accession}.sra.lock")
    
    if lock_file.exists():
        logger.warning(f"Stale lock file found for {accession}, removing it...")
        lock_file.unlink()
        
    cache_path = Path.home() / "ncbi" / "public" / "sra" / f"{accession}.sra"
    
    try:
        logger.info(f"Fetching {accession} from SRA to cache...")
        command = ["prefetch", accession]
        if verbose:
            command.append("--progress")
        run_command(command)
    except SystemExit:
        if cache_path.exists():
            logger.warning(f"Prefetch interrupted, cleaning up partial cache for {accession}")
            cache_path.unlink()
        raise
    
    logger.info(f"Converting {accession} to FASTQ...")
    run_sra_download(accession, output_dir, split_files, verbose)
    logger.info(f"Done: {accession}")

def concat_fasta(input_path: str, output_path: str):
    run_command(f'cat`{input_path}*.fasta`>`{output_path}.fasta', sep = "`")