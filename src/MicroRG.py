from ftplib import FTP
from pathlib import Path
from tqdm import tqdm
import subprocess
import pandas as pd

class MicroRG():

    def __init__(self, ref="RefSeq"):
        self.query = None
        self.ref = ref.lower()
        self.cwd = Path(".")
        
        # check assmebly report
        self.has_report = False
        self.dir_report = None
        self.check_assembly_report()

        self.results = None

    
    def check_assembly_report(self):

        # check existence
        if not (self.ref == "refseq" or self.ref == "genebank"):
            ValueError("Please provide valid reference database: RefSeq or Genebank")

        self.dir_report = self.cwd / f"assembly_summary_{self.ref}.txt"

        if not self.dir_report.exists():
            print("No assembly summary report.")
            print(f"Download the {self.ref} assembly summary report!")

            subprocess.call(["rsync", "-t", "-v",
            f"rsync://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_{self.ref}.txt", "./"])
        
        self.dir_report = self.cwd / f"assembly_summary_{self.ref}.txt"


    def search(self, query, full_genome=True, strain_specific=False):
        self.query = [item.lower() for item in query]
        n_query = len(query)
        n_found = 0
        n_strain_found = 0
        results = {}
        unfound = []

        with open(self.dir_report, 'r') as report:
            # skip two lines of header
            next(report)
            next(report)

            for line in report:
                line_split = line.split('\t')
                sp, strain, is_full, ftp_link = line_split[7].lower(), line_split[8].lower(), \
                    line_split[13], line_split[19]
                
                # check whether the reference genome is complete
                if full_genome:
                    if is_full != "Full":
                        continue

                # strain=ASM1234
                if "strain" in strain:
                    strain = strain.split('=')[1]

                # if strain specific
                if strain_specific:
                    if strain != "":
                        sp = sp + " " + strain

                if sp in self.query:
                    if sp not in results.keys():
                        n_found += 1
                        n_strain_found += 1
                        results[sp] = [(strain, ftp_link)]
                    else:
                        n_strain_found += 1
                        results[sp].append((strain, ftp_link))

        unfound = list(set(self.query) - set(results.keys()))
        self.results = results
        print(f"Search results: {n_found} / {n_query}")
        print(f"Number of strain: {n_strain_found}")
        print(f"Number of unfound species: {len(unfound)}")
        return results, unfound



    def download_ref_genome(self, dir_ref_output, dir_shell_output=None, as_shell_file=False, query=None):
        if query is not None:
            self.results, _ = self.search(query, strain_specific=False)
        
        # output directory
        dir_ref_output = Path(dir_ref_output)
        dir_fna = dir_ref_output / "fna"
        dir_gff = dir_ref_output / "gff"
        dir_fna.mkdir(exist_ok=True)
        dir_gff.mkdir(exist_ok=True)

        fna_links, gff_links = [], []

        for sp in self.results.keys():
                for strain_info in self.results[sp]:
                    _, file_link = strain_info[0], strain_info[1]
                    file_link = file_link.replace("https", "rsync")
                    folder_name = file_link.split('/')[-1]
                    fna_links.append(file_link + f"/{folder_name}_genomic.fna.gz")
                    gff_links.append(file_link + f"/{folder_name}_genomic.gff.gz")

        # output as shell file
        if as_shell_file:
            if dir_shell_output is None:
                ValueError("Requried a folder directory for the output shell files.")

            # fna files
            with open(dir_shell_output + "download_fna.sh", 'w') as file:
                for link in fna_links:
                    file.write(f"rsync -q {link} {dir_fna}\n")

            # gff
            with open(dir_shell_output + "download_gff.sh", 'w') as file:
                for link in gff_links:
                    file.write(f"rsync -q {link} {dir_gff}\n")
        
        else:
            # download via python which provides the progress bar
            if dir_ref_output.is_dir:
                print("Download fna files")
                for i in tqdm(range(len(fna_links))):
                    subprocess.call(["rsync", "-q", fna_links[i], dir_fna])
            
                print("Download gff files")
                for i in tqdm(range(len(gff_links))):
                    subprocess.call(["rsync", "-q", gff_links[i], dir_gff])

                print("Donload - Complete")
            else:
                FileExistsError("The output directory doesn't exist or isn't a folder.")

    def save_query_results(self, dir_output):
        with open(dir_output + "/query_results.tsv", 'w') as file:
            file.write(f"Species\tStrain\tFolder\tLink\n")
            for sp in self.results.keys():
                for sp in self.results.keys():
                    for strain_info in self.results[sp]:
                        strain, file_link = strain_info[0], strain_info[1]
                        file_link = file_link.replace("https", "rsync")
                        folder_name = file_link.split('/')[-1]
                        file.write(f"{sp}\t{strain}\t{folder_name}\t{file_link}\n")



    def _download_ref_genome(self, dir_ref_output, dir_shell_output=None, as_shell_file=False, query=None):
            if query is not None:
                self.results, _ = self.search(query, strain_specific=False)
            
            # output directory
            dir_ref_output = Path(dir_ref_output)
            dir_fna = dir_ref_output / "fna"
            dir_gff = dir_ref_output / "gff"
            dir_fna.mkdir(exist_ok=True)
            dir_gff.mkdir(exist_ok=True)

            fna_links, gff_links = [], []


            for sp in self.results.keys():
                fna_link = ""
                gff_link = ""

                for strain_info in self.results[sp]:
                    strain, file_link = strain_info[0], strain_info[1]
                    file_link = file_link.replace("https", "rsync")
                    folder_name = file_link.split('/')[-1]
                    fna_link = file_link + f"/{folder_name}_genomic.fna.gz"
                    gff_link = file_link + f"/{folder_name}_genomic.gff.gz"

                    if strain == "":
                        break

                fna_links.append(fna_link)
                gff_links.append(gff_link)
                    

            # output as shell file
            if as_shell_file:
                if dir_shell_output is None:
                    ValueError("Requried a folder directory for the output shell files.")

                # fna files
                with open(dir_shell_output + "download_fna.sh", 'w') as file:
                    for link in fna_links:
                        file.write(f"rsync -q {link} {dir_fna}\n")

                # gff
                with open(dir_shell_output + "download_gff.sh", 'w') as file:
                    for link in gff_links:
                        file.write(f"rsync -q {link} {dir_gff}\n")
            
            else:
                # download via python which provides the progress bar
                if dir_ref_output.is_dir:
                    print("Download fna files")
                    for i in tqdm(range(len(fna_links))):
                        subprocess.call(["rsync", "-q",fna_links[i], dir_fna])
                
                    print("Download gff files")
                    for i in tqdm(range(len(gff_links))):
                        subprocess.call(["rsync", "-q", gff_links[i], dir_gff])

                    print("Donload - Complete")
                else:
                    FileExistsError("The output directory doesn't exist or isn't a folder.")