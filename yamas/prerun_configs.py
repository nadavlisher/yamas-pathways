import os

from .utilities import run_cmd

command_based_on_os = {'Ubuntu': ['wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz'],
                       'CentOS': ['wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz'],
                       }

# Define required packages
pip_packages = ['rpy2', 'bs4', 'wget', 'ete3']
conda_packages = ['bioconda::bioconductor-rsubread', 'bioconda::bioconductor-genomicranges', 'bioconda::bioconductor-genomicalignments', 'bioconda::bioconductor-rtracklayer', 'r::r-jsonlite', 'bioconda::bioconductor-clusterprofiler=3.16']  # Modify as needed


def install_pip_packages():
    for package in pip_packages:
        run_cmd(['pip', 'install', package])


def install_conda_packages():
    for package in conda_packages:
        run_cmd(['conda', 'install', '-y', package])


def set_environment(operating_system_type):

    #defualt. it will update soon based on your os
    default_path= 'export PATH=$PATH:$PWD/sratoolkit.3.0.7-centos_linux64/bin'
    version="sratoolkit.3.0.7-centos_linux64"

    command = command_based_on_os[operating_system_type]
    print("######   Downloading SRA Toolkit")
    run_cmd(command)

    print("######   Extracting the contents of the tar file   ")
    run_cmd(['tar -vxzf sratoolkit.tar.gz > sra_toolkit_download_log.txt'])
    with open('sra_toolkit_download_log.txt', 'r') as file:
        last_line = file.readlines()[-1]
        version = last_line.split("/")[0]

    print("######   Installing pip packages   ")
    install_pip_packages()

    print("######   Installing conda packages   ")
    install_conda_packages()

    print(f"######   export PATH   ")
    # run_cmd([f'export PATH=$PATH:$PWD/{version}/bin'])
    os.environ['PATH'] = f"{os.environ['PATH']}:{os.getcwd()}/{version}/bin"

    print("######   Checking if we are ready to go...   ")
    run_cmd(['which fastq-dump > check_fastq-dump.txt'])
    with open('check_fastq-dump.txt', 'r') as check_file:
        content= check_file.read()
        if f'{version}/bin/fastq-dump' in content:
            print(f"\n######   Please run the following command in your terminal:   ######\n######   export PATH=$PATH:$PWD/{version}/bin   ######\n######   YAMAS is ready to run!   ###### ")
        else:
            print("######   YAMAS is NOT ready! Try again or you can set the environment by yourself.\n  Follow the instructions (from step 2) that can be found in git: https://github.com/YarinBekor/YaMAS.   ######\n    ")
