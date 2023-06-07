#! /bin/bash
#SBATCH -J scheduler
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1  
#SBATCH -o outErr/scheduler.o
#SBATCH --mem=1gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user@domain.edu
#SBATCH -A buyin user 

mkdir -p outErr
rm -r .snakemake/

src/dimorphite.py

snakemake -s Snakefile_site_screen --touch -j1
snakemake -s Snakefile_site_screen --cluster-config cluster.yaml --cluster 'sbatch -t {cluster.time} --mem={cluster.mem} -c {cluster.cpus} {cluster.other} -A merzjrke -o {cluster.out} -e {cluster.err}' -j 499 --restart-times 3 --rerun-incomplete

snakemake -s Snakefile_conf_gen --touch -j1
snakemake -s Snakefile_conf_gen --cluster-config cluster.yaml --cluster 'sbatch -t {cluster.time} --mem={cluster.mem} -c {cluster.cpus} {cluster.other} -A merzjrke -o {cluster.out} -e {cluster.err}' -j 499 --restart-times 3 --rerun-incomplete

#snakemake -s Snakefile_balloon --touch -j1
#snakemake -s Snakefile_balloon --cluster-config cluster.yaml --cluster 'sbatch -t {cluster.time} --mem={cluster.mem} -c {cluster.cpus} {cluster.other} -A merzjrke -o {cluster.out} -e {cluster.err}' -j 499 --restart-times 3 --rerun-incomplete

snakemake -s Snakefile_AG --touch -j1
snakemake -s Snakefile_AG --cluster-config cluster.yaml --cluster 'sbatch -t {cluster.time} --mem={cluster.mem} -c {cluster.cpus} {cluster.other} -A merzjrke -o {cluster.out} -e {cluster.err}' -j 499 --restart-times 3 --scheduler greedy

#snakemake -s Snakefile_fast --touch -j1 
#snakemake -s Snakefile_fast --cluster-config cluster.yaml --cluster 'sbatch -t {cluster.time} --mem={cluster.mem} -c {cluster.cpus} {cluster.other} -A merzjrke -o {cluster.out} -e {cluster.err}' -j 499 --restart-times 3 --rerun-incomplete

snakemake -s Snakefile_geom_opt --touch -j1
snakemake -s Snakefile_geom_opt --cluster-config cluster.yaml --cluster 'sbatch --time {cluster.time} --mem={cluster.mem} -c {cluster.cpus} {cluster.other} -A merzjrke -o {cluster.out} -e {cluster.err}' -j 499 --restart-times 3 --rerun-incomplete
