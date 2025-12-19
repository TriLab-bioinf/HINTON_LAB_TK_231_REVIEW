# ***Escherichia coli*** **Pathobionts and Crohn’s Disease:  Varied genetic paths leading to similar phenotypes**

*Melissa Arroyo-Mendoza, Hernan Lorenzi, Greg Phillips, and Deborah M. Hinton*

[![](https://zenodo.org/badge/DOI/10.5281/zenodo.17992151.svg)](https://doi.org/10.5281/zenodo.17992151)

## Generate protein-conservation plot for review

### Step-1: Identify Best Reciprocal Matches between LF82 proteome and the other 6 proteomes of interest

**Goal:** To identify conserved genes between *E. coli* LF82 and other *E. coli* strains.

#### Approach:

1- Download protein sequences for all genomes of interest from NCBI and place them in the `./ORTHOLOGS/proteins/` folder. The files should be named as follows:

-   LF82.faa
-   541-1.faa
-   541-15.faa
-   576-1.faa
-   MG1655.faa
-   NRG857c.faa
-   T75.faa

``` bash
mkdir -p ./ORTHOLOGS/proteins

for i in $(cat protein_sources.csv); do
    sname=$(echo ${i} | awk -F',' '{print $1}')
    url=$(echo ${i} | awk -F',' '{print $2}')
    echo Downloading proteins for ${sname}
    wget -O - ${url} | gunzip -c > ./ORTHOLOGS/proteins/${sname}.faa
done
```

2- Download annotation files in GFF format for all genomes of interest from NCBI and place them in the `./ORTHOLOGS/annotations/` folder. The files should be named as follows:

-   LF82.gff
-   541-1.gff
-   541-15.gff
-   576-1.gff
-   MG1655.gff
-   NRG857c.gff
-   T75.gff

``` bash
mkdir -p ./ORTHOLOGS/annotations

for i in $(cat gff_sources.csv); do
    sname=$(echo ${i} | awk -F',' '{print $1}')
    url=$(echo ${i} | awk -F',' '{print $2}')
    echo Downloading gff for ${sname}
    wget -O - ${url} | gunzip -c > ./ORTHOLOGS/annotations/${sname}.gff
done
```

3- Make BLASTP library for LF82 proteins (subject).

``` bash
module load blast/2.15.0+
makeblastdb -in ./ORTHOLOGS/proteins/LF82.faa -dbtype prot -out ./ORTHOLOGS/proteins/LF82.faa
```

4- Run BLASTP searched of the LF82 proteome (subject) vs all other proteomes (queries).

``` bash
# Use blast version 2.15.0+
module load blast/2.15.0+

mkdir -p ./ORTHOLOGS/blastp

for i in 541-1 541-15 576-1 MG1655 NRG857c T75; do
    echo Running blastp  of LF82 vs ${i}
    blastp -db ./ORTHOLOGS/proteins/LF82.faa \
        -query ./ORTHOLOGS/proteins/${i}.faa \
        -evalue 1e-10 \
        -out ./ORTHOLOGS/blastp/${i}_vs_LF82.bp \
        -outfmt "7 std qlen slen"
done
```

5- Select best reciprocal matches with `keep_best_hit.py`, filtering out hits with less than 50% coverage and 35% identity at the protein level.

``` bash
mkdir -p ./ORTHOLOGS/best_reciprocal_matches

for i in 541-1 541-15 576-1 MG1655 NRG857c T75; do
    echo Processing ${i}
    python ./keep_best_hit.py -i ./ORTHOLOGS/blastp/${i}_vs_LF82.bp \
    -o ./ORTHOLOGS/best_reciprocal_matches/${i}_vs_LF82.brm
done
```

6- Generate bar-plot for review in R using the following quarto script:

`./R/scripts/review.plot.qmd`
