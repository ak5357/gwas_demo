gwas_demo
================
ak5357
2024-09-27

Genome Wide Association Study. Tutorial located here:
<https://whussain2.github.io/Materials/Teaching/GWAS_R_2.html>

## Read the marker data

When we unzip the file downloaded from ricediversity.org, you will find
three files:

- .fam (contains names about the acessions/genotypes)
- .map (contains map information including chromosome, position and
  marker name)
- .ped (contains the marker allele data in 0, 1, 2, and 3 format; 2
  represents missing data). The files are plink converted files.

We recode the marker to 0 (homozygous, major allele), 1 (heterozygous),
and 2 (homozygous, minor allele). The missing are recoded to `NA`.

Finally, the marker data is converted to a matrix and transposed.

``` r
# Empty environment
rm(list = ls())

# Read the marker data files. Marker data
Geno = read_ped("data/sativas413.ped")
p = Geno$p #number of genotypes (samples)
n = Geno$n #number of markers (loci)
Geno = Geno$x #genotype matrix containing genetic data for those genotypes and markers

# Acession information
FAM = read.table("data/sativas413.fam")

# Map information
MAP = read.table("data/sativas413.map")

# Recode data in ped file
Geno[Geno == 2] = NA #Converting missing data to NA
Geno[Geno == 0] = 0  #Converting 0 data to 0
Geno[Geno == 1] = 1  #Converting 1 to 1
Geno[Geno == 3] = 2  #Converting 3 to 2

# Convert the marker data into matrix and transpose and check dimensions
Geno = matrix(Geno, nrow = p, ncol = n, byrow = TRUE) |> 
  t()

# See Dimensions of Geno
dim(Geno)
```

    ## [1]   413 36901

## Read the phenotypic data

In this section, we read the phenotypic data directly from
ricediversity.org and view the first 5 rows and columns of the data.

``` r
rice_pheno = read.table("http://www.ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt", 
    header = TRUE,
    stringsAsFactors = FALSE,
    sep = "\t")
# See first few columns and rows of the data
rice_pheno[1:5, 1:5]
```

    ##        HybID NSFTVID Flowering.time.at.Arkansas Flowering.time.at.Faridpur
    ## 1 081215-A05       1                   75.08333                         64
    ## 2 081215-A06       3                   89.50000                         66
    ## 3 081215-A07       4                   94.50000                         67
    ## 4 081215-A08       5                   87.50000                         70
    ## 5 090414-A09       6                   89.08333                         73
    ##   Flowering.time.at.Aberdeen
    ## 1                         81
    ## 2                         83
    ## 3                         93
    ## 4                        108
    ## 5                        101

Check the dimensions of the file.

``` r
dim(rice_pheno)
```

    ## [1] 413  38

Rename

``` r
rownames(Geno) = FAM$V2

#Check whether the operation worked
table(rownames(Geno) == rice_pheno$NSFTVID)
```

    ## 
    ## TRUE 
    ##  413
