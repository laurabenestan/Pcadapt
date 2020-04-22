# pcaAdapt

---

PCAadapt is a R package R package that performs genome scans to detect genes potentially under divergent selection based on a principal component analysis. To read more about the package see (Lu and Blum)[https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12592]

---
## 0. Prepare your data
---

** Prepare a .bed file** from your vcf file.
To do so, first use [VCFTOOLS](http://vcftools.sourceforge.net) in the terminal.
```{r, engine = 'bash', eval = FALSE}
vcftools --vcf nameofyourfile.vcf --plink-tped --out nameofyourfile
```

Then, use and [PLINK](http://zzz.bwh.harvard.edu/plink/).
```{r, engine = 'bash', eval = FALSE}
plink --tped nameofyourfile.tped --tfam nameofyourfile.tfam --make-bed --out nameofyourfile
```

The commend `--make-bed` will produce three files:
- a binary ped file (*.bed)
- the pedigree/phenotype information file (*.fam)
- an extended MAP file (*.bim) that contains information about the allele names, which would otherwise be lost in the .bed file

---
## 1. Download R package and input dataset
---

Download libraries. 
```{r}
library("pcadapt") 
library("vcfR")
library("plyr")
library('qvalue')
```

Import vcf. 
```{r}
data <- read.pcadapt("batch_1.bed", type = "bed")
# number of individuals detected:	44
# number of loci detected:		12735
```

---
## 2. Define the highest signal of genomic variation
---

Run pcadapt function by first perform it with a large enough number of principal components (e.g. K=20) in order to maximize your capacity of detecting any population structure in your dataset.
```{r}
data_pcadapt <- pcadapt(data, K = 20, ploidy=2) 
```

Visualize the genomic variation among your samples. 
Are you samples uniformly distributed or does some points cluster together?
```{r}
plot(data_pcadapt,option="scores",i=1,j=2)
```

Check the percent of variance explained by each principal component and select the minimum number of K.
```{r}
plot(data_pcadapt, option = "screeplot", K=5)
```

Here, we choose K=2 since, the highest signal of genomic variation is between two clusters. 
```{r}
pcaadapt_K2 <- pcadapt(data, K = 4)
```

---
## 3. Identify markers driving the highest genomic variation observed
---

After choosing the righ number of K to select, compute the test statistic based on a Principal Component Analysis.

Do a Mahattan plot on the P-values.
```{r}
plot(pcaadapt_K2, option = "manhattan", col, snp.info = NULL, plt.pkg = "ggplot")
```

Check the expected uniform distribution of the p-values using a Q-Q plot.
```{r}
plot(data_pcadapt, option = "qqplot")
```

Draw an histogram of the p-values: the excess of small p-values indicates the presence of outliers.
```{r}
hist(data_pcadapt$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
```

Choose a cutoff for outlier detection with the **Q-values**.
```{r}
qval <- qvalue(data_pcadapt$pvalues)$qvalues
alpha <- 0.01
outliers <- which(qval < alpha)
length(outliers)
```

Choose a cutoff for outlier detection with the **Bonferroni correction**.
```{r}
padj <- p.adjust(data_pcadapt$pvalues,method="bonferroni")
alpha <- 0.000000000000001
outliers <- which(padj < alpha)
length(outliers)
```

Choose a cutoff for outlier detection with the **Benjamini-Hochberg Procedure**.
```{r}
padj <- p.adjust(data_pcadapt$pvalues,method="BH")
alpha <- 0.01
outliers <- which(padj < alpha)
length(outliers)
```
Save the results.
```{r}
write.table(outliers, "386Outliers-bh.txt", sep="\t", quote=FALSE, row.names=FALSE) 
```
