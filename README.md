# pcaAdapt

---

PCAadapt is a R package R package that performs genome scans to detect genes potentially under divergent selection based on a principal component analysis. To read more about the package see (Lu and Blum)[https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12592]

---
## 0. Prepare your data
---

If you use a very recent the 4.3 version of vcftools, you need to **prepare a .bed file** from your vcf file to use pcadapt.
To do so, first use [VCFTOOLS](http://vcftools.sourceforge.net) in the terminal.
```{r, engine = 'bash', eval = FALSE}
vcftools --vcf nameofyourfile.vcf --plink-tped --out nameofyourfile
```

Then, use and [PLINK](http://zzz.bwh.harvard.edu/plink/).
```{r, engine = 'bash', eval = FALSE}
plink --tped nameofyourfile.tped --tfam nameofyourfile.tfam --make-bed --out nameofyourfile
```

The command `--make-bed` will produce three files:
- a binary ped file (*.bed)
- the pedigree/phenotype information file (*.fam)
- an extended MAP file (*.bim) that contains information about the allele names, which would otherwise be lost in the .bed file

Now extract the information of the order of each SNP in the vcf file. 
***Model organism***. If you have info on chromosome or scaffold position, you may want to use the first and second column of your vcf (i.e. CHROM POS).
***Non-model organism***. If you do not have this information, you may need to extract the third colum of your vcf (i.e. ID).
```{r, engine = 'bash', eval = FALSE}
grep -v "#" nameofyourfile.vcf | cut -f 3 > yournumberofsnps.txt
```
---
## 1. Download R package and input dataset
---

Download libraries. 
```{r}
library("pcadapt") 
library("vcfR")
library("plyr")
library('qvalue')
library("dplyr")
```

Import vcf or .bed (in this case check the info above).
```{r}
data <- read.pcadapt("batch_1.bed", type = "bed")
# number of individuals detected:	44
# number of loci detected:		12735
```

Import in which order individuals order are labelled in the .vcf using the .tfam file.
```{r}
ind <- read.table("44ind.tfam")
```

Import important biological information.
```{r}
ind_pop_mpa <- read.table("mpa_info_44ind.txt", header=TRUE, sep="\t")
```

Merge individuals and pop information
```{r}
ind_mpa <- merge(ind, ind_pop_mpa, by.x=c("V1"), by.y=c("IND"))
ind_mpa_info <- select(ind_mpa, V1,LAT,LON,DISTANCE,MPA,CATEGORY)
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

Create a dataframe gathering PCAadapt clustering information and MPA info
```{r}
pca_adapt_mpa <-cbind(data_pcadapt$scores, ind_mpa) 
```

Visualising the results according to MPA and Inside/Outside
```{r}
ggplot(pca_adapt_mpa, aes(x=pca_adapt_mpa$`1`, y=pca_adapt_mpa$`2`, shape = CATEGORY, fill= factor(MPA)))+
         geom_point(size=1.5)+
  scale_shape_manual(values=c(21, 24))+
  scale_fill_brewer(palette="Accent", guide=FALSE)+
  facet_wrap(~MPA)+
  theme_classic()+
  xlab("PC1")+
  ylab("PC2")+ 
  theme_bw()
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

Match the SNP with the Pvalues information
```{r}
snp_pvalue <- cbind(snp, data_pcadapt$pvalues) 
colnames(snp_pvalue) <- c("SNP","PVALUES")
```
By default, the parameter min.maf is set to 5%. P-values of SNPs with a minor allele frequency smaller than the threshold are not computed (NA is returned). Check you have no NA.
```{r}
sum(is.na(snp_pvalue$PVALUES))
```

Choose a cutoff for outlier detection with the **Q-values**.
```{r}
qval <- qvalue(snp_pvalue$PVALUES)$qvalues
alpha <- 0.01
outliers <- which(qval < alpha)
length(outliers)
```

Choose a cutoff for outlier detection with the **Bonferroni correction**.
```{r}
padj <- p.adjust(snp_pvalue$PVALUES,method="bonferroni")
alpha <- 0.00000004
outliers <- which(padj < alpha)
length(outliers)
```

Choose a cutoff for outlier detection with the **Benjamini-Hochberg Procedure**.
```{r}
padj <- p.adjust(snp_pvalue$PVALUES,method="BH")
alpha <- 0.01
outliers <- which(padj < alpha)
length(outliers)
```
Save the results.
```{r}
write.table(snp_pvalue, "677outliers-bh.txt", sep="\t", quote=FALSE, row.names=FALSE) 
```

## 4. Identify top outliers

Visualize the distribution of p-values
```{r}
quantile(snp_pvalue$PVALUES, probs = c(0.01, 0.99))
```

Get only the top 1% of the markers.
```{r}
top_1percent <- subset(snps_pvalue, snp_pvalue$PVALUES <= 4.823761e-39)
write.table(top_1percent, "Outliers.txt", sep="\t", quote=FALSE, row.names = FALSE)
```

