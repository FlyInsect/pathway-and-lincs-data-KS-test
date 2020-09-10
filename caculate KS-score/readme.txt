#KS test was performed on the differential gene expression value of each treatment (excluding DMSO) in the LINCS data and the gene expression value on the pathway, and the KS score value was calculated.  

1. Calculate the score of each processing (each column) of LINCS data  

2. Combine all columns.  

3. After normalization of all drug treatment scores.  

※ Sort by |score| from largest to smallest. Contains compounds that may treat and exacerbate the disease.  

※Only take score<0 data and sort by |score| from large to small. Contains only drugs that may be treated.  

For multiple processing data of a drug (drug_name), take the column with the smallest score value, and sort the drugs according to the score value from small to large.
