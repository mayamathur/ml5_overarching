

# Additional analyses for ML5

*Contact: Maya B. Mathur*

## Methods and results
For each study, we conducted analyses for three subsets: (1) all replications, regardless of which protocol they used; (2) only replications using the RP:P protocol; and (3) only replications using the Revised protocol. For each subset, we meta-analyzed the replications using methods that do not use any distributional assumptions and that provide valid inference even for small sample sizes (Hedges et al., 2010; Tipton, 2015). In the two cases in which the model failed to converge (noted in the column `robu.error` in the results file `results_table_full.csv`), we instead used parametric meta-analysis fit via REML (Viechtbauer, 2010). We conducted all analyses on the Fisher's $z$ scale and transformed all results back to the Pearson's $r$ scale for reporting.

We used robust methoods to estimate $P_{\text{orig}}$, $\widehat{P}_{>0}$, $\widehat{P}_{>0.10}$, and $\widehat{P}_{>0.20}$, where the thresholds for the latter three are on the Pearson's $r$ scale (Mathur & VanderWeele, in press; Mathur & VanderWeele, 2020). We also estimated the probability that the pooled replication estimate would be "statistically significant" and positive in sign, accounting for heterogeneity, if in fact it were drawn from the same distribution as the original study (`Psignif.agree`). The file `results_aggregated_by_subset.csv` aggregates these results by subset (all replications, RP:P, Revised), showing the mean heterogeneity estimate $\widehat{\tau}$, the percentage of studies with $\widehat{\tau}>0$, the mean values of $P_{\text{orig}}$, $\widehat{P}_{>0}$, $\widehat{P}_{>0.10}$, and $\widehat{P}_{>0.20}$, and . For replication sets that had a heterogeneity estimate $\widehat{\tau}=0$ exactly or which had only 1 replication, we simply report the $\widehat{P}_{>q}$ as 100% or 0% depending on whether the single point estimate was above or below $q$. Confidence intervals for $\widehat{P}_{>q}$ are omitted when they could not be estimated with BCa bootstrapping, when $\widehat{\tau}=0$ exactly, or when there was only 1 replication. 

Results of these analyses for all subsets are in `results_table_pretty.csv` and its more verbose version, `results_table_full.csv`. Results stratified by subset are in `results_table_pretty_combined.csv`, `results_table_pretty_rpp.csv`, and `results_table_pretty_revised.csv`. 

## Discussion
Interestingly, despite our close standardization of protocols across sites, 40% of replication sets within each of the 3 subsets had $\widehat{\tau}>0$. However, the above analyses accounting for heterogeneity still yield a strikingly clear story of non-replication. For example, even in the Revised subset, we estimated that only 30% of true population effects were above the modest effect size of $r=0.10$, and only 10% were above $r=0.20$. 

