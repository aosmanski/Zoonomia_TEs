# Zoonomia_TEs
Scripts and pipelines associated with the submitted publication: Osmanski et al. 2021. Insights into mammalian TE diversity via the curation of 200+ mammalian genome assemblies. Science: submitted. 



### template_RM_slurmC8.sh
Pipeline used to run RepeatMasker on all zoonomia assemblies


### template_rm2bed.sh
Pipleline used to convert repeatmasker output to .bed format. This pipeline merges overlapping hits based using 'lower_divergence' criterion.


### Flagship_Paper_Figure.R
Workflow used to create the Figure S1 in Christmas et al. 2021. Evolutionary constraint and innovation across hundreds of placental mammals. Science. 


### stacked_bar_circular_phylogeny.R
Workflow used to create figure 2. A phylogeny with stacked bar charts as tree tips depicting the proportions of various TE types found within the respective assembly. 


### TE_Assembly.R
Creates ultrametric phylogeny used for subsequent hierarchical Bayesian analyses. 


#### bat_tree.R
Workflow used to create Figure 4. A phylogeny with stacked bar chart tree tips depicting the TE content of Chiropterans. 


### TE_diversity.R
Runs the hierarchical Bayesian analyses behind figure 5.


### plot_TE_div_ass.R
Workflow used to plot figure 5. The recent mamammlian TE diversity using both Shannon's H & Pielou's J


### TE_diet.R
Runs the hierarchical Bayesian analyses behind figure 6.


### plot_diet.R
Workflow used to plot figure 6. Plot depicting recent DNA transposon accumulation among three dietary phenotypes: carnivore, herbivore, & omnivore. 
