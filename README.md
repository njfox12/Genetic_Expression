# Genetic_Expression

### Background
This project was a consulting project for my consulting course in my Master's study. This was a project in which I worked with two other individuals from my class to provide statistical consultation for a professor at Loyola University Chicago. We worked with Dr. Heather Wheeler, who works in the computer science department and runs a computational biology lab at the university. We were continuing a previous study she had done, where an elastic net was used to predict gene expression. The paper being referenced was "Survey of Heritability and Sparse Architecture of Gene Expression Traits Across Human Tissues" by Heather Wheeler et. al.

### Data
The data for this project was from a previous study, "Patterns of Cis Regulatory Variation in Diverse Human Populations" by Barbara Stranger et.al. Specifically, the data included SNP (single-nucleotide polymorphism) position and frequency with the corresponding level of expression. In this project, we focused on 107 individuals using chromosome 22, which has a total of 498 genes to be analyzed. 

### Project Goals
The goal of this project was to determine if any non-parametric models demonstrated better prediction than the elastic net model that was created in the previous study. We developed three different models, which built on each other. The first model was a principal components natural splines model (PCA-NS), the second was a k-nearest neighbor principal components natural splines model (KNN-PCA-NS), then the model which I was responsible for was the support vector machine k-nearest neighbor principal component analysis natural splines model (SVM-KNN-PCA-NS).

### References
Stranger, Barbara and Montogomery, Stephen and Dimas, Antigone, et al. (2012) "Patterns of Cis Regulatory Variation in Diverse Human Populations." PLOS Genetics 8(4): e1002639. doi: 10.1371/journal.pgen.1002639

Wheeler HE, Shah KP, Brenner J, Garcia T, Aquino-Michaels K, et al. (2016) "Survey of the Heritability and Sparse Architecture of Gene Expression Traits across Human Tissues." PLOS Genetics 12(11): e1006423. doi: 10.1371/journal.pgen.1006423
