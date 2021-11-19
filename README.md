# parMix : software tool for parental ancestry inference and genotype calling from multiple children
Software accompanies for "Joint Inference of Ancestry and Genotypes of Parents from a Small Number of Children", by Yiming Zhang and Yufeng Wu, manuscript, 2021. This paper is currently under review.

parMix is a software tool for jointly inferring parental ancestry and calling parental genotypes from data of a small number of children. The inputs for parMix are the phased genotypes of diploid individuals. 

parMix is an extended version of our previous tool PedMix [[1]](#1) (which can estimate the admixture proportion of recent ancestors from a single child). And different from PedMix, parMix can provide fine-scale inference, namely parental ancestry and genotypes at each single nucleotide polymorphism site. In the current version, parMix considers a single diploid family with more than one child, and each child belongs to an admixed population with two ancestral populations A and B. 

The key idea of parMix is constructing some HMMs (Hidden Markov Model) and then using the posterior decoding algorithm to call the ancestry and genotype of parents based on the genotypes of children. 

# Usage
parMix allows two types of inference. 1) Joint inference of parental ancestry and genotypes under the model with phasing errors. 2) Inference of parental ancestry and genotypes under the model without phasing errors. 

The first model infers the ancestry and genotypes for both parents at the same time, but the second infers the ancestry and genotypes for the single parent and runs independently for both parents.

Please read the tutorial of parMix for more details. (https://github.com/biotoolscoders/parmix/blob/main/Infer1Recom.py)

# Contact
If you have any questions, please email: yiming.zhang.cse@uconn.edu or yufeng.wu@uconn.edu.

# Cite of the paper
To be announced.

# Reference
<a id="1">[1]</a> 
Pei, J., Zhang, Y., Nielsen, R., & Wu, Y. (2020). 
Inferring the ancestry of parents and grandparents from genetic data. 
PLoS computational biology, 16(8), e1008065.

