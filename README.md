[![build](https://github.com/jvparidon/lmerMultiMember/actions/workflows/r.yml/badge.svg)](https://github.com/jvparidon/lmerMultiMember/actions)
[![codecov](https://codecov.io/gh/jvparidon/lmerMultiMember/branch/main/graph/badge.svg)](https://codecov.io/gh/jvparidon/lmerMultiMember)
## `lmerMultiMember`
Wrapper around `lme4::lmer` and `lme4::glmer` to add support for multimembership random effects.  
This package was developed by Jeroen van Paridon and Phillip Alday, based on an initial worked example of a multiple membership linear mixed model authored by Ben Bolker.

### What is a multimembership random effect?
Let's take authorship as an example: If we want to model some aspect of published journal articles, e.g. `citations ~ word_count`, we might want to account for variability caused by the authors of papers by using random intercepts. However, many papers don't have just a single author. How do we account for author variability then?  

One method could be to model each unique grouping of authors as a separate level of the random effect (i.e. "John A", "Mary B", "John A & Mary B", and "Mary B & John A" are four different levels of the random effect) but this has the downside of not correctly attributing the variance to the individual authors in the groups.  

Another option is to have separate random effects for each position in the author list (i.e. `(1|author1) + (1|author2) + (1|author3)` etc.) but while this does treat each author as a separate entity, it treats the same author as different people depending on where they occur in the author list, which is also a little strange.  

The most natural option would be to associate each observation (journal article) with multiple levels (authors) of a single random effect. This is a multimembership random effect.  

### How to specify a multimembership random effect
Specifying a multimembership model in `lmerMultiMember` works just like specifying any other mixed effects model in `lme4`, with the addition of a membership matrix (or weight matrix). This sparse matrix contains rows for all the observations in your dataset and columns for each unique group member. For each observation, the matrix contains 1s for all the group members associated with it, with 0s everywhere else. (If group members should not have equal association strength, the 1s in the membership matrix can, in principle, be replaces with fractions to represent relative contributions, etc.) The model syntax could then be e.g. `lmerMultiMember::lmer(citations ~ word_count + (1|author), memberships=list(author=membership_matrix))`.  

Creating a membership matrix might seem a little daunting, but we provide a helper function `lmerMultiMember::weights_from_vector` that creates a membership matrix from a vector of group memberships in comma-separated strings (e.g. `c("A,B", "A,B,C", "B", "C,A")`) so if you have group memberships documented in your dataset it should be easy enough to create a membership matrix for you model.  

### Contact
If you have any issues fitting models with `lmerMultiMember`, feel free to contact Jeroen van Paridon at [vanparidon@wisc.edu](mailto:vanparidon@wisc.edu).
