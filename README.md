MoffittFunctions
=============


MoffittFunctions is a collection of useful functions designed to assist in analysis and creation of professional reports. The current MoffittFunctions functions can be broken down to the following sections:

- Testing Functions
  - two_samp_bin_test
  - two_samp_cont_test
  - cor_test
- Fancy Output Functions
  - pretty_pvalues
  - pretty_pvalues_compareGroups
  - stat_paste
  - paste_tbl_grp
  - pretty_model_output
  - run_pretty_model_output
  - pretty_km_output
  - run_pretty_km_output
- Utility Functions
  - round_away_0
  - get_session_info
  - get_full_name
- Example Dataset
  - Bladder_Cancer


### Getting Started
You can install packages directly from GitLab  with a couple R lines of code.

```r

#Installing from Gitlab if you have ssh key set up
cred = git2r::cred_ssh_key(
	publickey = "MYPATH/.ssh/id_rsa.pub", 
	privatekey = "MYPATH/.ssh/id_rsa")

devtools::install_git(
    "git@gitlab.moffitt.usf.edu:ReproducibleResearch/MoffittFunctions.git", 
    credentials = cred, 
    build_opts = NULL)

# Installing from GitHub
remotes::install_git("https://github.com/z2thet/MoffittFunctions", build_opts = NULL)


# Loading MoffittFunctions and Example Dataset
library(MoffittFunctions)
data("Bladder_Cancer")
```

#### MoffittTemplates Package

The [MoffittTemplates](https://gitlab.moffitt.usf.edu:8000/ReproducibleResearch/R_Markdown_Templates) package makes extensive use of the `MoffittFunctions` package, and is a great way get started making professional statistical reports.
