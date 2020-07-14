---
title: "Connecting R and Python for CWTS Scientometrics Summer School"
author: "Alejandro Espinosa"
date: "`r Sys.Date()`"
output: html_document
---

## Summer School

The following script is my notes from the CWTS Scientometrics Summer School. These codes are available in the official GitHub of [CWTS](https://github.com/CWTSLeiden/CSSS). Also, in my particular case, I am using a Mac OS/Linux operation system. In Windows, please launch the "Anaconda prompt". 

In the following, I will only present Mac OS/Linux operation system.

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
#knitr::knit_engines$set(python = reticulate::eng_python)

```

## Conda environment

To begin, I will first create an environment that would contain all the python modules of the entire summer school using [conda](https://docs.conda.io/). After installing conda, I open my terminal and add the following line of codes: `conda create --name cwts python=3.8 anaconda` that would create a new environment call `cwts` and would also install python version 3.8. If you are willing to activate the environment from the terminal, then you would have to use `conda activate cwts`. Notice that I give the name `cwts` that could be modified for whatever you prefer.

We can install some modules from python before starting the summer school, such as 

- `python -m pip install numpy`
- `python -m pip install panda`
- `python -m pip install igraph`
- `python -m pip install matplotlib`
- `python -m pip install leidenalg`
- `python -m pip install networkx`

In my personal computer, I had some issues with the visualizations of `igraph`. Hence, I install some extra modules:

- `python -m pip install cairocffi`
- `python -m pip install pycairo`

Also, I install `cairo` using the following code suggested in the webpage of [pycairo](https://pycairo.readthedocs.io/en/latest/getting_started.html), that use [Homebrew](https://brew.sh):

- `brew install cairo pkg-config`

After installing some packages, you should now deactivate the session, and you would have to put the following line of code: `conda deactivate` in your terminal. If you are willing to check your environments, just press: `conda info --envs` in your terminal. If for some reasons you prefer to remove your entire environment, you would have to use the following code: `conda remove --name cwts --all`

Is often common for `python` users to use the `Jupyter Notebook`, that might be a better option if you are just going to use `python`. In which case, after activating your conda environment, you could set your working directory through `cd [DIRECTORY]` in the terminal. Then launch the Jupyter notebook by `jupyter notebook`. From there, you can start the [available script](https://github.com/CWTSLeiden/CSSS) of the summer school `01-basics.ipynb` and `02-advanced.ipynb`. We encourage you to check the Github! If you prefer a * more straightforward solution* to run the summer school, you could run online the python notebook available in the GitHub. 


## Connecting R and Python

We can run `Python` codes directly from `R`. The reason to use R environment rather than python is that most of the modules used in the summer school are programmed in python. Nonetheless, due that in my personal workflow, I tend to use statistical models for social networks analysis, I prefer connecting the **best of both worlds**. Some of these other statistical models that use social networks are the [exponential random graph model](http://statnet.org) or the [stochastic actor-oriented model](https://www.stats.ox.ac.uk/~snijders/siena/) coded in R.

Now that I have an environment called `cwts`, I could now start using R and Python together. We would need to install the `reticulate` package.
```{r, message=FALSE}
rm(list=ls())
library(here)
library(reticulate)

use_condaenv(condaenv = "cwts", conda = "auto", required = TRUE)
main <- import_main()

```

We can also install a package from R using anaconda of Python, and run any Python code in R using `py_run_string`
```{r}
#reticulate::conda_install(c("leidenalg", "igraph"), 
#                            envname = "cwts", pip = TRUE)

py_run_string("import numpy as np")
py_run_string("my_python_array = np.array([2,4,6,8])")
py_run_string("print(my_python_array)")

```

Also, we can use python directly from R adding the chunk `{python}` instead of `{r}` in the Rmarkdown. For example, we can replicate the classic model of Holland, Laskey & Leinhardt (1983) on Stochastic Block Model using [stochastic_block_model](https://networkx.github.io/documentation/stable/reference/generated/networkx.generators.community.stochastic_block_model.html)
```{python}
import networkx as nx
sizes = [75, 75, 300]
probs = [[0.25, 0.05, 0.02],
        [0.05, 0.35, 0.07],
        [0.02, 0.07, 0.40]]
        
g = nx.stochastic_block_model(sizes, probs, seed=0)
len(g)

H = nx.quotient_graph(g, g.graph['partition'], relabel=True)

for v in H.nodes(data=True):
      print(round(v[1]['density'], 3))
      
for v in H.edges(data=True):
      print(round(1.0 * v[2]['weight'] / (sizes[v[0]] * sizes[v[1]]), 3))

```

There are neat packages out there in `Python`, some of my favourites are:

1. [`python`](https://igraph.org/python/):  for general network-oriented implementations

2. [`NetworkX`](https://networkx.github.io):  for general network-oriented implementations

3. [`leiden_lag`](https://leidenalg.readthedocs.io/en/stable/): for community detection, with a special focus in the `leiden algorithm`

4. [`graph-tool`](https://graph-tool.skewed.de): for stochastic block models and other general network-oriented implementations
