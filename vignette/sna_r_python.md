# "Connecting R and Python for CWTS Scientometrics Summer School"

 [Alejandro Espinosa-Rada](https://www.research.manchester.ac.uk/portal/en/researchers/alejandro-espinosa(4ed72800-e02b-47a8-a958-640b6a07f563).html)

## Summer School

<div class="alert alert-success">
The following script is my notes from the CWTS Scientometrics Summer School. These codes are available in the official GitHub of [CWTS](https://github.com/CWTSLeiden/CSSS). Also, in my particular case, I am using a Mac OS/Linux operation system. In Windows, please launch the "Anaconda prompt". 
</div>

In the following, I will only present Mac OS/Linux operation system.

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
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

After installing some packages, you should now deactivate the session, and you would have to put the following line of code in your terminal: `conda deactivate`. If you are willing to check your environments, just press in your terminal: `conda info --envs`. If for some reasons you prefer to remove your entire environment, you would have to use the following code: `conda remove --name cwts --all`

Is often common for `Python` users to use the [`Jupyter Notebook`](https://jupyter.org), that might be a better option if you are just going to use `python`. In which case, after activating your conda environment, you could set your working directory through `cd [DIRECTORY]` in the terminal. Then after installying, launch the Jupyter notebook by `jupyter notebook`. From there, you can start the [available script](https://github.com/CWTSLeiden/CSSS) of the summer school `01-basics.ipynb` and `02-advanced.ipynb`. We encourage you to check the Github! If you prefer a *more straightforward solution* to run the summer school, you could run online the python notebook available in the GitHub. 


## Connecting R and Python

<div class="alert alert-success">
We can run `Python` codes directly from `R`. The reason to use `R` environment rather than `Python` is that most of the modules used in the summer school are programmed in `Python`. Nonetheless, due that in my personal workflow, I tend to use statistical models for social networks analysis, I prefer connecting the <span style="color:blue">**best of both worlds**</span>. Some of these other statistical models that use social networks such as:  

1. Different models available in [Statnet](http://statnet.org) (e.g. exponential random graph models, epidemiological models, relational event models, or the latent position and cluster models for statistical networks)

2. [Stochastic actor-oriented model](https://www.stats.ox.ac.uk/~snijders/siena/)

3. [Dynamic Network Actor-Oriented Model](https://github.com/snlab-ch/goldfish)
</div>

Now that I have an environment called `cwts`, I could now start using R and Python together. We would need to install the `reticulate` package.
```{r, message=FALSE}
# install.package("here")
# install.package("reticulate")
library(here)
library(reticulate)

use_condaenv(condaenv = "cwts", conda = "auto", required = TRUE)
main <- import_main()

```

We can also install a module from `Python` using `R`, and run any `Python` code in `R` using `py_run_string`
```{r}
# install a module from Python using R codes:
# conda_install(c("leidenalg", "igraph"), envname = "cwts", pip = TRUE)

# run a Python code from R:
py_run_string("import numpy as np")
py_run_string("my_python_array = np.array([2,4,6,8])")
py_run_string("print(my_python_array)")

```

Also, we can use `Python` directly from `R` adding the chunk `{python}` instead of `{r}` in the [`Rmarkdown`](https://rmarkdown.rstudio.com). For example, we can replicate the classic model of [Holland, Laskey & Leinhardt (1983)](https://www.sciencedirect.com/science/article/pii/0378873383900217) on Stochastic Block Model using [`stochastic_block_model`](https://networkx.github.io/documentation/stable/reference/generated/networkx.generators.community.stochastic_block_model.html)
```{python results="hide"}
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

1. [`python`](https://igraph.org/python/):  general purpose network analysis

2. [`NetworkX`](https://networkx.github.io):  general purpose network analysis

3. [`snap`](http://snap.stanford.edu): general purpose network analysis

4. [`metaknowledge`](https://metaknowledge.readthedocs.io/en/latest/): computational research in bibliometrics, scientometrics, and network analysis

5. [`leiden_lag`](https://leidenalg.readthedocs.io/en/stable/): for community detection, with a special focus in the `leiden algorithm`

6. [`graph-tool`](https://graph-tool.skewed.de): for stochastic block models and other general purpose network analysis

7. [`EstimNetDirected`](https://github.com/stivalaa/EstimNetDirected): Exponential random graph models for big networks