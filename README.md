## PARNAS ##
PARNAS identifies taxa that best represent diversity on a phylogenetic tree
and solves urgent needs in virology:
- Finding representative strains for detailed analysis (phenotypic characterization, Bayesian inference, etc.)
- Vaccine strain selection

PARNAS can take into account previous selections and user's constraints.
Additionally, PARNAS is flexible in allowing
arbitrary weighing of taxa, e.g., based on predicted fitness/antigenic drift.
More broadly PARNAS can be used to
- Downsample a large phylogeny while optimally preserving the underlying diversity
- Reduce redundancy among genetic/genomic sequences
- Identify key diversity groups on a phylogeny

PARNAS is faster and broader than [ADCL](https://matsen.github.io/pplacer/generated_rst/rppr_min_adcl_tree.html#rppr-min-adcl-tree) by Matsen et al. (Systematic Biology 2013), which solves a subset of PARNAS-enabled problems.

#### Installation ####
To install PARNAS, clone or download this project and run
`python setup.py install`. Note that PARNAS requires Python 3.7 or higher.

PARNAS depends on dendropy and Biopython for phylogenetic and MSA manipulations, numpy and numba for just-in-time compilation of the critical algorithms into machine code, and (optionally) phylo-treetime to infer ancestral substitutions along tree edges. These dependencies will be installed automatically.

## Tutorial ##

We use a human H1N1 (pdm09) dataset with HA sequences collected in 2020, downloaded from [IRD](fludb.org), for this tutorial.
The alignment and an inferred rooted tree can be found in the tutorial [folder](https://github.com/flu-crew/parnas/tutorial/H1N1_pdm_2020.zip).

The basic usage of PARNAS is to find a fixed number of representative taxa, as follows:<br>
`parnas -t H1N1_human_2020_IRD_CDS.rooted.tre -n 3 --color "H1N1_parnas_n3.tre"`<br>
PARNAS will identify 3 best representative strains and save a colored tree to H1N1_parnas_n3.tre.
Opening this tree in FigTree, will show the representatives and their respective clusters of strains with different colors. Below is `H1N1_parnas_n3.tre` output tree, when opened in FigTree. Each color corresponds to one PARNAS-selected representative.

<center>
<img src="tutorial/figures/H1N1_parnas_n3.png">
</center>

Additionally, in the output PARNAS specifies the amount of overall diversity covered by the chosen representatives, which is 43% for our dataset.

#### Determining the right number of representatives ####
The "diversity covered" metric calculated by PARNAS is a useful tool to determine how many representatives is sufficient.
To find the amount of diversity covered by different numbers of representatives, you can choose a large enough n, e.g., n=20, and run<br>
`parnas -t H1N1_human_2020_IRD_CDS.rooted.tre -n 20 --diversity "diversity_scores.csv"`<br>
This command will save a CSV file with optimal diversity scores for n between 2 and 20. For our dataset, it shows that only 6 representative strains are needed to cover over 60% of the overall diversity. Opening the CSV file in Excel/Numbers, one can then construct the following graph:

<center>
<img src="tutorial/figures/diversity_covered.png" width="450px">
</center>

#### Using prior representatives ####
It often may be useful to find new representatives and specify the previous strains/taxa that were treated as representatives. In our H1N1 dataset, we included two latest H1N1 vaccine strains, A/Wisconsin/588/2019 and A/Hawaii/70/2019. Using PARNAS we can fix these two strains as 'prior' representatives and find new representatives, so that old and new representatives work optimally together. We use `--prior-regex` option to notify PARNAS of the two vaccine strains:<br>
`parnas -t H1N1_human_2020_IRD_CDS.rooted.tre -n 3 --prior-regex "Vaccine.*" --color "H1N1_parnas_n3_vaccines.tre"`<br>
The result is shown in a tree below. PARNAS colors the prior representatives and the respective parts of the tree in red.

<center>
<img src="tutorial/figures/H1N1_parnas_n3_vaccines.png">
</center>

#### Specifying a coverage radius ####
You can notify PARNAS that a single taxon 'covers' other taxa within a fixed radius (using `--radius` option). PARNAS will then find representatives, which together cover as much diversity as possible. Alternatively, PARNAS can find the minimal number of representatives that jointly cover *all* diversity on the tree (`--cover` option). This feature is motivated by applications in virology, where evolutionary distance often correlates with antigenic drift. For example, in swine influenza A virus research, 5% (amino acid) divergence between HA sequences is considered a correlate of loss in cross-reactivity between strains. Therefore, specifying a radius on a tree corresponding to ~4% sequence divergence, can help identify good vaccine candidates.

To further facilitate this process, PARNAS can take amino acid sequence alignment and pass it to TreeTime to find ancestral AA substitutions along the tree edges. It will then re-scale the tree, so that edge lengths will reflect the number of substitutions along that edge, which will help specify a 4% divergence radius in terms of the # of amino acid substitutions.

`parnas -t H1N1_human_2020_IRD_CDS.rooted.tre --cover --threshold 96 --aa H1N1_human_2020_IRD_CDS.ha1.aln --color "parnas_96coverage.tre"`

In the above command we pass the AA alignment of HA1 sequences from our dataset with the `--aa` option, and specify the 96% threshold (i.e., 4% sequence divergence). Note that `--threshold` works as a substitute for `--radius`, when you would like PARNAS to use the alignment information and re-scale the tree. Running this command will show that a single strain (A/Wisconsin/32/2020) is sufficient to cover all diversity on our tree.

Next, we use a more restrictive threshold of 97% and also add the vaccine strains as the prior representatives. This way PARNAS can indicate the areas of the tree, which are not covered by the vaccine strains and suggest new representatives to solve this issue.

`parnas -t H1N1_human_2020_IRD_CDS.rooted.tre --cover --threshold 97 --aa H1N1_human_2020_IRD_CDS.ha1.aln --color "parnas_97coverage_vaccines.tre" --prior-regex "Vaccine.*"`

Opening `parnas_97coverage_vaccines.tre` in FigTree will show us that there are two clades in the tree, which are not covered by vaccine strains (the green and blue clades).

<center>
<img src="tutorial/figures/H1N1_parnas_97coverage_vaccines.png">
</center>

## PARNAS usage##

`parnas -t TREE [-n SAMPLES] [other options]`

*General options*

| Option | Description |
| --- | --- |
|--radius | Specify a radius (distance on a tree) so that every representative covers all diversity within that radius. PARNAS will then choose representatives that optimally cover as much diversity as possible |
| --prior-regex | Specify prior representatives with a regex. PARNAS will then identify representatives of 'new' diversity |
| --weights | Add a CSV file specifying weights for some or all taxa/strains |
| --cover | Instead of specifying the number of representatives, specify the radius and PARNAS will find representative that cover all diversity on the tree |

*Output options*

| Option | Description |
| --- | --- |

*Options to exclude/constrain taxa*

| Option | Description |
| --- | --- |

*Options to control sequence divergence*

| Option | Description |
| --- | --- |
