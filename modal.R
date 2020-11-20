library(bsplus)
modal_bootstrapping <-
  bs_modal(
    id = "modal_bootstrapping",
    title = h3("Bootstrapping"),
    body = div(
      p("Bootstraping is a statistical resampling method to estimate the accuracy of sample estimates. 
        Such resampling methods repeatedly draw (re)samples from a given sample, and evaluate the same statistics for each of them. 
        Thereby, the distributions of these statistics can be estimated including their mean and variance."),
      p("Bootstrapping is a resample method with replacement. Thus, a resample may contain the same data point of the original sample multiple times."),
      h4("Bootstrapping in HouseKeepR"),
      p("HouseKeepR uses bootstrapping to identify the house-keeping genes with the highest and at the same time most stable expression in condition versus control samples.
          Here, each bootstrap resample is a set of your selected data sets that can possibly contain the same selected data set multiple times.
          For each such bootstrap resample (set of data sets) genes are ranked by high average expression and low expression variance."),
      p("Finally, a final ranking of genes is produced by sorting genes by smallest average rank across resamples. In case of ties, smallest rank variance across resamples decides."),
      h4("Bootstrap sample size"),
      p("The number of data sets contained in each bootstrap resample can be specified."),
      strong("We strongly recommend leaving this value at the default (the number of selected data sets)."),
      h4("Bootstrap replications"),
      p("The more bootstrap replications are performed, i.e. the more bootstrap resamples are drawn, the more accurately distributions of sample statistics can be estimated. 
          Choosing lower values, will greatly reduce running time, but lower statistical accuracy of results."),
      strong("We strongly recommend not setting this value lower than the default (100).")
    ),
    size = "large"
  )

modal_ensembl <-
  bs_modal(
    id = "modal_ensembl",
    title = "ENSEMBL",
    body = div(
      p("HouseKeepR makes extensive use of the annotation data of the Ensembl database. In particular, we use Ensembl to map between various gene identifiers and symbols. 
        We also use it to identify paralogues genes and match genes of different organisms."),
      strong("We strongly recommend to use the most recent available Ensembl release for each organism."),
      p(strong("However, to ensure long-term reproducible results one must specify a particular Ensembl release for the analysis."))
    ),
    size = "large"
  )

modal_seed <-
  bs_modal(
    id = "modal_seed",
    title = "Random Seed",
    body = div(
      p("In order to ensure reproducible results a random seed can be specified that is used to initialize the random generator for the bootstrap samples."),
      strong("In order to achieve the same results across multiple runs set the seed to the same numerical value")
    ),
    size = "large"
  )
