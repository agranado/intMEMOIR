# Infering the history of cell divisions from stochastic recording events.  

intMEMOIR is a synthetic biology experimental technology that enables recording mutations in DNA. 
Here we perform inference of lineage relationships from the intMEMOIR reaout. 


## Usage

```R

source("memoirReconstruction.R")

reconstruction = reconstructLineage(ground_phylo = ground_truth,mu = mu, alpha = alpha , return_tree = T  )

```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.
