## Notes for using Docker with *pipeRline*

Need to: 
- `COPY` all files needed for running the pipeline into the Docker image, including:
    - _targets.R
    - _targets_packages.R
    - functions.R
    - themes.R

The R scripts can be sourced within the image, although it might be better to source them at runtime by the user, by running an R script already present within the image. 

The user, in the past, would clone the Github repo for [pipeRline](https://github.com/jackscanlan/piperline) and then run the scripts from within the cloned directory. 

You could run the Docker container within the same cloned directory and then source the files as-is. Nothing would need to change. User just adds their sample data to `./data` and adds a samplesheet too.

In the context of BASC, you clone the directory into your run folder, shifter the Docker image and then run locally, using a volume to persist the data (into `results` or something).