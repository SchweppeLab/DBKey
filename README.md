# DBKey
These instructions will help users to build the DBKey docker image, execute it, and connect to the shiny webserver.

## Clone the repo:
First, users will need to clone a local copy of the repository.  

``` git clone https://github.com/SchweppeLab/DBKey.git ```

Next, change directory into the project.
```
cd DBKey
```

## Building the docker image:
First, ensure that docker is installed.  You may do so by checking if the ```docker --version``` command returns a sensical result.
If it is not installed, please visit the docker website for the relevant directions to install docker on your host operating system.

Once docker is installed, we may build a copy of the image like so:

``` docker-compose build ```
or
``` docker compose build ```

This process may take several minutes for the relevant base images to be pulled, and dependencies to be installed.

## Running a docker container:

Finally, one may launch a container using the following command:

```docker-compose up -d dbkey```
or
```docker compose up -d```

In this case, we have opened the local host http port 3838 (i.e., 127.0.0.1:3838) for the Shiny server to communicate through.


Simply open your favorite web browser and navigate to ```127.0.0.1:3838```

## Usage 
 - DBKey accepts Prosit-TMT .msp files and SpectraST .sptxt files. 
 - For Prosit .msp files, metadata needs to be manually entered (fragmentation energy, mass analyzer, CID/HCD)
 - Library entries lacking metadata will not be used in RTLS methods. Library entries are only considered RTLS candidates if FT/IT and CID/HCD matches acquired scan
 - If no ion annotations are provided, RTLS auto-generates them  prior to running.
 - Precursor m/zs are recalculated based on amino acid sequence and specifed modifications.
 - Top N and intensity cutoff will remove fragment peaks that don't pass the filters.
 - MassOffset take a .csv input with two columns, "Sequence" and "massOffset". Sequence is a list of library entries to apply the corresponding massOffset to.
 - To verify compatibility with RTLS, .db files can be examined in mzVault software. (https://thermo.flexnetoperations.com/)
 - Example input files are provided in [examples](https://github.com/SchweppeLab/DBKey/tree/main/examples)

## Reference
Real-time spectral library matching for sample multiplexed quantitative proteomics.
Christopher D McGann, William D Barshop, Jesse Canterbury, Chuwei Lin, Wassim Gabriel, Mathias Wilhelm, Graeme McAlister, Devin K Schweppe
bioRxiv 2023.02.08.527705; doi: https://doi.org/10.1101/2023.02.08.527705 
