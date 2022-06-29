# MSPtoDB
These instructions will help users to build the MSPtoDB docker image, execute it, and connect to the shiny webserver.

## Clone the repo:
First, users will need to clone a local copy of the repository.  

``` git clone https://github.com/SchweppeLab/MSPtoDB.git ```

Next, change directory into the project and switch to the 'docker' branch of the SchweppeLab/MSPtoDB project.
```
cd MSPtoDB
git checkout chunkedDockerToMain
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

```docker-compose up -d msptodb```
or
```docker compose up -d```

In this case, we have opened the local host http port 3838 (i.e., 127.0.0.1:3838) for the Shiny server to communicate through.

Simply open your favorite web browser and navigate to ```127.0.0.1:3838```
