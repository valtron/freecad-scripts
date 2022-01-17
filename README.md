# FreeCAD scripts

## Usage with a docker container

- Download and run the container with `docker run --pull always --rm -p 8888:8888 symbench/freecad-scripts`.
- Connect to the jupyter server within the container using your browser at `https://127.0.0.1:8888/?token...`
- Look at the `models/pv_example.ipynb` example on how to use the `freecad_scripts` library.

## Usage on your local machine

- Make sure that you have the `freecad-python3` (version 0.19) and `gmsh` packages installed. 
  On ubuntu 21.04 you would run the command `sudo apt-get install freecad-python3 gmsh`.
- Clone the `freecad_scripts` repository and install it with `pip3 install -e .`. 
  This will allow you to update the code (e.g. with `git pull`) without reinstallation.
- Use `freecad-scripts --help` to list the supported subcommands.

## Creating the docker image (only for maintainers)

- Clone or update the `freecad_scripts` repository from https://github.com/symbench/freecad-scripts
- Run `docker build -t symbench/freecad-scripts .` to create the latest version
- Run `docker push symbench/freecad-scripts:latest` to publish it on docker hub.
