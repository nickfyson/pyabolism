# Pyabolism Jupyter Notebook (Docker Container)

## What it Gives You

* An unlicensed installation of the Gurobi optimization library, with bindings installed for the Conda Python 2.7.x environments. (See below for licensing details.)
* The Pyabolism python module installed into the Conda Python 2.7.x enviroment.
* A folder with example scripts, located in the default Jupyter working directory.
* Docker container is based on jupyter/scipy-notebook, so also includes...
    * Jupyter Notebook server (v4.0.x or v3.2.x, see tag)
    * Conda Python 3.4.x and Python 2.7.x environments
    * pandas, matplotlib, scipy, seaborn, scikit-learn, scikit-image, sympy, cython, patsy, statsmodel, cloudpickle, dill, numba, bokeh pre-installed
    * Unprivileged user `jovyan` (uid=1000, configurable, see options) in group `users` (gid=100) with ownership over `/home/jovyan` and `/opt/conda`
    * [tini](https://github.com/krallin/tini) as the container entrypoint and [start-notebook.sh](../minimal-notebook/start-notebook.sh) as the default command
    * Options for HTTPS, password auth, and passwordless `sudo`

### Python 2.7.x only

Note that while the Jupyter installation includes a Python 3 kernel, the Pyabolism module currently only works with Python 2.7.x.

## Usage Instructions

This image comes with an installation of the Gurobi optimization library, but before running Pyabolism code you must obtain a valid license. Gurobi is free for academic use, and the necessary activation code can be easily obtained through their website, [www.gurobi.com](http://www.gurobi.com).

For purposes of licensing Gurobi when spinning up this Docker container you should include the additional option --net=host, dictating that docker shouldn't containerize networking. See the [Docker networking documentation](https://docs.docker.com/articles/networking/) for details, but in short this means that the Gurobi license will be bound to the docker host, rather than a specific instance of the container. This is vital for reusing the Pyabolism image without obtaining a fresh license every time. Using this networking option means we define the port Jupyter connects through in an alternative way. To simplify the running of this container we include a Makefile, for which you must open a terminal inside the folder 'jupyter-notebook'.

The first stage is to rename env_make_example as env_make and customise the port number as required. Then the container should be built and run with the default options...
```bash
    make build run
```

The Jupyter notebook should be viewable on the relevant port of your docker machine IP address. For example, if you stuck with the default port of 8888...

```bash
    http://123.123.123.123:8888
```

To license this instance of the gurobi installation, open a terminal from within the Jupyter interface by selecting 'Terminal' from the 'New' menu. Type the following command using the license code optained from Gurobi, and sticking with default options for save location...

```bash
    grbgetkey xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
```

Returning to the home of the Jupyter notebook you should see a new file 'gurobi.lic'. Copy and paste the contents of this file into a new file named 'gurobi.lic' placed in the jupyter-notebook folder. Building the Docker container will automatically include all files with the .lic suffix in the root of the the container, and this can be used to ensure that all new containers are licensed for Gurobi (so long as the host computer stays the same). Now rebuild and run the container.
```bash
make build run
```

You should have a fully licensed Gurobi installation in the resulting container, which will work so long as the container is instantiated with the correct settings and on the same host computer. If sharing your docker snapshot with a 3rd party, or transferring to another host system, you will need to rerun the licensing command with a new license code.

In general use, the Pyabolism Docker container is best run in the background, using...
```bash
make start
```

To use the Pyabolism container to do your own work, you need to map a folder from your local machine to the 'work' folder in the container. This can be achieved following the example in the example env_make file. Simply replace the example path with the full path to a folder on your computer, and uncomment the line. Then re-run the container with...
```bash
make rm start
```
