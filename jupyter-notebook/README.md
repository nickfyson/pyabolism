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

For purposes of licensing Gurobi when spinning up this Docker container you should include the additional option ```--net=host```, dictating that docker shouldn't containerize networking. See the [Docker networking documentation](https://docs.docker.com/articles/networking/) for details, but in short this means that the Gurobi license will be bound to the docker host, rather than a specific instance of the container. This is vital for reusing the Pyabolism image without obtaining a fresh license every time. Using this networking option means we define the port Jupyter connects through in an alternative way (see below).

For convenience we name the container 'pyabolism', and the full command you need to run is...

```bash
docker run -d --net=host -e PORT=8888 --name=pyabolism nickfyson/pyabolism
```

Once the command completes, the Jupyter notebook should be viewable on port 8888 of your docker machine IP address. For example...

```bash
    http://123.123.123.123:8888
```

To license your gurobi installation, open a terminal from within the Jupyter interface by selecting 'Terminal' from the 'New' menu. Type the following command using the license code optained from Gurobi, and sticking with default options for save location...

```bash
grbgetkey xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
```

Finally, returning to the command prompt on your host computer, take a snapshot of the running container. This preserves the licensing of Gurobi for later use...

```bash
docker commit pyabolism pyabolism_lic
```

Now in future a fully licensed Pyabolism container can be created using the following command...

```bash
docker run -d --net=host -e PORT=8888 --name=pyabolism pyabolism_lic
```

If sharing your docker snapshot with a 3rd party, or transferring to another host system, you will need to rerun the licensing command with a new license code.
