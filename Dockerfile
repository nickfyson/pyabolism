
FROM jupyter/scipy-notebook

# File Author / Maintainer
MAINTAINER Nick Fyson

# Update the repository sources list
RUN apt-get -y update

RUN apt-get -y install wget
RUN apt-get -y install libxml2
RUN apt-get -y install libxml2-dev
RUN apt-get -y install zlib1g zlib1g-dev 
RUN apt-get -y install bzip2 libbz2-dev 


# install libsbml into python2 CONDA environment
# early to avoid caching issues!
RUN ($CONDA_DIR/envs/python2/bin/pip install python-libsbml) 

RUN ($CONDA_DIR/envs/python2/bin/pip install xlrd) 

# install gurobi library
RUN wget -q -O - http://packages.gurobi.com/gurobi-apt.sh | sudo sh

# install Pyabolism into python2 CONDA version
COPY . /tmp/pyabolism
RUN (cd /tmp/pyabolism && $CONDA_DIR/envs/python2/bin/python setup.py install && rm -rf /tmp/pyabolism)

# install python bindings into python2 CONDA version
RUN (cd /opt/gurobi*/linux64 && $CONDA_DIR/envs/python2/bin/python setup.py install)

# copy examples directory
COPY examples /home/jovyan/work/examples/
RUN chown -R $NB_USER:users /home/$NB_USER/work/examples

