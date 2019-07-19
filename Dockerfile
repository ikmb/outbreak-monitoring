FROM continuumio/anaconda
LABEL authors="Marc P. Hoeppner" description="Docker image containing software dependencies for the IKMB outbreak pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /opt/conda/envs/outbreak-monitoring-1.0/bin:$PATH

RUN apt-get update && apt-get -y install procps make gcc  git build-essential autotools-dev automake libsparsehash-dev libboost-all-dev

RUN cd /opt && \
	git clone --recurse-submodules https://github.com/bcgsc/biobloom.git build_bloom && \
	cd build_bloom && \
	./autogen.sh && \
	./configure --prefix=/opt/biobloom && make && \
	cd /opt && rm -Rf build_bloom

ENV PATH /opt/biobloom/bin


