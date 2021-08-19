FROM debian:11

RUN apt-get update
RUN apt-get install -y wget \
                       make \
                       ocl-icd-opencl-dev \
                       vim \
                       g++ libpython3-dev \
                       pkg-config
RUN apt-get install -y python3-dev
RUN apt-get install -y freeglut3-dev libglu1-mesa-dev

COPY . sibernetic

RUN cd sibernetic/intel-opencl-tmp && \
  cp -R intel-opencl/* / && \
  ldconfig && \
  cd .. && \
  rm -r intel-opencl-tmp

RUN cp -R /opt/intel/opencl/include/CL /usr/include/

COPY makefile sibernetic/makefile


ENV PYTHON=python3.9
RUN cd sibernetic && \
    make clean && make all  # Use python 3 libs
