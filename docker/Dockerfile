FROM debian:11

ENV PYTHON_VERSION 3.9
RUN apt-get update && \
    apt-get install -y wget \
                       make \
                       ocl-icd-opencl-dev \
                       vim \
                       g++ libpython3-dev \
                       "^python${PYTHON_VERSION}-dev$" \
                       freeglut3-dev libglu1-mesa-dev

COPY intel-opencl-tmp intel-opencl-tmp

RUN cd intel-opencl-tmp && \
  cp -R intel-opencl/* / && \
  ldconfig && \
  cd .. && \
  rm -r intel-opencl-tmp

RUN cp -R /opt/intel/opencl/include/CL /usr/include/

COPY sibernetic/inc inc
COPY sibernetic/src src
COPY sibernetic/makefile makefile


ENV PYTHON_CONFIG /usr/bin/python${PYTHON_VERSION}-config
RUN make clean && make all  # Use python 3 libs
