FROM debian:11

ENV PYTHON_VERSION 3.9
# libdpkg-perl is needed because, somehow, the CFLAGS in the python config script
# adds a flag that depends on a SPECS file from libdpkg-perl
RUN apt-get update && \
    apt-get install -y wget \
                       make \
                       ocl-icd-opencl-dev \
                       g++ \
                       libdpkg-perl \
                       "^python${PYTHON_VERSION}-dev$" \
                       python3-pip \
                       freeglut3-dev libglu1-mesa-dev

COPY intel-opencl intel-opencl

RUN cp -R intel-opencl/* / && \
  ldconfig && \
  rm -r intel-opencl

RUN cp -R /opt/intel/opencl/include/CL /usr/include/

ARG FFMPEG
ENV FFMPEG=${FFMPEG:-true}

RUN if [ "$FFMPEG" = true ] ; then \
        apt-get update && apt-get install -y libswscale-dev libavutil-dev libavcodec-dev libavformat-dev ; \
    fi

COPY sibernetic/configuration configuration

RUN pip3 install numpy

COPY sibernetic/inc inc
COPY sibernetic/src src
COPY sibernetic/makefile makefile
COPY sibernetic/main_sim.py main_sim.py

ENV PYTHON_CONFIG /usr/bin/python${PYTHON_VERSION}-config
ENV PYTHONPATH $PYTHONPATH:.
RUN make clean && make debug  # Use python 3 libs
ENTRYPOINT ["Release/Sibernetic"]
