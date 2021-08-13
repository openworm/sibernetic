#!/bin/bash -e

if [ "$DISPLAY" ] ; then
    XSOCK=/tmp/.X11-unix
    XAUTH=/run/user/$(id -u)/.docker.xauth
    if [ ! -e $XAUTH ] ; then
        xauth -f "$XAUTH" generate "$DISPLAY" . untrusted
    fi
    xauth nlist $DISPLAY | sed -e 's/^..../ffff/' | xauth -f $XAUTH nmerge -
    XPART="-e DISPLAY=$DISPLAY -v $XSOCK:$XSOCK -v $XAUTH:$XAUTH -e XAUTHORITY=$XAUTH --net=host"
fi

docker run -u ${RUN_UID:-1000} \
    -v $(readlink -f simulations):/simulations \
    -v $(readlink -f video):/tmp/video \
    $XPART \
    -it \
    --device=/dev/dri:/dev/dri \
    sibernetic:latest \
    "$@"
