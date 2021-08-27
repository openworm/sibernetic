main: sibernetic intel-opencl Dockerfile
	docker build -t sibernetic .
	date > main

intel-opencl:
	mkdir -p intel-opencl-tmp/intel-opencl
	cd intel-opencl-tmp && \
		wget http://registrationcenter-download.intel.com/akdlm/irc_nas/11396/SRB5.0_linux64.zip && \
		unzip -of SRB5.0_linux64.zip && \
		tar -C intel-opencl -Jxf intel-opencl-r5.0-63503.x86_64.tar.xz && \
		tar -C intel-opencl -Jxf intel-opencl-devel-r5.0-63503.x86_64.tar.xz && \
		tar -C intel-opencl -Jxf intel-opencl-cpu-r5.0-63503.x86_64.tar.xz && \
		mv intel-opencl ../
	rm -rf intel-opencl-tmp

sibernetic: ../src ../makefile ../inc 
	rm -rf sibernetic
	mkdir sibernetic
	cp -r ../src ../makefile ../inc ../main_sim.py ../configuration sibernetic/

clean-docker:
	docker rmi sibernetic

clean: clean-docker
	rm -rf sibernetic
	docker ps -aq | xargs docker rm
	docker rmi sibernetic
	docker image prune -f
	rm main

.PHONY: clean clean-docker
