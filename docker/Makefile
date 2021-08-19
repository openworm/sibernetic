main: sibernetic Dockerfile
	cp Dockerfile sibernetic/
	docker build -t sibernetic sibernetic
	date > main

sibernetic: ../src intel-opencl-tmp ../makefile ../inc
	rm -rf sibernetic
	mkdir sibernetic
	cp -r ../src sibernetic/
	cp -r intel-opencl-tmp sibernetic/
	cp ../makefile sibernetic/
	cp -r ../inc sibernetic/

clean-docker:
	docker rmi sibernetic

clean: clean-docker
	rm -rf sibernetic
	docker ps -aq | xargs docker rm
	docker rmi sibernetic
	docker image prune -f
	rm main

.PHONY: clean clean-docker
