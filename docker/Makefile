main: sibernetic intel-opencl-tmp Dockerfile
	docker build -t sibernetic .
	date > main

sibernetic: ../src ../makefile ../inc 
	rm -rf sibernetic
	mkdir sibernetic
	cp -r ../src sibernetic/
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
