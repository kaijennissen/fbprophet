IMAGE = cassandra:1.0
CONTAINER = cassandra-container

THIS_FILE := $(realpath $(lastword $(MAKEFILE_LIST)))
THIS_FILE_DIR := $(shell dirname $(THIS_FILE))

build: 
	DOCKER_BUILDKIT=1 docker build --tag $(IMAGE)  \
	    -f Dockerfile \
        .

run:
	docker run --rm \
		-it \
		-p 8888:8888 \
	    -v $(THIS_FILE_DIR):/home/jovyan/work/ \
	    --name $(CONTAINER) \
	    $(IMAGE)
