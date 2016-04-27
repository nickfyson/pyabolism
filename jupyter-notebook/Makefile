include env_make
NS = nickfyson
VERSION ?= latest

REPO = pyabolism
NAME = pyabolism
INSTANCE = default

.PHONY: build push shell run start stop rm release

build:
	docker build -t $(NS)/$(REPO):$(VERSION) .

# push:
# 	docker push $(NS)/$(REPO):$(VERSION)

shell:
	docker exec -it --user=jovyan $(NAME)-$(INSTANCE) bash

run:	rm
	docker run --name $(NAME)-$(INSTANCE) $(VOLUMES) $(PORTS) $(ENV) $(NS)/$(REPO):$(VERSION)
					# || docker start $(NAME)-$(INSTANCE) & docker attach $(NAME)-$(INSTANCE)

start:
	if docker ps --all | grep -q $(NAME)-$(INSTANCE) ; then \
		docker start $(NAME)-$(INSTANCE); \
	else \
		docker run -d --name $(NAME)-$(INSTANCE) $(VOLUMES) $(PORTS) $(ENV) $(NS)/$(REPO):$(VERSION); \
	fi;

stop:
	if docker ps | grep -q $(NAME)-$(INSTANCE) ; then docker stop $(NAME)-$(INSTANCE); fi

rm:
	if docker ps --all | grep -q $(NAME)-$(INSTANCE) ; then docker rm -f $(NAME)-$(INSTANCE); fi

dev:	rm build start

# release: build
# 	make push -e VERSION=$(VERSION)

default: build
