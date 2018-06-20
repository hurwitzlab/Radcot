SING242 = /usr/local/bin/singularity #2.4.2-dist
SING231 = /home/vagrant/bin/singularity #2.3.1-dist

SINGULARITY = $(SING231)
IMG_NAME = radcot.img
DEF_NAME = image.def
SIZE = 4444

help:
	sudo $(SINGULARITY) exec $(IMG_NAME) run.sh -h

run:
	sudo $(SINGULARITY) run $(IMG_NAME)

clean:
	rm -f $(IMG_NAME)

img: clean 
#	sudo $(SINGULARITY) build $(IMG_NAME) $(DEF_NAME) #actually the 2.4.2 way
	sudo rm -rf /tmp/* #since /tmp is bound to /tmp on host
	sudo $(SINGULARITY) create --size $(SIZE) $(IMG_NAME)
	sudo $(SINGULARITY) bootstrap $(IMG_NAME) $(DEF_NAME)

shell:
	sudo $(SINGULARITY) shell --writable $(IMG_NAME)
#	sudo $(SINGULARITY) shell $(IMG_NAME)

scp:
	scp $(IMG_NAME) stampede2:/work/03859/sdaniel/stampede2/radcot/stampede/
