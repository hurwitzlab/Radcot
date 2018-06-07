APP = radcot-0.0.1
EMAIL = scottdaniel@email.arizona.edu

clean:
	find . \( -name \*.out -o -name \*.param \) -exec rm {} \;

up:
	copy_from_manifest.py

app:
	apps-addupdate -F app.json

test: clean
	sbatch test.sh

job:
	jobs-submit -F job.json

jobs-template:
	jobs-template -A $(APP) > job.json

jobs-submit:
	jobs-submit -F job.json

list:
	files-list -S data.iplantcollaborative.org scottdaniel/applications/$(APP)/stampede
	
shell:
	sudo $(SINGULARITY) shell --writable $(IMG_NAME)
#	sudo $(SINGULARITY) shell $(IMG_NAME)

scp:
	scp $(IMG_NAME) stampede2:/work/03859/sdaniel/stampede2/radcot/stampede/
