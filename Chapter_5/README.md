# AMR_web
Developing a website to run a user friendly tool of our AMR phenoptype prediciton models

To pull the singularity image used for the decision tree smk:
```
singularity pull --arch amd64 library://ldillon/amrwebsite/amrmlpipeline:6
```

To pull the singularity image used for the CNN smk:
```
singularity pull --arch amd64 library://lucyd/machinelearning/mlpackages:1
```
