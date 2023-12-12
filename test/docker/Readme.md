# What is this?

Dockerfile used to generate `ghcr.io/hsu-hpc/testmamico`, a container used for the github workflow that executes the CI job `test`.

# How to update the container image?

- `cd` here
- Modify the Dockerfile if necessary, e.g. if new packages are needed 
- Execute `docker build -t testmamico .`
- Get a [personal access token](https://docs.github.com/de/packages/working-with-a-github-packages-registry/working-with-the-container-registry)
- `export CR_PAT=YOUR_TOKEN`
- `echo $CR_PAT | docker login ghcr.io -u [username] --password-stdin`
- `docker tag testmamico ghcr.io/hsu-hpc/testmamico`
- `docker push ghcr.io/hsu-hpc/testmamico:latest`

# How to test the container locally

`docker run --rm -it --entrypoint bash testmamico`