# Docker Image for Running Ubiquity in Azure Machine Learning

## Build and Push to DockerHub
```
docker build -t ubiquity_aml .

docker run --name ubiquity_aml --rm -p 8787:8787 nlmixr_aml

docker image tag ubiquity_aml cford38/ubiquity_aml:latest
docker push cford38/ubiquity_aml:latest
```