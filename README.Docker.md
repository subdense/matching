
### Requirements

Install docker and docker-compose

(on ubuntu for docker compose: `sudo apt install docker-compose-v2`)

If you are working behind a http/https proxy, add the corresponding lines to your `compose.yaml` file:
```
...
   build:
      context: .
      args:
        HTTP_PROXY: ...
        HTTPS_PROXY: ...
...
```

### Building and running your application

Create the output path 'output_data' with write privileges.

When you're ready, start your application by running:
`sudo docker compose up --build`.

Your application will run and put the output files in 'output_data'.

To run a specific command other than the default (which is `python run.py`), you can use `docker run`, for example:
`sudo docker compose run --build --rm matching python run.py -layer1 data/bati/neudorf_2012.shp -layer2 data/bati/neudorf_2022.shp -output_prefix test_docker`


### Deploying your application to the cloud

First, build your image, e.g.: `docker build -t myapp .`.
If your cloud uses a different CPU architecture than your development
machine (e.g., you are on a Mac M1 and your cloud provider is amd64),
you'll want to build the image for that platform, e.g.:
`docker build --platform=linux/amd64 -t myapp .`.

Then, push it to your registry, e.g. `docker push myregistry.com/myapp`.

Consult Docker's [getting started](https://docs.docker.com/go/get-started-sharing/)
docs for more detail on building and pushing.

### References
* [Docker's Python guide](https://docs.docker.com/language/python/)
