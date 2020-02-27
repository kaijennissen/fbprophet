FROM jupyter/scipy-notebook:2fd856878b83
USER root

RUN apt-get update &&\
    apt-get install -y \
	apt-utils \
    gcc \ 
    fortune \
    cowsay \
    vim

COPY requirements.txt /tmp/
RUN pip install -r /tmp/requirements.txt
COPY startup.sh /
EXPOSE 8888

#CMD ["/bin/sh", "startup.sh"]
CMD ["jupyter","lab","--allow-root","--ip","0.0.0.0"]
