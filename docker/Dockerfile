FROM tensorflow/tensorflow:nightly-py3-jupyter

RUN apt-get update &&\
    apt-get install -y \
	apt-utils \
    	gcc \ 
    	fortune \
    	cowsay \
   	 vim

COPY requirements.txt /tmp/
COPY requirements_git.txt /tmp/

RUN pip install -r /tmp/requirements.txt
RUN pip install -r /tmp/requirements_git.txt

EXPOSE 8888
#CMD ["/bin/sh", "startup.sh"]
CMD ["jupyter","lab","--allow-root","--ip","0.0.0.0"]
