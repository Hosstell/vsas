FROM python:3

COPY ./main.py /
COPY ./requirements.txt /

WORKDIR /in

RUN pip3 install -r /requirements.txt

ENTRYPOINT ["python3", "/main.py"]