FROM python:3

COPY ./main.py /

WORKDIR /in

RUN pip3 install matplotlib
RUN pip3 install numpy

ENTRYPOINT ["python3", "/main.py"]