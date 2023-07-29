FROM python:3

COPY ./main.py /

WORKDIR /in

RUN pip3 install -r requirements.txt

ENTRYPOINT ["python3", "/main.py"]