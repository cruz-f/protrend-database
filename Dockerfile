FROM python:3.7.7-slim

MAINTAINER BYOSYSTEMS group <biosystems.um@gmail.com>

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

WORKDIR /home/protrend-database

COPY requirements.txt /home/protrend-database
COPY run.sh /home/protrend-database

RUN pip install -r requirements.txt

COPY src/protrend /home/protrend-database/protrend


CMD ["/bin/bash", "run.sh"]