from flask import Flask
from celery import Celery, current_app
from dotenv import load_dotenv, find_dotenv
import os

basedir = os.path.abspath(os.path.dirname(__file__))
load_dotenv(find_dotenv())


def make_celery(app):
    celery = Celery(
        app.import_name,
        backend=os.environ.get('result_backend'),
        broker=os.environ.get('broker_url')
    )

    celery.conf.update(app.config)

    class ContextTask(celery.Task):
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return self.run(*args, **kwargs)

    celery.Task = ContextTask
    return celery



