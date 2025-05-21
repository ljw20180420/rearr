from celery import Celery
import os

try:
    celery_broker = os.environ["CELERY_BROKER"]
except KeyError:
    celery_broker = "amqp://localhost:5672"
try:
    celery_backend = os.environ["CELERY_BACKEND"]
except KeyError:
    celery_backend = "redis://localhost:6379/0"

celeryApp = Celery(
    broker=celery_broker,
    backend=celery_backend,
    include=['celery_project.tasks'],
    broker_connection_retry_on_startup=True,
    timezone='Asia/Shanghai'
)

# set_default make default_app = celeryApp. This is necessary because in flask function, _tls.current_app = None, which causes celery.current_app to be default_app. For more details, see https://stackoverflow.com/questions/26527214/why-celery-current-app-refers-the-default-instance-inside-flask-view-functions.
celeryApp.set_default()