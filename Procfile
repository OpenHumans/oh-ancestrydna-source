release: python manage.py migrate
web: gunicorn oh_data_uploader.wsgi --log-file=-
worker: celery -A main worker --concurrency=1
