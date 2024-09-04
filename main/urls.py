from django.urls import re_path

from . import views

urlpatterns = [
    re_path(r'^$', views.index, name='index'),
    re_path(r'^complete/?$', views.complete, name='complete'),
    re_path(r'^logout/?$', views.logout_user, name='logout'),
    re_path(r'^overview/?$', views.overview, name='overview'),
    re_path(r'^upload_simple/?$', views.upload_old, name='upload_old'),
    re_path(r'^about/?$', views.about, name='about'),
    re_path(r'^list/?$', views.list_files, name='list'),
    re_path(r'^delete/(?P<file_id>\w+)/?$', views.delete_file, name='delete'),
    re_path(r'^trigger_processing/?$', views.trigger, name='trigger'),
]
