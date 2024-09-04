from django.urls import re_path

from . import views

app_name = 'project-admin'

urlpatterns = [
    re_path(r'^$', views.home, name='home'),
    re_path(r'^config-general-settings/?$',
        views.config_general_settings, name='config-general-settings'),
    re_path(r'^config-oh-settings/?$',
        views.config_oh_settings, name='config-oh-settings'),
    re_path(r'^config-file-settings/?$',
        views.config_file_settings, name='config-file-settings'),
    re_path(r'^config-homepage-text/?$',
        views.config_homepage_text, name='config-homepage-text'),
    re_path(r'^login/?$', views.admin_login, name='login'),
    re_path(r'^add-file/?$', views.add_file, name='add-file'),
    re_path(r'^delete-file/(?P<file_id>\w+)/?$', views.delete_file,
        name='delete-file'),
]
