import os
basedir = os.path.abspath(os.path.dirname(__file__))

username = os.getlogin()

class Config(object):
    DEBUG = False
    TESTING = False
    CSRF_ENABLED = True
    SECRET_KEY = 'UofACH!@3'

class MainConfig(Config):
    MAIN_ENV_PATH = '/home/' + username + '/anaconda3/envs/PRIMITI'
    SUB_ENV_PATH = '/home/' + username + '/anaconda3/envs/PRIMITI_iLearn'


