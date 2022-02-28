### Pull latest files from github

`sudo git clone https://github.com/fulminatingmoat/imb-web-util.git /srv/flask-web-util`

### Enter directory

`cd /srv/flask-web-util`

### Create virtual environment

`python3 -m venv venv`

### Install python dependencies

`source ./venv/bin/activate`

`pip install -r requirements.txt`

### Systemd service file

/etc/systemd/system/flask-web-util.service

```
[Unit]
Description=Gunicorn instance to serve the web-utils run by flask
After=network.target

[Service]
User=root
Group=www-data
WorkingDirectory=/srv/flask-web-util/
Environment="PATH=/srv/flask-web-util/venv/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin"
ExecStart=/srv/flask-web-util/venv/bin/gunicorn --workers 3 --bind unix:flask-web-util.sock -m 007 'app:app'

[Install]
WantedBy=multi-user.target
```

`systemctl daemon reload`

`systemctl enable --now flask-web-util.service`

### NGINX site file

Site file needs domain name website is hosted on under instead of [[SERVER_NAME]]

/etc/nginx/sites-available/flask-web-util

`ln -s /etc/nginx/sites-available/flask-web-util /etc/nginx/sites-enabled`

```
server {
        server_name [[SERVER_NAME]];

        root /srv/flask-web-util/static;

        #add_header Content-Security-Policy "default-src  'self';";

        location / {
                gzip_static     on;
                sendfile        on;
                try_files $uri @proxy;
        }

        location @proxy {
                include proxy_params;
                proxy_pass http://unix:/srv/flask-web-util/flask-web-util.sock;
        }
}
```

### Upload flow

#### Client Side: templates/index.html

File is split into chunks, uploaded individually to static/chunks, and sent to the server

#### Server Side: app.py

combined into completed file in static/completed, completed file opened by pandas and saved as pickle in
static/compressed, reducing file size

### Result flow

Results are saved in static/task_results by celery task id, before calculating hash and moving to static/results under
the file's hash to prevent duplicates

### DISH implementation

The python version of DISH is implemented in lib/DISH.py, when analysis is started, the impute function will be run with
files referenced by hash given by the user

### File uploads

.csv files expect , delimiter

.tsv files expect \t delimiter

.txt files can be , or \t or space delimited

.feather and .parquet or .pq files are accepted