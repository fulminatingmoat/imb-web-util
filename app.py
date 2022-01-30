import hashlib
import os

import pandas as pd
from flask import Flask, render_template, request, redirect, url_for, flash, jsonify, Response, send_from_directory
from ext import make_celery
from lib import impute

CHUNK_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static', 'chunks')
COMPLETED_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static', 'completed')
RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static', 'results')
COMPRESSED_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static', 'compressed')

app = Flask(__name__)

app.config.update(
    broker_url='amqp://fulminatingmoat:initialpassword@172.16.11.7:49156/imb',
    result_backend='redis://172.16.11.7:49153/0'
)


celery = make_celery(app)


@celery.task(bind=True)
def impute_task(*args, **kwargs):
    return impute(*args, **kwargs)


def convert_to_pickle(file_id, file_name):
    ext = file_name.split('.')[-1]
    filename = f"{file_id}.{ext}"
    filepath = os.path.join(COMPLETED_DIR, filename)
    pickle_filepath = os.path.join(COMPRESSED_DIR, f"{file_id}.pkl")
    if os.path.exists(filepath) and not os.path.exists(pickle_filepath):
        match ext:
            case 'csv':
                pd.read_csv(filepath, sep=None).to_pickle(os.path.join(COMPRESSED_DIR, f"{file_id}.pkl"))
            case 'tsv':
                pd.read_csv(filepath, sep='\t').to_pickle(os.path.join(COMPRESSED_DIR, f"{file_id}.pkl"))
            case 'feather':
                pd.read_feather(filepath).to_pickle(os.path.join(COMPRESSED_DIR, f"{file_id}.pkl"))
            case 'parquet' | 'pq':
                pd.read_parquet(filepath).to_pickle(os.path.join(COMPRESSED_DIR, f"{file_id}.pkl"))


@app.route('/')
def hello_world():  # put application's code here
    #result = impute_task.delay()
    return render_template('index.html')


def get_chunk_name(file_id, chunk_num):
    return f"{file_id}.{chunk_num}"


def validate_chunk(flow_dict):
    chunk_num = int(flow_dict.get('flowChunkNumber'))
    chunk_size = int(flow_dict.get('flowChunkSize'))
    total_size = int(flow_dict.get('flowTotalSize'))
    file_name = flow_dict.get('flowFilename')
    file_id = flow_dict.get('flowIdentifier')
    current_chunk_size = int(flow_dict.get('flowCurrentChunkSize'))
    total_chunks = int(flow_dict.get('flowTotalChunks'))
    upload_token = flow_dict.get('flowUploadToken')
    if chunk_num > total_chunks:
        return True
    # f"{file_id}.{file_name.split('.')[-1]}"
    chunk = os.path.join(CHUNK_DIR, file_id, f"{file_id}.{chunk_num}")
    completed_file = os.path.join(COMPLETED_DIR, "{file_id}.{file_name.split('.')[-1]}")
    compressed_file = os.path.join(COMPRESSED_DIR, f"{file_id}.pkl")
    if os.path.exists(compressed_file):
        return True
    try:
        if f"{file_id}.{chunk_num}" not in os.listdir(os.path.join(CHUNK_DIR, file_id)):
            if completed_file not in os.listdir(COMPLETED_DIR):
                return False
            if os.path.getsize(completed_file) != total_size:
                return False
    except FileNotFoundError:
        return False
    if os.path.getsize(chunk) != current_chunk_size:
        return False

    return True


@app.route('/upload', methods=['GET', 'POST'])
def results():
    if request.method == 'POST':
        file = request.files['file']
        print(request.form['flowChunkNumber'])
        if not validate_chunk(request.form):
            chunk_name = get_chunk_name(request.form['flowIdentifier'], request.form['flowChunkNumber'])
            try:
                os.mkdir(os.path.join(CHUNK_DIR, request.form['flowIdentifier']))
            except FileExistsError:
                pass
            filename = os.path.join(CHUNK_DIR, request.form['flowIdentifier'], chunk_name)
            file.save(filename)
            if chunk_name in os.listdir(os.path.join(CHUNK_DIR, request.form['flowIdentifier'])):
                return Response(status=200)
            else:
                return Response(status=203)
        else:
            return Response(status=203)
    else:
        if validate_chunk(request.args):
            return Response(status=200)
        else:
            return Response(status=204)


@app.route('/finished', methods=['POST'])
def combine():
    print(request.get_json())
    json = request.get_json()
    filename = f"{json['hash']}.{json['file'].split('.')[-1]}"
    with open(os.path.join(COMPLETED_DIR, filename), "wb") as f:
        for chunk in range(1, len(os.listdir(os.path.join(CHUNK_DIR, json['hash'])))+1):
            with open(os.path.join(CHUNK_DIR, json['hash'], f"{json['hash']}.{chunk}"), "rb") as c:
                f.write(c.read())
    print(os.path.getsize(os.path.join(COMPLETED_DIR, filename)), json['size'])
    md5_hash = hashlib.md5()
    with open(os.path.join(COMPLETED_DIR, filename), "rb") as f:
        # Read and update hash in chunks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            md5_hash.update(byte_block)
        print(md5_hash.hexdigest(), json['hash'])
    if md5_hash.hexdigest() == json['hash']:
        convert_to_pickle(md5_hash.hexdigest(), filename)
        return Response(status=200)
    else:
        return Response(status=529)


@app.route('/header/<file_hash>', methods=['GET'])
def get_header(file_hash):
    file = pd.read_pickle(os.path.join(COMPRESSED_DIR, f"{file_hash}.pkl"))
    head = file.head(3).iloc[:, :6]
    #return jsonify(head.to_dict())
    data = head.to_dict(orient='split')
    data['shape'] = file.shape
    return jsonify(data)


@app.route('/analyze' , methods=['POST'])
def analyze():
    json = request.get_json()
    reference_file = json['reference_file_hash']
    map_file = json['reference_map_hash']
    summary_file = json['summary_stats_hash']
    genome_build = json['genome_build']
    custom = json['custom']
    print(type(custom))
    result = impute_task.delay(reference_file=reference_file, map_file=map_file, summary_file=summary_file,
                               genome_build=genome_build, custom=custom)
    return jsonify({'result': url_for('taskstatus', task_id=result.id)}), 202, {'Location': url_for('taskstatus', task_id=result.id)}


@app.route('/status/<task_id>')
def taskstatus(task_id):
    task = impute_task.AsyncResult(task_id)
    if task.state == 'PENDING':
        # job did not start yet
        response = {
            'state': task.state,
            'current': 0,
            'total': 1,
            'status': 'Pending...'
        }
    elif task.state != 'FAILURE':
        response = {
            'state': task.state,
            'progress': task.info.get('progress', 0),
            'total': task.info.get('total', 1),
            'status': task.info.get('status', '')
        }
        if 'result' in task.info:
            response['result'] = task.info['result']
    else:
        # something went wrong in the background job
        response = {
            'state': task.state,
            'current': 1,
            'total': 1,
            'status': str(task.info),  # this is the exception raised
        }
    return jsonify(response)


@app.route('/result/<result_hash>')
def serve_result(result_hash):
    return send_from_directory(os.path.join(RESULT_DIR), f"{result_hash}.tsv")


if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=True)
