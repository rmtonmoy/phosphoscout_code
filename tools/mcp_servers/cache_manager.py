
from pathlib import Path
import json
import os

PROJECT_ROOT = Path(__file__).resolve().parents[2]
cache_path = str(PROJECT_ROOT / 'data' / 'cache' / 'uniprot')




def create_if_not_exists(lock, entity):
    with lock:
        path = f"{cache_path}/{entity}.json"
        if os.path.isfile(path) == False:
            #print("oopss")
            with open(f'{cache_path}/{entity}.json', 'w') as f:
                json.dump({}, f)

def write_whole_json(lock, entity, upd_data):
    raise AssertionError("Tried to write!")

    with lock:
        with open(f'{cache_path}/{entity}.json', 'w') as f:
            json.dump(upd_data, f)


def write_json(lock, entity, key, value):
    raise AssertionError("Tried to write!")

    with lock:
        if key == None or value == None:
            print(key)
            print(value)
            raise AssertionError("Tried to write None")

        with open(f'{cache_path}/{entity}.json', 'r') as f:
            upd_data = json.load(f)

        upd_data[key] = value
        with open(f'{cache_path}/{entity}.json', 'w') as f:
            json.dump(upd_data, f)



def read_whole_json(lock, entity):
    with lock:
        with open(f'{cache_path}/{entity}.json', 'r') as f:
            return json.load(f)


def read_json(lock, entity, key):
    with lock:
        with open(f'{cache_path}/{entity}.json', 'r') as f:
            data = json.load(f)
            if key in data.keys():
                return data[key]
            else:
                return None
