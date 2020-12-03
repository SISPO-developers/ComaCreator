import json
import time

def wait_processes(pids, sleep_time, logger):
    while len(pids) > 0:
        logger.info("wait 5 seconds before next polling")
        time.sleep(sleep_time)
        for i in range(0, len(pids)):
            poll = pids[i].is_alive()
            if poll is False:
                pids[i].join()
                del pids[i]
                break


def print_text(fname, txt, logger):
    if txt is not None:
        text = txt.format(fname)
        if logger is not None:
            logger.info(text)
        else:
            print(text)


def read_json(fname, txt=None, logger=None):
    print_text(fname, txt, logger)
    with open(fname, "r") as json_file:
        return json.load(json_file)
    return None


def save_json(fname, task, txt=None, logger=None):
    print_text(fname, txt, logger)
    with open(fname, 'w') as outfile:
        json.dump(task, outfile)


