import argparse
import asyncio
import sys
import time
import xml.etree.ElementTree as ET
from threading import Lock

import aiohttp


def progress(count, total, suffix=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s  %s\r' % (bar, percents, '%', suffix))
    sys.stdout.flush()  # As suggested by Rom Ruben


argParser = argparse.ArgumentParser()
argParser.add_argument("-s", "--search", type=str, nargs='+', required=True, help="Search query")
argParser.add_argument("-o", "--output", type=str, required=True, help="Output file name")
args = argParser.parse_args()


SEARCH = ' '.join(args.search).strip("'").strip('"')
OUTPUT_FILENAME = args.output

print("SEARCH QUERY:", SEARCH)
print("OUTPUT FILENAME:", OUTPUT_FILENAME)


async def get_all_ids_of_seq():
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "nucleotide",
        "term": SEARCH,
        "RetMax": 100,
        "RetStart": 0
    }

    res = []
    async with aiohttp.ClientSession() as session:

        while True:
            try:
                async with session.get(url, params=params) as result:
                    content = await result.content.read()
                    data = ET.fromstring(content)
                    ids = [id_.text for id_ in data[3]]

                    if len(ids) == 0:
                        break

                    res.extend(ids)
                    params["RetStart"] += params["RetMax"]
            except:
                pass

    return res

TOTAL = 0
DOWNLOAD = 0
DOWNLOAD_LOCK = Lock()


async def get_fasta_by_id(ids):
    global TOTAL, DOWNLOAD, DOWNLOAD_LOCK
    res = []

    for id_ in ids:

        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            "db": "nucleotide",
            "id": id_,
            "rettype": "fasta",
            "retmode": "text"
        }

        fasta = None
        async with aiohttp.ClientSession() as session:
            tries = 3
            while tries > 0:
                time.sleep(1)
                try:
                    async with session.get(url, params=params) as result:
                        if result.status != 200:
                            tries -= 1
                            continue

                        fasta = await result.content.read()
                        break
                except:
                    pass

        with DOWNLOAD_LOCK:
            DOWNLOAD += 1
            progress(DOWNLOAD, TOTAL)

        if fasta is None:
            print(f"Ошибка скачивания fasta. id сущности = {id_}")

        res.append(fasta)

    return res


async def main():
    global TOTAL, DOWNLOAD, DOWNLOAD_LOCK
    fastas_ids = await get_all_ids_of_seq()
    fastas_ids = fastas_ids
    TOTAL = len(fastas_ids)

    get_fasta_by_id_tasks_1 = asyncio.create_task(get_fasta_by_id(fastas_ids[0::3]))
    get_fasta_by_id_tasks_2 = asyncio.create_task(get_fasta_by_id(fastas_ids[1::3]))
    get_fasta_by_id_tasks_3 = asyncio.create_task(get_fasta_by_id(fastas_ids[2::3]))

    fastas = await asyncio.gather(
        get_fasta_by_id_tasks_1,
        get_fasta_by_id_tasks_2,
        get_fasta_by_id_tasks_3,
        return_exceptions=True
    )
    decode_b = lambda x: x.decode()
    fastas = ["".join(map(decode_b, x)) for x in fastas]
    fastas = "".join(fastas)

    file = open(OUTPUT_FILENAME, "w")
    file.write(fastas)
    file.close()

asyncio.run(main())
