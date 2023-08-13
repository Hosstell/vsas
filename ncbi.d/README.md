NCBI.D: Скачивание fasta данных с сайта ncbi
Источник: https://hub.docker.com/r/vsas/ncbi.d
Новая версия: **1.1.0**

**--- Запуск ---**
**Команда запуска:**
``` docker run --rm -v "$(pwd):/in" vsas/ncbi.d '("Potato virus Y"[Organism] OR ("Potato virus Y"[Organism] OR potato virus Y[All Fields])) AND ("2018/01/01"[PDAT] : "2023/12/31"[PDAT]) AND ("9000"[SLEN] : "9100"[SLEN])' -o potato_virus_y.txt```

**--- PS ---**  
**Чтобы обновить локальный образ:**   
``` docker pull vsas/ncbi.d:latest ```