Расчет количества определенных нуклеотидов в начале ридов  
Источник: https://hub.docker.com/r/vsas/freads  
Новая версия: **1.6.0**

**--- Запуск ---**  
**Команда запуска:**  
``` docker run --rm -v "$(pwd):/in" vsas/freads -f ./miRNA_S7872Nr2.1.fastq.gz.sam -c 3 ```

-f - Входной .sam файл
-с - Количество первых нуклеотидов к риде, которые мы учитываем в результате. [Не обязательный. Значение по умолчанию - 1]
-o - Префикс названий выходных файлов. Он нужен при обработке нескольких файлов. 

После выполнения команды, в папке появится файлы:
- miRNA_S7872Nr2.1.fastq.gz.sam.c3.0.5.png (3 первых нуклеотида для последовательностей с флагом 0 в sam-файле)
- miRNA_S7872Nr2.1.fastq.gz.sam.c3.0.3.png (3 последних нуклеотида для последовательностей с флагом 0 в sam-файле)
- miRNA_S7872Nr2.1.fastq.gz.sam.c3.16.5.png (3 первых нуклеотида для последовательностей с флагом 16 в sam-файле)
- miRNA_S7872Nr2.1.fastq.gz.sam.c3.16.3.png (3 последних нуклеотида для последовательностей с флагом 16 в sam-файле)

**--- Изменение ---**  
**Параметр -o**  
Команда: 
``` 
docker run --rm -v "$(pwd):/in" vsas/freads -f ./miRNA_S7872Nr2.1.fastq.gz.sam 
```  
Выходные файлы:
```
- miRNA_S7872Nr2.1.fastq.gz.sam.c3.0.5.png
- miRNA_S7872Nr2.1.fastq.gz.sam.c3.0.3.png
- miRNA_S7872Nr2.1.fastq.gz.sam.c3.16.5.png
- miRNA_S7872Nr2.1.fastq.gz.sam.c3.16.3.png
```

Команда:   
``` 
docker run --rm -v "$(pwd):/in" vsas/freads -f ./miRNA_S7872Nr2.1.fastq.gz.sam -o newname
```  
Выходные файлы:
```
- newname.c3.0.5.png
- newname.c3.0.3.png
- newname.c3.16.5.png
- newname.c3.16.3.png
```

**Возможность обрабатываться несколько файлов в один результат**  
Пример:
```
docker run --rm -v "$(pwd):/in" vsas/freads -f miRNA_S7872Nr2.1.fastq.gz.sam,miRNA_S7872Nr3.1.fastq.gz.sam -o -o newname
```

**--- PS ---**  
**Чтобы обновить локальный образ:**   
``` docker pull vsas/freads:latest ```