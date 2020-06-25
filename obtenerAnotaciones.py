import pandas as pd
from pybiomart import Dataset

nameFicETo = 'D:\\CosasTFG\\sScTomato.tsv'
nameFicETa = 'D:\\CosasTFG\\sScThaliana.tsv'
nameFicATo = 'D:\\CosasTFG\\anotacionesTomato.tsv'
nameFicATa = 'D:\\CosasTFG\\anotacionesThaliana.tsv'

# Snapdragon-Tomate
dfSalidaTomato = pd.read_csv(nameFicETo, sep='\t',
                             names=['line', 'query acc.ver', 'subject acc.ver', '% identity',
                                    'alignment length', 'mismatches', 'gap opens',
                                    'q. start', 'q. end', 's. start',
                                    's. end', 'evalue', 'bit score'],
                             usecols=['query acc.ver', 'subject acc.ver', '% identity',
                                      'alignment length', 'mismatches', 'gap opens',
                                      'q. start', 'q. end', 's. start',
                                      's. end', 'evalue', 'bit score'])

# Convertimos en un mapa los resultados
matchesTomato = {}

for i, row in dfSalidaTomato.iterrows():
    genSnapdragon = row['query acc.ver']
    genTomato = row['subject acc.ver']
    evalue = row['evalue']

    matchesTomato[genSnapdragon] = (genTomato, evalue)

# Snapdragon-Thaliana
dfSalidaThaliana = pd.read_csv(nameFicETa, sep='\t',
                               names=['line', 'query acc.ver', 'subject acc.ver', '% identity',
                                      'alignment length', 'mismatches', 'gap opens',
                                      'q. start', 'q. end', 's. start',
                                      's. end', 'evalue', 'bit score'],
                               usecols=['query acc.ver', 'subject acc.ver', '% identity',
                                        'alignment length', 'mismatches', 'gap opens',
                                        'q. start', 'q. end', 's. start',
                                        's. end', 'evalue', 'bit score'])

# Convertimos en un mapa los resultados
matchesThaliana = {}

for i, row in dfSalidaThaliana.iterrows():
    genSnapdragon = row['query acc.ver']
    genThaliana = row['subject acc.ver']
    evalue = row['evalue']

    matchesThaliana[genSnapdragon] = (genThaliana, evalue)

# Ahora sacamos las anotaciones y las guardamos en un mapa

# Tomato
anotacionesTomato = {}

datasetTomato = Dataset(name='slycopersicum_eg_gene',
                        virtual_schema='plants_mart',
                        host='http://plants.ensembl.org')

# Para cada gen, lanzamos una consulta a Biomart
for genSnapdragon in matchesTomato:
    print('tomato')
    genTomato, evalue = matchesTomato[genSnapdragon]
    resultTomato = datasetTomato.query(
        attributes=['ensembl_gene_id', 'ensembl_transcript_id', 'go_id', 'go_linkage_type', 'namespace_1003'],
        filters={'link_ensembl_transcript_stable_id': genTomato})

    # Si encuentra match:
    if len(resultTomato) > 0:
        # Comprobamos que no haya valores nan
        lr = resultTomato['GO term accession'].tolist()
        lr = [(e, evalue) for e in lr if str(e) != 'nan']
        # Si sigue habiendo anotaciones, las guardamos
        if len(lr) > 0:
            anotacionesTomato[genSnapdragon] = lr

# Thaliana
anotacionesThaliana = {}

datasetThaliana = Dataset(name='athaliana_eg_gene',
                          virtual_schema='plants_mart',
                          host='http://plants.ensembl.org')

# Para cada gen, lanzamos una consulta a Biomart
for genSnapdragon in matchesThaliana:
    print('h')
    genThaliana, evalue = matchesThaliana[genSnapdragon]
    resultThaliana = datasetThaliana.query(
        attributes=['ensembl_gene_id', 'ensembl_transcript_id', 'go_id', 'go_linkage_type', 'namespace_1003'],
        filters={'link_ensembl_transcript_stable_id': genThaliana})

    # Si encuentra match:
    if len(resultThaliana) > 0:
        # Comprobamos que no haya valores nan
        lr = resultThaliana['GO term accession'].tolist()
        lr = [(e, evalue) for e in lr if str(e) != 'nan']
        # Si sigue habiendo anotaciones, las guardamos
        if len(lr) > 0:
            anotacionesThaliana[genSnapdragon] = lr

# Abrimos el fichero de salida y escribimos las anotaciones en formato
# genSnapdragon \t anotacion \t e-valor (del gen del que se saca la anotación)
fATo = open(nameFicATo, 'w+', buffering=1)

for gen in anotacionesTomato:
    l = anotacionesTomato[gen]
    for anotacion, evalue in l:
        fATo.write(gen + '\t' + anotacion + '\t' + str(evalue) + '\n')

# Abrimos el fichero de salida y escribimos las anotaciones en formato
# genSnapdragon \t anotacion \t e-valor (del gen del que se saca la anotación)
fATa = open(nameFicATa, 'w+', buffering=1)

for gen in anotacionesThaliana:
    l = anotacionesThaliana[gen]
    for anotacion, evalue in l:
        fATa.write(gen + '\t' + anotacion + '\t' + str(evalue) + '\n')
