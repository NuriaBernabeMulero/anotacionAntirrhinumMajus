# -*- coding: utf-8 -*-

import pandas as pd

# Valor máximo de e-value permitido
MAX_EV = 0.01

# Nombres de los ficheros de entrada y salida para los BLAST
nameFic_E_SDTo_tbx = 'salidaSDTomato_tbx'
nameFic_S_SDTo_tbx = 'salidaSDTomatoFiltrada_tbx.tsv'
nameFic_E_SDTa_tbx = 'salidaSDThaliana_tbx'
nameFic_S_SDTa_tbx = 'salidaSDThalianaFiltrada_tbx.tsv'
nameFic_E_To_tbx = 'salidaTomato_tbx'
nameFic_S_To_tbx = 'salidaTomatoFiltrada_tbx.tsv'
nameFic_E_Ta_tbx = 'salidaThaliana_tbx'
nameFic_S_Ta_tbx = 'salidaThalianaFiltrada_tbx.tsv'
nameFic_E_SDTo_dcmb = 'salidaSDTomato_dcmb'
nameFic_S_SDTo_dcmb = 'salidaSDTomatoFiltrada_dcmb.tsv'
nameFic_E_SDTa_dcmb = 'salidaSDThaliana_dcmb'
nameFic_S_SDTa_dcmb = 'salidaSDThalianaFiltrada_dcmb.tsv'
nameFic_E_To_dcmb = 'salidaTomato_dcmb'
nameFic_S_To_dcmb = 'salidaTomatoFiltrada_dcmb.tsv'
nameFic_E_Ta_dcmb = 'salidaThaliana_dcmb'
nameFic_S_Ta_dcmb = 'salidaThalianaFiltrada_dcmb.tsv'

# TBLASTX
# Tomato-Snapdragon
df_chunk_SDTomato_tbx = pd.read_csv(nameFic_E_SDTo_tbx, sep='\t', chunksize=1000000,
                                    names=['query acc.ver', 'subject acc.ver', '% identity',
                                           'alignment length', 'mismatches', 'gap opens',
                                           'q. start', 'q. end', 's. start',
                                           's. end', 'evalue', 'bit score'])

# Conservamos solo aquellos resultados que no superen el valor máximo de e-value
chunk_list_SDTomato_tbx = []

for chunk in df_chunk_SDTomato_tbx:
    print('tbx: tomato-snapdragon')
    chunk_filter = chunk[chunk['evalue'] <= MAX_EV]

    chunk_list_SDTomato_tbx.append(chunk_filter)

dfSDTomato_tbx = pd.concat(chunk_list_SDTomato_tbx)

# Escribimos los resultados que nos hemos quedado en el fichero de salida
dfSDTomato_tbx.to_csv(nameFic_S_SDTo_tbx, sep='\t', header=False)

# Thaliana-Snapdragon
df_chunk_SDThaliana_tbx = pd.read_csv(nameFic_E_SDTa_tbx, sep='\t', chunksize=1000000,
                                      names=['query acc.ver', 'subject acc.ver', '% identity',
                                             'alignment length', 'mismatches', 'gap opens',
                                             'q. start', 'q. end', 's. start',
                                             's. end', 'evalue', 'bit score'])

# Conservamos solo aquellos resultados que no superen el valor máximo de e-value
chunk_list_SDThaliana_tbx = []

for chunk in df_chunk_SDThaliana_tbx:
    print('tbx: thaliana-snapdragon')
    chunk_filter = chunk[chunk['evalue'] <= MAX_EV]

    chunk_list_SDThaliana_tbx.append(chunk_filter)

dfSDThaliana_tbx = pd.concat(chunk_list_SDThaliana_tbx)

# Escribimos los resultados que nos hemos quedado en el fichero de salida
dfSDThaliana_tbx.to_csv(nameFic_S_SDTa_tbx, sep='\t', header=False)

# Snapdragon-Tomato
df_chunk_Tomato_tbx = pd.read_csv(nameFic_E_To_tbx, sep='\t', chunksize=1000000,
                                  names=['query acc.ver', 'subject acc.ver', '% identity',
                                         'alignment length', 'mismatches', 'gap opens',
                                         'q. start', 'q. end', 's. start',
                                         's. end', 'evalue', 'bit score'])

# Conservamos solo aquellos resultados que no superen el valor máximo de e-value
chunk_list_Tomato_tbx = []

for chunk in df_chunk_Tomato_tbx:
    print('tbx: snapdragon-tomato')
    chunk_filter = chunk[chunk['evalue'] <= MAX_EV]

    chunk_list_Tomato_tbx.append(chunk_filter)

dfTomato_tbx = pd.concat(chunk_list_Tomato_tbx)

# Escribimos los resultados que nos hemos quedado en el fichero de salida
dfTomato_tbx.to_csv(nameFic_S_To_tbx, sep='\t', header=False)

# Snapdragon-Thaliana
df_chunk_Thaliana_tbx = pd.read_csv(nameFic_E_Ta_tbx, sep='\t', chunksize=1000000,
                                    names=['query acc.ver', 'subject acc.ver', '% identity',
                                           'alignment length', 'mismatches', 'gap opens',
                                           'q. start', 'q. end', 's. start',
                                           's. end', 'evalue', 'bit score'])

# Conservamos solo aquellos resultados que no superen el valor máximo de e-value
chunk_list_Thaliana_tbx = []

for chunk in df_chunk_Thaliana_tbx:
    print('tbx: snapdragon-thaliana')
    chunk_filter = chunk[chunk['evalue'] <= MAX_EV]

    chunk_list_Thaliana_tbx.append(chunk_filter)

dfThaliana_tbx = pd.concat(chunk_list_Thaliana_tbx)

# Escribimos los resultados que nos hemos quedado en el fichero de salida
dfThaliana_tbx.to_csv(nameFic_S_Ta_tbx, sep='\t', header=False)

# DC-MEGABLAST
# Tomato-Snapdragon
df_chunk_SDTomato_dcmb = pd.read_csv(nameFic_E_SDTo_dcmb, sep='\t', chunksize=1000000,
                                     names=['query acc.ver', 'subject acc.ver', '% identity',
                                            'alignment length', 'mismatches', 'gap opens',
                                            'q. start', 'q. end', 's. start',
                                            's. end', 'evalue', 'bit score'])

# Conservamos solo aquellos resultados que no superen el valor máximo de e-value
chunk_list_SDTomato_dcmb = []

for chunk in df_chunk_SDTomato_dcmb:
    print('dcmb: tomato-snapdragon')
    chunk_filter = chunk[chunk['evalue'] <= MAX_EV]

    chunk_list_SDTomato_dcmb.append(chunk_filter)

dfSDTomato_dcmb = pd.concat(chunk_list_SDTomato_dcmb)

# Escribimos los resultados que nos hemos quedado en el fichero de salida
dfSDTomato_dcmb.to_csv(nameFic_S_SDTo_dcmb, sep='\t', header=False)

# Thaliana-Snapdragon
df_chunk_SDThaliana_dcmb = pd.read_csv(nameFic_E_SDTa_dcmb, sep='\t', chunksize=1000000,
                                       names=['query acc.ver', 'subject acc.ver', '% identity',
                                              'alignment length', 'mismatches', 'gap opens',
                                              'q. start', 'q. end', 's. start',
                                              's. end', 'evalue', 'bit score'])

# Conservamos solo aquellos resultados que no superen el valor máximo de e-value
chunk_list_SDThaliana_dcmb = []

for chunk in df_chunk_SDThaliana_dcmb:
    print('dcmb: thaliana-snapdragon')
    chunk_filter = chunk[chunk['evalue'] <= MAX_EV]

    chunk_list_SDThaliana_dcmb.append(chunk_filter)

dfSDThaliana_dcmb = pd.concat(chunk_list_SDThaliana_dcmb)

# Escribimos los resultados que nos hemos quedado en el fichero de salida
dfSDThaliana_dcmb.to_csv(nameFic_S_SDTa_dcmb, sep='\t', header=False)

# Snapdragon-Tomato
df_chunk_Tomato_dcmb = pd.read_csv(nameFic_E_To_dcmb, sep='\t', chunksize=1000000,
                                   names=['query acc.ver', 'subject acc.ver', '% identity',
                                          'alignment length', 'mismatches', 'gap opens',
                                          'q. start', 'q. end', 's. start',
                                          's. end', 'evalue', 'bit score'])

# Conservamos solo aquellos resultados que no superen el valor máximo de e-value
chunk_list_Tomato_dcmb = []

for chunk in df_chunk_Tomato_dcmb:
    print('dcmb: snapdragon-tomato')
    chunk_filter = chunk[chunk['evalue'] <= MAX_EV]

    chunk_list_Tomato_dcmb.append(chunk_filter)

dfTomato_dcmb = pd.concat(chunk_list_Tomato_dcmb)

# Escribimos los resultados que nos hemos quedado en el fichero de salida
dfTomato_dcmb.to_csv(nameFic_S_To_dcmb, sep='\t', header=False)

# Snapdragon-Thaliana
df_chunk_Thaliana_dcmb = pd.read_csv(nameFic_E_Ta_dcmb, sep='\t', chunksize=1000000,
                                     names=['query acc.ver', 'subject acc.ver', '% identity',
                                            'alignment length', 'mismatches', 'gap opens',
                                            'q. start', 'q. end', 's. start',
                                            's. end', 'evalue', 'bit score'])

# Conservamos solo aquellos resultados que no superen el valor máximo de e-value
chunk_list_Thaliana_dcmb = []

for chunk in df_chunk_Thaliana_dcmb:
    print('dcmb: thaliana-snapdragon')
    chunk_filter = chunk[chunk['evalue'] <= MAX_EV]

    chunk_list_Thaliana_dcmb.append(chunk_filter)

dfThaliana_dcmb = pd.concat(chunk_list_Thaliana_dcmb)

# Escribimos los resultados que nos hemos quedado en el fichero de salida
dfThaliana_dcmb.to_csv(nameFic_S_Ta_dcmb, sep='\t', header=False)

# Mostramos por pantalla el número de resultados que quedan en cada dataframe
print('sd to tbx: ' + str(len(dfSDTomato_tbx)))
print('sd ta tbx: ' + str(len(dfSDThaliana_tbx)))
print(' to tbx: ' + str(len(dfTomato_tbx)))
print(' ta tbx: ' + str(len(dfThaliana_tbx)))
print('sd to dcmb: ' + str(len(dfSDTomato_dcmb)))
print('sd ta dcmb: ' + str(len(dfSDThaliana_dcmb)))
print(' to dcmb: ' + str(len(dfTomato_dcmb)))
print(' ta dcmb: ' + str(len(dfThaliana_dcmb)))
