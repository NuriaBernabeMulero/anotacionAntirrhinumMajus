# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

# Valor mínimo para el porcentaje de identidad
MIN_PID = 70

# Nombre de los ficheros en los que almacenaremos la salida
nameFicSTo = 'sScTomato_70.tsv'
nameFicSTa = 'sScThaliana_70.tsv'

# Nombre de los ficheros de entrada
nameFicGenesNoAnotados = 'genes_no_anotados.tsv'
nameFic_E_SDTo_dcmb = 'salidaSDTomatoFiltrada_dcmb.tsv'
nameFic_E_To_dcmb = 'salidaTomatoFiltrada_dcmb.tsv'
nameFic_E_SDTa_dcmb = 'salidaSDThalianaFiltrada_dcmb.tsv'
nameFic_E_Ta_dcmb = 'salidaThalianaFiltrada_dcmb.tsv'
nameFic_E_SDTo_tbx = 'salidaSDTomatoFiltrada_tbx.tsv'
nameFic_E_To_tbx = 'salidaTomatoFiltrada_tbx.tsv'
nameFic_E_SDTa_tbx = 'salidaSDThalianaFiltrada_tbx.tsv'
nameFic_E_Ta_tbx = 'salidaThalianaFiltrada_tbx.tsv'

# Leemos el fichero de genes_no_anotados y lo convertimos en una lista
genesNoAnotados = [g.strip() for g in open(nameFicGenesNoAnotados).readlines()]

# DC-MEGABLAST
# Tomato-Snapdragon
dfSalidaSDTomato_dcmb = pd.read_csv(nameFic_E_SDTo_dcmb, sep='\t',
                                    names=['line', 'query acc.ver', 'subject acc.ver', '% identity',
                                           'alignment length', 'mismatches', 'gap opens',
                                           'q. start', 'q. end', 's. start',
                                           's. end', 'evalue', 'bit score'],
                                    usecols=['query acc.ver', 'subject acc.ver', '% identity',
                                             'alignment length', 'mismatches', 'gap opens',
                                             'q. start', 'q. end', 's. start',
                                             's. end', 'evalue', 'bit score'])

# Convertimos los resultados tomato-snapdragon que pasen el umbral en un mapa
mSDTomato_dcmb = {}

for i, row in dfSalidaSDTomato_dcmb.iterrows():
    print('dcmb: tomato-snapdragon = ' + str(i))
    genTomato = row['query acc.ver']
    genSnapdragon = row['subject acc.ver']
    pID = row['% identity']
    if pID >= MIN_PID:
        if genSnapdragon in mSDTomato_dcmb:
            mSDTomato_dcmb[genSnapdragon].append(genTomato)
        else:
            mSDTomato_dcmb[genSnapdragon] = [genTomato]

# Snapdragon-Tomato
dfSalidaTomato_dcmb = pd.read_csv(nameFic_E_To_dcmb, sep='\t',
                                  names=['line', 'query acc.ver', 'subject acc.ver', '% identity',
                                         'alignment length', 'mismatches', 'gap opens',
                                         'q. start', 'q. end', 's. start',
                                         's. end', 'evalue', 'bit score'],
                                  usecols=['query acc.ver', 'subject acc.ver', '% identity',
                                           'alignment length', 'mismatches', 'gap opens',
                                           'q. start', 'q. end', 's. start',
                                           's. end', 'evalue', 'bit score'])

# Creamos un mapa donde guardaremos los resultados snapdragon-tomato
matchesTomato_dcmb = {}

# Creamos un nuevo dataframe donde introduciremos aquellos resultados que pasen los filtros
dfSalidaTomatoFiltrada_dcmb = pd.DataFrame(columns=['query acc.ver', 'subject acc.ver', '% identity',
                                                    'alignment length', 'mismatches', 'gap opens',
                                                    'q. start', 'q. end', 's. start',
                                                    's. end', 'evalue', 'bit score'])

# Extraemos una lista de todos los genes de snapdragon para los que se encuentra resultado
# Para cada uno de ellos se cogen las parejas asociadas, se ordenan por e-value y se elige la mejor cuyo porcentaje de
# identidad sea igual o mayor al umbral y sea recíproca
l = list(np.unique(list(dfSalidaTomato_dcmb['query acc.ver'])))

for genSnapdragon in l:
    print('dcmb: snapdragon-tomato')
    df = dfSalidaTomato_dcmb.loc[dfSalidaTomato_dcmb['query acc.ver'] == genSnapdragon]
    df = df.sort_values('evalue')
    for i, row in df.iterrows():
        print('dcmb: snapdragon-tomato = ' + str(i))
        genTomato = row['subject acc.ver']
        evalue = row['evalue']
        pID = row['% identity']
        if genSnapdragon in mSDTomato_dcmb and genTomato in mSDTomato_dcmb[genSnapdragon] and pID >= MIN_PID:
            print('dcmb: snapdragon-tomato: found ' + genSnapdragon)
            matchesTomato_dcmb[genSnapdragon] = (genTomato, evalue)
            dfSalidaTomatoFiltrada_dcmb = dfSalidaTomatoFiltrada_dcmb.append(row)
            break

# Thaliana-Snapdragon
dfSalidaSDThaliana_dcmb = pd.read_csv(nameFic_E_SDTa_dcmb, sep='\t',
                                      names=['line', 'query acc.ver', 'subject acc.ver', '% identity',
                                             'alignment length', 'mismatches', 'gap opens',
                                             'q. start', 'q. end', 's. start',
                                             's. end', 'evalue', 'bit score'],
                                      usecols=['query acc.ver', 'subject acc.ver', '% identity',
                                               'alignment length', 'mismatches', 'gap opens',
                                               'q. start', 'q. end', 's. start',
                                               's. end', 'evalue', 'bit score'])

# Convertimos los resultados thaliana-snapdragon que pasen el umbral en un mapa
mSDThaliana_dcmb = {}

for i, row in dfSalidaSDThaliana_dcmb.iterrows():
    print('dcmb: thaliana-snapdragon = ' + str(i))
    genThaliana = row['query acc.ver']
    genSnapdragon = row['subject acc.ver']
    pID = row['% identity']
    if pID >= MIN_PID:
        if genSnapdragon in mSDThaliana_dcmb:
            mSDThaliana_dcmb[genSnapdragon].append(genThaliana)
        else:
            mSDThaliana_dcmb[genSnapdragon] = [genThaliana]

# Snapdragon-Thaliana
dfSalidaThaliana_dcmb = pd.read_csv(nameFic_E_Ta_dcmb, sep='\t',
                                    names=['line', 'query acc.ver', 'subject acc.ver', '% identity',
                                           'alignment length', 'mismatches', 'gap opens',
                                           'q. start', 'q. end', 's. start',
                                           's. end', 'evalue', 'bit score'],
                                    usecols=['query acc.ver', 'subject acc.ver', '% identity',
                                             'alignment length', 'mismatches', 'gap opens',
                                             'q. start', 'q. end', 's. start',
                                             's. end', 'evalue', 'bit score'])

# Creamos un mapa donde guardaremos los resultados snapdragon-thaliana
matchesThaliana_dcmb = {}

# Creamos un nuevo dataframe donde introduciremos aquellos resultados que pasen los filtros
dfSalidaThalianaFiltrada_dcmb = pd.DataFrame(columns=['query acc.ver', 'subject acc.ver', '% identity',
                                                      'alignment length', 'mismatches', 'gap opens',
                                                      'q. start', 'q. end', 's. start',
                                                      's. end', 'evalue', 'bit score'])

# Extraemos una lista de todos los genes de snapdragon para los que se encuentra resultado
# Para cada uno de ellos se cogen las parejas asociadas, se ordenan por e-value y se elige la mejor cuyo porcentaje de
# identidad sea igual o mayor al umbral y sea recíproca
l = list(np.unique(list(dfSalidaThaliana_dcmb['query acc.ver'])))

for genSnapdragon in l:
    print('dcmb: snapdragon-thaliana')
    df = dfSalidaThaliana_dcmb.loc[dfSalidaThaliana_dcmb['query acc.ver'] == genSnapdragon]
    df = df.sort_values('evalue')
    for i, row in df.iterrows():
        print('dcmb: snapdragon-thaliana = ' + str(i))
        genThaliana = row['subject acc.ver']
        evalue = row['evalue']
        pID = row['% identity']
        if genSnapdragon in mSDThaliana_dcmb and genThaliana in mSDThaliana_dcmb[genSnapdragon] and pID >= MIN_PID:
            print('dcmb: snapdragon-thaliana: found ' + genSnapdragon)
            matchesThaliana_dcmb[genSnapdragon] = (genThaliana, evalue)
            dfSalidaThalianaFiltrada_dcmb = dfSalidaThalianaFiltrada_dcmb.append(row)
            break

# TBLASTX
# Tomato-Snapdragon
df_chunk_SalidaSDTomato_tbx = pd.read_csv(nameFic_E_SDTo_tbx, sep='\t', chunksize=1000000,
                                          names=['line', 'query acc.ver', 'subject acc.ver', '% identity',
                                                 'alignment length', 'mismatches', 'gap opens',
                                                 'q. start', 'q. end', 's. start',
                                                 's. end', 'evalue', 'bit score'],
                                          usecols=['query acc.ver', 'subject acc.ver', '% identity',
                                                   'alignment length', 'mismatches', 'gap opens',
                                                   'q. start', 'q. end', 's. start',
                                                   's. end', 'evalue', 'bit score'])

# Convertimos los resultados tomato-snapdragon que pasen el umbral en un mapa
mSDTomato_tbx = {}

for chunk in df_chunk_SalidaSDTomato_tbx:
    print('tbx: tomato-snapdragon')
    for i, row in chunk.iterrows():
        print('tbx: tomato-snapdragon = ' + str(i))
        genTomato = row['query acc.ver']
        genSnapdragon = row['subject acc.ver']
        pID = row['% identity']
        if pID >= MIN_PID:
            if genSnapdragon in mSDTomato_tbx:
                mSDTomato_tbx[genSnapdragon].append(genTomato)
            else:
                mSDTomato_tbx[genSnapdragon] = [genTomato]

# Snapdragon-Tomato
df_chunk_SalidaTomato_tbx = pd.read_csv(nameFic_E_To_tbx, sep='\t', chunksize=1000000,
                                        names=['line', 'query acc.ver', 'subject acc.ver', '% identity',
                                               'alignment length', 'mismatches', 'gap opens',
                                               'q. start', 'q. end', 's. start',
                                               's. end', 'evalue', 'bit score'],
                                        usecols=['query acc.ver', 'subject acc.ver', '% identity',
                                                 'alignment length', 'mismatches', 'gap opens',
                                                 'q. start', 'q. end', 's. start',
                                                 's. end', 'evalue', 'bit score'])

# Creamos un mapa donde guardaremos los resultados snapdragon-tomato
matchesTomato_tbx = {}

# Creamos un nuevo dataframe donde introduciremos aquellos resultados que pasen los filtros
dfSalidaTomatoFiltrada_tbx = pd.DataFrame(columns=['query acc.ver', 'subject acc.ver', '% identity',
                                                   'alignment length', 'mismatches', 'gap opens',
                                                   'q. start', 'q. end', 's. start',
                                                   's. end', 'evalue', 'bit score'])

for chunk in df_chunk_SalidaTomato_tbx:
    # Para cada chunk, extraemos una lista de todos los genes de snapdragon para los que se encuentra resultado
    # Para cada uno de ellos se cogen las parejas asociadas, se ordenan por e-value y se elige la mejor cuyo porcentaje
    # de identidad sea igual o mayor al umbral y sea recíproca
    l = list(np.unique(list(chunk['query acc.ver'])))

    for genSnapdragon in l:
        print('tbx: snapdragon-tomato')
        df = chunk.loc[chunk['query acc.ver'] == genSnapdragon]
        df = df.sort_values('evalue')
        for i, row in df.iterrows():
            print('tbx: snapdragon-tomato = ' + str(i))
            genTomato = row['subject acc.ver']
            evalue = row['evalue']
            pID = row['% identity']
            if genSnapdragon in mSDTomato_tbx and genTomato in mSDTomato_tbx[genSnapdragon] and pID >= MIN_PID:
                print('tbx: snpadragon-tomato: found ' + genSnapdragon)
                if genSnapdragon not in matchesTomato_tbx:
                    matchesTomato_tbx[genSnapdragon] = (genTomato, evalue)
                    dfSalidaTomatoFiltrada_tbx = dfSalidaTomatoFiltrada_tbx.append(row)
                # Ahora es necesario tener en cuenta que podría haberse encontrado un resultado mejor en otro chunk
                else:
                    m, e = matchesTomato_tbx[genSnapdragon]
                    if e > evalue:
                        matchesTomato_tbx[genSnapdragon] = (genTomato, evalue)

                        dfSalidaTomatoFiltrada_tbx = dfSalidaTomatoFiltrada_tbx.loc[
                            dfSalidaTomatoFiltrada_tbx['query acc.ver'] != genSnapdragon]

                        dfSalidaTomatoFiltrada_tbx = dfSalidaTomatoFiltrada_tbx.append(row)
                break

# Thaliana-Snapdragon
df_chunk_SalidaSDThaliana_tbx = pd.read_csv(nameFic_E_SDTa_tbx, sep='\t', chunksize=1000000,
                                            names=['line', 'query acc.ver', 'subject acc.ver', '% identity',
                                                   'alignment length', 'mismatches', 'gap opens',
                                                   'q. start', 'q. end', 's. start',
                                                   's. end', 'evalue', 'bit score'],
                                            usecols=['query acc.ver', 'subject acc.ver', '% identity',
                                                     'alignment length', 'mismatches', 'gap opens',
                                                     'q. start', 'q. end', 's. start',
                                                     's. end', 'evalue', 'bit score'])

# Convertimos los resultados thaliana-snapdragon que pasen el umbral en un mapa
mSDThaliana_tbx = {}

for chunk in df_chunk_SalidaSDThaliana_tbx:
    print('tbx: thaliana-snapdragon')
    for i, row in chunk.iterrows():
        print('tbx: thaliana-snapdragon = ' + str(i))
        genThaliana = row['query acc.ver']
        genSnapdragon = row['subject acc.ver']
        pID = row['% identity']
        if pID >= MIN_PID:
            if genSnapdragon in mSDThaliana_tbx:
                mSDThaliana_tbx[genSnapdragon].append(genThaliana)
            else:
                mSDThaliana_tbx[genSnapdragon] = [genThaliana]

# Snapdragon-Thaliana
df_chunk_SalidaThaliana_tbx = pd.read_csv(nameFic_E_Ta_tbx, sep='\t', chunksize=1000000,
                                          names=['line', 'query acc.ver', 'subject acc.ver', '% identity',
                                                 'alignment length', 'mismatches', 'gap opens',
                                                 'q. start', 'q. end', 's. start',
                                                 's. end', 'evalue', 'bit score'],
                                          usecols=['query acc.ver', 'subject acc.ver', '% identity',
                                                   'alignment length', 'mismatches', 'gap opens',
                                                   'q. start', 'q. end', 's. start',
                                                   's. end', 'evalue', 'bit score'])

# Creamos un mapa donde guardaremos los resultados snapdragon-thaliana
matchesThaliana_tbx = {}

# Creamos un nuevo dataframe donde introduciremos aquellos resultados que pasen los filtros
dfSalidaThalianaFiltrada_tbx = pd.DataFrame(columns=['query acc.ver', 'subject acc.ver', '% identity',
                                                     'alignment length', 'mismatches', 'gap opens',
                                                     'q. start', 'q. end', 's. start',
                                                     's. end', 'evalue', 'bit score'])

for chunk in df_chunk_SalidaThaliana_tbx:
    # Para cada chunk, extraemos una lista de todos los genes de snapdragon para los que se encuentra resultado
    # Para cada uno de ellos se cogen las parejas asociadas, se ordenan por e-value y se elige la mejor cuyo porcentaje
    # de identidad sea igual o mayor al umbral y sea recíproca
    l = list(np.unique(list(chunk['query acc.ver'])))

    for genSnapdragon in l:
        print('tbx: snapdragon-thaliana')
        df = chunk.loc[chunk['query acc.ver'] == genSnapdragon]
        df = df.sort_values('evalue')
        for i, row in df.iterrows():
            print('tbx: snapdragon-thaliana = ' + str(i))
            genThaliana = row['subject acc.ver']
            evalue = row['evalue']
            pID = row['% identity']
            if genSnapdragon in mSDThaliana_tbx and genThaliana in mSDThaliana_tbx[genSnapdragon] and pID >= MIN_PID:
                print('tbx: snapdragon-thaliana: found ' + genSnapdragon)
                if genSnapdragon not in matchesThaliana_tbx:
                    matchesThaliana_tbx[genSnapdragon] = (genThaliana, evalue)
                    dfSalidaThalianaFiltrada_tbx = dfSalidaThalianaFiltrada_tbx.append(row)
                # Ahora es necesario tener en cuenta que podría haberse encontrado un resultado mejor en otro chunk
                else:
                    m, e = matchesThaliana_tbx[genSnapdragon]
                    if e > evalue:
                        matchesThaliana_tbx[genSnapdragon] = (genThaliana, evalue)

                        dfSalidaThalianaFiltrada_tbx = dfSalidaThalianaFiltrada_tbx.loc[
                            dfSalidaThalianaFiltrada_tbx['query acc.ver'] != genSnapdragon]

                        dfSalidaThalianaFiltrada_tbx = dfSalidaThalianaFiltrada_tbx.append(row)
                break

# Ahora vamos a seleccionar el mejor resultado de todos los propuestos por cada algoritmo para cada gen

# Variables para contabilizar cuántos genes provienen de cada programa
nToTBX = 0
nToDCMB = 0
nTaTBX = 0
nTaDCMB = 0

# Tomate
dfSalidaTomato = pd.DataFrame(columns=['query acc.ver', 'subject acc.ver', '% identity',
                                       'alignment length', 'mismatches', 'gap opens',
                                       'q. start', 'q. end', 's. start',
                                       's. end', 'evalue', 'bit score'])

# Iteramos sobre los resultados del tblastx, ya que son más
# Para cada uno comprobamos si también tiene resultado con dc-megablast y de ser así nos quedamos el de mejor e-value
for i, row in dfSalidaTomatoFiltrada_tbx.iterrows():
    print('tomato: ' + str(i))
    genSnapdragon = row['query acc.ver']
    genTomato = row['subject acc.ver']
    evalue = row['evalue']

    if genSnapdragon in matchesTomato_dcmb:
        m, e = matchesTomato_dcmb[genSnapdragon]
        if e < evalue:
            dfSalidaTomato = dfSalidaTomato.append(
                dfSalidaTomatoFiltrada_dcmb.loc[dfSalidaTomatoFiltrada_dcmb['query acc.ver'] == genSnapdragon])
            nToDCMB += 1
        else:
            dfSalidaTomato = dfSalidaTomato.append(row)
            nToTBX += 1
    else:
        dfSalidaTomato = dfSalidaTomato.append(row)
        nToTBX += 1

# Thaliana
dfSalidaThaliana = pd.DataFrame(columns=['query acc.ver', 'subject acc.ver', '% identity',
                                         'alignment length', 'mismatches', 'gap opens',
                                         'q. start', 'q. end', 's. start',
                                         's. end', 'evalue', 'bit score'])

# Iteramos sobre los resultados del tblastx, ya que son más
# Para cada uno comprobamos si también tiene resultado con dc-megablast y de ser así nos quedamos el de mejor e-value
for i, row in dfSalidaThalianaFiltrada_tbx.iterrows():
    print('thaliana: ' + str(i))
    genSnapdragon = row['query acc.ver']
    genThaliana = row['subject acc.ver']
    evalue = row['evalue']

    if genSnapdragon in matchesThaliana_dcmb:
        m, e = matchesThaliana_dcmb[genSnapdragon]
        if e < evalue:
            dfSalidaThaliana = dfSalidaThaliana.append(
                dfSalidaThalianaFiltrada_dcmb.loc[dfSalidaThalianaFiltrada_dcmb['query acc.ver'] == genSnapdragon])
            nTaDCMB += 1
        else:
            dfSalidaThaliana = dfSalidaThaliana.append(row)
            nTaTBX += 1
    else:
        dfSalidaThaliana = dfSalidaThaliana.append(row)
        nTaTBX += 1

# Comprobamos si hay genes de tomate que solo tienen resultado con dc-megablast y de ser así los añadimos
l = [gen for gen in matchesTomato_dcmb if gen not in matchesTomato_tbx]

if len(l) > 0:
    for genSnapdragon in l:
        dfSalidaTomato.append(
            dfSalidaTomatoFiltrada_dcmb.loc[dfSalidaTomatoFiltrada_dcmb['query acc.ver'] == genSnapdragon])
        nToDCMB += 1

# Comprobamos si hay genes de thaliana que solo tienen resultado con dc-megablast y de ser así los añadimos
l = [gen for gen in matchesThaliana_dcmb if gen not in matchesThaliana_tbx]

if len(l) > 0:
    for genSnapdragon in l:
        dfSalidaThaliana.append(
            dfSalidaThalianaFiltrada_dcmb.loc[dfSalidaThalianaFiltrada_dcmb['query acc.ver'] == genSnapdragon])
        nTaDCMB += 1

# Escribimos los ficheros con los resultados finales
dfSalidaTomato.to_csv(nameFicSTo, sep='\t', header=False)
dfSalidaThaliana.to_csv(nameFicSTa, sep='\t', header=False)

# Mostramos por pantalla algunos datos
print('genes no anotados que quedan sin anotar: ' +
      str(len([gen for gen in genesNoAnotados if (gen not in matchesTomato_tbx and gen not in matchesTomato_dcmb) and
               (gen not in matchesThaliana_tbx and gen not in matchesThaliana_dcmb)])))

print('tomato tbx: ' + str(nToTBX))
print('tomato dcmb: ' + str(nToDCMB))
print('thaliana tbx: ' + str(nTaTBX))
print('thaliana dcmb: ' + str(nTaDCMB))
