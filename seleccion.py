import numpy as np
import pandas as pd
from goatools.obo_parser import GODag

# Umbral máximo para el e-value
MAX_EV = 0.01

# Nombre de los ficheros de entrada y salida
nameFicGenesAnotados = 'D:\\CosasTFG\\anotacionGO.tsv'
nameFicGenesNoAnotados = 'D:\\CosasTFG\\genes_no_anotados.tsv'
nameFicATo = 'D:\\CosasTFG\\anotacionesTomato.tsv'
nameFicATa = 'D:\\CosasTFG\\anotacionesThaliana.tsv'
nameFicS = 'D:\\CosasTFG\\salidaScript.tsv'

# Leemos el fichero de anotacionesGO y lo convertimos en un diccionario
anotacionesSnapdragon = {}

with open(nameFicGenesAnotados) as f:
    for line in f:
        row = line.strip().split('\t')
        gen = row[0]
        goTerm = row[1]
        if gen in anotacionesSnapdragon:
            anotacionesSnapdragon[gen].append(goTerm)
        else:
            anotacionesSnapdragon[gen] = [goTerm]

# Leemos el fichero de genes_no_anotados y lo convertimos en una lista
genesNoAnotados = [g.strip() for g in open(nameFicGenesNoAnotados).readlines()]

# Leemos las anotaciones y las convertimos en mapas
# Podríamos modificar el parámetro de e-value máximo para solo aceptar anotaciones provenientes de parejas muy similares
# Para el resto del proceso no utilizamos el e-value

# Tomato
dfTomato = pd.read_csv(nameFicATo, sep='\t',
                       names=['gen', 'anotacion', 'evalue'])

anotacionesTomato = {}

for i, row in dfTomato.iterrows():
    gen = row['gen']
    anotacion = row['anotacion']
    ev = row['evalue']
    
    if ev <= MAX_EV:
        if gen in anotacionesTomato:
            anotacionesTomato[gen].append(anotacion)
        else:
            anotacionesTomato[gen] = [anotacion]

# Thaliana
dfThaliana = pd.read_csv(nameFicATa, sep='\t',
                         names=['gen', 'anotacion', 'evalue'])

anotacionesThaliana = {}

for i, row in dfThaliana.iterrows():
    gen = row['gen']
    anotacion = row['anotacion']
    ev = row['evalue']

    if ev <= MAX_EV:
        if gen in anotacionesThaliana:
            anotacionesThaliana[gen].append(anotacion)
        else:
            anotacionesThaliana[gen] = [anotacion]

# Ahora elegimos las posibles anotaciones
godag = GODag("D:\\CosasTFG\\go-basic.obo")

# Recomendación: dejar descomentada solo una de las dos partes (evaluación o anotación de no anotados)

######################## Para comparar con los ya anotados ###########################################
# Variables para contabilizar anotaciones no encontradas, FN, FP, TP, TN
'''SDnnP = 0
SDnP = 0
PnSD = 0
SDsP = 0
nSDnL = 0

# Utilizamos solo los genes anotados
# Para cada gen extraemos sus anotaciones
for gen in anotacionesSnapdragon:
    print('evaluación: ' + gen)

    aSD = anotacionesSnapdragon[gen]
    aTo = []
    aTa = []

    if gen in anotacionesTomato:
        aTo = anotacionesTomato[gen]

    if gen in anotacionesThaliana:
        aTa = anotacionesThaliana[gen]

    aP = list(np.unique(aTo + aTa))

    # Comprobamos cuáles de las que tiene asignadas no están entre las posibles
    SDnnP += len([gen for gen in aSD if gen not in aP])

    l = []

    if len(aP) > 0:
        # Si solo encuentra anotaciones en un genoma, las copiamos
        if len(aTo) == 0 or len(aTa) == 0:
            l = aP
        else:
            # Si encuentra en ambos genomas, copiamos las que estén en los dos
            lb = [a for a in aP if a in aTo and a in aTa]
            aPs1 = [a for a in aP if (a not in aTo or a not in aTa) and a in godag]

            ######### Sección para evaluar el método Ayllón ############
            # Comentar para evaluar el método conservador
            # En esta sección comprobamos si una anotación tiene descendientes entre las propuestas y de no ser así
            # la copia también
            for a in aPs1:
                for b in aPs1:
                    if a != b and len(common_parent_go_ids([a, b], godag)) != 0:
                        go_root = deepest_common_ancestor([a, b], godag)
                        n = semantic_distance(a, b, godag)

                        if go_root == a or go_root == b:
                            if go_root == a:
                                if b not in l:
                                    l.append(b)
                                    
                                if a in l:
                                    l.remove(a)
                            else:
                                if a not in l:
                                    l.append(a)
                                
                                if b in l:
                                    l.remove(b)
                        else:
                            if a not in l:
                                l.append(a)
                                
                            if b not in l:
                                l.append(b)
            ############################################################

            l = list(np.unique(l + lb))

    # Contabilizamos FN, FP, TP, TN
    SDnP += len([a for a in aSD if a not in l])
    PnSD += len([a for a in l if a not in aSD])
    SDsP += len([a for a in aSD if a in l])
    nSDnL += len([a for a in aP if a not in aSD and a not in l])

# Imprimimos por pantalla los datos
print('no encontradas: ' + str(SDnnP))
print('FN: ' + str(SDnP))
print('FP: ' + str(PnSD))
print('TP: ' + str(SDsP))
print('TN: ' + str(nSDnL))

print('no anotados con anotación: ' +
      str(len([g for g in genesNoAnotados if g in anotacionesTomato or g in anotacionesThaliana])))'''
#######################################################################################################

########################## Para escribir los no anotados ##############################################
ficheroSalidaScript = open(nameFicS, 'w+', buffering=1)

# Utilizamos solo los genes no anotados
# Para cada gen extraemos sus anotaciones
for gen in genesNoAnotados:
    aTo = []
    aTa = []

    if gen in anotacionesTomato:
        aTo = anotacionesTomato[gen]

    if gen in anotacionesThaliana:
        aTa = anotacionesThaliana[gen]

    aP = list(np.unique(aTo + aTa))

    l = []

    if len(aP) > 0:
        # Si solo encuentra anotaciones en un genoma, las copiamos
        if len(aTo) == 0 or len(aTa) == 0:
            l = aP
        # Si encuentra en ambos genomas, copiamos las que estén en los dos
        else:
            lb = [a for a in aP if a in aTo and a in aTa]
            aPb = [a for a in aP if (a not in aTo or a not in aTa) and a in godag]

            ######### Sección para obtener anotaciones mediante el método Ayllón ############
            # Comentar para evaluar el método conservador
            # En esta sección comprobamos si una anotación tiene descendientes entre las propuestas y de no ser así
            # la copia también
            '''for a in aPb:
                for b in aPb:
                    if a != b and len(common_parent_go_ids([a, b], godag)) != 0:
                        go_root = deepest_common_ancestor([a, b], godag)
                        n = semantic_distance(a, b, godag)

                        if go_root == a or go_root == b:
                            if go_root == a:
                                if b not in l:
                                    l.append(b)

                                if a in l:
                                    l.remove(a)
                            else:
                                if a not in l:
                                    l.append(a)

                                if b in l:
                                    l.remove(b)
                        else:
                            if a not in l:
                                l.append(a)

                            if b not in l:
                                l.append(b)'''
            #################################################################################

            l = list(np.unique(l + lb))

        # Cuando acabamos de seleccionar las anotaciones para un gen, las escribimos
        for e in l:
            ficheroSalidaScript.write(gen + '\t' + e + '\n')

# Escribimos también las que ya teníamos para obtener un fichero completo
# Se puede comentar esta parte para obtener un fichero con solo los no anotados para los que obtenemos anotación
for gen in anotacionesThaliana:
    l = anotacionesThaliana[gen]

    for a in l:
        ficheroSalidaScript.write(gen + '\t' + a + '\n')
#######################################################################################################
