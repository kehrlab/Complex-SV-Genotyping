from Bio import SeqIO
from Bio.Seq import Seq
import json
import sys
import gc
import copy

def createEndNode(referenceIdx, junction, direction, margin):
    node = {}
    if direction < 0:
        chrom = junction['rNameLeft']
        start = junction['xLeft'] - margin 
        end = junction['xLeft'] 
    elif direction > 0:
        chrom = junction['rNameRight']
        start = junction['xRight'] 
        end = junction['xRight'] + margin 
    node['reference'] = str(chrom) + ':' + str(start) + '-' + str(end)
    node['name'] = 'ref-' + node['reference']
    return node

def createMainRefNode(referenceIdx, firstJunction, lastJunction):
    node = {}
    if lastJunction['rNameRight'] != firstJunction['rNameLeft']:
        return None

    if firstJunction['xLeft'] + 1 > lastJunction['xRight'] - 1:
        return None

    chrom = firstJunction['rNameLeft']
    start = firstJunction['xLeft'] + 1
    end = lastJunction['xRight'] - 1
    
    node['reference'] = str(chrom) + ':' + str(start) + '-' + str(end)
    node['name'] = 'ref-' + node['reference']
    return node

def createMiddleNode(referenceIdx, leftJunction, rightJunction):
    node = {}
    chrom = leftJunction['rNameRight']
    start = leftJunction['xRight']
    end = rightJunction['xLeft']

    if end == (start + 1):
        return None
    
    if start > end:
        node = createMiddleNodeSeq(referenceIdx, leftJunction, rightJunction)
    else:
        node = createMiddleNodeRef(referenceIdx, leftJunction, rightJunction)

    return node

def createMiddleNodeRef(referenceIdx, leftJunction, rightJunction):
    node = {}
    chrom = leftJunction['rNameRight']
    start = leftJunction['xRight']
    end = rightJunction['xLeft']

    if end == (start + 1):
        return None
    #sequence = referenceIdx[chrom][start:(end + 1)]

    node['reference'] = str(chrom) + ':' + str(start) + '-' + str(end)
    node['name'] = 'ref-' + node['reference']

    return node

def createMiddleNodeSeq(referenceIdx, leftJunction, rightJunction):
    node = {}
    chrom = leftJunction['rNameRight']
    start = min(leftJunction['xRight'], rightJunction['xLeft']) - 1
    end = max(leftJunction['xRight'], rightJunction['xLeft']) - 1
    name = 'var-' + str(chrom) + ':' + str(start + 1) + '-' + str(end + 1)

    sequence = referenceIdx[chrom][start:(end + 1)]
    if leftJunction['directionRight'] == 'left' and rightJunction['directionLeft'] == 'right':
        sequence = sequence.reverse_complement()
        name = 'var-' + str(chrom) + ':' + str(end + 1) + '-' + str(start + 1)
    sequence = str(sequence.seq)

    node['name'] = name
    node['sequence'] = sequence

    return node

def adjustJunctions(structure, chromosome):
    if len(structure) != 1:
        return structure
    k = list(structure.keys())[0]

    if (structure[k]["rNameLeft"] != chromosome or structure[k]["rNameRight"] != chromosome):
        return None
    if (structure[k]["xLeft"] <= structure[k]["xRight"]):
        return structure
    
    structure[str(int(k) + 1)] = copy.deepcopy(structure[k])
    structure[str(int(k) + 2)] = copy.deepcopy(structure[k])
    structure[k]["xLeft"] = structure[k]["xRight"] - 1
    structure[str(int(k) + 2)]["xRight"] = structure[str(int(k) + 1)]["xLeft"] + 1
    return structure

if len(sys.argv) != 4:
    print("Usage: python convert_to_graphs.py <variantFile> <referenceGenome> <outputDir>")
    exit(1)

variantFile = sys.argv[1]
referenceFile = sys.argv[2]
outputDir = sys.argv[3]

print("Create reference index...")
referenceIndex = SeqIO.index(referenceFile, "fasta")

with open(variantFile) as jsonFile:
    variantData = json.load(jsonFile)

counter = 0
for variant in variantData:
    variantDescription = variantData[variant]['VAR']
    if len(variantDescription) > 1:
        print("Skipping ", variant, " due to to more than one involved chromosome")
        continue
    
    # for the progress bar...
    counter += 1
    percent = counter / len(variantData) * 100

    nodes = []
    edges = []
    paths = []
    target_regions = []

    sequences = ['ALT', 'REF']
    model_name = variant
    
    
    chrom = list(variantDescription.keys())[0]
    junctions = variantDescription[chrom]
    junctions = adjustJunctions(junctions, chrom)
    if junctions == None:
        continue

    junctionKeys = [int(k) for k in junctions.keys()]
    junctionKeys.sort()

    ################################## nodes ##################################
    
    nodes.append({
        'name' : 'source', 
        'sequence' : 'NNNNNNNNNN'
        })
    nodes.append(createEndNode(referenceIndex, junctions[str(junctionKeys[0])], -1, 150))

    for i in range(len(junctionKeys) - 1):
        leftJunction = junctions[str(junctionKeys[i])]
        rightJunction = junctions[str(junctionKeys[i+1])]
        nodes.append(createMiddleNode(referenceIndex, leftJunction, rightJunction))

    tempNode = createMainRefNode(referenceIndex, junctions[str(junctionKeys[0])], junctions[str(junctionKeys[-1])])
    if tempNode != None:
        nodes.append(tempNode)

    nodes.append(createEndNode(referenceIndex, junctions[str(junctionKeys[-1])], 1, 150))
    nodes.append({
        'name' : 'sink', 
        'sequence' : 'NNNNNNNNNN'
        })

    # make sure node names are unique
    for n in range(1, len(nodes) - 1):
        nodes[n]['name'] += ('_' + str(n))

    ################################## edges ##################################

    # add source edge
    edges.append({
        'from' : 'source', 
        'to' : nodes[1]['name'], 
        'name' : 'source_' + nodes[1]['name']
        })

    # add reference edges
    if tempNode != None:
        edges.append({
            'from' : nodes[1]['name'], 
            'to' : nodes[-3]['name'], 
            'name' : nodes[1]['name'] + '_' + nodes[-3]['name'], 
            'sequences' : ['REF']
            })
        edges.append({
            'from' : nodes[-3]['name'], 
            'to' : nodes[-2]['name'], 
            'name' : nodes[-3]['name'] + '_' + nodes[-2]['name'], 
            'sequences' : ['REF']
            })
    else:
        edges.append({
            'from' : nodes[1]['name'], 
            'to' : nodes[-2]['name'], 
            'name' : nodes[1]['name'] + '_' + nodes[-2]['name'], 
            'sequences' : ['REF']
            })


    # add alt edges
    if tempNode != None:
        maxIdx = len(nodes) - 3
    else:
        maxIdx = len(nodes) - 2
    if len(junctionKeys) > 1:
        for i in range(1, maxIdx - 1):
            edges.append({
                'from' : nodes[i]['name'], 
                'to' : nodes[i+1]['name'], 
                'name' : nodes[i]['name'] + '_' + nodes[i+1]['name'], 
                'sequences' : ['ALT']
                })
    edges.append({
    'from' : nodes[maxIdx - 1]['name'], 
    'to' : nodes[-2]['name'], 
    'name' : nodes[maxIdx - 1]['name'] + '_' + nodes[-2]['name'], 
    'sequences' : ['ALT']
    })

    # add sink edge
    edges.append({
        'to' : 'sink', 
        'from' : nodes[-2]['name'], 
        'name' : nodes[-2]['name'] + '_sink'
        })

    ##########################################################################

    # add paths
    if tempNode != None:
        paths.append({
            'nodes' : [nodes[1]['name'], nodes[-3]['name'], nodes[-2]['name']],
            'path_id' : 'REF|1',
            'sequence' : 'REF'
        })
    else:
        paths.append({
            'nodes' : [nodes[1]['name'], nodes[-2]['name']],
            'path_id' : 'REF|1',
            'sequence' : 'REF'
        })

    altPath = [nodes[i]['name'] for i in range(1, maxIdx)]
    altPath.append(nodes[-2]['name'])
    paths.append({
       'nodes' : altPath,
       'path_id' : 'ALT|1',
       'sequence' : 'ALT'
    })
    
    # add target region
    r = junctions[str(junctionKeys[0])]['rNameLeft'] + ':'
    r += str(junctions[str(junctionKeys[0])]['xLeft'] - 150) + '-'
    r += str(junctions[str(junctionKeys[-1])]['xRight'] + 150)
    target_regions.append(r)

    ################## merge identical nodes #################################
    # i = 0
    # while i < len(nodes):
    #     match = False
    #     for j in range(i+1, len(nodes)):
    #         if nodes[j]['name'] == nodes[i]['name']:
    #             if 'sequence' in nodes[i].keys() and 'sequence' in nodes[j].keys() and nodes[i]['sequence'] == nodes[j]['sequence']:
    #                 match = True
    #                 nodes.pop(j)
    #                 break
    #     if not match:
    #         i += 1

    ################################ write ###################################
    graph = {
        'edges' : edges,
        'model_name' : model_name,
        'nodes' : nodes,
        'paths' : paths,
        'sequencenames' : sequences,
        'target_regions' : target_regions
    }
    filename = outputDir + '/' + variant + '_graph.json'
    with open(filename, 'w') as jsonOut:
        json.dump(graph, jsonOut, indent = 4)

    print("Progress: ", int(percent), "%\t(",variant,")               ", end = "\r", sep = "")
    gc.collect()
