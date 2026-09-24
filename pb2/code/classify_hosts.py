#!/usr/bin/env python3
"""Classify hosts: Mammalia (natural) / lab-mammal (sensitivity only) / Aves / exclude.
Dictionary-first with strain-name fallback; emits audit TSV of every distinct host string."""
import csv, re, sys
from collections import Counter

MAMMAL_NATURAL = [r'homo sapiens', r'\bbos\b', r'\bbos taurus\b', r'cattle', r'\bcow\b', r'dairy', r'bovine',
    r'felis', r'\bcat\b', r'feline', r'lynx', r'puma', r'panthera', r'leopard', r'tiger', r'\blion\b', r'acinonyx',
    r'canis', r'\bdog\b', r'canine', r'vulpes', r'\bfox\b', r'alopex', r'nyctereutes', r'raccoon dog', r'coyote', r'wolf',
    r'neovison', r'\bmink\b', r'mustela(?! putorius furo)', r'marten', r'martes', r'fisher', r'gulo', r'wolverine',
    r'seal', r'phoca', r'halichoerus', r'mirounga', r'pusa', r'erignathus', r'hydrurga', r'leptonychotes', r'lobodon',
    r'otaria', r'sea lion', r'zalophus', r'eumetopias', r'arctocephalus', r'fur seal', r'callorhinus', r'walrus', r'odobenus',
    r'tursiops', r'dolphin', r'phocoena', r'porpoise', r'whale', r'balaenoptera', r'physeter', r'delphin', r'globicephala',
    r'procyon', r'raccoon\b', r'mephitis', r'skunk', r'spilogale', r'conepatus', r'nasua', r'coati',
    r'ursus', r'\bbear\b', r'cougar', r'mountain lion', r'puma concolor', r'bovidae', r'bobcat', r'pinniped', r'serval', r'caracal', r'ocelot', r'margay', r'jaguarundi', r'herpestes', r'mongoose', r'paradoxurus', r'civet', r'genetta', r'viverr',
    r'sus scrofa', r'\bswine\b', r'\bpig\b', r'capra', r'\bgoat\b', r'ovis', r'sheep', r'equus', r'horse',
    r'odocoileus', r'\bdeer\b', r'capreolus', r'rangifer', r'reindeer', r'alces', r'moose', r'antilocapra', r'bison',
    r'lama', r'alpaca', r'vicugna', r'camel', r'oryctolagus', r'rabbit', r'lepus', r'hare',
    r'pteropus', r'\bbat\b', r'myotis', r'rhinolophus', r'rousettus', r'eptesicus',
    r'macaca', r'primate', r'gorilla', r'pteronura', r'enhydra', r'lontra', r'lutra', r'otter', r'aonyx',
    r'meles\b', r'badger', r'melogale', r'taxidea', r'eira', r'galictis', r'poecilogale', r'ictonyx']
MAMMAL_LAB = [r'mus musculus', r'\bmouse\b', r'\bmice\b', r'mesocricetus', r'hamster',
    r'mustela putorius furo', r'\bferret\b', r'cavia', r'guinea pig', r'rattus', r'\brat\b']
AVES = [r'chicken', r'gallus', r'turkey', r'meleagris', r'\bduck\b', r'anas', r'mallard', r'pintail', r'teal', r'wigeon',
    r'shoveler', r'gadwall', r'pochard', r'scaup', r'eider', r'merganser', r'shelduck', r'goose', r'\bgeese\b', r'anser',
    r'branta', r'swan', r'cygnus', r'gull', r'larus', r'tern', r'sterna', r'skua', r'crow', r'corvus', r'raven', r'magpie',
    r'starling', r'sturnus', r'sparrow', r'passer', r'finch', r'robin', r'thrush', r'turdus', r'warbler', r'swift',
    r'swallow', r'pigeon', r'columba', r'dove', r'quail', r'pheasant', r'partridge', r'guinea fowl', r'numida',
    r'eagle', r'hawk', r'falcon', r'kestrel', r'vulture', r'condor', r'buzzard', r'kite', r'harrier', r'osprey', r'accipiter',
    r'owl', r'strix', r'bubo', r'tyto', r'penguin', r'spheniscus', r'pelican', r'cormorant', r'phalacrocorax', r'heron',
    r'egret', r'ardea', r'stork', r'ibis', r'spoonbill', r'flamingo', r'crane', r'\bgrus\b', r'rail', r'coot', r'fulica',
    r'sandpiper', r'plover', r'snipe', r'curlew', r'godwit', r'shintail', r'avocet', r'stilt', r'oystercatcher', r'lapwing',
    r'parrot', r'cockatoo', r'macaw', r'budgerigar', r'ostrich', r'emu', r'rhea', r'pheasant', r'toucan', r'hornbill',
    r'woodpecker', r'kingfisher', r'bee-eater', r'roller', r'hoopoe', r'myna', r'jay', r'nutcracker', r'waxwing',
    r'bird', r'avian', r'wild bird', r'poultry', r'broiler', r'layer', r'breeder', r'backyard', r'chukar', r'guillemot',
    r'razorbill', r'puffin', r'auk', r'petrel', r'shearwater', r'albatross', r'gannet', r'morus', r'sula', r'frigatebird',
    r'tropicbird', r'grebe', r'loon', r'gavia', r'ptarmigan', r'grouse', r'turkey vulture', r'cathartes', r'caracara',
    r'coragyps', r'vulture', r'sanderling', r'murre', r'uria', r'haliaeetus', r'bufflehead', r'bucephala', r'pelecan',
    r'spatula', r'chroicocephalus', r'grackle', r'quiscalus', r'calidris', r'dunlin', r'kittiwake', r'rissa', r'aythya',
    r'somateria', r'redhead', r'podiceps', r'coscoroba', r'chenonetta', r'stercorarius', r'leucophaeus', r'leucocarbo',
    r'widgeon', r'wigeon', r'aix ', r'brant', r'scoter', r'melanitta', r'thalasseus', r'buteo', r'goshawk', r'gypaetus',
    r'circus', r'elanus', r'ictinia', r'leucophaeus', r'onychoprion', r'gelochelidon', r'hydroprogne', r'chroicocephalus', r'goldeneye', r'cairina', r'muscovy', r'fulmar', r'mareca', r'gadwell', r'peacock', r'pavo', r'willet', r'tringa', r'procellaria', r'merlin', r'falco', r'mergus', r'puffinus', r'ardenna', r'oceanodroma', r'ceryle', r'megaceryle', r'sturnella', r'agelaius', r'euphagus', r'molothrus', r'icterus', r'piranga', r'cardinalis', r'pheucticus', r'passerina', r'spiza', r'zenaida', r'streptopelia', r'geopelia', r'ptychoramphus', r'cerorhinca', r'fratercula', r'cephphus', r'brachyramphus', r'synthliboramphus', r'aethia', r'alle\b', r'larosterna', r'chilidonias', r'phaetusa', r'rhodostethia', r'pagophila', r'xema', r'creagrus', r'leucophaeus', r'dromaius', r'\bemu\b', r'skimmer', r'rynchops', r'jaeger', r'raptor']
EXCLUDE = [r'house fly', r'musca domestica', r'insect', r'pet food', r'^env$', r'^u$', r'^environment$', r'environment', r'water', r'soil', r'sediment', r'swab', r'feces', r'faeces', r'milk(?!.*cow)', r'unknown',
    r'\bn/?a\b', r'laboratory', r'vaccine', r'reference']

def classify(host, strain):
    h = (host or '').lower().strip()
    s = (strain or '').lower()
    # strain fallback: A/<host>/...
    st = ''
    m = re.match(r'^a/([^/]+)/', s)
    if m: st = m.group(1)
    text = h if h else st
    for p in EXCLUDE:
        if re.search(p, text): return 'exclude'
    for p in MAMMAL_LAB:
        if re.search(p, text): return 'mammal_lab'
    for p in MAMMAL_NATURAL:
        if re.search(p, text): return 'mammal'
    for p in AVES:
        if re.search(p, text): return 'avian'
    if h: return 'unresolved_host'
    if st: return 'unresolved_strain'
    return 'exclude'

if __name__ == '__main__':
    files = sys.argv[1:-1]
    out_audit = sys.argv[-1]
    ctr = Counter(); examples = {}
    for f in files:
        for row in csv.DictReader(open(f), delimiter='\t'):
            key = (row['host'] or f"strain:{row['strain'].split('/')[1] if '/' in row['strain'] else '?'}")
            cls = classify(row['host'], row['strain'])
            ctr[(key, cls)] += 1
            examples.setdefault((key, cls), (row['host'], row['strain'], row['accession']))
    w = csv.writer(open(out_audit, 'w'), delimiter='\t')
    w.writerow(['host_key', 'class', 'n_cds', 'example_host', 'example_strain', 'example_acc'])
    for (k, c), n in sorted(ctr.items(), key=lambda x: (x[0][1], -x[1])):
        eh, es, ea = examples[(k, c)]
        w.writerow([k, c, n, eh, es, ea])
    summ = Counter()
    for (k, c), n in ctr.items(): summ[c] += n
    print(dict(summ))
