const validTargets = {
    'rawData': ['removeDuplicates'],
    'removeDuplicates': ['noDuplicates'],
    'noDuplicates': ['demultiplex'],
    'demultiplexAuxiliary': ['demultiplex'],
    'markers': ['buildMarker'],
    'buildMarker': ['demultiplexAuxiliary'],
    'demultiplex': ['noMix'],
    'noMix': ['sxPostProcess'],
    'sxPostProcess': ['toAlign'],
    'toAlign': ['rearrange'],
    'reference': ['defaultDirection', 'rearrange'],
    'defaultDirection': ['direction_file'],
    'direction_file': ['rearrange'],
    'rearrange': ['alignments'],
    'alignments': [],
    'genome': ['indexGenome', 'sxGetReference'],
    'indexGenome': ['genomeIndex'],
    'genomeIndex': ['sxGetReference'],
    'plasmid': ['sxGetMarkers', 'sxGetReference'],
    'sxGetReference': ['reference'],
    'sxGetMarkers': ['markers']
}

export {
    validTargets
}
