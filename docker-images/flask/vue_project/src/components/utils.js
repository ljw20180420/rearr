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
    'reference': ['defaultCorrect', 'rearrange'],
    'defaultCorrect': ['correct'],
    'correct': ['rearrange'],
    'rearrange': ['alignments'],
    'alignments': [],
    'genome': ['indexGenome', 'sxGetReference'],
    'indexGenome': ['genomeIndex'],
    'genomeIndex': ['sxGetReference'],
    'csvfile': ['sxGetMarkers', 'sxGetReference'],
    'sxGetReference': ['reference'],
    'sxGetMarkers': ['markers']
}

export {
    validTargets
}
