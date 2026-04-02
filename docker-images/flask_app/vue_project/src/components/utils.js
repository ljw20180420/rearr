const validTargets = {
    'rawData': ['removeDuplicates'],
    'removeDuplicates': ['noDuplicates'],
    'noDuplicates': ['demultiplex'],
    'demultiplexAuxiliary': ['demultiplex'],
    'markers': ['buildMarker'],
    'buildMarker': ['demultiplexAuxiliary'],
    'demultiplex': ['noMix'],
    'noMix': ['postProcess'],
    'postProcess': ['toAlign'],
    'toAlign': ['rearrange'],
    'reference': ['defaultDirection', 'rearrange'],
    'defaultDirection': ['direction_file'],
    'direction_file': ['rearrange'],
    'rearrange': ['alignments'],
    'alignments': [],
}

export {
    validTargets
}
