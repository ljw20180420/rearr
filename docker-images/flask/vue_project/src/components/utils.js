const validTargets = {
    'rawData': ['removeDuplicates'],
    'removeDuplicates': ['noDuplicates'],
    'noDuplicates': ['demultiplex'],
    'demultiplexAuxiliary': ['demultiplex'],
    'spliters': ['buildSpliter'],
    'buildSpliter': ['demultiplexAuxiliary'],
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
    'csvfile': ['sxGetSpliters', 'sxGetReference'],
    'sxGetReference': ['reference'],
    'sxGetSpliters': ['spliters']
}

export {
    validTargets
}
