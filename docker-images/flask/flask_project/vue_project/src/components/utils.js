function validTargets(source) {
    switch(source) {
        case 'target file':
        case 'pair file':
            return ['removeDuplicates'];
        case 'removeDuplicates':
            return ['file without duplicates'];
        case 'file without duplicates':
        case 'target spliter index':
        case 'pair spliter index':
            return ['demultiplex'];
        case 'target spliter':
            return ['buildSpliter/target spliter'];
        case 'buildSpliter/target spliter':
            return ['target spliter index'];
        case 'pair spliter':
            return ['buildSpliter/pair spliter'];
        case 'buildSpliter/pair spliter':
            return ['pair spliter index'];
        case 'demultiplex':
            return ['demultiplex file'];
        case 'demultiplex file':
            return ['sxPostProcess'];
        case 'sxPostProcess':
            return ['file of reads to align'];
        case 'file of reads to align':
        case 'file of reference':
            return ['rearrange'];
        case 'rearrange':
            return ['alignments'];
        case 'alignments':
            return [];
        case 'genome':
            return ['indexGenome', 'getReference'];
        case 'indexGenome':
            return ['genome index'];
        case 'genome index':
            return ['getReference'];
        case 'csvfile':
            return ['getSpliters', 'getReference'];
        case 'getReference':
            return ['file of reference'];
        case 'getSpliters':
            return ['target spliter', 'pair spliter'];
    }
}

function isValidConnection(connection) {
    return validTargets(connection.source).includes(connection.target);
}

function initRunJobNode(id, x, y, title = "run job node") {
    return {
        id: id,
        type: 'runJobNode',
        active: false,
        position: {
            x: x,
            y: y
        },
        title: title
    }
}

function initDataTankNode(id, x, y, data, title = "data tank node") {
    return {
        id: id,
        type: 'dataTankNode',
        active: false,
        position: {
            x: x,
            y: y
        },
        data: data.map(initDataTankNodeData),
        title: title
    }
}

function initDataTankNodeData({type, name, value, options}) {
    switch(type) {
        case 'value':
            return {
                type: type,
                name: name,
                value: value
            };
        case 'select':
            return {
                type: type,
                name: name,
                value: value,
                options: options
            };
        case 'file':
            return {
                type: type,
                taskId: null,
                name: name,
                value: [null]
            };
        case 'files':
            return {
                type: type,
                taskId: null,
                name: name,
                value: new Array(value).fill(null)
            };
    }
}

export {
    validTargets,
    isValidConnection,
    initRunJobNode,
    initDataTankNode
}
