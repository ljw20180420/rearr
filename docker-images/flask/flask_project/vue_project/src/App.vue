<script setup>
import { ref, markRaw } from 'vue';
import { VueFlow, useVueFlow } from '@vue-flow/core';
import { MiniMap } from '@vue-flow/minimap'
import { Controls } from '@vue-flow/controls'
import { Background } from '@vue-flow/background'
import dataTank from './components/dataTank.vue';
import runJob from './components/runJob.vue';
import dataTunnel from './components/dataTunnel.vue';
import { validTargets } from './components/utils.js';

function initRunJobNode(id, x, y, title) {
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

function initDataTankNode(id, x, y, title, data) {
    return {
        id: id,
        data: data,
        type: 'dataTankNode',
        active: false,
        position: {
            x: x,
            y: y
        },
        title: title
    }
}

const base_url = import.meta.env.BASE_URL;

const nodeTypes = {
  dataTankNode: markRaw(dataTank),
  runJobNode: markRaw(runJob)
};

const edgeTypes = {
  dataTunnelEdge: markRaw(dataTunnel)
};

const nodes = ref([
  initDataTankNode(
    'csvfile', 0, 0,
    "Shi Xing's .csv file contains hints to extract reference from genome. See sx/sxExtractSpliter.md.",
    {
      '.csv file': {
        type: 'file'
      }
    }
  ),
  initRunJobNode(
    'sxGetReference', 2100, 0,
    "Extract reference from genome based on .csv file in the same format as Shi Xing. See sx/getSxCsvFileRef."
  ),
  initDataTankNode(
    'reference', 2400, 0,
    'Reads in .post file are aligned to these references. See core/Rearrangement/rearr.md.',
    {
      '.ref file': {
        type: 'file'
      }
    }
  ),
  initRunJobNode(
    'rearrange', 3300, 500,
    'Align reads in .post file to references.'
  ),
  initDataTankNode(
    'alignments', 3600, 500,
    'Chimeric alignment results.',
    {
      '.alg file': {
        type: 'file'
      }
    }
  ),
  initRunJobNode(
    'defaultCorrect', 2700, 0,
    'Generate default .correct file with all up.'
  ),
  initDataTankNode(
    'correct', 3000, 0,
    'Correct chimeric alignments. For each junction between two adjacent references, the alignment try to generate blunt end for either upstream or downstream. A default .correct file with all up can be generated from .ref file for convenience. See core/Rearrangement/correct_micro_homology.awk.',
    {
      '.correct file': {
        type: 'file'
      }
    }
  ),
  initDataTankNode(
    'genome', 1200, 50,
    'Genome file to extract references from.',
    {
      'genome file': {
        type: 'file'
      }
    }
  ),
  initRunJobNode(
    'indexGenome', 1500, 50,
    'Build bowtie2 index for genome.'
  ),
  initDataTankNode(
    'genomeIndex', 1800, 50,
    'Bowtie2 index of genome. Extension base numbers are feed to sxGetReference together with genome index to extract references. Rearr need extended reference to catch templated insertion.',
    {
      'extensions': {
        'cut1 upstream': {
          type: 'value',
          value: 50
        },
        'cut1 downstream': {
          type: 'value',
          value: 0
        },
        'cut2 upstream': {
          type: 'value',
          value: 10
        },
        'cut1 downstream': {
          type: 'value',
          value: 100
        }
      },
      'genome index': {
        '1.bt2': {
          type: 'file'
        },
        '2.bt2': {
          type: 'file'
        },
        '3.bt2': {
          type: 'file'
        },
        '4.bt2': {
          type: 'file'
        },
        'rev.1.bt2': {
          type: 'file'
        },
        'rev.2.bt2': {
          type: 'file'
        }
      }
    }
  ),
  initDataTankNode(
    'rawData', 600, 1000,
    'Files containing rawdata. See core/removeDuplicates.md.',
    {
      '.fastq files': [
        {
          type: 'file'
        },
        {
          type: 'file'
        }
      ]
    }
  ),
  initRunJobNode(
    'removeDuplicates', 900, 1000,
    'Remove duplicated reads in rawData fastq files. See core/removeDuplicates.md.'
  ),
  initDataTankNode(
    'noDuplicates', 1200, 1000,
    'File of side-by-side reads in rawData fastq files without duplicates. See core/removeDuplicates.md.',
    {
      '.noDup file': {
        type: 'file'
      }
    }
  ),
  initRunJobNode(
    'sxGetSpliters', 300, 1500,
    "Extract spliter from .csv file in the same format as Xing Shi. See sx/sxExtractSpliter.md."
  ),
  initDataTankNode(
    'spliters', 600, 1500,
    'Reads for demultiplexing .noDup file. Spliters maps to rawData fastq files bijectively. See core/demultiplex.md.',
    {
      '.fasta files': [
        {
          type: 'file'
        },
        {
          type: 'file'
        }
      ]
    }
  ),
  initRunJobNode(
    'buildSpliter', 900, 1500,
    'Build bowtie2 index for spliters.'
  ),
  initDataTankNode(
    'demultiplexAuxiliary', 1200, 1500,
    'Bowtie2 index of spliters. minScores are used to filter low-quality maps in demultiplex step. See core/demultiplex.md.',
    {
      'auxiliaries': [
        {
          'minScore': {
            type: "value",
            value: 30
          },
          'spliterIndex': {
            '1.bt2': {
              type: 'file'
            },
            '2.bt2': {
              type: 'file'
            },
            '3.bt2': {
              type: 'file'
            },
            '4.bt2': {
              type: 'file'
            },
            'rev.1.bt2': {
              type: 'file'
            },
            'rev.2.bt2': {
              type: 'file'
            }
          }
        },
        {
          'minScore': {
            type: "value",
            value: 100
          },
          'spliterIndex': {
            '1.bt2': {
              type: 'file'
            },
            '2.bt2': {
              type: 'file'
            },
            '3.bt2': {
              type: 'file'
            },
            '4.bt2': {
              type: 'file'
            },
            'rev.1.bt2': {
              type: 'file'
            },
            'rev.2.bt2': {
              type: 'file'
            }
          }
        }
      ]
    }
  ),
  initRunJobNode(
    'demultiplex', 1500, 1500,
    'Demultiplex .noDup reads by spliters. See core/demultiplex.md.'
  ),
  initDataTankNode(
    'noMix', 1800, 1500,
    "Demultiplexed .noDup reads. Minimal base number is feed to sxPostProcess together with .demultiplex file to filter too short reads after remove 5' spliter and 3' adapter. See core/demultiplex.md.",
    {
      'minimal base number': {
        type: 'value',
        value: 30
      },
      '.demultiplex file': {
        type: 'file'
      }
    }
  ),
  initRunJobNode(
    'sxPostProcess', 2100, 1500,
    "Remove 5' spliter and 3' adapter, and filter reads which are too short for alignment. sxPostProcess only supports data in the same format as Xing Shi. See sx/sxCutR2AdapterFilterCumulate."
  ),
  initDataTankNode(
    'toAlign', 2400, 1500,
    'Clean target reads ready for alignment. Alignment score settings are based on affine gap penalty. See rearrangment --help for details about s0, s1, s2, u, v, ru, rv, qu, qv. See core/Rearrangement/rearr.md.',
    {
      'align scores': {
        's0': {
          type: 'value',
          value: -6
        },
        's1': {
          type: 'value',
          value: 4
        },
        's2': {
          type: 'value',
          value: 2
        },
        'u': {
          type: 'value',
          value: -3
        },
        'v': {
          type: 'value',
          value: -9
        },
        'ru': {
          type: 'value',
          value: 0
        },
        'rv': {
          type: 'value',
          value: 0
        },
        'qu': {
          type: 'value',
          value: 0
        },
        'qv': {
          type: 'value',
          value: -5
        }
      },
      '.post file': {
        type: 'file'
      }
    }
  ),
]);

let source_target_pairs = [];
for (const node of nodes.value) {
  const source = node.id;
  for (const target of validTargets[source]) {
    source_target_pairs.push(
      {
        id: `${source}-${target}`,
        source: source,
        target: target,
        type: 'dataTunnelEdge'
      }
    )
  }
};
const edges = ref(source_target_pairs);

const { onConnect, addEdges } = useVueFlow();
onConnect(param => {
  debugger;
  param.id = `${param.source}-${param.target}`;
  param.type = 'dataTunnelEdge';
  addEdges(param);
})

window.onbeforeunload = () => fetch(`${base_url}/stop`);
</script>

<template>
  <div style="height: 900px">
    <VueFlow
      v-model:nodes="nodes"
      v-model:edges="edges"
      :node-types="nodeTypes"
      :edge-types="edgeTypes"
    >
      <MiniMap pannable zoomable nodeColor="black" maskColor="rgb(0, 0, 0, 0.7)" />
      <Controls />
      <Background variant="dots" gap=100 size=4 patternColor="black" />
    </VueFlow>
  </div>
</template>

<style>
/* import the necessary styles for Vue Flow to work */
@import '@vue-flow/core/dist/style.css';
/* import the default theme, this is optional but generally recommended */
@import '@vue-flow/core/dist/theme-default.css';
/* import minimap style */
@import '@vue-flow/minimap/dist/style.css';
/* import controls style */
@import '@vue-flow/controls/dist/style.css';

.vue-flow__node {
    background: black;
    color: white;
    border: 2px solid rgb(64, 64, 64);
    border-radius: 2px;
    box-shadow: 0 0 0 2px rgb(128, 128, 128);
    padding: 0px;
    white-space: nowrap;
}

a:link {
  color: red;
  background-color: transparent;
  text-decoration: none;
}

a:visited {
  color: orange;
  background-color: transparent;
  text-decoration: none;
}
</style>
